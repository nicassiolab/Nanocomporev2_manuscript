using CSV
using SQLite
using DataFrames
using Base.Threads
using StatsBase
using Peaks
using HypothesisTests
using Random
using ArgParse


function parse_commandline()
    args = ArgParseSettings()
    @add_arg_table! args begin
        "--db", "-d"
            help = "Nanocompore database with read level results."
			required = true
        "--results_tsv", "-r"
            help = "Nanocompore result TSV file."
			required = true
        "--mod_samples", "-m"
            help = "The labels of the modified samples (as comma separated values)."
			required = true
        "--threshold_prob", "-t"
            help = "The minimum probability to consider a position modified [0, 100]."
    		arg_type = Int
    		default = 50
        "--permutations", "-p"
            help = "Number of permutations to test."
    		arg_type = Int
    		default = 10_000
        "--transcripts"
            help = "A file with one transcript name per line. Will only process the listed transcripts."
    end
	parse_args(args)
end


function contingency(a::Vector{Int64}, b::Vector{Int64})::Matrix{Integer}
	[sum(a .== 0 .&& b .== 0) sum(a .== 0 .&& b .== 1)
	 sum(a .== 1 .&& b .== 0) sum(a .== 1 .&& b .== 1)]
end


function main()
	args = parse_commandline()

	sample_list = join(map(s -> "\"$s\"", split(args["mod_samples"], ",")), ", ")

    results = DataFrame(CSV.File(args["results_tsv"], delim = '\t', select = ["ref_id", "pos", "GMM_chi2_qvalue", "GMM_LOR"]))
    conn = SQLite.DB(args["db"])
    
    
    sig_results = results[results.GMM_chi2_qvalue .!== missing .&& results.GMM_chi2_qvalue .< 0.01 .&& abs.(results.GMM_LOR) .> 0.8, :]
    
    multimod_refs = filter(pair -> pair[2] > 1, countmap(sig_results.ref_id)) |> keys |> Set
    
    transcripts = DBInterface.execute(conn, "SELECT id, name FROM transcripts") |> DataFrame
    transcripts = transcripts[in.(transcripts.name, Ref(multimod_refs)), :]

    if !isnothing(args["transcripts"])
        transcript_set = readlines(args["transcripts"]) |> Set
        transcripts = transcripts[in.(transcripts.name, Ref(transcript_set)), :]
    end

    print_lock = ReentrantLock()
    
    
    println(join(["reference", "pos1", "pos2", "chi2", "pvalue", "phi", "emp_pvalue", "residual00", "residual01", "residual10", "residual11", "min_expected", "fisher_pvalue", "a", "b", "c", "d"], "\t"))
    
    # transcript = first(eachrow(transcripts))
    Threads.@threads for transcript in eachrow(transcripts)
    	reads = DBInterface.execute(conn, "SELECT mod_probs FROM read_results WHERE transcript_id = $(transcript.id) AND sample IN ($sample_list)") |> DataFrame
    	mod_probs = hcat(map(probs -> reinterpret.(Int8, probs), reads.mod_probs)...) |> transpose
    	valid = mod_probs .>= 0
    	mods = Int.(mod_probs .> args["threshold_prob"])
    
    	ref_sig_results = sig_results[sig_results.ref_id .== transcript.name, :]
    
    	# show(ref_sig_results)
    
    	significance = zeros(maximum(ref_sig_results.pos) .+ 1)
    	significance[ref_sig_results.pos .+ 1] = -log10.(ref_sig_results.GMM_chi2_qvalue)
    
    	peaks = argmaxima(significance, 5, strict = false)
    
		output = []
    	for (i, pos1) in enumerate(peaks)	
    		for pos2 in peaks[i+1:end]
    
    			valid_pairs = valid[:, pos1] .&& valid[:, pos2]
    			
    			mod1 = mods[valid_pairs, pos1]
    			mod2 = mods[valid_pairs, pos2]
    
				# We add one to make sure we don't get counts of zero
				# which breaks the test
    			observed = contingency(mod1, mod2) .+ 1
    
    			res = ChisqTest(observed)
    			pval = pvalue(res)
    			chi2 = res.stat
    			phi = sqrt(chi2 / sum(valid_pairs))
				min_expected = minimum(res.expected)
				
				# fisher_pval = missing
				# if min_expected < 5
				fisher_pval = pvalue(FisherExactTest(observed'...))
				# end
    
    			emp_stats = []
    			for p in 1:args["permutations"]
    				shuffle!(mod2)			
    				
    				perm_res = ChisqTest(contingency(mod1, mod2) .+ 1)
    				push!(emp_stats, perm_res.stat)
    			end
    
    			emp_pval = (sum(emp_stats .>= chi2) + 1) / (args["permutations"] + 1)

				# We subtract one from the positions to revert them back to 0-based
				push!(output, [transcript.name, pos1-1, pos2-1, chi2, pval, phi, emp_pval, res.residuals[1, 1], res.residuals[1, 2], res.residuals[2, 1], res.residuals[2, 2], min_expected, fisher_pval, observed[1, 1], observed[1, 2], observed[2, 1], observed[2, 2]])
    			
    		end
    	end
		lock(print_lock) do
			for line in output
		    	println(join(string.(line), "\t"))
			end
		end
    end
end


main()

