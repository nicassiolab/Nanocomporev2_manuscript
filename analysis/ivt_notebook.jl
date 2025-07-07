### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 41818e42-3cb8-4d48-93d9-b53e0eea7873
using Pkg

# ╔═╡ 82dba028-3fbe-44e2-9e8f-377f738363a7
Pkg.add("Peaks")

# ╔═╡ 9fbe17fb-5b2a-4cda-a70d-1598de0a6ab0
Pkg.add("Clustering")

# ╔═╡ 692af519-1ac2-483c-86c5-f1801d8426bd
Pkg.add("AlgebraOfGraphics")

# ╔═╡ 9ba8c503-f59e-4d6a-b97c-9826279d744a
Pkg.add("GenomicAnnotations")

# ╔═╡ bcd1df63-e0b6-467b-9858-de5a15a61765
Pkg.add("StatsPlots")

# ╔═╡ 93db1a5d-4de2-482b-a4e2-540db2059689
begin
	using DataFrames
	using CSV
	using MLBase
	using CairoMakie
	using Random
	using StatsBase
	using Printf
	# using MultipleTesting
end

# ╔═╡ c88509a9-5a1d-4895-9ae7-5174cc199772
using Peaks

# ╔═╡ 76d3d1c4-0e01-4360-bc3d-bc274ff1c39c
using AlgebraOfGraphics

# ╔═╡ 164f8adc-3bba-11f0-3c64-19ee3bf9097e
include("lib.jl")

# ╔═╡ bc77ac69-225a-4f81-b219-15683bb5eb4a
# ╠═╡ disabled = true
#=╠═╡
function simple_peaks(df, radius)
	N = size(df, 1)
	keep = repeat([true], N)
	for i in 1:N
		for j in max(1, i - radius):min(N, i + radius)
			if df[i, :ref_id] == df[j, :ref_id] &&
			   abs(df[i, :pos] - df[j, :pos]) <= radius &&
			   df[i, :predicted] > df[j, :predicted]
				keep[j] = false
			end
		end
	end
	return df[keep, :]
end
  ╠═╡ =#

# ╔═╡ dedfd9b6-9b9e-4149-834a-547cbac99873
# ╠═╡ disabled = true
#=╠═╡
# function peaks(df, radius)
# 	all_positions = DataFrames.combine(groupby(df, :ref_id),
# 									   :pos => (p -> minimum(p):maximum(p)) => :pos)
# 	all_positions = leftjoin(all_positions, df,
# 			 				 on = [:ref_id, :pos])
# 	all_positions = sort(all_positions, [:ref_id, :pos])
# 	# all_positions[:, :predicted] = coalesce.(all_positions.predicted, 0)
# 	peak_indices = argmaxima(all_positions.predicted, radius, strict = false)
# 	all_positions[peak_indices, :]
# end
  ╠═╡ =#

# ╔═╡ 1795bc42-bdc8-4685-bb52-39271581bc06
function peaks(df, radius)
	all_positions = DataFrames.combine(groupby(df, :ref_id),
									   :pos => (p -> 1:maximum(p)) => :pos)
	all_positions = leftjoin(all_positions, df,
			 				 on = [:ref_id, :pos])
	all_positions = sort(all_positions, [:ref_id, :pos])

	peaks = DataFrames.combine(groupby(all_positions, :ref_id),
				   	 	       :predicted_raw => (p -> argmaxima(p, radius, strict = false)) => :pos)
	sort(leftjoin(peaks, df,
			 	  on = [:ref_id, :pos]),
		 [:ref_id, :pos])
end

# ╔═╡ 88dac2fb-7062-4e75-aaaf-c2596d90b299
function annotate_peaks(peaks, mod_ref, writer_motifs, writer_mods)
	local df = copy(peaks)
	df[:, :peak_id] = 1:nrow(df)
	df[:, :start] = max.(df.genomicPos .- 4, 0)
	# end is non-inclusive in BED files
	df[:, :end] = df.genomicPos .+ 5
	local ranges = df[:, [:peak_id, :chr, :strand]]
	ranges[:, :pos] = range.(df.start, df.end)
	local mod_overlaps = innerjoin(flatten(ranges, :pos), mod_ref,
			  					   on = [:chr, :strand, :pos])
	rename!(mod_overlaps, :pos => :mod_pos)
	annot_peaks = leftjoin(df, mod_overlaps,
			 			   on = [:peak_id, :chr, :strand])
	for (writer, motif) in writer_motifs
		matched = annot_peaks.mod .=== missing .&&
				  occursin.(motif, annot_peaks.ref_kmer)
		annot_peaks[matched, :mod] .= writer_mods[writer]
		annot_peaks[matched, :source] .= "motif"
	end
	annot_peaks
end

# ╔═╡ 5a2ad3e8-9f14-49ad-a78b-6fdd4be05220
# ╠═╡ disabled = true
#=╠═╡
function expand_pos(df, pos_col, left, right)
	N = left + right + 1
	expanded_df = repeat(df, N)
	pos_offsets = sort(repeat(-left:right, size(df, 1)))
	expanded_df[:, pos_col] += pos_offsets
	expanded_df
end
  ╠═╡ =#

# ╔═╡ 8445b562-cb03-4a7e-911b-5dc44778ee5c
# ╠═╡ disabled = true
#=╠═╡
begin
	glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/high_confidence_GLORI_sites.bed"))
	glori[!, "pos"] = glori.start .+ 2
	glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
	glori[!, "modified"] .= 1
	glori
end
  ╠═╡ =#

# ╔═╡ 8b930a4c-cef4-49e6-8a76-d59d4df8a349
	begin
	glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k.bed"))
	# glori[!, "pos"] = glori.start .+ 2
	glori[:, :pos] = glori.start
	glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
	glori[!, "modified"] .= 1
	glori
end

# ╔═╡ a6d02ee4-2e77-4520-b140-36fb2150b1b0
binned_glori = unique(glori[!, ["chrom", "strand", "bin", "modified"]])

# ╔═╡ b89d56f0-5104-4b0d-a1ae-6090d262db25
writer_motifs = Dict(
	# m1A
	"TRMT6-61A" => r"GTTC[AG]A",
	
	# pseudouridine
	"PUS4" => r"[AG][AG]TTC[AGCT]A", # https://doi.org/10.1038/s41587-021-00915-6
	"PUS7" => r"T[ACTG]TA[AG]",
	
	# m5C
	"NSUN6" => r"CTCC[AT]",
	
	# m3C
	"METTL8_1" => r"GCCTCC[AG]",
	"METTL8_2" => r"AGGC[AT]GG",
	
	# m7G
	"METTL1" => r"[AG]AGGT",
	
	# Nm
	"FTSJ3" => r"[ACTG]AGATC[AG][AG]", #[AG][AG][AG]",
	
	# ac4C
	# "NAT10" =>  r"C[ACTG][ACTG]C"
	
	# m6A
	"METTL3" => r"[GAT][AG]AC[ATC]"
)

# ╔═╡ 566437e3-1abe-49bf-a262-e9a4c702432e
writer_mods = Dict(
	"TRMT6-61A" => "m1A",
	"PUS4" => "Y",
	"PUS7" => "Y",
	"NSUN6" => "m5C",
	"METTL8_1" => "m3C",
	"METTL8_2" => "m3C",
	"METTL1" => "m7G",
	"FTSJ3" => "Nm",
	"NAT10" => "ac4C",
	"METTL3" => "m6A"
)

# ╔═╡ b5d8eafa-9ea4-4c59-bef8-3070126bb095
rmbase_psu = DataFrame(CSV.File("/home/mzdravkov/rmbase_mods/human.hg38.Pseudo.result.col29.bed",
				   header = ["chr", "start", "end", "mod_id", "score", "strand", "mod_type", "support_num", "support_list", "support_list_sub", "pub_list", "cell_list", "seq_type_list", "gene_id", "transcript_id", "gene_name", "gene_type", "region", "seq", "motif_score", "conserved_sites_list", "snoRNA_detail_info", "snoRNA_guide_site", "snoRNA_name_list", "snoRNA_database", "writer_loc", "writer_id", "writer_name_list", "source"]))

# ╔═╡ c9fabd63-2131-41fc-a50a-bda27319bc6d
begin
	atoi = DataFrame(CSV.File(
		"/home/mzdravkov/rmbase_mods/human.hg38.RNA-editing.result.col29.bed",
		delim='\t',
		header = vcat(["chr", "start", "end", "mod_id", "score", "strand", "mod_type"], repeat(["."], 22))))[:, 1:7]
	for col in ["support_num", "support_list", "support_list_sub", "pub_list", "cell_list"]
		atoi[:, col] .= "na"
	end
	atoi
end

# ╔═╡ 43a5d6b6-652d-403c-b4fd-945ba2b80ead
begin
	rmbase_all = copy(rmbase_psu[:, 1:12])
	for file in ["human.hg38.m1A.result.col29.bed",
				 "human.hg38.m5C.result.col29.bed",
				 "human.hg38.m7G.result.col29.bed",
				 "human.hg38.Nm.result.col29.bed",
				 "human.hg38.otherMod.result.col29.bed"]
				 # "human.hg38.RNA-editing.result.col29.bed"]
		local mods = DataFrame(CSV.read("/home/mzdravkov/rmbase_mods/$file",
										DataFrame,
										select = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
									    header = ["chr", "start", "end", "mod_id", "score", "strand", "mod_type", "support_num", "support_list", "support_list_sub", "pub_list", "cell_list"],
									    ntasks = 1))[:, 1:12]
		rmbase_all = vcat(rmbase_all, mods)
	end
	# rmbase_all = rmbase_all[occursin.("HEK293T", rmbase_all.cell_list), :]
	rmbase_all
end

# ╔═╡ 5abecd8c-e06c-477c-8938-545f64591385
hela_ac4C = DataFrame(CSV.File("/home/mzdravkov/GSE162043_RedaCTseq_ac4C_sites_GSE162043_RedaCTseq_ac4C_sites_hg38_lift.bed", header = [:chr, :start, :end, :name, :score, :strand]))

# ╔═╡ b0a1b932-7f0b-40ac-a2d6-e4867d2b54da
begin
	m7GDB_ngs = sort(DataFrame(CSV.File("/home/mzdravkov/base_human_filtered.csv")),
					 [:seqnames, :start])
	m7GDB_ngs.start .-= 1
end

# ╔═╡ 7ca021f1-532f-4a3c-98b6-0ff08c27d01e
m5C_atlas = DataFrame(CSV.File("/home/mzdravkov/hglift_m5C_to_grch38.bed",
							   header = [:chr, :start, :end, :name, :score, :strand]))

# ╔═╡ 33a21a05-dac1-40e0-97a4-8dae37fa8bfc
begin
	local glori_sel = glori[:, [:chrom, :strand, :pos]]
	rename!(glori_sel, [:chr, :strand, :pos])
	glori_sel[:, :mod] .= "m6A"
	glori_sel[:, :source] .= "HEK293T"
	
	local rmbase_sel = rmbase_all[:, [:chr, :strand, :start, :mod_type, :cell_list]]
	rename!(rmbase_sel, :mod_type => :mod, :start => :pos, :cell_list => :source)

	local atoi_sel = atoi[:, [:chr, :strand, :start, :mod_type, :cell_list]]
	rename!(atoi_sel, :start => :pos, :mod_type => :mod, :cell_list => :source)

	local m7g_sel = m7GDB_ngs[:, [:seqnames, :strand, :start, :source]]
	m7g_sel[:, :mod] .= "m7G"
	rename!(m7g_sel, :seqnames => :chr, :start => :pos)

	local m5c_sel = m5C_atlas[:, [:chr, :strand, :start, :name]]
	m5c_sel[:, :mod] .= "m5C"
	rename!(m5c_sel, :start => :pos, :name => :source)

	local ac4c_sel = hela_ac4C[:, [:chr, :strand, :start]]
	ac4c_sel[:, :mod] .= "ac4C"
	ac4c_sel[:, :source] .= "HeLa"
	rename!(ac4c_sel, :start => :pos)

	mod_ref = vcat(glori_sel, rmbase_sel, atoi_sel, m7g_sel, m5c_sel, ac4c_sel)

	mod_ref[:, :mod] = ifelse.(map(m -> m in Set(["Um", "Am", "Cm", "Gm"]),  mod_ref.mod), "Nm", mod_ref.mod)
	mod_ref
end

# ╔═╡ 9232e724-d8df-4970-931c-f1c11f0ce859
mod_ref |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/combined_mod_ref.tsv", delim = '\t')

# ╔═╡ 82a9954d-2c23-46bc-a0bb-c2fc078eabb1
begin
	local counts = DataFrame([(k, v) for (k, v) in countmap(mod_ref.mod)],
							 [:mod, :count])
	counts = counts[counts.count .> 5, :]
	sort!(counts, :count, rev = true)
	
	local f = Figure()
	local ax = Axis(f[1, 1],
				    xlabel = "Modification type",
				    ylabel = "Log count",
				    xticks = (1:nrow(counts), counts.mod),
				    xticklabelrotation = 45,
				    yscale = log10)
	
	barplot!(ax, 1:nrow(counts), counts.count)
	
	f
end

# ╔═╡ 57a4e1e9-59cb-4180-b735-9da13feb64a6
begin
	local df = copy(mod_ref)
	rename!(df, :pos => :start)
	df[:, :end] = df.start .+ 1
	df[:, :score] .= 0
	local mod_ref_bed = df[:, [:chr, :start, :end, :mod, :score, :strand, :source]]
	CSV.write("/home/mzdravkov/mod_ref.bed", mod_ref_bed, delim = '\t')
end

# ╔═╡ fbd57196-5c5b-49a9-a5f5-a39d25a33c87
ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_v2_0_0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
		LOR_threshold = 0.8)

# ╔═╡ d41911f6-ab79-4a7f-b245-76ef3cfcd1e6
begin
	peaks_ivt = peaks(ivt, 4)
	sig_peaks_ivt = peaks_ivt[peaks_ivt.predicted .!== missing .&&
						 	  peaks_ivt.predicted .>= -log10(0.01), :]
end

# ╔═╡ 5aaee16e-1121-4bf7-b4c6-60941b62414c
begin
	local df = copy(sig_peaks_ivt)
	df[:, :start] = max.(df.genomicPos .- 4, 0)
	# end is non-inclusive in BED files
	df[:, :end] = df.genomicPos .+ 5
	local peaks_bed = sort(df[:, [:chr, :start, :end, :ref_id, :predicted, :strand]],
		 				   [:chr, :start])
	CSV.write("/home/mzdravkov/ivt_sig_peaks.bed", peaks_bed, delim = '\t')
end

# ╔═╡ 4016b836-3669-4098-9725-d3708755db40
an_peaks_ivt = annotate_peaks(peaks_ivt, mod_ref, writer_motifs, writer_mods)

# ╔═╡ 390466c1-9b31-460b-82cd-676d85ca41da
an_sig_peaks_ivt = annotate_peaks(sig_peaks_ivt, mod_ref, writer_motifs, writer_mods)

# ╔═╡ 081adcae-8bfe-4c6f-9a17-1ddd8e9dd87e
begin
	local bed = an_sig_peaks_ivt[:, [:chr, :genomicPos, :mod, :predicted, :strand]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, :mod] = coalesce.(bed.mod, ".")
	bed[:, [:chr, :genomicPos, :end, :mod, :predicted, :strand]] |> CSV.write("/home/mzdravkov/ivt_mods.bed", delim = '\t', header = false)
end

# ╔═╡ 3a91b67e-8efa-4e84-a519-6426594ea30e
begin
	local bed = an_sig_peaks_ivt[:, [:chr, :genomicPos, :predicted]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, [:chr, :genomicPos, :end, :predicted]] |> CSV.write("/home/mzdravkov/ivt_mods.bedGraph", delim = '\t', header = false)
end

# ╔═╡ a4d697a0-71b6-4e01-879e-6a91d9060b8d
function format_power(base, pow)
	superscript_digits = Dict(
	    '0' => '⁰',
	    '1' => '¹',
	    '2' => '²',
	    '3' => '³',
	    '4' => '⁴',
	    '5' => '⁵',
	    '6' => '⁶',
	    '7' => '⁷',
	    '8' => '⁸',
	    '9' => '⁹'
	)
	
	superscript_number = join([superscript_digits[digit]
							   for digit in string(pow)])
	"$(base)$(superscript_number)"
end

# ╔═╡ 16109b94-06e5-4285-9333-5b8a98761ef7
begin
	local df = copy(mod_ref)
	df[:, :category] = map(r -> if r.source !== missing &&
								occursin("HEK293T", r.source)
									1
								elseif r.source !== missing && occursin("HEK293", r.source)
									2
								else
							        3
								end,
								eachrow(df))
	local counts = DataFrames.combine(groupby(df, [:mod, :category]), nrow => :count)
	local totals = DataFrames.combine(groupby(counts, :mod), :count => sum => :count)
	local common_mods = totals[totals.count .>= 10, :mod] |> Set
	counts = counts[map(m -> m in common_mods, counts.mod), :]
	totals = Dict(zip(totals.mod, totals.count))
	counts[:, :total] = map(m -> totals[m], counts.mod)
	sort!(counts, :total, rev = true)
	local mods = unique(counts.mod)
	local mods_ids = Dict(collect(zip(mods, 1:length(mods))))
	local colormap = Dict(1 => :green,
						  2 => :orange,
						  3 => :lightgrey,
						  4 => :lightyellow)
	local f = Figure()
	local max_pow = Int(ceil(log10.(max(values(totals)...))))
	local ax = Axis(f[1, 1],
					title = "Modifications in the reference",
				    xlabel = "Modification type",
				    ylabel = "log₁₀(count)",
				    xticks = (1:length(unique(counts.mod)), unique(counts.mod)),
				    xticklabelrotation = 45,
				    # yscale = log10)
				    yticks = (0:max_pow, [format_power(10, i) for i in 0:max_pow]))
	barplot!(ax,
			 map(m -> mods_ids[m], counts.mod),
			 map(m -> log10(totals[m]), counts.mod) .* counts.count ./ map(m -> totals[m], counts.mod),
	         stack = counts.category,
	         color = map(c -> colormap[c], counts.category),
			 strokecolor = :black,
			 strokewidth = 0.25)

	Legend(f[2, 1],
	       [PolyElement(color = :green),
			PolyElement(color = :orange),
		    PolyElement(color = :lightgrey)],
	       ["HEK293T", "HEK293", "Other/Unspecified"],
		   orientation = :horizontal)
	
	f
end

# ╔═╡ 23ce8398-5f79-47d4-afe9-c46818011b9a
begin
	local df = innerjoin(mod_ref, ivt,
		  				 on = [:chr, :strand, :pos => :genomicPos],
		  				 makeunique = true)
	println(nrow(df))
	df[:, :category] = map(r -> if r.source !== missing &&
								occursin("HEK293T", r.source)
									1
								elseif r.source !== missing && occursin("HEK293", r.source)
									2
								else
							        3
								end,
								eachrow(df))
	local counts = DataFrames.combine(groupby(df, [:mod, :category]), nrow => :count)
	local totals = DataFrames.combine(groupby(counts, :mod), :count => sum => :count)
	local common_mods = totals[totals.count .>= 10, :mod] |> Set
	counts = counts[map(m -> m in common_mods, counts.mod), :]
	totals = Dict(zip(totals.mod, totals.count))
	counts[:, :total] = map(m -> totals[m], counts.mod)
	sort!(counts, :total, rev = true)
	local mods = unique(counts.mod)
	local mods_ids = Dict(collect(zip(mods, 1:length(mods))))
	local colormap = Dict(1 => :green,
						  2 => :orange,
						  3 => :lightgrey,
						  4 => :lightyellow)
	local f = Figure()
	local max_pow = Int(ceil(log10.(max(values(totals)...))))
	local ax = Axis(f[1, 1],
					title = "Reference modifications that are covered by Nanocompore",
				    xlabel = "Modification type",
				    ylabel = "log₁₀(count)",
				    xticks = (1:length(unique(counts.mod)), unique(counts.mod)),
				    xticklabelrotation = 45,
				    # yscale = log10)
				    yticks = (0:max_pow, [format_power(10, i) for i in 0:max_pow]))
	ylims!(ax, (0, max_pow + 0.5))
	local max_categories = DataFrames.combine(groupby(counts, :mod),
											  :category => maximum => :max_cat)
	local max_categories = Dict(zip(max_categories.mod, max_categories.max_cat))
	local bar_labels = ifelse.(counts.category .== map(m -> max_categories[m], counts.mod),
							   map(t -> "$(Int(t))", counts.total),
							   "")
	barplot!(ax,
			 map(m -> mods_ids[m], counts.mod),
			 map(m -> log10(totals[m]), counts.mod) .* counts.count ./ map(m -> totals[m], counts.mod),
	         stack = counts.category,
	         color = map(c -> colormap[c], counts.category),
			 bar_labels = bar_labels,
			 strokecolor = :black,
			 strokewidth = 0.25)

	Legend(f[2, 1],
	       [PolyElement(color = :green),
			PolyElement(color = :orange),
		    PolyElement(color = :lightgrey)],
	       ["HEK293T", "HEK293", "Other/Unspecified"],
		   orientation = :horizontal)
	
	f
end

# ╔═╡ 25f512fd-e5f8-40d7-a7eb-2d04e18b884b
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT/WT peaks",
				    xlabel = "|LOR| > 0.5",
				    ylabel = "-log₁₀(q-value)",
				    xticks = 0:0.5:9,
				    yticks = 0:25:300)
	xlims!(ax, (0, 9))
	local df = an_peaks_ivt[an_peaks_ivt.GMM_chi2_pvalue .!== missing .&& an_peaks_ivt.GMM_chi2_qvalue .!= 1, :]
	scatter!(ax,
			 abs.(df.GMM_LOR),
			 -log10.(df.GMM_chi2_qvalue),
			 color = map(r -> if r.mod === missing
				 				:grey
							  elseif occursin("motif", r.source)
								:orange
							  else
								:red
							  end, eachrow(df)),
			 markersize = 5)
	vlines!(ax, [0.8], linestyle = :dash, color = :black, linewidth = 0.8)
	hlines!(ax, [2], linestyle = :dash, color = :black, linewidth = 0.8)
	Legend(f[2, 1],
	       [PolyElement(color = :grey),
			PolyElement(color = :orange),
			PolyElement(color = :red)],
	       ["No reference overlaps", "Motif match", "Overlaps reference mods"],
		   orientation = :horizontal)
	f
end

# ╔═╡ 6669794a-cebf-46f7-8d1b-5bb91ccb7280
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT/WT MA plot - peaks",
				    xlabel = "Total coverage in both conditions",
					yticks = -9:1:9,
				    ylabel = "LOR",
				    xscale = log10)
	xlims!(ax, (10, 12000))
	ylims!(ax, (-9, 9))
	local df = an_peaks_ivt
	scatter!(ax, df.mean_cov, df.GMM_LOR,
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 16fff5ab-601e-42d7-9f74-1670690733bb
begin
	local g = DataFrames.combine(
		groupby(an_sig_peaks_ivt, [:peak_id, :mod]),
		:source => (
		  source -> min(map(s -> if occursin("HEK293T", s)
				  					1
								 elseif occursin("HEK293", s)
									2
		  						 elseif occursin("motif", s)
									4
								 else
									3
								 end,
							ifelse.(source .=== missing, "", source))...)
		  # s -> any(occursin.("HEK293T",
				#  			 ifelse.(s .=== missing, "", s)))
		) => :category)
	local counts = DataFrames.combine(groupby(g, [:mod, :category]),
									  nrow => :count)
	counts[:, :mod] = ifelse.(counts.mod .=== missing,
							  "Not matched",
							  counts.mod)
	local totals = DataFrames.combine(groupby(counts, :mod), :count => sum => :count)
	local totals = Dict(zip(totals.mod, totals.count))
	counts[:, :total] = map(m -> totals[m], counts.mod)
	sort!(counts, :total, rev = true)
	local mods = unique(counts.mod)
	local mods_ids = Dict(collect(zip(mods, 1:length(mods))))
	local f = Figure()
	local max_pow = Int(ceil(log10.(max(values(totals)...))))
	local colormap = Dict(1 => :green,
						  2 => :orange,
						  3 => :lightgrey,
						  4 => :lightyellow,
						  5 => :black)
	local ax = Axis(f[1, 1],
					title = "Putative modifications in significant peaks",
				    xlabel = "Modification type",
				    ylabel = "log₁₀(count)",
				    xticks = (1:length(unique(counts.mod)), unique(counts.mod)),
				    xticklabelrotation = 45,
				    # yscale = log10)
				    yticks = (0:max_pow, [format_power(10, i) for i in 0:max_pow]))
	ylims!(ax, (0, max_pow + 1))
	local max_categories = DataFrames.combine(groupby(counts, :mod),
											  :category => maximum => :max_cat)
	local max_categories = Dict(zip(max_categories.mod, max_categories.max_cat))
	local bar_labels = ifelse.(counts.category .== map(m -> max_categories[m], counts.mod),
							   map(t -> "$(Int(t))", counts.total),
							   "")
	counts[counts.mod .== "Not matched", :category] .= 5
	barplot!(ax,
			 map(m -> mods_ids[m], counts.mod),
			 map(m -> log10(totals[m]), counts.mod) .* counts.count ./ map(m -> totals[m], counts.mod),
	         stack = counts.category,
	         color = map(c -> colormap[c], counts.category),
			 bar_labels = bar_labels,
			 strokecolor = :black,
			 strokewidth = 0.25)

	Legend(f[2, 1],
	       [PolyElement(color = :green),
			PolyElement(color = :orange),
		    PolyElement(color = :lightgrey),
		    PolyElement(color = :lightyellow)],
	       ["HEK293T", "HEK293", "Other/Unspecified", "Motif"],
		   "Attribution source",
		   orientation = :horizontal)
	# local f = Figure()
	# local ax = Axis(f[1, 1],
	# 			    title = "Detected modifications (WT / IVT)",
	# 			    xticks = (1:length(mods), unique(counts.mod)),
	# 			    yticks = (0:500:12000, [if mod(n, 1000) == 0 "$n" else "" end
	# 										for n in 0:500:12000]))
	# 			    # ytickformat = "{:d}")
	# local colormap = Dict(true => :green,
	# 					  # 2 => :orange,
	# 					  false => :lightgrey,
	# 					  4 => :lightyellow)
	# barplot!(ax,
	# 		 map(m -> mods_ids[m], counts.mod),
	# 		 counts.count,
	#          stack = counts.category,
	#          color = map(c -> colormap[c], counts.category),
	# 		 strokecolor = :black,
	# 		 strokewidth = 0.25)
	# Legend(f[2, 1],
	# 	   [PolyElement(color = :green),
	# 		PolyElement(color = :orange),
	# 		PolyElement(color = :lightgrey),
	# 		PolyElement(color = :lightyellow)],
	#        ["HEK293T", "HEK293", "Unspecified/Other cell line", "Motif"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	f
end

# ╔═╡ 86179fb3-d3d0-4f5b-a8d3-14fb43dac834
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
				    xlabel = "Coverage",
				    ylabel = "# positions",
				    # xscale = log10)
				    xticks = 0:2000:10000)
	hist!(ax, ivt.mean_cov, bins = 100)
	vlines!(ax, [200], color = :red, linestyle = :dash, linewidth = 1)
	f
end

# ╔═╡ 5fde98ae-bdf2-4dcb-90a0-f5017b97262a


# ╔═╡ d90664a5-8e7d-469b-8b05-b103d829f672
# ╠═╡ disabled = true
#=╠═╡
tx_ref["ENST00000001008.6|ENSG00000004478.8|OTTHUMG00000090429.4|OTTHUMT00000206861.3|FKBP4-201|FKBP4|3715|protein_coding|"][(1763 + 1 + 4):(1763 + 1 + 8 + 4)]
  ╠═╡ =#

# ╔═╡ ab03ca5a-3186-493e-9e3e-13367a85e624
# ╠═╡ disabled = true
#=╠═╡
Pkg.add("FASTX")
  ╠═╡ =#

# ╔═╡ cefdac37-7b63-4c8c-ad92-41ca18b78452
# ╠═╡ disabled = true
#=╠═╡
using FASTX
  ╠═╡ =#

# ╔═╡ 140e4f5f-3c19-4c94-8829-4a46a95eb6cd
# ╠═╡ disabled = true
#=╠═╡
begin
	tx_ref = Dict()
	FASTAReader(open("/work/mzdravkov/gencode.v41.transcripts.fa")) do reader
	   for record in reader
		   tx_ref[identifier(record)] = sequence(record)
	   end
	end
end
  ╠═╡ =#

# ╔═╡ a3f63803-5803-4942-8dce-b1b90b86f6a1
function between(x, s, e)
	s <= x < e
end

# ╔═╡ 9513c89b-dd9e-40dd-91be-c480438e3f07
function overlaps(mod_ref, row)
	matches = mod_ref[mod_ref.chr .== row.chr .&&
		    		  mod_ref.strand .== row.strand .&&
		    		  between.(mod_ref.pos, row.start, row.end), :]
	Dict(zip(matches.mod, matches.pos))
end

# ╔═╡ 56ed77ee-fcfb-498a-ae5c-f9e77ceb971a
begin
	local df = copy(mod_ref)
	df = df[map(m -> m in Set(["m6A", "m5C", "Nm", "m7G", "m1A", "ac4C", "A-I", "Y", "m3C"]), df.mod), :]
	df[:, :source] = map(r -> if occursin("HEK293T", r.source)
							   		 "HEK293T"
								   elseif occursin("HEK293", r.source)
									 "HEK293"
								   else
								   	 "Other"
								   end, eachrow(df))
		
	local table = combine(
			groupby(df, :mod),
			nrow => :count,
		    :source => (s -> 100*sum(s .== "HEK293T")/length(s)) => :HEK293T,
		    :source => (s -> 100*sum(s .== "HEK293")/length(s)) => :HEK293,
		    :source => (s -> 100*sum(s .== "Other")/length(s)) => :Other)

	local covered_ref_mods = combine(
		groupby(innerjoin(mod_ref, ivt,
						  on = [:chr, :strand, :pos => :genomicPos],
						  makeunique = true),
				:mod),
		nrow => :ivt_covered)

	table = leftjoin(table, covered_ref_mods, on = :mod)

	local mod_peaks = combine(groupby(an_sig_peaks_ivt, :mod),
							  nrow => :peaks)

	table = outerjoin(table, mod_peaks,
					  on = :mod,
			 		  matchmissing = :equal)
	sort!(table, :count)
	ivt_mod_table = table
end

# ╔═╡ e429bac0-0add-4e15-8a7a-6ab22a5c9943
begin
	local ref = copy(mod_ref[:, [:mod, :source]])
	ref[:, :source] = map(
		s -> if occursin("HEK293T", s)
				"HEK293T"
			 elseif occursin("HEK293", s)
				"HEK293"
			 else
				"Other"
			 end,
		mod_ref.source)
	local ref_counts = combine(groupby(ref, [:mod, :source]), nrow => :count)


	local covered_ref = innerjoin(mod_ref, ivt,
						  on = [:chr, :strand, :pos => :genomicPos],
						  makeunique = true)[:, [:mod, :source]]
	covered_ref[:, :source] = map(
		s -> if occursin("HEK293T", s)
				"HEK293T"
			 elseif occursin("HEK293", s)
				"HEK293"
			 else
				"Other"
			 end,
		covered_ref.source)
	local covered_ref_counts = combine(groupby(covered_ref, [:mod, :source]),
									   nrow => :count)

	local an_sig_peaks_ivt_cp = copy(an_sig_peaks_ivt[:, [:mod, :source]])
	an_sig_peaks_ivt_cp[:, :source] = map(
		s -> if s === missing
			    "Unclassified"
			 elseif s === "motif"
				"Motif"
			 elseif occursin("HEK293T", s)
				"HEK293T"
			 elseif occursin("HEK293", s)
				"HEK293"
			 else
				"Other"
			 end,
		an_sig_peaks_ivt_cp.source)
	local detected_counts = combine(groupby(an_sig_peaks_ivt_cp[an_sig_peaks_ivt_cp.mod .!== missing, :], [:mod, :source]),
								    nrow => :count)

	local colors = Dict(
		"HEK293T" => :green,
		"HEK293" => :orange,
		"Motif" => :lightyellow,
		"Other" => :lightgrey,
		"Unclassified" => :black
	)
	
	local unclassified = ivt_mod_table[ivt_mod_table.mod .=== missing, :peaks][1]
	local table = ivt_mod_table[ivt_mod_table.mod .!== missing, :]
	
	local f = Figure(size = (800, 600))
	local ax = Axis(f[1, 1])
				    # aspect = 1)
	# xlims!(ax, (6, 44))
	local rows = 0:10:(10+10*nrow(table))
	for (h, x) in zip(["Modification", "Reference", "IVT covered", "In peaks"], [0, 25, 55, 85])
		if h == "Modification"
			text!(ax, x, (10*nrow(table));
				  text = h,
				  font = :bold,
				  align = (:left, :center),
				  fontsize = 16)
		else
			text!(ax, x, (10*nrow(table));
				  text = h,
				  font = :bold,
				  align = (:left, :center),
				  fontsize = 16)
		end
	end
	for (mod, y) in zip(table.mod, rows)
		text!(ax, 0, y;
			  text = coalesce(mod, "Unclassified"),
			  align = (:left, :center),
			  fontsize = 15)
		hlines!(ax, [y - 5], color = :black, linewidth = 0.5)
		hlines!(ax, [y + 5], color = :black, linewidth = 0.5)
	end
	for (mod, c, y) in zip(table.mod, table.count, rows)
		text!(ax, 25, y;
			  text = ifelse(c === missing, "", string(c)),
			  align = (:right, :center),
			  fontsize = 15)
		rc = ref_counts[ref_counts.mod .=== mod, :]
		pie!(ax, 33, y, rc.count,
			 color = map(s -> colors[s], rc.source),
		 	 radius = 4)
	end
	for (mod, c, y) in zip(table.mod, table.ivt_covered, rows)
		text!(ax, 55, y;
			  text = ifelse(c === missing, "", string(c)),
			  align = (:right, :center),
			  fontsize = 15)
		rc = covered_ref_counts[covered_ref_counts.mod .=== mod, :]
		pie!(ax, 63, y, rc.count,
			 color = map(s -> colors[s], rc.source),
		 	 radius = 4)
	end
	for (mod, c, y) in zip(table.mod, table.peaks, rows)
		text!(ax, 85, y;
			  text = ifelse(c === missing, "", string(c)),
			  align = (:right, :center),
			  fontsize = 15)
		rc = detected_counts[detected_counts.mod .=== mod, :]
		pie!(ax, 93, y, rc.count,
			 color = map(s -> colors[s], rc.source),
		 	 radius = 4)
	end

	text!(ax, 0, -12,
		  text = "Total",
		  fontsize = 15,
		  align = (:left, :center))
	text!(ax, 25, -12,
		  text = string(sum(ref_counts.count)),
		  fontsize = 15,
		  align = (:right, :center))
	text!(ax, 55, -12,
		  text = string(sum(covered_ref_counts.count)),
		  fontsize = 15,
		  align = (:right, :center))
	text!(ax, 85, -12,
		  text = string(sum(detected_counts.count)),
		  fontsize = 15,
		  align = (:right, :center))
	
	text!(ax, 0, -22,
		  text = "Unclassified peaks",
		  fontsize = 15,
		  align = (:left, :center))
	text!(ax, 85, -22,
		  text = string(unclassified),
		  fontsize = 15,
		  align = (:right, :center))
	
	Legend(f[2, 1],
		   [PolyElement(color = :green),
			PolyElement(color = :orange),
			PolyElement(color = :lightgrey),
			PolyElement(color = :lightyellow)],
	       ["HEK293T", "HEK293", "Unspecified/Other cell line", "Motif"],
		   "Attribution source",
		   orientation = :horizontal,
		   framevisible = false)
	
	hidedecorations!(ax)
	hidespines!(ax)
	f
end

# ╔═╡ 2b618c14-b478-4d5b-abfc-19bbeb073fde
md"""
## STORM/WT
"""

# ╔═╡ ea86f47c-ed99-4108-bb5c-c9ee225dc7ef
storm = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_test_v2_release/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
	    LOR_threshold = 0.8)

# ╔═╡ 7fcfeb9f-5e35-4408-8cfb-b24c5f52c6a2
begin
	peaks_storm = peaks(storm, 4)
	sig_peaks_storm = peaks_storm[peaks_storm.predicted .!== missing .&&
						  	      peaks_storm.predicted .>= -log10(0.01), :]
end

# ╔═╡ eb2f9c7c-b5bb-493d-8cc1-abdc452e646d
an_peaks_storm = annotate_peaks(peaks_storm, mod_ref, writer_motifs, writer_mods)

# ╔═╡ ed261b2a-459a-440e-99ad-91afdeda2952
an_sig_peaks_storm = annotate_peaks(sig_peaks_storm, mod_ref, writer_motifs, writer_mods)

# ╔═╡ dc85be9a-f0af-4b8b-924e-758b3ce45a7a
begin
	local bed = an_sig_peaks_storm[:, [:chr, :genomicPos, :mod, :predicted, :strand]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, :mod] = coalesce.(bed.mod, ".")
	bed[:, [:chr, :genomicPos, :end, :mod, :predicted, :strand]] |> CSV.write("/home/mzdravkov/storm_mods.bed", delim = '\t', header = false)
end

# ╔═╡ bfe35e43-60fd-4ba9-9f05-f5a7f868e011
begin
	local bed = an_sig_peaks_storm[:, [:chr, :genomicPos, :predicted]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, [:chr, :genomicPos, :end, :predicted]] |> CSV.write("/home/mzdravkov/storm_mods.bedGraph", delim = '\t', header = false)
end

# ╔═╡ a81dc14a-fd54-417b-bcd7-84615a62aa37
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "STORM/WT peaks",
				    xlabel = "|LOR| > 0.5",
				    ylabel = "-log₁₀(q-value)",
					xticks = 0:0.5:9,
				    yticks = 0:25:300)
	xlims!(ax, (0, 9))
	local df = an_peaks_storm[an_peaks_storm.GMM_chi2_pvalue .!== missing .&& an_peaks_storm.GMM_chi2_qvalue .!= 1, :]
	scatter!(ax,
			 abs.(df.GMM_LOR),
			 -log10.(df.GMM_chi2_qvalue),
			 color = map(r -> if r.mod === missing
				 				:grey
							  elseif occursin("motif", r.source)
								:orange
							  else
								:red
							  end, eachrow(df)),
			 markersize = 5)
	vlines!(ax, [0.8], linestyle = :dash, color = :black, linewidth = 0.8)
	hlines!(ax, [2], linestyle = :dash, color = :black, linewidth = 0.8)
	Legend(f[2, 1],
	       [PolyElement(color = :grey),
			PolyElement(color = :orange),
			PolyElement(color = :red)],
	       ["No reference overlaps", "Motif match", "Overlaps reference mods"],
		   orientation = :horizontal)
	f
end

# ╔═╡ a2dbec3c-961f-4904-961d-be2d9d538720
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "STORM/WT MA plot - peaks",
				    xlabel = "Total coverage in both conditions",
					yticks = -9:1:9,
				    ylabel = "LOR",
				    xscale = log10)
	xlims!(ax, (10, 12000))
	ylims!(ax, (-9, 9))
	local df = an_peaks_storm
	scatter!(ax, df.mean_cov, df.GMM_LOR,
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ bf514f88-b2f4-4a35-8c4a-08fec382e74c
md"""
## Single molecule IVT/WT
"""

# ╔═╡ c1002495-1c24-450c-9585-1aeb7d588bb5
import MultipleTesting

# ╔═╡ 823b1914-ab79-499a-98ec-348f24b0107a
begin
	comod = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/single_molecule_comodification2.tsv"))
	comod = comod[.! isnan.(comod.pvalue) .&& comod.expected .!= 0 .&& comod.observed .>= 5, :]
	comod[:, :qvalue] = MultipleTesting.adjust(comod.pvalue, MultipleTesting.BenjaminiHochberg())
	comod
end

# ╔═╡ 80836ca5-66c5-4f9c-a815-9c57230f6948
sig_comod = comod[comod.qvalue .<= 0.01, :]

# ╔═╡ 66983a7b-8a6a-4b9f-91ba-bb40b5b5ebce
combine(groupby(sig_comod, [:reference, :pos1, :pos2]), nrow => :count)

# ╔═╡ 6c7f7c71-ed77-46b1-97a7-0d8b9c8d758e
begin
	local f = Figure(size = (1200, 800))
	local ax1 = Axis(f[1, 1],
					 title = "Co-unmodified (0-0)",
					 xlabel = "log₂(Obs/Exp)",
					 ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)")
	local ax2 = Axis(f[1, 2],
					 title = "Mutually exclusive (0-1/1-0)",
					 xlabel = "log₂(Obs/Exp)",
					 ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)")
	local ax3 = Axis(f[1, 3],
					 title = "Co-modified (1-1)",
					 xlabel = "log₂(Obs/Exp)",
					 ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)")

	local s00 = comod.state1 .== 0 .&& comod.state2 .== 0
	local s01 = comod.state1 .== 0 .&& comod.state2 .== 1
	local s10 = comod.state1 .== 1 .&& comod.state2 .== 0
	local s11 = comod.state1 .== 1 .&& comod.state2 .== 1

	local effsize = log2.(comod.observed ./ comod.expected)
	
	scatter!(ax1,
			 effsize[s00],
			 -log10.(comod[s00, :pvalue]),
			 color = ifelse.(comod[s00, :qvalue] .<= 0.01, :red, :grey),
			 markersize = 5)
	
	local selected = s00 .&& -log10.(comod.qvalue) .> 60
	text!(ax1,
		  [Point2f(x, y) for (x, y) in zip(effsize[selected], -log10.(comod.pvalue[selected]))];
		  text = map(r -> split(r, "|")[5], comod.reference[selected]),
		  align = (:left, :bottom),
		  fontsize = 9)

	println("Significant 0-0: ", nrow(comod[s00 .&& comod.qvalue .<= 0.01, :]))

	scatter!(ax2,
			 effsize[s01 .|| s10],
			 -log10.(comod[s01 .|| s10, :pvalue]),
			 color = ifelse.(comod[s01 .|| s10, :qvalue] .<= 0.01, :red, :grey),
			 markersize = 5)
	local selected = (s01 .|| s10) .&& -log10.(comod.qvalue) .> 110
	local offsets = Dict("SERBP1-203" => (0.3, 3),
						 "YBX3-201" => (0, 2),
						 "TUBA1B-201" => (0.01, -4),
						 "PGAM1-201" => (0.01, -4),
						 "HSPD1-202" => (-0.04, 3),
						 "EEF2-201" => (0.08, -2))
	local offset = m -> m in keys(offsets) ? offsets[m] : (0, 0)
	text!(ax2,
		  [Point2f(x + 0.01 + offset(m)[1],
				   y + 1 + offset(m)[2])
		   for (x, y, m) in zip(effsize[selected],
							    -log10.(comod.pvalue[selected]),
							    map(r -> String(split(r, "|")[5]), comod.reference[selected]))];
		  text = map(r -> split(r, "|")[5], comod.reference[selected]),
		  align = (:left, :bottom),
		  fontsize = 9)
	println("Significant 0-1: ", nrow(comod[s01 .&& comod.qvalue .<= 0.01, :]))
	println("Significant 1-0: ", nrow(comod[s10 .&& comod.qvalue .<= 0.01, :]))

	scatter!(ax3,
			 effsize[s11],
			 -log10.(comod[s11, :pvalue]),
			 color = ifelse.(comod[s11, :qvalue] .<= 0.01, :red, :grey),
			 markersize = 5)
	local selected = s11 .&& -log10.(comod.qvalue) .> 75
	text!(ax3,
		  [Point2f(x, y) for (x, y) in zip(effsize[selected], -log10.(comod.pvalue[selected]))];
		  text = map(r -> split(r, "|")[5], comod.reference[selected]),
		  align = (:left, :bottom),
		  fontsize = 9)

	println("Significant 1-1: ", nrow(comod[s11 .&& comod.qvalue .<= 0.01, :]))

	linkaxes!(ax1, ax2, ax3)

	Legend(f[2, 1:3],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 29c659f1-4e5a-4c65-aa53-33dea7d730f3
sig_comod[-log10.(sig_comod.qvalue) .> 120, :]

# ╔═╡ 1d820491-8699-4324-b113-b4d7ac0a819a
begin
	local f = Figure(size = (1200, 800))
	local ax = Axis(f[1, 1],
					title = "Distance between significant associated sites",
				    xlabel = "Distance (nt)",
				    ylabel = "Count",
				    xticks = 0:50:1400)

	df = sig_comod[:, [:pos1, :pos2]] |> unique

	hist!(ax, abs.(df.pos2 .- df.pos1), bins=nrow(df), color = :black)
	f
end

# ╔═╡ bc42854a-ecbe-46c2-b289-114b12cf6fe9
begin
	local p1 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	rename!(p1, :mod => :mod1)
	local p2 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	rename!(p2, :mod => :mod2)

	an_sig_comod = leftjoin(leftjoin(sig_comod, p1,
		 	 		        		 on = [:reference => :ref_id, :pos1 => :pos]),
			     		    p2,
			 				on = [:reference => :ref_id, :pos2 => :pos])
end

# ╔═╡ 26d0dc17-cd92-47d4-bec0-33ad54c002cd
begin
	local f = Figure()


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)

	local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
	# mods = ifelse.(mods .=== missing, "Unknown", mods)
	local M = length(mods)

	local ax = Axis(f[1, 1],
					title = "Number of associated modification pair sites by modification types",
				    xticks = (1:M, mods),
				    yticks = (1:M, mods))

	df = unique(df[:, [:pos1, :pos2, :mod1, :mod2]])

	local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)

	local t = 0
	local counts = zeros((length(mods), length(mods)))
	for (i, m1) in enumerate(mods)
		for (j, m2) in enumerate(mods)
			if j < i
				continue
			end
			c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
						(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
			counts[i, j] = c
			t += c
		end
	end
	println(t)

	heatmap!(ax, 1:M, 1:M, counts, colormap = :heat)
	text!(ax,
		  [string(Int(counts[i, j])) for i in 1:M for j in i:M],
		  position = [Point2f(x, y) for x in 1:M for y in x:M],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[i, j] > 2000 ? :white : :black for i in 1:M for j in i:M])

	f
end

# ╔═╡ c2cf87ff-d11c-4bb3-b44d-6a16f4040c88
begin
	local mods = ("m5C", "Y")
	local df = an_sig_comod[(an_sig_comod.mod1 .=== mods[1] .&& an_sig_comod.mod2 .=== mods[2]) .||
		(an_sig_comod.mod1 .=== mods[2] .&& an_sig_comod.mod2 .=== mods[1]), :]
	df[:, :dir] = df.observed .> df.expected
	combine(groupby(df, [:state1, :state2, :dir]), nrow => :count)
end

# ╔═╡ fd549de3-8c1f-4358-8ec0-48fdd85d0bd1
function non_trivial(row)
	if row.mod1 !== row.mod2
		true
	elseif (row.state1, row.state2) == (0, 0) && row.observed > row.expected
		false
	elseif (row.state1, row.state2) == (0, 1) && row.observed < row.expected
		false
	elseif (row.state1, row.state2) == (1, 0) && row.observed < row.expected
		false
	elseif (row.state1, row.state2) == (1, 1) && row.observed > row.expected
		false
	else
		true
	end
end

# ╔═╡ b37f9ca0-5a53-4716-9589-9b57268f72b2
begin
	local f = Figure()


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)

	df = filter(non_trivial, df)

	local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
	# mods = ifelse.(mods .=== missing, "Unknown", mods)
	local M = length(mods)

	local ax = Axis(f[1, 1],
					title = "Number of non-trivial associated modification pair sites by modification types",
				    xticks = (1:M, mods),
				    yticks = (1:M, mods))

	df = unique(df[:, [:pos1, :pos2, :mod1, :mod2]])
	
	local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)

	local t = 0
	local counts = zeros((length(mods), length(mods)))
	for (i, m1) in enumerate(mods)
		for (j, m2) in enumerate(mods)
			if j < i
				continue
			end
			c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
						(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
			counts[i, j] = c
			t += c
		end
	end
	println(t)

	heatmap!(ax, 1:M, 1:M, counts, colormap = :heat)
	text!(ax,
		  [string(Int(counts[i, j])) for i in 1:M for j in i:M],
		  position = [Point2f(x, y) for x in 1:M for y in x:M],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[i, j] > 1500 ? :white : :black for i in 1:M for j in i:M])

	f
end

# ╔═╡ 67638f93-2a44-4772-9c35-2e728cebc277


# ╔═╡ c99dbb5d-77da-4c2e-9b7a-874f5ead6574
map(println, sort(unique(map(r -> r[6], split.(filter(r -> non_trivial(r) && r.mod1 !== missing && r.mod2 !== missing, an_sig_comod).reference, "|")))))

# ╔═╡ 409ec9ef-22c6-4192-a53c-382d059ff299
begin
	local f = Figure(size = (1200, 800))
	local ax = Axis(f[1, 1],
					title = "Distance between significant non-trivially associated sites",
				    xlabel = "Distance (nt)",
				    ylabel = "Count",
				    xticks = 0:50:1400)
	
	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)
	df = filter(non_trivial, df)

	hist!(ax, abs.(df.pos2 .- df.pos1), bins=1400)
	f
end

# ╔═╡ 664f5ea5-378a-45cc-a6cf-cab779b7f933
begin
	# local df = copy(an_sig_comod)
	
	# sort(combine(groupby(df, [:pattern, :direction]), nrow => :count), [:pattern, :direction])
	
	local f = Figure()


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)
	df[:, :pattern] = map(r -> Dict((0, 0) => "00", (0, 1) => "01", (1, 0) => "01", (1, 1) => "11")[(r.state1, r.state2)], eachrow(df))
	df[:, :direction] = sign.(log2.(df.observed ./ df.expected))
	df = df[df.pattern .== "01" .&& df.direction .== 1, :]

	local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
	# mods = ifelse.(mods .=== missing, "Unknown", mods)
	local M = length(mods)

	local ax = Axis(f[1, 1],
					title = "Number of associated modification pair sites by modification types",
				    xticks = (1:M, mods),
				    yticks = (1:M, mods))

	local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)

	local t = 0
	local counts = zeros((length(mods), length(mods)))
	for (i, m1) in enumerate(mods)
		for (j, m2) in enumerate(mods)
			if j < i
				continue
			end
			c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
						(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
			counts[i, j] = c
			t += c
		end
	end
	println(t)

	heatmap!(ax, 1:M, 1:M, counts, colormap = :heat)
	text!(ax,
		  [string(Int(counts[i, j])) for i in 1:M for j in i:M],
		  position = [Point2f(x, y) for x in 1:M for y in x:M],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[i, j] > 4000 ? :white : :black for i in 1:M for j in i:M])

	f
end

# ╔═╡ d1c2e993-4ba4-46a5-840c-38bfbf397858
begin
	# local df = copy(an_sig_comod)
	
	# sort(combine(groupby(df, [:pattern, :direction]), nrow => :count), [:pattern, :direction])
	
	local f = Figure(size = (1600, 950))

	for (c, pattern) in enumerate(["00", "01", "11"]), (r, direction) in enumerate([1, -1])

		local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)
		df[:, :pattern] = map(r -> Dict((0, 0) => "00", (0, 1) => "01", (1, 0) => "01", (1, 1) => "11")[(r.state1, r.state2)], eachrow(df))
		df[:, :direction] = sign.(log2.(df.observed ./ df.expected))
		
	
		local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
		# mods = ifelse.(mods .=== missing, "Unknown", mods)
		local M = length(mods)

		dir = direction == 1 ? "common" : "rare"
		local ax = Axis(f[r, c],
						title = "# mod-associated pairs | pattern: $pattern direction: $dir",
					    xticks = (1:M, mods),
					    yticks = (1:M, mods))
	
		local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)
	
		local t = 0
		local counts = zeros((length(mods), length(mods)))
		for (i, m1) in enumerate(mods)
			for (j, m2) in enumerate(mods)
				if j < i
					continue
				end
				c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
							(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
				counts[i, j] = c
				t += c
			end
		end
		println(t)
	
		df = df[df.pattern .== pattern .&& df.direction .== direction, :]
		local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)
	
		local t = 0
		local subcounts = zeros((length(mods), length(mods)))
		for (i, m1) in enumerate(mods)
			for (j, m2) in enumerate(mods)
				if j < i
					continue
				end
				c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
							(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
				subcounts[i, j] = c
				t += c
			end
		end
	
		heatmap!(ax, 1:M, 1:M, subcounts./(counts.+1), colormap = :heat, colorrange = [0, 1])
		text!(ax,
			  [string(Int(subcounts[i, j])) for i in 1:M for j in i:M],
			  position = [Point2f(x, y) for x in 1:M for y in x:M],
			  align = (:center, :center),
			  fontsize = 14,
			  color = [subcounts[i, j]/(counts[i, j] + 1) > 0.7 ? :white : :black for i in 1:M for j in i:M])
	end

	Colorbar(f[1:2, 4], colorrange = (0, 1), colormap = :heat)

	f
end

# ╔═╡ cb81e3ba-5eea-4169-a4ff-5b7433e8b171
groupby(an_sig_comod[map(v -> v == ["m6A", "Y"], Vector.(eachrow(coalesce.(an_sig_comod[:, [:mod1, :mod2]], "Unknown")))), :], [:reference, :pos1, :pos2]) |> println

# ╔═╡ 07f53334-5428-4890-b02a-9e3c5de7b22c
begin
	local f = Figure(size = (500, 1000))


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)

	axes = []

	for (i, (mod, mod2)) in enumerate([("m6A", "m5C"),
									   ("m6A", "Y"),
									   ("m6A", "m1A"),
									   ("m5C", "Y"),
									   ("m6A", "Nm")])
		# "m5C", "m1A", "Y", "Nm"])
		ax = Axis(f[i, 1],
				  title = "Distance between associated $mod-$mod2 pairs",
				  xticks = 0:100:1200)
		instances = df[(df.mod1 .== mod2 .&& df.mod2 .== mod) .||
		   			   (df.mod1 .== mod .&& df.mod2 .== mod2), [:pos1, :pos2]]

		instances = unique(instances)


		hist!(ax, abs.(instances.pos2 .- instances.pos1), bins = nrow(instances))
		push!(axes, ax)
	end

	linkxaxes!(axes...)
	
	f
end

# ╔═╡ 40bbff56-4264-4d9a-af0e-59485c923f92
begin
	local f = Figure(size = (400, 1200))

	local df = sig_comod

	# local axes = []
	

	local colors = [:orange, :green, :purple, :grey]
	for (i, pattern) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)])
		ax = Axis(f[i, 1],
		  title = "Distance between associated pairs with pattern $pattern")

		instances = df[df.state1 .== pattern[1] .&&
		   			   df.state2 .== pattern[2], :]

		instances = instances[abs.(instances.pos2 .- instances.pos1) .> 10, :]

		if sum(instances.observed .> instances.expected) > nrow(instances)/2
			hist!(ax, abs.(instances[instances.observed .> instances.expected, :].pos2 .- instances[instances.observed .> instances.expected, :].pos1),
					 label = "obs>exp",
					 color = (:orange, 0.5),
				   	 bins = nrow(instances[instances.observed .> instances.expected, :]))
			hist!(ax, abs.(instances[instances.observed .< instances.expected, :].pos2 .- instances[instances.observed .< instances.expected, :].pos1),
					 label = "obs<exp",
					 color = (:blue, 0.5),
				     bins = nrow(instances[instances.observed .< instances.expected, :]))
		else
			hist!(ax, abs.(instances[instances.observed .< instances.expected, :].pos2 .- instances[instances.observed .< instances.expected, :].pos1),
				 label = "obs<exp",
				 color = (:blue, 0.5),
			     bins = nrow(instances[instances.observed .< instances.expected, :]))
			hist!(ax, abs.(instances[instances.observed .> instances.expected, :].pos2 .- instances[instances.observed .> instances.expected, :].pos1),
				 label = "obs>exp",
				 color = (:orange, 0.5),
			   	 bins = nrow(instances[instances.observed .> instances.expected, :]))
		end
		# hist!(ax, abs.(instances.pos2 .- instances.pos1),
		# 	  	 # label = ifelse.(instances.observed .> instances.expected,
		# 					# 	 "obs>exp",
		# 					# 	 "obs<exp"),
		# 		 color = ifelse.(instances.observed .> instances.expected,
		# 						 :orange,
		# 						 :blue),

		# 	     bins = nrow(instances))

		axislegend(ax)
		# xlims!(ax, (10, 500))

		# push!(axes, ax)
	end

	# linkxaxes!(axes...)
	# linkyaxes!(axes...)

	f
end

# ╔═╡ 0ee4e9af-4b1e-41fb-8746-12ff51d20e27
an_sig_comod[occursin.("lnc", an_sig_comod.reference), :]

# ╔═╡ 9503a9e1-8c7e-4dfe-aeed-015d7ec1095c
import SQLite

# ╔═╡ 58e326c7-e5b5-4f97-884e-7321a8084139
import Clustering

# ╔═╡ e30b68c9-9aaa-44f7-86d9-c8e1ab9254dd
map(r -> println(split(r, "|")[6]), unique(sig_comod.reference))

# ╔═╡ 3dc430dc-a218-42bf-90da-8ba8ba6dd56a
begin
	local db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))


	local reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	local probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	probs = ifelse.(probs .== -1, 0, probs)
	probs ./= 100

	# local c = Clustering.kmeans(probs, 4)
	# local idx =  sortperm(assignments(c))
	# # M = M[:,idx]
	# println(idx)

	local N = size(probs, 2)
	local D = size(probs, 1)
	
	dicer_probs = probs
	dicer_distances = zeros((N, N))
	for i in 1:N
		for j in 1:N
			dicer_distances[i, j] = sum((probs[:, i] .- probs[:, j]) .^ 2)/D
		end
	end
end

# ╔═╡ 46b1d7d3-2e9e-4d19-91d8-07150e981ef0
begin
	local f = Figure()

	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))

	local ax = Axis(f[1, 1],
					title = "Modification probability for significant DANCR sites",
					xlabel = "Significant positions",
					ylabel = "Reads",
				    xticks = (1:length(positions), map(string, positions)),
				    xticklabelrotation = 45)

	local clust = Clustering.hclust(dicer_distances, branchorder = :barjoseph)
	heatmap!(ax, dicer_probs[:, clust.order])

	Colorbar(f[1, 2], colorrange = [0, 1])

	f
end

# ╔═╡ 809b0c5c-1538-4bbe-bb4d-293f6058b2a6
function arc_plot2(ref; range = nothing)
	name = split(ref, "|")[5]

	df = an_sig_comod[an_sig_comod.reference .== ref, :]
	if range !== nothing
		df = df[between.(df.pos1, range[1], range[2]) .&& between.(df.pos2, range[1], range[2]), :]
	end
	positions = unique(union(df.pos1, df.pos2))
	
	f = Figure(size = (1000, 700))
	ax = Axis(f[1, 1],
			  title = "$name modification site associations",
			  xlabel = "Reference position",
			  # xticks = 700:20:960,
			  yticksvisible = false,
			  yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	# local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end
	allmods = comod[comod.reference .== ref, :]
	if range !== nothing
		allmods = allmods[between.(allmods.pos1, range[1], range[2]) .&& between.(allmods.pos2, range[1], range[2]), :]
	end
	for row in eachrow(allmods)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for row in eachrow(df)
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			colormap = :reds
			colorrange = [-200, 10]
		elseif row.state1 == 0 && row.state2 == 1
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [0, 10]
		elseif row.state1 == 1 && row.state2 == 0
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [10, 0]
		else
			colormap = :greens
			colorrange = [-200, 10]
		end
		direction = row.observed > row.expected ? 1 : -1
		y = row.observed > row.expected ? 2.5 : -2.5
		style = row.observed > row.expected ? :solid : :dash
		arcsize = π/10
		for i in 0:9
			arcstart = direction*i*arcsize
			arcend = direction*(i+1)*arcsize
			arc!(ax, Point2f(center, y), radius, arcstart, arcend,
				 color = i,
				 colorrange = colorrange,
				 colormap = colormap,
				 linestyle = :solid,
				 linewidth = 2)
				 # linewidth = -log10(row.qvalue))
		end
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :red),
			PolyElement(color = :green),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["Non-modified", "Modified", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ beb418fd-d39b-4c5e-bf9d-98e936b9c819
arc_plot2("ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|")

# ╔═╡ de7c2cc7-7931-4290-97be-39b24634d755
arc_plot2("ENST00000370990.5|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000095785.1|SERBP1-202|SERBP1|1621|protein_coding|"; range = (1100, 1450))

# ╔═╡ 088e365f-05f2-429c-9fa2-763ddf30b82b
arc_plot2("ENST00000605788.6|ENSG00000166197.17|OTTHUMG00000018944.5|OTTHUMT00000050012.2|NOLC1-209|NOLC1|3723|protein_coding|"; range = (3650, Inf))

# ╔═╡ ea219179-f224-4aa3-bb1d-fdf494742363
begin
	local df = an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
	local positions = unique(union(df.pos1, df.pos2))
	
	local f = Figure(size = (1000, 700))
	local ax = Axis(f[1, 1],
				    title = "DANCR (lncRNA) modification site associations",
				    xlabel = "Reference position",
					xticks = 700:20:960,
				    yticksvisible = false,
				    yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	local mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end
	for row in eachrow(comod[comod.reference .== ref, :])
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	hlines!(ax, [0, 0], [minimum(keys(mods)), maximum(keys(mods))], color = :lightgrey, alpha = 0.5, linewidth = 16)

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for row in eachrow(df)
		println(row.pos1, " ", row.pos2, " ", row.state1, " ", row.state2, " ", row.expected, " ", row.observed)
		# angle = row.state1 == row.state2 == 0 ? -π : π
		# y = row.state1 == row.state2 == 0 ? -2.5 : 2.5
		# color = row.observed > row.expected ? :green : :red
		# color, cmap = colormap[(row.state1, row.state2)]
		# style = row.state1 == row.state2 == 0 ? :dash : :solid
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			colormap = :reds
			colorrange = [-200, 10]
		elseif row.state1 == 0 && row.state2 == 1
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [0, 10]
		elseif row.state1 == 1 && row.state2 == 0
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [10, 0]
		else
			colormap = :greens
			colorrange = [-200, 10]
		end
		direction = row.observed > row.expected ? 1 : -1
		y = row.observed > row.expected ? 2.5 : -2.5
		style = row.observed > row.expected ? :solid : :dash
		arcsize = π/10
		for i in 0:9
			arcstart = direction*i*arcsize
			arcend = direction*(i+1)*arcsize
			arc!(ax, Point2f(center, y), radius, arcstart, arcend,
				 color = i,
				 colorrange = colorrange,
				 colormap = colormap,
				 linestyle = :solid,
				 linewidth = -log10(row.qvalue))
		end
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :red),
			PolyElement(color = :green),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["Non-modified", "Modified", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ bcc0fec3-9925-4c52-8411-a59583f9171e
begin
	local df = an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
	local positions = unique(union(df.pos1, df.pos2))
	
	local f = Figure(size = (1000, 700))
	local ax = Axis(f[1, 1],
				    title = "DANCR (lncRNA) modification site associations",
				    xlabel = "Reference position",
					xticks = 700:20:960,
				    yticksvisible = false,
				    yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	local mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end

	for row in eachrow(comod[comod.reference .== ref, :])
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for row in eachrow(df)
		println(-log10(row.qvalue))
		println(row.pos1, " ", row.pos2, " ", row.state1, " ", row.state2, " ", row.expected, " ", row.observed)
		# angle = row.state1 == row.state2 == 0 ? -π : π
		# y = row.state1 == row.state2 == 0 ? -2.5 : 2.5
		# color = row.observed > row.expected ? :green : :red
		# color, cmap = colormap[(row.state1, row.state2)]
		# style = row.state1 == row.state2 == 0 ? :dash : :solid
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			color = :blue
		elseif row.state1 == 0 && row.state2 == 1
			color = :orange
		elseif row.state1 == 1 && row.state2 == 0
			color = :orange
		else
			color = :purple
		end
		angle = row.observed > row.expected ? π : -π
		y = row.observed > row.expected ? 2.5 : -2.5
		style = row.observed > row.expected ? :solid : :dash
		arc!(ax, Point2f(center, y), radius, 0, angle,
			 color = color,
			 # linewidth = 2,
			 linewidth = -log10(row.qvalue))
		
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :blue),
			PolyElement(color = :orange),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["State 0-0", "State 0-1/1-0", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ b117d0ed-b44b-46a6-b9ea-204e1d5195e9
function arc_plot(ref; range = nothing)
	name = split(ref, "|")[5]
	# df = an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
	df = an_sig_comod[an_sig_comod.reference .== ref, :]
	if range !== nothing
		df = df[between.(df.pos1, range[1], range[2]) .&& between.(df.pos2, range[1], range[2]), :]
	end
	positions = unique(union(df.pos1, df.pos2))
	
	f = Figure(size = (1000, 700))
	ax = Axis(f[1, 1],
		  	  title = "$name modification site associations",
			  xlabel = "Reference position",
			  # xticks = 700:20:960,
			  yticksvisible = false,
			  yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	# local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	local mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end
	allmods = comod[comod.reference .== ref, :]
	if range !== nothing
		allmods = allmods[between.(allmods.pos1, range[1], range[2]) .&& between.(allmods.pos2, range[1], range[2]), :]
	end
	for row in eachrow(allmods)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	max_radius = maximum(df.pos2 .- df.pos1)/2
	arc_offset = log10(max_radius)

	if max_radius < 16
		ylims!(ax, (-16, 16))
	end

	hlines!(ax, [0, 0], [minimum(keys(mods)), maximum(keys(mods))], color = :lightgrey, alpha = 0.5, linewidth = 16)

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for (i, row) in enumerate(eachrow(df))
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			color = :blue
		elseif row.state1 == 0 && row.state2 == 1
			color = :orange
		elseif row.state1 == 1 && row.state2 == 0
			color = :orange
		else
			color = :purple
		end
		angle = row.observed > row.expected ? π : -π
		y = row.observed > row.expected ? arc_offset : -arc_offset
		Vector.(eachrow(df[1:i, [:pos1, :pos2]]))
		prev = df[1:i, :]
		prevcount = sum(prev.pos1 .== row.pos1 .&& prev.pos2 .== row.pos2 .&& (prev.observed .> prev.expected) .== (row.observed > row.expected))
		y *= prevcount * 1.5*arc_offset
		if y > 0
			y = max(y, max_radius/16)
		end
		if y < 0
			y = min(y, -max_radius/16)
		end

		style = row.observed > row.expected ? :solid : :dash
		arc!(ax, Point2f(center, y), radius, 0, angle,
			 color = color,
			 linewidth = 2)
			 # linewidth = -log10(row.qvalue))
		
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :blue),
			PolyElement(color = :orange),
			PolyElement(color = :purple),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["Co-unmodified (0-0)", "Mutually exclusive (0-1/1-0)", "Co-modified (1-1)", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 5408e59e-53ae-4ef9-b9ae-e86be3119575
arc_plot("ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|")

# ╔═╡ 72faa038-d71c-46f9-8110-d68b7633f2e6
arc_plot("ENST00000370990.5|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000095785.1|SERBP1-202|SERBP1|1621|protein_coding|"; range = (1100, 1500))

# ╔═╡ 9b4cf8f2-84ae-4aba-981a-7ef0e669a453
an_sig_comod[occursin.("SERBP1-202", an_sig_comod.reference), :]

# ╔═╡ d7fb857c-0871-4be2-8322-e36ab95e6f19
arc_plot("ENST00000370994.8|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000025983.2|SERBP1-203|SERBP1|6676|protein_coding|"; range = (1120, 1650))

# ╔═╡ 05ff34ad-5cfc-48a4-bfc9-4d4c7cfda4d1
arc_plot("ENST00000552461.5|ENSG00000089157.16|OTTHUMG00000169317.10|OTTHUMT00000403466.1|RPLP0-226|RPLP0|2591|retained_intron|"; range = (2500, Inf))

# ╔═╡ 171a33c2-e532-4167-82a6-27d830402104
arc_plot("ENST00000419869.7|ENSG00000023734.11|OTTHUMG00000168791.2|OTTHUMT00000401114.2|STRAP-202|STRAP|1874|protein_coding|", range = (1800, Inf))

# ╔═╡ 657e7156-57ef-4f64-a4ba-13af4bb2784e
arc_plot("ENST00000264720.8|ENSG00000115207.15|OTTHUMG00000097785.9|-|GTF3C2-201|GTF3C2|3607|protein_coding|", range = (2000, Inf))

# ╔═╡ 89114c37-3d83-4e4c-a46f-32cd04e85bf0
sort(an_sig_comod[an_sig_comod.mod1 .!== an_sig_comod.mod2 .&& an_sig_comod.mod1 .!== missing .&& an_sig_comod.mod2 .!== missing, :], :qvalue)

# ╔═╡ d6b27bbe-824e-492b-9aea-dd5041f4917f
sort(innerjoin(innerjoin(an_sig_comod[an_sig_comod.mod1 .!== an_sig_comod.mod2 .&& an_sig_comod.mod1 .!== missing .&& an_sig_comod.mod2 .!== missing .&& (.! occursin.("GTF3C2", an_sig_comod.reference))
, :],
		  	an_sig_peaks_ivt[an_sig_peaks_ivt.source .!== "motif", [:ref_id, :pos]],
		  	on = [:reference => :ref_id, :pos1 => :pos]),
		  an_sig_peaks_ivt[an_sig_peaks_ivt.source .!== "motif", [:ref_id, :pos]],
		  on = [:reference => :ref_id, :pos2 => :pos]), :qvalue)

# ╔═╡ e0bfd78c-25fa-4752-a4d7-660aa4b6cd66
begin
	local db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	local ref = "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]


	local reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	local orig_probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	local probs = ifelse.(orig_probs .== -1, 0, orig_probs)
	probs ./= 100

	# local c = Clustering.kmeans(probs, 4)
	# local idx =  sortperm(assignments(c))
	# # M = M[:,idx]
	# println(idx)

	local N = size(probs, 2)
	local D = size(probs, 1)
	
	slc5a6_probs = probs
	slc5a6_distances = zeros((N, N))
	for i in 1:N
		for j in 1:N
			# slc5a6_distances[i, j] = sqrt(sum((probs[:, i] .- probs[:, j]) .^ 2))/D
			slc5a6_distances[i, j] = sum(abs.(probs[:, i] .- probs[:, j]))/D
		end
	end
	slc5a6_probs = ifelse.(orig_probs .== -1, missing, probs)
end

# ╔═╡ b6a72db1-ec4d-4cc1-8174-c6179316dec7
function get_probs(ref)
	
	db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]


	reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	orig_probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	probs = ifelse.(orig_probs .== -1, 0, orig_probs)
	probs ./= 100

	# local c = Clustering.kmeans(probs, 4)
	# local idx =  sortperm(assignments(c))
	# # M = M[:,idx]
	# println(idx)

	N = size(probs, 2)
	D = size(probs, 1)
	
	distances = zeros((N, N))
	for i in 1:N
		for j in 1:N
			# slc5a6_distances[i, j] = sqrt(sum((probs[:, i] .- probs[:, j]) .^ 2))/D
			distances[i, j] = sum(abs.(probs[:, i] .- probs[:, j]))/D
		end
	end
	probs = ifelse.(orig_probs .== -1, missing, probs)
end

# ╔═╡ f2829ecb-cbdf-4e0d-953c-0072f03546ac
begin
	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"
	local dancr_probs = get_probs(ref)
	DataFrame(transpose(dancr_probs), map(string, unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2])))) |> CSV.write("/home/mzdravkov/dancr_probs.tsv", delim = '\t')
end

# ╔═╡ 74279cfc-ab7a-4e9a-b048-ad17a378c20f
begin
	local ref = "ENST00000253024.10|ENSG00000130726.12|OTTHUMG00000183546.3|OTTHUMT00000467074.2|TRIM28-201|TRIM28|3364|protein_coding|"
	local trim28_probs = get_probs(ref)
	DataFrame(transpose(trim28_probs), map(string, unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2])))) |> CSV.write("/home/mzdravkov/trim28_probs.tsv", delim = '\t')
end

# ╔═╡ 5b5358c8-f35c-473f-9e92-56bdb4dbc6a7
arc_plot("ENST00000370990.5|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000095785.1|SERBP1-202|SERBP1|1621|protein_coding|")#
, range = (3200, Inf))

# ╔═╡ 9031f497-26d2-4b3b-ba34-3438df409fd9
an_sig_comod[an_sig_comod.reference .== "ENST00000316448.10|ENSG00000179218.15|OTTHUMG00000180574.3|OTTHUMT00000451952.2|CALR-201|CALR|1901|protein_coding|"


, :]

# ╔═╡ ad11f743-1eee-4298-87bf-620a0f1cd1ec
an_sig_peaks_ivt[an_sig_peaks_ivt.ref_id .== "ENST00000264720.8|ENSG00000115207.15|OTTHUMG00000097785.9|-|GTF3C2-201|GTF3C2|3607|protein_coding|"


, :]

# ╔═╡ f793c68f-731e-426f-ad6e-ef5419f16370
an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]

# ╔═╡ cc022b63-4a04-4142-976c-02756047076d
begin
	local f = Figure()

	local ref = "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]

	local ax = Axis(f[1, 1],
					title = "Modification probability for significant SLC5A6 sites",
					xlabel = "Significant positions",
					ylabel = "Reads",
				    xticks = (1:length(positions), map(string, positions)),
				    xticklabelrotation = 45)

	local clust = Clustering.hclust(slc5a6_distances, branchorder = :barjoseph)
	heatmap!(ax, slc5a6_probs[:, clust.order])

	Colorbar(f[1, 2], colorrange = [0, 1])

	f
	# clust
end

# ╔═╡ a6b2a25c-6fb2-4dd0-a7f1-ea16ceb63a97
DataFrame(transpose(slc5a6_probs), map(string, unique(union(comod[comod.reference .== "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|", :pos1], comod[comod.reference .== "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|", :pos2])))) |> CSV.write("/home/mzdravkov/slc5a6_probs.tsv", delim = '\t')

# ╔═╡ 69798f9b-f9a8-41f7-a2ff-8f333a6ed728
begin
	local df = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/datasets/miR-eCLIP/sample_comparisons/HEK293T_DMSO_vs_HEK293T_STM.miRNA_target_sample_comparisons.final.tsv"))
	rename!(df, [:chr, :start, :end, :stran, :gene, :gene_id, :feature, :mirna, :avg_ctrl, :avg_storm, :lfc, :pscore, :qscore])
	mireclip = df
end

# ╔═╡ 741a21f9-ea94-4620-83a6-a43dd5fc300d
begin
	local df = copy(an_sig_comod)
	df[:, :gene_id] = map(rx -> rx[2], split.(df.reference, "|"))
	df = innerjoin(
			innerjoin(df, rename(ivt[:, [:ref_id, :pos, :genomicPos]], :genomicPos => :genomicPos1),
				   	  on = [:reference => :ref_id, :pos1 => :pos]),
			rename(ivt[:, [:ref_id, :pos, :genomicPos]], :genomicPos => :genomicPos2),
			on =  [:reference => :ref_id, :pos2 => :pos])
	df = innerjoin(df, mireclip,
			  	   on = [:gene_id])

	# df[between.(df.genomicPos1, df.start, df.end) .|| between.(df.genomicPos2, df.start, df.end), :]

	df[df.qscore .>= 2, :]
end

# ╔═╡ e3f624cc-1c45-461a-ac19-bed551b1ba4e
md"""
## STORM resistant mods
"""

# ╔═╡ 5e8c4ba1-9f8d-40df-b02e-cfeec4bf4d57
# ╠═╡ disabled = true
#=╠═╡
begin
	# local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A" .&& an_sig_peaks_ivt.source .=== "motif", :]
	local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A", :]
	println(length(unique(ivt_m6a.peak_id)))
	local refjoin = innerjoin(ivt_m6a, storm,
			  				  on = [:ref_id],
			  				  makeunique = true)
	refjoin = refjoin[between.(refjoin.genomicPos_1, refjoin.start, refjoin.end), :]

	local agg_ivt_peaks = combine(groupby(refjoin, :peak_id),
							   nrow => :overlaps,
							   :GMM_chi2_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_pval)
	storm_res_peaks = agg_ivt_peaks[agg_ivt_peaks.overlaps .== 9 .&& agg_ivt_peaks.min_storm_pval .> 0.2, :peak_id] |> unique
end
  ╠═╡ =#

# ╔═╡ e8cee36d-6886-4e27-806e-b4df1dadff7d
begin
	# local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A" .&& an_sig_peaks_ivt.source .=== "motif", :]
	local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A", :]
	println(length(unique(ivt_m6a.peak_id)))
	local refjoin = innerjoin(ivt_m6a, storm,
			  				  on = [:ref_id],
			  				  makeunique = true)
	refjoin = refjoin[between.(refjoin.genomicPos_1, refjoin.start, refjoin.end), :]

	local agg_ivt_peaks = combine(groupby(refjoin, :peak_id),
							   nrow => :overlaps,
							   :GMM_chi2_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_pval,
							   :KS_intensity_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_ksi_pval,
							   :KS_dwell_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_ksd_pval)
	storm_res_peaks = agg_ivt_peaks[agg_ivt_peaks.overlaps .== 9 .&& agg_ivt_peaks.min_storm_pval .> 0.2 .&& agg_ivt_peaks.min_storm_ksi_pval .> 0.1 .&& agg_ivt_peaks.min_storm_ksd_pval .> 0.1, :peak_id] |> unique
end

# ╔═╡ 59a166ff-402f-44cd-91d4-1f58ca676e33
storm[storm.GMM_chi2_pvalue .=== missing, :]

# ╔═╡ 25c85ec0-af9e-4f74-bb9f-df7325604100
an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :source] |> countmap

# ╔═╡ 4cea1683-7f74-42dd-a0d5-5d66eca41253
map(r -> split(r, "|")[6], an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&& an_sig_peaks_ivt.source .== "HEK293T", :ref_id]) |> unique |> length

# ╔═╡ 58a03eb1-4847-4f48-a86e-9da1233023a3
map(r -> println(r[6]), split.(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&& an_sig_peaks_ivt.source .== "HEK293T", :ref_id] |> unique, "|")) 

# ╔═╡ 56b54a7c-fb57-4aa0-9e84-45ade0d06b5c
storm[storm.ref_id .== "ENST00000372586.4|ENSG00000175283.8|OTTHUMG00000020765.2|OTTHUMT00000054515.2|DOLK-201|DOLK|2074|protein_coding|" .&& storm.pos .== 1709, :]

# ╔═╡ 52a1db73-e957-41fd-a6dc-f483b874e9f2
begin
	local df = an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :]

	

	# for (i, sample) in enumerate(["WT_1", "IVT_1", "WT_2", "IVT_2", "WT_SPK", "IVT_3"])
	# 	violin!(ax, i, df[:, "$(sample)_mod"], label = sample)
	# end
	
	# f

	# df = stack(df, [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod])
	df = stack(df, [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod])

	df[:, :cluster] = ifelse.(occursin.("_mod", df.variable), "modified", "unmodified")
	df[:, :sample] = replace.(df.variable, "_mod" => "", "_unmod" => "")

	local samples = df.sample |> unique
	local sample_ids = Dict(zip(samples, 1:length(samples)))

	local f = Figure()
	local ax = Axis(f[1, 1],
				    xticks = (1:length(samples), samples))

	ylims!(ax, (0, 300))
	
	violin!(ax,
		    map(s -> sample_ids[s], df.sample),
		    df.value,
		    side = ifelse.(df.cluster .== "unmodified", :left, :right),
		    color = ifelse.(df.cluster .== "unmodified", :blue, :orange))
	
	f
end

# ╔═╡ fbd6502e-efa4-4e52-81db-3bb998314b14
an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&& an_sig_peaks_ivt.ref_id .== "ENST00000372586.4|ENSG00000175283.8|OTTHUMG00000020765.2|OTTHUMT00000054515.2|DOLK-201|DOLK|2074|protein_coding|" .&& an_sig_peaks_ivt.pos .== 1709, :]

# ╔═╡ bf177ece-0e74-45f2-8bb0-2b5a2ac2ddbe
begin
	local df = an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :]
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, [:WT_1_mod, :WT_1_unmod]] ./= sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, [:WT_2_mod, :WT_2_unmod]] ./= sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, [:WT_SPK_mod, :WT_SPK_unmod]] ./= sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	counts[:, [:IVT_1_mod, :IVT_1_unmod]] ./= sum.(eachrow(counts[:, [:IVT_1_mod, :IVT_1_unmod]]))
	counts[:, [:IVT_2_mod, :IVT_2_unmod]] ./= sum.(eachrow(counts[:, [:IVT_2_mod, :IVT_2_unmod]]))
	counts[:, [:IVT_3_mod, :IVT_3_unmod]] ./= sum.(eachrow(counts[:, [:IVT_3_mod, :IVT_3_unmod]]))

	# local sums = sum.(eachrow(counts))

	# counts ./ sums

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Per sample proportion of reads")

	f
	# transpose(Array(counts))
end

# ╔═╡ 8d2aafde-d745-4a0d-bcdb-d82a7877292f
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :])

	df[:, :mean_wt_mod] = sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	df[:, :glori] = df.source .== "HEK293T"
	df = sort(df, [:glori, :mean_wt_mod])

	local non_glori_count = sum(.! df.glori)
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, :WT_1] = counts.WT_1_mod ./ sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, :WT_2] = counts.WT_2_mod ./ sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, :WT_SPK] = counts.WT_SPK_mod ./ sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	counts[:, :IVT_1] = counts.IVT_1_mod ./ sum.(eachrow(counts[:, [:IVT_1_mod, :IVT_1_unmod]]))
	counts[:, :IVT_2] = counts.IVT_2_mod ./ sum.(eachrow(counts[:, [:IVT_2_mod, :IVT_2_unmod]]))
	counts[:, :IVT_3] = counts.IVT_3_mod ./ sum.(eachrow(counts[:, [:IVT_3_mod, :IVT_3_unmod]]))

	cols = [:WT_1, :WT_2, :WT_SPK, :IVT_1, :IVT_2, :IVT_3]
	counts = counts[:, cols]
	# local counts = df[:, cols]

	# local sums = sum.(eachrow(counts))


	# counts ./ sums

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
					yticks = ([900, non_glori_count, 2200, nrow(counts)],
							  ["DRACH", "$non_glori_count", "GLORI+", "$(nrow(counts))"]),
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	lines!(ax, [0.5, 6.5], [non_glori_count, non_glori_count], color = :red)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Proportion of reads modified")

	f
	# transpose(Array(counts))
end

# ╔═╡ 3d10cafd-f87c-40a7-ba15-f64e924cc686
md"""
## STORM/IVT
"""

# ╔═╡ 28779e98-1a0e-481f-bb10-c6973bd98976
storm_ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/STORM_IVT_eventalign_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
		LOR_threshold = 0.8)

# ╔═╡ 40173831-c379-4e62-90fc-f7f2e55ff8c6
begin
	peaks_storm_ivt = peaks(storm_ivt, 4)
	sig_peaks_storm_ivt = peaks_storm_ivt[peaks_storm_ivt.predicted .!== missing .&&
						 	  peaks_storm_ivt.predicted .>= -log10(0.01), :]
end

# ╔═╡ 95a3d0c3-09c8-490d-ae5c-993d8ec009e8
an_peaks_storm_ivt = annotate_peaks(peaks_storm_ivt, mod_ref, writer_motifs, writer_mods)

# ╔═╡ 2dbf183b-159e-44bf-87e5-ae9d6517eeb3
an_peaks_storm_ivt

# ╔═╡ 4444b683-1571-4df8-9841-0e28df365dec
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :])

	local r = copy(an_peaks_storm_ivt)
	r[:, :range] = range.(r.pos .- 4, r.pos .+ 4)
	r = flatten(r, :range)

	local j = leftjoin(df, r,
			 		   on = [:ref_id, :pos => :range],
			 		   makeunique = true)

	(j[j.source .=== "HEK293T", :].GMM_chi2_qvalue_1 .<= 0.01) |> countmap
end

# ╔═╡ 757fa776-c741-4eda-9e85-c1110231fcd4
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :])
	
	local refjoin = innerjoin(df, storm_ivt,
			  				  on = [:ref_id],
			  				  makeunique = true)
	refjoin = refjoin[between.(refjoin.genomicPos_1, refjoin.start, refjoin.end), :]

	local df_with_storm_ivt = combine(groupby(refjoin, :peak_id),
							   nrow => :overlaps,
							   :GMM_chi2_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_ivt_pval)

	local j = leftjoin(df, df_with_storm_ivt,
		     		   on = [:peak_id])
	storm_resistant_storm_ivt = j
	storm_resistant_storm_ivt_sig = j[j.min_storm_ivt_pval .< 0.01, :peak_id] |> unique

	j[:, :glori] = j.source .=== "HEK293T"
	j[:, :storm_ivt_sig] = j.min_storm_ivt_pval .< 0.01
	combine(groupby(j, [:glori, :storm_ivt_sig]), nrow => :count)
	# (leftjoin(df[df.source .=== "HEK293T", :], df_with_storm_ivt,
	# 	      on = [:peak_id]).min_storm_ivt_pval .< 0.01) |> countmap
	# df_with_storm_ivt[df_with_storm_ivt.overlaps .== 9 .&& df_with_storm_ivt.min_storm_pval .> 0.2, :peak_id] |> unique |> length
end

# ╔═╡ b8437ca7-b885-4dad-a72d-dc0b9c9a6be5


# ╔═╡ 90c789e2-999c-4e6c-8d49-0926ca080e71
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&&
		map(p -> p in storm_resistant_storm_ivt_sig, an_sig_peaks_ivt.peak_id), :])

	df[:, :mean_wt_mod] = sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	df[:, :glori] = df.source .== "HEK293T"
	df = sort(df, [:glori, :mean_wt_mod])

	local non_glori_count = sum(.! df.glori)
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, :WT_1] = counts.WT_1_mod ./ sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, :WT_2] = counts.WT_2_mod ./ sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, :WT_SPK] = counts.WT_SPK_mod ./ sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	counts[:, :IVT_1] = counts.IVT_1_mod ./ sum.(eachrow(counts[:, [:IVT_1_mod, :IVT_1_unmod]]))
	counts[:, :IVT_2] = counts.IVT_2_mod ./ sum.(eachrow(counts[:, [:IVT_2_mod, :IVT_2_unmod]]))
	counts[:, :IVT_3] = counts.IVT_3_mod ./ sum.(eachrow(counts[:, [:IVT_3_mod, :IVT_3_unmod]]))

	cols = [:WT_1, :WT_2, :WT_SPK, :IVT_1, :IVT_2, :IVT_3]
	counts = counts[:, cols]
	# local counts = df[:, cols]

	# local sums = sum.(eachrow(counts))


	# counts ./ sums

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
					yticks = ([900, non_glori_count, 2200, nrow(counts)],
							  ["DRACH", "$non_glori_count", "GLORI+", "$(nrow(counts))"]),
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	lines!(ax, [0.5, 6.5], [non_glori_count, non_glori_count], color = :red)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Proportion of reads modified")

	f
	# transpose(Array(counts))
end

# ╔═╡ f1795d88-1e1f-4aa7-b245-b2f2cef7db17
begin
	local df = an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :]
end

# ╔═╡ 967bcd01-03a6-4cc5-ac07-b997b65d8a8f
storm_resistant_storm_ivt

# ╔═╡ 322240e8-a405-4ccf-8d94-3c52aa60eb26
glori_stm1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432635_GLORI-293T-STM20-1_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ e698c85d-ab27-4dda-9349-cf679b22a0f5
glori_stm2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432636_GLORI-293T-STM20-2_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ 4415b732-d854-4eb9-9ae6-ed872b087c3e
glori_stm = innerjoin(glori_stm1, glori_stm2,
		  		      on = [:Chr, :Strand, :Sites],
		  			  makeunique = true)

# ╔═╡ 8dd89e2e-78e3-4d2a-b2f8-c9d29e949752
glori_raw1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ c817df72-596c-41be-9d55-80dcbe73c485
glori_raw2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ 76fbfe22-1223-454d-8f6c-0505cdd6c8c4
glori_raw = innerjoin(glori_raw1, glori_raw2,
					  on = [:Chr, :Strand, :Sites],
					  makeunique = true)

# ╔═╡ e33af351-2a7f-45fc-be7e-b36467b3e324
begin
	local t = innerjoin(storm_resistant_storm_ivt, glori_stm,
	          			on = [:chr => :Chr, :strand => :Strand])

	resistant_multi_glori = innerjoin(t[between.(t.Sites, t.start, t.end), :],
			  glori_raw,
			  on = [:chr => :Chr, :strand => :Strand, :Sites],
			  renamecols = "" => "_ctrl")
	resistant_multi_glori[:, :stm_mean_ratio] = (resistant_multi_glori.Ratio .+ resistant_multi_glori.Ratio_1)./2
	resistant_multi_glori[:, :ctrl_mean_ratio] = (resistant_multi_glori.Ratio_ctrl .+ resistant_multi_glori.Ratio_1_ctrl)./2
	resistant_multi_glori
end

# ╔═╡ b2fabd16-f406-4b70-b5ab-343ff87b6a05
resistant_multi_glori[log2.(resistant_multi_glori.stm_mean_ratio ./ resistant_multi_glori.ctrl_mean_ratio) .> -0.2, :]

# ╔═╡ 7b353ba1-cc33-405b-962e-5b7f82d1323b
storm_resistant_storm_ivt.glori |> countmap

# ╔═╡ cb6f8399-addd-4095-98bf-915f1b38ee4d
(data(resistant_multi_glori) * mapping(:ctrl_mean_ratio, :stm_mean_ratio) * visual(Scatter) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; xlabel = "CTRL mean mod ratio", ylabel = "STM mean mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

# ╔═╡ 3cbff0f7-eb89-49e5-88e4-43ec57ab87b3
begin
	local t = innerjoin(glori_raw, glori_stm,
		  				on = [:Chr, :Strand, :Sites],
		  				renamecols = "" => "_stm")
	t[:, :ctrl_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :stm_mean_ratio] = (t.Ratio_stm .+ t.Ratio_1_stm)./2
	t[:, :resistant] = log2.(t.stm_mean_ratio ./ t.ctrl_mean_ratio) .> -0.1
	glori_storm_resistant = t[t.resistant, :]
	println(sum(t.resistant))
	data(t) * mapping(:ctrl_mean_ratio, :stm_mean_ratio, color = :resistant) * visual(Scatter; markersize = 5) |> draw(; axis = (; xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 46ce6ac5-c128-46df-8d38-90776f42861b
innerjoin(glori_storm_resistant, ivt,
		  on = [:Chr => :chr, :Strand => :strand, :Sites => :genomicPos])

# ╔═╡ ca061644-6cf0-4547-bb0d-c21642f76f3f
begin
	local t = innerjoin(storm[storm.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))

	c = round(cor(t.glori_mean_ratio, t.pred_ratio), digits = 2)
	println(c)
	
	(data(t) * mapping(:glori_mean_ratio, :pred_ratio) * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Nanocompore (hard assignment) and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ d636baaf-b65f-4c5e-a9ed-775a7f1c9909
storm_soft = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_soft_clustering/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
	    LOR_threshold = 0.8)

# ╔═╡ 134e34a4-e467-49e6-9c39-12123401d1b7
begin
	local t = innerjoin(storm_soft[storm_soft.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))

	local c = round(cor(t.glori_mean_ratio, t.pred_ratio), digits = 2)
	println(c)
	
	(data(t) * mapping(:glori_mean_ratio, :pred_ratio) * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Nanocompore (soft assignment) and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ d62aed96-41e6-4a01-a614-8f7c2d440b38
begin
	local t = storm_soft
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	t
end

# ╔═╡ d9be422c-3630-4c86-859a-3df1cd725405
function calc_mod_ratio(df, samples)
	mod_cols = filter(n -> any(occursin.(samples, n)) && occursin("_mod", n), names(df))
	all_cols = filter(n -> any(occursin.(samples, n)) && (occursin("_mod", n) || occursin("_unmod", n)), names(df)) 
	sum.(eachrow(df[:, mod_cols])) ./ sum.(eachrow(df[:, all_cols]))
end

# ╔═╡ 1bf15795-87e8-4cf7-92b1-9db7b102450a
begin
	local t = copy(storm_resistant_storm_ivt)
	# t[:, :wt_mod_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	
	t = innerjoin(t, storm_ivt,
		  	      on = [:ref_id, :pos],
			      renamecols = "" => "_stm")
	t[:, :wt_mod_ratio] = calc_mod_ratio(t, ["WT_1", "WT_2", "WT_SPK"])
	t[:, :stm_mod_ratio] = calc_mod_ratio(t, ["STORM_1", "STORM_2", "STORM_SPK"])
	t[(t.wt_mod_ratio .+ t.stm_mod_ratio)./2 .> 0.7 .&& log2.(t.stm_mod_ratio ./ t.wt_mod_ratio) .> -0.2, :]
end

# ╔═╡ 10a0495d-b861-4959-a5c2-1d74b6a0f5c6
# ╠═╡ disabled = true
#=╠═╡
begin
	local t = copy(an_sig_peaks_ivt)
	
	
	t = innerjoin(t, storm_ivt,
			      on = [:ref_id, :pos],
			      renamecols = "" => "_stm")

	t[:, :mod_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	t[:, :mod_ratio_stm] = sum.(eachrow(t[:, [:STORM_1_mod_stm, :STORM_2_mod_stm, :STORM_SPK_mod_stm]])) ./ sum.(eachrow(t[:, [:STORM_1_mod_stm, :STORM_2_mod_stm, :STORM_SPK_mod_stm, :STORM_1_unmod_stm, :STORM_2_unmod_stm, :STORM_SPK_unmod_stm]]))

	t
	# t = innerjoin(t, glori_raw,
	#           	  on = [:chr => :Chr, :strand => :Strand])

	# t = t[between.(t.Sites, t.start, t.end), :]
end
  ╠═╡ =#

# ╔═╡ 0e7cfe6c-a700-4ec1-b8dc-a67f246782b3
begin
	local wtstm = copy(storm[storm.predicted .> 2, :])
	wtstm[:, :wt_mod_ratio] = calc_mod_ratio(wtstm, ["WT_1", "WT_2", "WT_SPK"])
	local wtivt = copy(ivt[ivt.predicted .> 2, :])
	wtivt[:, :wt_mod_ratio] = calc_mod_ratio(wtivt, ["WT_1", "WT_2", "WT_SPK"])
	local t = innerjoin(wtstm, wtivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")
	
	# t = innerjoin(t, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
				  # on = [:ref_id, :pos])
	
	# data(t) * mapping(:wt_mod_ratio => "WT mod. ratio (WT/STORM)", :wt_mod_ratio_ivt => "WT mod. ratio (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
	data(t) * mapping(:GMM_LOR => abs => "|LOR| (WT/STORM)", :GMM_LOR_ivt => abs => "|LOR| (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
end

# ╔═╡ 3b6c3b40-d54e-41c2-afac-6bf8234ae648
begin
	local wtstm = innerjoin(storm, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
						    on = [:ref_id, :pos])
	wtstm[:, :wt_mod_ratio] = calc_mod_ratio(wtstm, ["WT_1", "WT_2", "WT_SPK"])
	local wtivt = copy(ivt[ivt.predicted .>= 0, :])
	wtivt[:, :wt_mod_ratio] = calc_mod_ratio(wtivt, ["WT_1", "WT_2", "WT_SPK"])
	local t = innerjoin(wtstm, wtivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")

	data(t) * mapping(:wt_mod_ratio => "WT mod. ratio (WT/STORM)", :wt_mod_ratio_ivt => "WT mod. ratio (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
end

# ╔═╡ f97950ff-8c2b-4fc4-99e2-51c29bf7b489
begin
	local wtstm = innerjoin(storm, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
						    on = [:ref_id, :pos])
	wtstm[:, :wt_mod_ratio] = calc_mod_ratio(wtstm, ["WT_1", "WT_2", "WT_SPK"])
	local wtivt = copy(ivt[ivt.predicted .>= 0, :])
	wtivt[:, :wt_mod_ratio] = calc_mod_ratio(wtivt, ["WT_1", "WT_2", "WT_SPK"])
	local t = innerjoin(wtstm, wtivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")

	data(t) * mapping(:GMM_LOR => abs => "|LOR| (WT/STORM)", :GMM_LOR_ivt => abs => "|LOR| (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
end

# ╔═╡ b9ee5aa2-ced4-4f71-b43a-bbfd6d5afef9
begin
	local wtstm = innerjoin(storm, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
						    on = [:ref_id, :pos])
	wtstm[:, :storm_mod_ratio] = calc_mod_ratio(wtstm, ["STORM_1", "STORM_2", "STORM_SPK"])
	local stmivt = copy(storm_ivt[storm_ivt.predicted .>= 0, :])
	wtivt[:, :storm_mod_ratio] = calc_mod_ratio(stmivt, ["STORM_1", "STORM_2", "STORM_SPK"])
	local t = innerjoin(wtstm, stmivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")

	data(t) * mapping(:storm_mod_ratio, :storm_mod_ratio_ivt) * visual(Scatter, markersize = 5) |> draw
	
end

# ╔═╡ bd94a0c6-0e84-4699-b449-36e45b2c2c23
storm_resistant_storm_ivt

# ╔═╡ b79d6323-c6be-4b18-ae93-1a1ce48a1db3
function get_read_probs(ref)
	
	db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]


	reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	orig_probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	probs = ifelse.(orig_probs .== -1, missing, orig_probs)
	probs ./= 100
end

# ╔═╡ ca38417c-1d5c-45c6-816c-d2d6898af7db
glori_storm_resistant

# ╔═╡ 3bfb00e8-f963-4c03-825f-71c4343074fa
begin
	eclip_files = ["GSM4970401_meCLIP_HEK293_rep1.m6aList_highConfidence_3054.bed",
				   "GSM4970403_meCLIP_HEK293_rep2.m6aList_highConfidence_1775.bed",
				   "GSM4970405_meCLIP_HEK293_rep3.m6aList_highConfidence_2911.bed"]

	eclip_reps = []
	for fn in eclip_files
		df = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data//GSE147440_RAW/" * fn, header = ["chr", "start", "end", "name", "data", "strand"]))
		df = DataFrames.transform(df, :data => ByRow(x -> map(v -> parse(Float64, v), split(x, "|"))) => [:mod_ratio, :nmod, :nreads])
		df = select!(df, Not(:data))
		push!(eclip_reps, df)
	end
	eclip_reps
end

# ╔═╡ 1c622f87-b763-4f8b-a253-fd38a1fd2830
begin
	local c = copy(sig_peaks_storm)
	c[:, :start] = c.genomicPos .- 4
	c[:, :end] = c.genomicPos .+ 4
	local t = innerjoin(c, eclip_reps[1],
		  				on = [:chr, :strand],
		  				makeunique = true)
	t[between.(t.start_1, t.start, t.end), :]
end

# ╔═╡ 45a5a8f7-a1bb-434c-b882-3b24d3efa412
begin
	local t = innerjoin(storm_resistant_storm_ivt, eclip_reps[1],
		  				on = [:chr, :strand],
		  				makeunique = true)
	t[between.(t.start_1, t.start, t.end), :]
end

# ╔═╡ 96482819-8af4-4ce7-bfda-2e20ac05a66c
begin
	local encore_experiments = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/ENCORE_m6A/experiment_report_2025_7_1_10h_33m.tsv", header = 2))
	encore_experiments = innerjoin(encore_experiments, DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/ENCORE_m6A/ENCORE_m6a_associated_RBPs_filtered.tsv")),
								   on = [:Accession => :accession])
	encore_experiments[:, :bed] = [string(last(split(u, "/"))) 
								   for u in encore_experiments.url]
	local encore_m6a_rbps = []
	for experiment in eachrow(encore_experiments)
		bed = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/ENCORE_m6A/$(experiment.bed)", header = [:chr, :start, :end, :name, :x, :strand, :a, :b, :c, :d], select = [:chr, :start, :end, :name, :strand]))
		bed[:, :rbp] .= experiment.rbp
		bed[:, :cell] .= experiment["Biosample term name"]
		push!(encore_m6a_rbps, bed)
	end
	encore_m6a_rbp_peaks = reduce(vcat, encore_m6a_rbps)
end

# ╔═╡ 6a70233a-4754-49e8-bd81-292d364b4031
begin
	local t = innerjoin(storm_resistant_storm_ivt, encore_m6a_rbp_peaks,
		  				on = [:chr, :strand],
		  				makeunique = true)
	combine(groupby(sort(t[between.(t.genomicPos, t.start_1, t.end_1), :], [:ref_id, :pos]), [:chr, :strand, :genomicPos, :glori]), nrow => :count)
end

# ╔═╡ 445fe109-2361-4bea-bd03-74f3f12dca65
begin
	local t = innerjoin(storm_resistant_storm_ivt,
			  encore_m6a_rbp_peaks,
		  	  on = [:chr, :strand],
		  	  makeunique = true)
	eclip_overlaps = t[between.(t.genomicPos, t.start_1, t.end_1), :peak_id] |> unique

	conf_storm_resistant = copy(storm_resistant_storm_ivt)
	conf_storm_resistant[:, :overlaps_eclip] = map(p -> p in eclip_overlaps, storm_resistant_storm_ivt.peak_id)

	local conf = r -> r.storm_ivt_sig && (r.glori || r.overlaps_eclip)
	conf_storm_resistant = filter(conf, conf_storm_resistant)

	conf_storm_resistant = innerjoin(conf_storm_resistant, storm_ivt,
		 				 			 on = [:ref_id, :pos],
		  						     renamecols = "" => "_stmivt")
	conf_storm_resistant[:, :wt_mod_ratio] = calc_mod_ratio(conf_storm_resistant, ["WT_1", "WT_2", "WT_SPK"])
	conf_storm_resistant[:, :stm_mod_ratio] = calc_mod_ratio(conf_storm_resistant, ["STORM_1", "STORM_2", "STORM_SPK"])
	conf_storm_resistant
end

# ╔═╡ f7a4a1b1-5970-4e1c-a37a-5db901c1115c
combine(groupby(conf_storm_resistant, [:glori, :overlaps_eclip]), nrow => :count)

# ╔═╡ 4605495f-6df1-44fb-8a33-6e84bebe83df
data(conf_storm_resistant) * mapping(:wt_mod_ratio => "WT mod. ratio", :stm_mod_ratio => "STORM mod. ratio") * visual(Scatter) |> draw(; axis = (; xticks = 0:0.1:1, yticks = 0:0.1:1))

# ╔═╡ 38f8fbed-20a9-46ab-92e9-edf3cb693bfa
conf_storm_resistant

# ╔═╡ f84bfc2a-c037-4627-b0c6-4f47d4f45a1b
(
	data(DataFrame(x = [0, 5000], y = [0, 0])) * mapping(:x, :y) * visual(Lines; color = :red) +
	data(conf_storm_resistant) * mapping((:mean_cov, :mean_cov_stmivt) => ((a, b) -> (a+b)/4), (:stm_mod_ratio, :wt_mod_ratio) => ((s, w) -> 100*(s - w))) * visual(Scatter)
	
) |> draw(; axis = (xticks = (0:500:5000, [t % 1000 == 0 ? "$t" : "" for t in 0:500:5000]), yticks = -80:10:80, ytickformat = "{}%", xlabel = "Mean coverage (per condition)", ylabel = "STORM stoichiometry - WT stoichiometry"))

# ╔═╡ a0798ff2-c72f-4de5-b31a-3c8a3bacf117
begin
	local breast_cancer_genes = Set(split("MDC1;SETD2;HSP90AB1;CITED2;EHMT2;FHL1;OLA1;HDLBP;HMGB2;CASC3;RPL8;PHB2;PTPRF;LAPTM4B;CHAF1A;SDF4;DHX16;GLUL;RBM5;PELP1;SRRM2;BRD2;ACTR2;HSP90AA1;SEMA6A;CSNK2A1;HLA-C;RPL13A;DDX54;LSS;SREBF2;WDR77;NME1;RHOB;CDC25B;ILF3;TCP1;RAB35;MYH9;MCM5;PTMA;GAPDH;CREB5;NOTCH2;KDM3B;HM13;GMPR2;MAZ;RPL36A;SRI;HOXB13;HNRNPDL;DVL1;TP53BP2;ABL1;SSR1;FLNA;PRAME;CSK;MAP2K7;SF3B1;RPL19;PAK4;JUN;SLC35A4;CDKN2C;XRCC5;MRPL28;NR2F2;RRP1B;DEK;FKBP1A;SQLE;HNRNPM;SLC6A8;HNRNPK;KRT18;SCD;FASN;HNRNPD;ZYX;GNAS;MAP3K10;TPT1;EIF3A", ";"))

	local t = copy(conf_storm_resistant)
	t[:, :gene] = map(r -> string(split(r, "|")[6]), t.ref_id)
	t[:, :breast_cancer] = map(g -> g in breast_cancer_genes, t.gene)
	t

	(
	data(DataFrame(x = [0, 5000], y = [0, 0])) * mapping(:x, :y) * visual(Lines; color = :red) +
	data(t) * mapping((:mean_cov, :mean_cov_stmivt) => ((a, b) -> (a+b)/4), (:stm_mod_ratio, :wt_mod_ratio) => ((s, w) -> 100*(s - w)), color = :breast_cancer) * visual(Scatter)
	
) |> draw(; axis = (xticks = (0:500:5000, [t % 1000 == 0 ? "$t" : "" for t in 0:500:5000]), yticks = -80:10:80, ytickformat = "{}%", xlabel = "Mean coverage (per condition)", ylabel = "STORM stoichiometry - WT stoichiometry"))
end

# ╔═╡ 49252d74-9a99-4795-9527-3e7eccf20eed
begin
	local refs = Set(unique(conf_storm_resistant.ref_id))

	local resist = combine(groupby(conf_storm_resistant, :ref_id), nrow => :res_mods)
	resist = Dict(zip(resist.ref_id, resist.res_mods))

	local t = copy(sig_peaks_storm[map(r -> r in refs, sig_peaks_storm.ref_id), :])

	local sensitive = combine(groupby(t, :ref_id), nrow => :sensitive_mods)
	sensitive = Dict(zip(sensitive.ref_id, sensitive.sensitive_mods))

	refs = [r for r in refs]
	local counts = DataFrame(ref = refs,
		 		  			 nresist = [get(resist, r, 0) for r in refs],
			  	   			 nsensitive = [get(sensitive, r, 0) for r in refs])

	data(counts) * mapping(:nsensitive => "Number of STORM-sensitive sites") * frequency() |> draw(; axis = (; title = "Transcripts with STM2457-resistant m6A sites by number of sensitive sites", ylabel = "Number of transcripts"))
	
	# combine(groupby(counts, [:nresist, :nsensitive]), nrow => :count)
	# sort(counts, [:nresist, :nsensitive])
	
end

# ╔═╡ e9788e02-f17f-4bd8-b344-caa18863e92e
conf_storm_resistant_20perc = conf_storm_resistant[abs.(conf_storm_resistant.stm_mod_ratio .- conf_storm_resistant.wt_mod_ratio) .< 0.2, :]

# ╔═╡ 3266d3d6-64a2-4936-91ac-059a1ccf76c4
import GenomicAnnotations

# ╔═╡ da65e2fe-54a8-4ca6-b599-411ee2a3e9d4
conf_storm_resistant_gencode = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/gencode41_annotation_for_conf_storm_resistant.gtf", header = [:chr, :source, :feature, :start, :end, :score, :strand, :frame, :attr]))

# ╔═╡ aa699be6-0fa9-468d-847c-b7b3b7f91c72
function calculate_utr_cds_lengths(transcript_id::String, gtf::DataFrame)
	transcript_rows = gtf[occursin.(transcript_id, gtf.attr), :]

	rev = transcript_rows[1, :strand] == '-'

	sort!(transcript_rows, [:start], rev = rev)

    five_prime_utr_len = 0
    cds_len = 0
    three_prime_utr_len = 0

	utr_5p = true

    for row in eachrow(transcript_rows)
        feature_type = row.feature
        start = row.start
        end_ = row.end

        if feature_type == "UTR" && utr_5p
            five_prime_utr_len += end_ - start + 1
			utr_5p = false
        elseif feature_type == "CDS"
            cds_len += end_ - start + 1
        elseif feature_type == "UTR"
            three_prime_utr_len += end_ - start + 1
        end
    end

    return (five_prime_utr_len, cds_len, three_prime_utr_len)
end

# ╔═╡ 007e6893-0718-41c0-8ca0-03fd367f6cc7
begin
	local t = copy(conf_storm_resistant)
	t[:, :len] = map(r -> parse(Int64, split(r, "|")[7]), t.ref_id)
	local transcripts = map(r -> string(split(r, "|")[1]), t.ref_id)
	local all_lens = map(t -> calculate_utr_cds_lengths(t, conf_storm_resistant_gencode),
					 transcripts)
	# t[:, :utr5_len] = first.(lens)
	# t[:, :cds_len] = second.(lens)
	# t[:, :utr3_len] = last.(lens)

	t[:, :meta_pos] = [
		if pos <= utr5_len
			(pos / utr5_len)/4
		elseif pos <= utr5_len + cds_len
			1/4 + ((pos - utr5_len) / cds_len)/2
		else
			3/4 + ((pos - utr5_len - cds_len) / utr3_len)/4
		end
		for (pos, (utr5_len, cds_len, utr3_len)) in zip(t.pos, all_lens)
	]

	t = t[.! isinf.(t.meta_pos), :]
	# t

	# all_lens[5], t[5, :]

	(
		data(DataFrame(x = [1/4, 1/4], y = [0, 12])) * mapping(:x, :y) * visual(Lines, color = :red, linestyle = :dash) +
		data(DataFrame(x = [3/4, 3/4], y = [0, 12])) * mapping(:x, :y) * visual(Lines, color = :red, linestyle = :dash) +
		data(t) * mapping(:meta_pos) * AlgebraOfGraphics.histogram(bins = 50)
	)|> draw(; axis = (; xticks = ([0, 1/4, 3/4, 1], ["TSS", "Start codon", "Stop codon", "TTS"]), ylabel = "Count", yticks = 1:12))
end

# ╔═╡ 7212801b-ee05-4699-b4fb-c02377710e10
calculate_utr_cds_lengths("ENST00000470886.1", conf_storm_resistant_gencode)

# ╔═╡ 6a423afd-67f4-4152-9aa8-7d69c4262f63


# ╔═╡ c02b568d-e607-4729-9cef-1d189ae22bc1
begin
	local filepath = "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_sites.tsv"
	local df = copy(storm_resistant_storm_ivt)
	rename!(df,
			:mean_cov => :coverage,
		    :start => :peak_start,
		    :end => :peak_end,
		    :overlaps => :storm_ivt_overlaps,
		    :min_storm_ivt_pval => :storm_ivt_min_pvalue)
	sort!(df, [:source, :ref_id, :pos])
	

	local f = open(filepath, "w")
	# local data = read(f);
	# seekstart(f);
	local comment = """
	# A table with putative STM2457-resistant sites.
	# These are peaks that were significant when comparing the IVT and WT with Nanocompore, have been inferred as m6A (because they've been detected with GLORI or due to overlapping a DRACH motif), but were not significant when comparing WT to STORM.
	# The first columns are the results from Nanocompore for the IVT/WT comparison, a few columns have been added at the end. Specifically, whether the peak was also significant when comparing the STORM and IVT conditions.
	# More technically:
	# 1. We took the IVT/WT Nonocompore results.
	# 2. We peak called with min distance between peaks equal to 5nt.
	# 3. We selected only the significant peaks at q-value < 0.01 and |GMM_LOR| > 0.8.
	# 4. We selected only the sites where the whole peak (the 9-mer centered at the peak's position) is covered by results in the Nanocompore's comparison of STORM/WT and none of them is significant even at non-strict thresholds of the non-multiple-test-corrected pvalues (GMM_chi2_pvalue > 0.2 && KS_intensity_pvalue > 0.1 && KS_dwell_pvalue 0.1).
	# 5. We combined the resulting peaks with information from Nanocompore's results from the STORM/IVT comparison. We marked the peak as storm_ivt_sig if there was at least one position overlapping the peak that was present in the STORM/IVT results with a GMM_chi2_pvalue < 0.01.
	"""
	write(f, comment)
	close(f)
	df |> CSV.write(filepath, delim = '\t', append = true, writeheader = true)
end

# ╔═╡ 7aa12ba0-e6e7-4c36-9620-7ed443013d40
begin
	rna_stability = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/pelizzola_2016_rna_stability_hek293t.csv"))
	rna_stability[:, :halflife] = (log(2) ./ rna_stability.Deg) ./ 60
	rna_stability
end

# ╔═╡ a11ab145-e6d2-4452-8df8-6976f5f67453
begin
	local t = copy(conf_storm_resistant)
	local s = copy(rna_stability)
	s[:, :gene] = map(r -> string(split(r, ".")[1]), s.Gene)
	t[:, :gene] = map(r -> string(split(split(r, "|")[2], ".")[1]), t.ref_id)
	local j = innerjoin(unique(t[:, [:gene, :mod]]), s[:, [:gene, :halflife]],
			  			on = [:gene])
	println(quantile(j.halflife), mean(j.halflife))
	data(j) * mapping(:halflife => round => "Half-life (h)") * frequency() |> draw
end

# ╔═╡ 046608c2-5c79-4651-bcfe-1b8b965e63c2
begin
	local t = copy(conf_storm_resistant)
	local s = copy(rna_stability)
	s[:, :gene] = map(r -> string(split(r, ".")[1]), s.Gene)
	t[:, :gene] = map(r -> string(split(split(r, "|")[2], ".")[1]), t.ref_id)
	local j = unique(leftjoin(s, t[:, [:gene, :mod]], on = :gene))
	j[:, :mod] = ifelse.(j.mod .=== missing, "No STORM-insensitive m6As", "STORM-insensitive m6As")
	data(j) * mapping(:mod => "", :halflife => "Half-lif (h)") * visual(Violin) |> draw(; axis = (; yticks = 0:5:30))
end

# ╔═╡ e2bf25a2-b5c2-48be-95b3-9ee2905b7b35
begin
	local t = copy(conf_storm_resistant)
	local s = copy(rna_stability)
	s[:, :gene] = map(r -> string(split(r, ".")[1]), s.Gene)
	t[:, :gene] = map(r -> string(split(split(r, "|")[2], ".")[1]), t.ref_id)
	local j = innerjoin(unique(t[:, [:gene]]), s[:, [:gene, :Deg]],
			  			on = [:gene])
	local concentration = (k, t) -> exp(-1*k*t)
	local ngenes = nrow(j)
	local concentrations = combine(
			groupby(j, :gene),
			:gene => (_ -> 0:24) => :t_h,
			:Deg => (k -> concentration.(k, 0:60:24*60)) => :concentration)

	data(combine(groupby(concentrations, :t_h), :concentration => mean => :mean_concentration)) * mapping(:t_h => "Time (h)", :mean_concentration => "Mean concentration") * visual(Lines) |> draw(; axis = (; xticks = 0:4:24))
end

# ╔═╡ 5a1dfded-7f69-4fae-ad5c-9242686026d8
begin
	local t = innerjoin(conf_storm_resistant, encore_m6a_rbp_peaks,
						on = [:chr, :strand],
						renamecols = "" => "_rbp")
	t = t[between.(t.genomicPos, t.start_rbp, t.end_rbp), :]
	data(combine(groupby(t, [:peak_id, :rbp_rbp]), nrow => :count)) * mapping(:rbp_rbp => "RBP") * frequency() |> draw(; axis = (; xticklabelsize = 10))
end

# ╔═╡ 77dbbadb-2588-4db7-8e84-5b375dbec467
conf_storm_resistant

# ╔═╡ ce3be394-4d52-40c6-9df2-4cad1ca55fd7
begin
	local df = copy(conf_storm_resistant)

	df[:, :mean_stm_mod] = sum.(eachrow(df[:, [:STORM_1_mod_stmivt, :STORM_2_mod_stmivt, :STORM_SPK_mod_stmivt]])) ./ sum.(eachrow(df[:, [:STORM_1_mod_stmivt, :STORM_2_mod_stmivt, :STORM_SPK_mod_stmivt, :STORM_1_unmod_stmivt, :STORM_2_unmod_stmivt, :STORM_SPK_unmod_stmivt]]))
	df[:, :glori] = df.source .== "HEK293T"
	df = sort(df, [:glori, :mean_stm_mod])

	local non_glori_count = sum(.! df.glori)
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:STORM_1_mod_stmivt, :STORM_1_unmod_stmivt, :STORM_2_mod_stmivt, :STORM_2_unmod_stmivt, :STORM_SPK_mod_stmivt, :STORM_SPK_unmod_stmivt, :IVT_1_mod_stmivt, :IVT_1_unmod_stmivt, :IVT_2_mod_stmivt, :IVT_2_unmod_stmivt, :IVT_3_mod_stmivt, :IVT_3_unmod_stmivt]
	local counts = df[:, cols] .+ 0.0

	counts[:, :STORM_1] = counts.STORM_1_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_1_mod_stmivt, :STORM_1_unmod_stmivt]]))
	counts[:, :STORM_2] = counts.STORM_2_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_2_mod_stmivt, :STORM_2_unmod_stmivt]]))
	counts[:, :STORM_3] = counts.STORM_SPK_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_SPK_mod_stmivt, :STORM_SPK_unmod_stmivt]]))

	counts[:, :IVT_1] = counts.IVT_1_mod_stmivt ./ sum.(eachrow(counts[:, [:IVT_1_mod_stmivt, :IVT_1_unmod_stmivt]]))
	counts[:, :IVT_2] = counts.IVT_2_mod_stmivt ./ sum.(eachrow(counts[:, [:IVT_2_mod_stmivt, :IVT_2_unmod_stmivt]]))
	counts[:, :IVT_3] = counts.IVT_3_mod_stmivt ./ sum.(eachrow(counts[:, [:IVT_3_mod_stmivt, :IVT_3_unmod_stmivt]]))

	cols = [:STORM_1, :STORM_2, :STORM_3, :IVT_1, :IVT_2, :IVT_3]
	counts = counts[:, cols]
	# local counts = df[:, cols]

	# local sums = sum.(eachrow(counts))


	# counts ./ sums

	println(cor(counts.STORM_1, counts.STORM_2))
	println(cor(counts.STORM_1, counts.STORM_3))
	println(cor(counts.STORM_2, counts.STORM_3))

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Siteas",
					yticks = ([900, non_glori_count, 2200, nrow(counts)],
							  ["DRACH", "$non_glori_count", "GLORI+", "$(nrow(counts))"]),
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	lines!(ax, [0.5, 6.5], [non_glori_count, non_glori_count], color = :red)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Proportion of reads modified")

	f
	# transpose(Array(counts))
end

# ╔═╡ 9f184944-6093-487f-a095-98a8572af620


# ╔═╡ ff570f16-c8fb-4e0d-8851-a0fb93b23645
begin
	local mettl16_kd_peaks = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSE90914_m6a_peaks_mettl16_kd.csv"))
	local t = innerjoin(conf_storm_resistant, mettl16_kd_peaks,
					    on = [:chr, :strand],
					    renamecols = "" => "_m16")
	t = t[between.(t.genomicPos, t.start_m16, t.stop_m16), :]
end

# ╔═╡ b23610f5-240b-4385-955b-216b3221d1d8
conf_storm_resistant

# ╔═╡ 119b9407-c1cd-4eff-a0a6-993b9b05f50b
map(r -> println(split(split(r, "|")[2], ".")[1]), unique(conf_storm_resistant_20perc.ref_id)) |> length

# ╔═╡ f208fdf1-e949-468a-a987-d798baa939ca
map(r -> println(split(r, "|")[6]), unique(conf_storm_resistant.ref_id)) |> length

# ╔═╡ 434ef740-5386-4a64-9ebe-8211e40b0c2b
map(r -> println(">$(split(r.ref_id, "|")[5])-$(r.pos)\n$(r.ref_kmer)"), eachrow(conf_storm_resistant)) |> length

# ╔═╡ f9dbfebc-9c8e-451a-a46c-03b548b630de
md"""
 I took the storm-resistant glori+ sites (43) and ran motif enrichment using https://meme-suite.org/meme/. The DRACH motif was found in 28 sites
"""

# ╔═╡ 43603d8b-82dd-4670-a999-51802c6173df
stm_resist_bp = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_genes_all_GO_biological_process_enirchments.txt", delim = '\t'))

# ╔═╡ d5a93508-e559-4767-a078-54f9bb8fb3bf
stm_resist_cc = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_genes_all_GO_cellular_compartment_enirchments.txt", delim = '\t'))

# ╔═╡ 4ad4579a-4b64-460b-9ec4-a80f86ee391a
stm_resist_mf = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_genes_all_GO_molecular_function_enirchments.txt", delim = '\t'))

# ╔═╡ e86d7b59-449e-4e8f-8ad2-1c708ea15076
begin
	local cc = stm_resist_cc[stm_resist_cc.Benjamini .< 0.01, :]
	local mf = stm_resist_mf[stm_resist_mf.Benjamini .< 0.01, :]

	local Ncc = nrow(cc)
	local Nmf = nrow(mf)

	local counts = zeros((Ncc, Nmf))

	for (c, rcc) in enumerate(eachrow(cc)), (m, rmf) in enumerate(eachrow(mf))
		set_cc = Set(split(rcc.Genes, ", "))
		set_mf = Set(split(rmf.Genes, ", "))

		counts[c, m] = length(intersect(set_cc, set_mf))
	end
	local f = Figure(size = (800, 800))
	local ax = Axis(f[1, 1],
				    xticks = (1:Nmf, map(last, split.(mf.Term, "~"))),
					yticks = (1:Ncc, map(last, split.(cc.Term, "~"))),
				    xticklabelrotation = 0.7)

	heatmap!(ax, transpose(counts), colormap = :heat)

	text!(ax,
		  [string(Int(counts[j, i])) for i in 1:Nmf for j in 1:Ncc],
		  position = [Point2f(x, y) for x in 1:Nmf for y in 1:Ncc],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[j, i] > 280 ? :white : :black for i in 1:Nmf for j in 1:Ncc])
	
	
	f
end

# ╔═╡ 6412e0da-abe4-44fb-9e06-7cce1835d686
begin
	local cc = stm_resist_cc[stm_resist_cc.Benjamini .< 0.01, :]
	local bp = stm_resist_bp[stm_resist_bp.Benjamini .< 0.01, :]

	local Ncc = nrow(cc)
	local Nbp = nrow(bp)

	local counts = zeros((Ncc, Nbp))

	for (c, rcc) in enumerate(eachrow(cc)), (m, rbp) in enumerate(eachrow(bp))
		set_cc = Set(split(rcc.Genes, ", "))
		set_bp = Set(split(rbp.Genes, ", "))

		counts[c, m] = length(intersect(set_cc, set_bp))
	end
	local f = Figure(size = (800, 800))
	local ax = Axis(f[1, 1],
				    xticks = (1:Nbp, map(last, split.(bp.Term, "~"))),
					yticks = (1:Ncc, map(last, split.(cc.Term, "~"))),
				    xticklabelrotation = 0.7)

	heatmap!(ax, transpose(counts), colormap = :heat)

	text!(ax,
		  [string(Int(counts[j, i])) for i in 1:Nbp for j in 1:Ncc],
		  position = [Point2f(x, y) for x in 1:Nbp for y in 1:Ncc],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[j, i] > 280 ? :white : :black for i in 1:Nbp for j in 1:Ncc])
	
	
	f
end

# ╔═╡ 2d9caf6d-154d-408e-8e4d-aef4f1c3ff78
md"""
## SCRATCHPAD
"""

# ╔═╡ 1a8aac68-4308-4700-ad2d-8e6989757416
storm[storm.mean_cov .< 30, :]

# ╔═╡ 4b5089dd-74aa-4dbf-8e26-bdaecc0b5b36
begin
	local non_m6a_storm = an_sig_peaks_storm[an_sig_peaks_storm.mod .!== missing .&&
										 an_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	# local ivt_covered = with_ivt.predicted_1 .!== missing
	# with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	println(sum(with_ivt.predicted_1 .< 2))
	println(sum(with_ivt.predicted_1 .>= 2))

	with_ivt[:, :color] = map(
		s -> if occursin("HEK293T", s)
				:green
			 elseif occursin("HEK293", s)
				:orange
			 elseif s == "motif"
				:lightyellow
			 else
				:lightgrey
			 end,
		with_ivt.source)

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Non-m6A (based on GLORI) mods significant in STORM/WT",
				    xlabel = "STORM/WT -log10(q-value)",
				    ylabel = "IVT/WT -log10(q-value)")

	scatter!(ax,
			 with_ivt.predicted,
			 with_ivt.predicted_1,
			 color = with_ivt.color,
			 strokewidth = 1,
			 strokecolor = :black,
			 markersize = 7)
	text!(ax, with_ivt.predicted .+ rand(nrow(with_ivt)) .* 2, with_ivt.predicted_1 .+ rand(nrow(with_ivt)) .* 2, text = with_ivt.mod, fontsize = ifelse.(with_ivt.predicted .> 60 .|| with_ivt.predicted_1 .> 15, 8, 0))

	Legend(f[2, 1],
		   [PolyElement(color = :green),
			PolyElement(color = :orange),
			PolyElement(color = :lightgrey),
			PolyElement(color = :lightyellow)],
	       ["HEK293T", "HEK293", "Unspecified/Other cell line", "Motif"],
		   "Attribution source",
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 7361c452-d2c9-4793-855e-055911702d8a
begin
	local non_m6a_storm = an_sig_peaks_storm[an_sig_peaks_storm.mod .!== missing .&&
									  	     an_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	# local ivt_covered = with_ivt.predicted_1 .!== missing
	# with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	println(sum(with_ivt.predicted_1 .< 2))
	println(sum(with_ivt.predicted_1 .>= 2))


	local inhib = with_ivt[with_ivt.predicted_1 .< 2, :]

	inhib = combine(groupby(inhib, [:peak_id, :mod]),
					names(inhib) .=> first .=> names(inhib))
	
	# with_ivt[:, :color] = map(
	# 	s -> if occursin("HEK293T", s)
	# 			:green
	# 		 elseif occursin("HEK293", s)
	# 			:orange
	# 		 elseif s == "motif"
	# 			:lightyellow
	# 		 else
	# 			:lightgrey
	# 		 end,
	# 	with_ivt.source)
	local colormap = Dict(zip(unique(with_ivt.mod), 1:length(unique(with_ivt.mod))))

	local f = Figure(size = (1400, 1100))
	local ax = Axis(f[1, 1],
					title = "Putative non-m6A inhibited by STORM",
				    xlabel = "|LOR|",
				    ylabel = "WT/STORM -log10(q-value)")
	xlims!(ax, (0.7, 4.5))
	
	local non_drachs = .! occursin.(writer_motifs["METTL3"], inhib.ref_kmer)
	# println(inhib[non_drachs, :])

	scatter!(ax,
			 abs.(inhib.GMM_LOR),
			 inhib.predicted_raw,
			 color = map(m -> colormap[m], inhib.mod),
			 colormap = :tab10,
			 colorrange = (1, length(colormap)),
			 # strokewidth = 1,
			 # strokecolor = :black,
			 marker = ifelse.(non_drachs,
						 :circle,
						 :star5),
			 markersize = 15)

	
	local xoffsets = abs.(randn(sum(non_drachs))) ./ 3
	local yoffsets = randn(sum(non_drachs)) .* 3
	
	text!(ax,
		  abs.(inhib[non_drachs, :GMM_LOR]) .+ xoffsets,
		  disallowmissing(inhib[non_drachs, :predicted_raw]) .+ yoffsets,
		  text = inhib[non_drachs, :ref_kmer],
		  align = (:left, :center))
		  # fontsize = ifelse.(abs.(inhib.GMM_LOR) .> 1 .|| inhib.predicted_raw .> 15, 12, 0))

	for (x, y, x2, y2) in zip(abs.(inhib[non_drachs, :GMM_LOR]),
							  inhib[non_drachs, :predicted_raw],
							  abs.(inhib[non_drachs, :GMM_LOR]) .+ xoffsets,
							  inhib[non_drachs, :predicted_raw] .+ yoffsets)
		lines!(ax, [x, x2], [y, y2])
	end

	
	Legend(f[2, 1],
		   [PolyElement(color = c,
						colormap=:tab10,
						colorrange = (1, length(colormap)))
			for c in values(colormap)],
		   [m for m in keys(colormap)],
		   "Inferred modification identity",
		   orientation = :horizontal,
		   framevisible = false)
	
	f
	# inhib[non_drachs, :]
end

# ╔═╡ ae549f94-6e59-48f2-aedd-5703b7301202
groupby(storm[storm.chr .== "chr2" .&& storm.genomicPos .> 27200310 .&& storm.genomicPos .< 27200333, [:ref_id, :genomicPos, :pos, :GMM_chi2_pvalue, :GMM_LOR, :GMM_chi2_qvalue]], :ref_id)

# ╔═╡ 0a46f8ed-2a03-4b04-98be-5d85473343ea
groupby(ivt[ivt.chr .== "chr2" .&& ivt.genomicPos .> 27200305 .&& ivt.genomicPos .< 27200333, [:ref_id, :genomicPos, :pos, :GMM_chi2_qvalue, :predicted, :mean_cov, :predicted_raw, :GMM_LOR]], :ref_id)

# ╔═╡ 3c36cc71-db84-4f59-81fe-e4dd7761f72b
annot_sig_peaks[occursin.("UBAC2", annot_sig_peaks.ref_id), :]

# ╔═╡ e9d02b89-5a52-424f-a371-d0ca248bf316
annot_peaks[occursin.("UBAC2", annot_peaks.ref_id), :]

# ╔═╡ 913d87d1-6d0a-4db3-9e46-f8a10d247dd5
peaks_ivt[occursin.("UBAC2", peaks_ivt.ref_id), :]

# ╔═╡ faddf0dd-f931-4379-a13a-021bd5b23cae
peaks(ivt[occursin.("ENST00000403766.8", ivt.ref_id), :], 4)

# ╔═╡ c4a7d77d-68cd-40ef-a168-6dbf9dc410b8
ivt[occursin.("ENST00000403766.8", ivt.ref_id), :][1420:1530, [:pos, :GMM_chi2_pvalue, :GMM_chi2_qvalue, :GMM_LOR, :predicted]]

# ╔═╡ d3567161-530b-476a-b339-a0b2d415da7e
argmaxima(ivt[occursin.("ENST00000403766.8", ivt.ref_id), :][1420:1530, :predicted])

# ╔═╡ ff2cea1e-4bb4-4cd2-a1f4-a68cfeb6efde
begin
	local non_m6a_storm = an_sig_peaks_storm[an_sig_peaks_storm.mod .!== missing .&&
										 an_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	# local ivt_covered = with_ivt.predicted_1 .!== missing
	# with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	println(sum(with_ivt.predicted_1 .< 2))
	println(sum(with_ivt.predicted_1 .>= 2))


	local promote = with_ivt[with_ivt.predicted_1 .>= 2, :]

	promote = combine(groupby(promote, [:peak_id, :mod]),
					  names(promote) .=> first .=> names(promote))
	
	# with_ivt[:, :color] = map(
	# 	s -> if occursin("HEK293T", s)
	# 			:green
	# 		 elseif occursin("HEK293", s)
	# 			:orange
	# 		 elseif s == "motif"
	# 			:lightyellow
	# 		 else
	# 			:lightgrey
	# 		 end,
	# 	with_ivt.source)
	local colormap = Dict(zip(unique(with_ivt.mod), 1:length(unique(with_ivt.mod))))

	local f = Figure(size = (1400, 1100))
	local ax = Axis(f[1, 1],
					title = "Putative non-m6A promoted by STORM",
				    xlabel = "|LOR|",
				    ylabel = "WT/STORM -log10(q-value)")
	xlims!(ax, (0.7, 4.5))
	
	local non_drachs = .! occursin.(writer_motifs["METTL3"], promote.ref_kmer)

	scatter!(ax,
			 abs.(promote.GMM_LOR),
			 promote.predicted_raw,
			 color = map(m -> colormap[m], promote.mod),
			 colormap = :tab10,
			 colorrange = (1, length(colormap)),
			 # strokewidth = 1,
			 # strokecolor = :black,
			 marker = ifelse.(non_drachs,
						 :circle,
						 :star5),
			 markersize = 15)

	
	local xoffsets = abs.(randn(sum(non_drachs))) ./ 3
	local yoffsets = randn(sum(non_drachs)) .* 3
	
	text!(ax,
		  abs.(promote[non_drachs, :GMM_LOR]) .+ xoffsets,
		  disallowmissing(promote[non_drachs, :predicted_raw]) .+ yoffsets,
		  text = promote[non_drachs, :ref_kmer],
		  align = (:left, :center))

	for (x, y, x2, y2) in zip(abs.(promote[non_drachs, :GMM_LOR]),
							  promote[non_drachs, :predicted_raw],
							  abs.(promote[non_drachs, :GMM_LOR]) .+ xoffsets,
							  promote[non_drachs, :predicted_raw] .+ yoffsets)
		lines!(ax, [x, x2], [y, y2])
	end

	
	Legend(f[2, 1],
		   [PolyElement(color = c,
						colormap=:tab10,
						colorrange = (1, length(colormap)))
			for c in values(colormap)],
		   [m for m in keys(colormap)],
		   "Inferred modification identity",
		   orientation = :horizontal,
		   framevisible = false)
	
	f
	# promote[non_drachs, :]
end

# ╔═╡ 26299f99-3851-4130-84b6-eaac2c59fbd8
Set.(skipmissing.(zip([ifelse.(occursin.(motif, annot_sig_peaks.ref_kmer), writer_mods[writer], missing) for (writer, motif) in writer_motifs]...)))

# ╔═╡ da893d5b-1fb4-41f1-a8c9-16c09196739a
begin
	local motif_supported = zip(
	[ifelse.(occursin.(motif, annot_sig_peaks.ref_kmer),
			 writer_mods[writer],
			 missing)
	 for (writer, motif) in writer_motifs]...) |> 
	d -> skipmissing.(d) |>
	d -> Set.(d)
	# map(s -> annot_sig_peaks.mod in s, motif_supported)
	map(((m, s),) -> m in s, zip(annot_sig_peaks.mod, motif_supported)) |> countmap
end

# ╔═╡ f474351b-6b02-4402-a722-9bd03ac6b1b9
combine(groupby(annot_sig_peaks, [:peak_id, :mod]),
		names(annot_sig_peaks) .=> first)

# ╔═╡ 78ac5d01-2ded-41c5-9b6b-64dfc94a2141
begin
	local ambig_peaks = filter(r -> r.unique_mods > 1, eachrow(combine(groupby(annot_sig_peaks, :peak_id),
		:mod => length∘unique => :unique_mods))).peak_id |> Set
	annot_sig_peaks[map(p -> p in ambig_peaks, annot_sig_peaks.peak_id), :]
end

# ╔═╡ 69b7aa3a-7d0f-4af9-8c71-f1479d60225b
first(annot_sig_peaks)

# ╔═╡ 577ba306-497a-4cfc-8ea4-bdbc8afb2076
annot_sig_peaks.mod |> unique

# ╔═╡ f621f4d2-eb36-4733-84a1-9a77f3d88cb9
begin
	local non_m6a_storm = annot_sig_peaks_storm[annot_sig_peaks_storm.mod .!== missing .&&
										 annot_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)

	with_ivt[with_ivt.predicted_1 .> 75, :]
end

# ╔═╡ 31aa6d10-3ba3-4af0-9458-2c6a60aa01d1
mod_storm = annotate_mods(storm)

# ╔═╡ 552238b4-26f1-4a3c-9706-65410a8ea4d1
# mod_peak_storm = simple_peaks(mod_storm, 4)
# mod_peak_storm = peaks(mod_storm, 4)

# ╔═╡ 966e302e-5ab9-4806-ac72-244eb4be2aae
mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&& mod_peak_storm.mod_type .!== missing, :]

# ╔═╡ eb32157e-74aa-4e61-9b4c-3b57bfa1cb0a
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = leftjoin(non_m6a_storm, annot_peaks,
			 				  on = [:ref_id, :pos],
							  makeunique = true)
	local ivt_covered = with_ivt.predicted_1 .!== missing
	with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Non-m6A mods significant in STORM/WT",
				    xlabel = "STORM/WT -log10(q-value)",
				    ylabel = "IVT/WT -log10(q-value)")

	scatter!(ax, with_ivt.predicted, with_ivt.predicted_1, markersize = 5, color = ivt_covered)
	# text!(ax, with_ivt.predicted .+ rand(nrow(with_ivt)) .* 2, with_ivt.predicted_1 .+ rand(nrow(with_ivt)) .* 2, text = with_ivt.mod_type, fontsize = 2 .*log2.(with_ivt.predicted))
	f
end

# ╔═╡ 734aebf9-6da9-4e6e-bc79-2905bcd75e7a
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
end

# ╔═╡ a6b67788-7967-47fc-9d47-bc4e7c202f5a
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = innerjoin(storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	local refs = Set(with_ivt.ref_id)
	mod_peak_storm[map(ref -> ref in refs, mod_peak_storm.ref_id) .&& mod_peak_storm.mean_m6A_level .!== missing .&& occursin.("ENST00000621592.8", mod_peak_storm.ref_id), :]
end

# ╔═╡ 9a545506-f35d-48b8-aeff-ca24d8ab6442
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = leftjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	local refs = Set(with_ivt.ref_id)
	map(r -> println(split(r, "|")[6]), [r for r in refs])
	nothing
end

# ╔═╡ 3b9fd69d-ca22-4f05-a563-3492ddd66285
# ╠═╡ disabled = true
#=╠═╡
begin
	local a = glori[:, [:chrom, :start, :end, :strand]]
	rename!(a, [:chr, :start, :end, :strand])
	a[:, :name] .= "m6A"
	local b = rmbase_all[:, [:chr, :start, :end, :strand, :mod_type]]
	rename!(b, [:chr, :start, :end, :strand, :name])
	local c = hela_ac4C[:, [:chr, :start, :end, :strand]]
	c[:, :name] .= "ac4C"
	local d = atoi[:, [:chr, :start, :end, :strand]]
	d[:, :name] .= "A-I"
	local combined = vcat(a, b, c, d)
	combined[:, :score] .= 0
	combined[:, [:chr, :start, :end, :name, :score, :strand]] |> CSV.write("/home/mzdravkov/mods.bed", delim = '\t', header = false)
end
  ╠═╡ =#

# ╔═╡ 4167bc84-b881-4eb0-b140-aa35c0acd365
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(ivt_m6a_motif))
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT significant DRACH not in GLORI",
				    xlabel = "IVT/WT -log10(P-value)",
				    ylabel = "STORM/WT -log10(P-value)")
	println(sum([x >= 2 for x in t.predicted_1 if x !== missing]))
	scatter!(ax, t.predicted, t.predicted_1, markersize = 5)
	f
end

# ╔═╡ c92f80c5-6f47-47c3-92d8-8d1cc9fd6c42
# ╠═╡ disabled = true
#=╠═╡
storm_ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/STORM_IVT_eventalign_v2_0_0_mapq_gt_0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false)
  ╠═╡ =#

# ╔═╡ 50443a47-3600-4733-8b77-1d56717b00b4
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(ivt_m6a_motif), " ", nrow(t))
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT significant DRACH not in GLORI",
				    xlabel = "IVT/WT -log10(P-value)",
				    ylabel = "IVT/STORM -log10(P-value)")
	println(sum([x >= 2 for x in t.predicted_1 if x !== missing]))
	println(sum([x <= 2 for x in t.predicted_1 if x !== missing]))
	scatter!(ax, t.predicted, t.predicted_1, markersize = 5)
	f
end


# ╔═╡ 3060c2aa-078b-43e3-8451-b17642f470eb
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	local res = DataFrames.combine(groupby(t, [:wt_storm, :ivt_storm]), nrow => :count)
	res[:, :ivt_wt] .= true
	res[:, [:ivt_wt, :wt_storm, :ivt_storm, :count]]
end


# ╔═╡ a65930bd-58ca-4326-a04d-0416614b7abb
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	t[t.wt_storm .=== true .&& t.ivt_storm .=== true, :]
end


# ╔═╡ 5b380673-9c93-4e6f-a158-efe469ab8c52
mod_storm_ivt = annotate_mods(storm_ivt)

# ╔═╡ e265dcee-67b7-4238-987a-19242ffe9f03
mod_peak_storm_ivt = peaks(mod_storm_ivt, 4)

# ╔═╡ d3dd2610-f200-4b22-bfb3-6d28df96de13
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, sort(expand_labels(mod_peak_storm, :pos, 4, 4), [:ref_id, :pos]),
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, sort(expand_labels(mod_peak_storm_ivt, :pos, 4, 4), [:ref_id, :pos]),
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	DataFrames.combine(groupby(t, [:wt_storm, :ivt_storm]), nrow => :count)
end


# ╔═╡ e6e36c0d-f522-4959-a979-20573ffda8c2
begin
	# local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
	# 	  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(storm, ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted .>= 2
	t[:, :ivt_wt] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	local res = DataFrames.combine(groupby(t, [:ivt_wt, :wt_storm, :ivt_storm]), nrow => :count)
	# res[:, :ivt_wt] .= true
	res[:, [:ivt_wt, :wt_storm, :ivt_storm, :count]]
end


# ╔═╡ 55243de1-4f6e-4184-a1c3-40b3e24b1a1b
function genomic_peaks(df, radius)
	all_positions = DataFrames.combine(groupby(df, [:chr, :strand]), :genomicPos => (p -> minimum(p):maximum(p)) => :genomicPos)
	all_positions = leftjoin(all_positions, df,
			 				 on = [:chr, :strand, :genomicPos])
	all_positions = sort(all_positions, [:chr, :strand, :genomicPos])
	# all_positions[:, :predicted] = coalesce.(all_positions.predicted, 0)
	peak_indices = argmaxima(all_positions.predicted, radius, strict = false)
	all_positions[peak_indices, :]
end

# ╔═╡ ee927604-3622-4dd0-8f02-b6838e036fc1
mod_all_ivt = annotate_mods(ivt, only_significant = false)

# ╔═╡ 2254d304-dc0a-4dae-87d7-f3cd18ebe152
begin
	local hek293T_mods = mod_all_ivt[(mod_all_ivt.cell_list .!== missing .&& occursin.("HEK293T", mod_all_ivt.cell_list)) .|| (mod_all_ivt.mod_type .!== missing .&& mod_all_ivt.mod_type .== "m6A"), :]
	# sum(hek293T_mods.predicted .> -log10(0.01)) / nrow(hek293T_mods)
	for mod in unique(hek293T_mods.mod_type)
		t = hek293T_mods[hek293T_mods.mod_type .!== missing .&& hek293T_mods.mod_type .== mod, :]
		println(mod, " ", sum(t.predicted .> -log10(0.01)) / nrow(t))
	end
end

# ╔═╡ 4e1e1f43-de2e-4584-ac04-31bcfd7036ad
mod_all_ivt_p12 = annotate_mods(ivt_p12, only_significant = false)

# ╔═╡ 5fbffc42-e89e-47d8-b493-e297ad66ca0b
begin
	local hek293T_mods = mod_all_ivt_p12[(mod_all_ivt_p12.cell_list .!== missing .&& occursin.("HEK293T", mod_all_ivt_p12.cell_list)) .|| (mod_all_ivt_p12.mod_type .!== missing .&& mod_all_ivt_p12.mod_type .== "m6A"), :]
	for mod in unique(hek293T_mods.mod_type)
		t = hek293T_mods[hek293T_mods.mod_type .!== missing .&& hek293T_mods.mod_type .== mod, :]
		println(mod, " ", sum(t.predicted .> -log10(0.01)) / nrow(t))
	end
end

# ╔═╡ 8201e87f-19c7-4b61-bfb8-46ed75bfedcc
# genomic_peaks(ivt, 4)

# ╔═╡ 3dd015d7-3c68-4b7b-9331-d6bd51c3702d
rmbase_all.mod_type |> countmap

# ╔═╡ e9c15162-c3d4-449f-88db-0280b5874a69
md"""
## Indirect recall analysis via GLORI
"""

# ╔═╡ 70267a1f-34da-49c6-a09f-1c214d0dfecb
begin
	rna004_500 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_500_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=true)
	binned_rna004_500 = annotate_binned_results(bin(rna004_500, BIN_SIZE),
										        binned_glori)
	:ok
end

# ╔═╡ c15eca41-5fda-4869-a298-4991fa4ad72d
begin
	ivt_500 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_depth_500_v2_0_0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=true)
	binned_ivt_500 = annotate_binned_results(bin(ivt_500, BIN_SIZE),
										     binned_glori)
	:ok
end

# ╔═╡ a0e6818d-4d4d-4644-84ab-a4c45c9dcb26
binned_storm_ivt_common = innerjoin(binned_rna004_500, binned_ivt_500,
							  		on = [:chr, :strand, :bin],
							  		makeunique = true)


# ╔═╡ 36477658-b2f5-4ff4-b57a-72367ed61819
begin
	local f = Figure()
	local ax = Axis(f[1, 1], aspect = 1)
	xlims!(ax, (0, 200))
	ylims!(ax, (0, 200))
 	local t = binned_storm_ivt_common[abs.(binned_storm_ivt_common.mean_cov .- binned_storm_ivt_common.mean_cov_1) .< 50 .&& binned_storm_ivt_common.modified .== 1, :]
 	scatter!(ax, t.predicted, t.predicted_1, markersize = 5)
	f
end

# ╔═╡ 4593f1eb-c752-4043-98d9-1a00954d738e
scatter(binned_storm_ivt_common.mean_cov, binned_storm_ivt_common.mean_cov_1, markersize = 5)

# ╔═╡ 0f4380ab-a77f-4413-823f-c16044efd721
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT MA plot - positions",
				    xlabel = "Total coverage in both conditions",
					yticks = -9:1:9,
				    ylabel = "LOR",
				    xscale = log10)
	local df = ivt[abs.(ivt.GMM_LOR) .> 0.01, :]
	# local m = df.GMM_LOR
	# local a = log10.(df.mean_cov)
	scatter!(ax, df.mean_cov, df.GMM_LOR,
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 91d20be7-a7a4-4120-8756-e8e5656ec804
annot_peaks

# ╔═╡ 0e7d33ba-f19d-45b3-8606-0fe9eb1845c3
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1])
	local df = annot_peaks[annot_peaks.GMM_chi2_pvalue .!== missing, :]
	local expected = -log10.((1:nrow(df))./nrow(df))
	scatter!(ax,
			 expected,
			 -log10.(sort(df.GMM_chi2_pvalue .+ 1e-310)),
			 markersize = 5)
	lines!(ax, [0, 4], [0, 4], color = :red)
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	f
end
  ╠═╡ =#

# ╔═╡ 05558267-7c22-4480-99f7-847b02a273cd
import StatsPlots

# ╔═╡ a1052731-3d45-42e5-b851-9dc0a8cf8177
import Distributions

# ╔═╡ c71344c3-4faa-4849-a35a-6f511eaf4f34
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT/WT qq-plot",
				    xlabel = "Expected p-value (-log₁₀-scale)",
				    ylabel = "Observed p-value (-log₁₀-scale)")
	local df = peaks_ivt[peaks_ivt.GMM_chi2_pvalue .!== missing, :]
	local expected = -log10.((1:nrow(df))./nrow(df))
	local observed = -log10.(sort(df.GMM_chi2_pvalue .+ 1e-310))
	
	qqplot!(ax,
			expected,
			observed,
		    qqline = :identity,
		    markercolor = :blue,
			markersize = 5,
		    color = :red)
	
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	
	local median_observed_chisq = median(quantile.(Distributions.Chisq(n), df.GMM_chi2_pvalue))
	local expected_chisq = quantile(Distributions.Chisq(n), 0.5)
	local λ = median_observed_chisq/expected_chisq
	println("Lambda gc: ", median_observed_chisq/expected_chisq)

	text!(ax, 0.1, 280; text = "λgc = $(round(λ; digits = 4))")
	
	f
end

# ╔═╡ 7c1bce27-3c08-4b0a-9364-9e992a793a97
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "STORM/WT qq-plot",
				    xlabel = "Expected p-value (-log₁₀-scale)",
				    ylabel = "Observed p-value (-log₁₀-scale)")
	local df = peaks_storm[peaks_storm.GMM_chi2_pvalue .!== missing, :]
	local expected = -log10.((1:nrow(df))./nrow(df))
	local observed = -log10.(sort(df.GMM_chi2_pvalue .+ 1e-310))
	
	qqplot!(ax,
			expected,
			observed,
		    qqline = :identity,
		    markercolor = :blue,
			markersize = 5,
		    color = :red)
	
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	
	local median_observed_chisq = median(quantile.(Distributions.Chisq(n), df.GMM_chi2_pvalue))
	local expected_chisq = quantile(Distributions.Chisq(n), 0.5)
	local λ = median_observed_chisq/expected_chisq
	println("Lambda gc: ", median_observed_chisq/expected_chisq)

	text!(ax, 0.1, 280; text = "λgc = $(round(λ; digits = 4))")
	
	f
end

# ╔═╡ d3fc7047-56a6-4dce-96d3-61189556299e
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Mod. associations qq-plot",
				    xlabel = "Expected p-value (-log₁₀-scale)",
				    ylabel = "Observed p-value (-log₁₀-scale)")
	local df = comod
	local expected = -log10.((1:nrow(df))./nrow(df))
	local observed = -log10.(sort(df.pvalue))
	
	qqplot!(ax,
			expected,
			observed,
		    qqline = :identity,
		    markercolor = :blue,
			markersize = 5,
		    color = :red)
	
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	
	local median_observed_chisq = median(quantile.(Distributions.Chisq(n), df.pvalue))
	local expected_chisq = quantile(Distributions.Chisq(n), 0.5)
	local λ = median_observed_chisq/expected_chisq
	println("Lambda gc: ", median_observed_chisq/expected_chisq)

	text!(ax, 0.1, 280; text = "λgc = $(round(λ; digits = 4))")
	
	f
end

# ╔═╡ ce1658a1-67fd-4743-90f3-9dc914d857be
StatsPlots.qqplot(Distributions.Uniform(), -log10.(annot_peaks[annot_peaks.GMM_chi2_pvalue .!== missing, :GMM_chi2_pvalue]), qqline = :R)

# ╔═╡ 749d91db-cebd-43f2-a22b-248fd9b8587d
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT q-value / coverage",
				    xlabel = "Total coverage in both conditions",
					# yticks = -9:1:9,
				    ylabel = "-log10(q-value)",
				    xscale = log10)
	local df = annot_peaks
	# local m = df.GMM_LOR
	# local a = log10.(df.mean_cov)
	scatter!(ax, df.mean_cov, -log10.(df.GMM_chi2_qvalue),
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 730a04cc-b8db-406a-94ad-f951baf08adf
# ╠═╡ disabled = true
#=╠═╡
begin
	local glori1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv"))
	local glori2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv"))

	local common = innerjoin(glori1, glori2,
			  				 on = [:Chr, :Strand, :Sites],
			  				 renamecols = "_1" => "_2")
	# Convert to 0-based indexing to conform to the BED file format standard
	common[:, :Sites] .-= 1
	common[:, :ratio] = (common.Ratio_1 .+ common.Ratio_2)./2
	common[:, :end] = common.Sites .+ 1

	rename!(common,
		    :Chr => :chr,
		    :Sites => :start,
		    :Strand => :strand,
		    :Gene_1 => :gene)
	
	common = common[:, [:chr, :start, :end, :gene, :ratio, :strand]]
	common |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ 8805231b-7e9f-4dec-a999-ac8e750a49b9
# ╠═╡ disabled = true
#=╠═╡
begin
	local glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
	glori[:, :modified] .= 1
	glori = glori[:, [:chr, :strand, :start, :modified, :ratio]]
	local ref002 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx.tsv"))
	local annot_ref = leftjoin(ref002, glori,
							   on = [:chr, :strand, :genomicPos => :start])
	annot_ref[:, :modified] = ifelse.(annot_ref.modified .=== missing, 0, 1)
	annot_ref[:, :ratio] = ifelse.(annot_ref.ratio .=== missing, 0, annot_ref.ratio)
	sort(annot_ref, [:ref_id, :pos]) |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ a6b2b943-c350-4d99-9124-cfaaed1ddd85
# ╠═╡ disabled = true
#=╠═╡
begin
	local glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
	glori[:, :modified] .= 1
	glori = glori[:, [:chr, :strand, :start, :modified, :ratio]]
	local ref004 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx.tsv"))
	local annot_ref = leftjoin(ref004, glori,
							   on = [:chr, :strand, :genomicPos => :start])
	annot_ref[:, :modified] = ifelse.(annot_ref.modified .=== missing, 0, 1)
	annot_ref[:, :ratio] = ifelse.(annot_ref.ratio .=== missing, 0, annot_ref.ratio)
	sort(annot_ref, [:ref_id, :pos]) |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═41818e42-3cb8-4d48-93d9-b53e0eea7873
# ╠═164f8adc-3bba-11f0-3c64-19ee3bf9097e
# ╠═93db1a5d-4de2-482b-a4e2-540db2059689
# ╠═bc77ac69-225a-4f81-b219-15683bb5eb4a
# ╠═82dba028-3fbe-44e2-9e8f-377f738363a7
# ╠═c88509a9-5a1d-4895-9ae7-5174cc199772
# ╠═dedfd9b6-9b9e-4149-834a-547cbac99873
# ╠═1795bc42-bdc8-4685-bb52-39271581bc06
# ╠═88dac2fb-7062-4e75-aaaf-c2596d90b299
# ╠═5a2ad3e8-9f14-49ad-a78b-6fdd4be05220
# ╠═8445b562-cb03-4a7e-911b-5dc44778ee5c
# ╠═8b930a4c-cef4-49e6-8a76-d59d4df8a349
# ╠═a6d02ee4-2e77-4520-b140-36fb2150b1b0
# ╠═b89d56f0-5104-4b0d-a1ae-6090d262db25
# ╠═566437e3-1abe-49bf-a262-e9a4c702432e
# ╠═b5d8eafa-9ea4-4c59-bef8-3070126bb095
# ╠═c9fabd63-2131-41fc-a50a-bda27319bc6d
# ╠═43a5d6b6-652d-403c-b4fd-945ba2b80ead
# ╠═5abecd8c-e06c-477c-8938-545f64591385
# ╠═b0a1b932-7f0b-40ac-a2d6-e4867d2b54da
# ╠═7ca021f1-532f-4a3c-98b6-0ff08c27d01e
# ╠═33a21a05-dac1-40e0-97a4-8dae37fa8bfc
# ╠═9232e724-d8df-4970-931c-f1c11f0ce859
# ╠═82a9954d-2c23-46bc-a0bb-c2fc078eabb1
# ╠═16109b94-06e5-4285-9333-5b8a98761ef7
# ╠═57a4e1e9-59cb-4180-b735-9da13feb64a6
# ╠═fbd57196-5c5b-49a9-a5f5-a39d25a33c87
# ╠═d41911f6-ab79-4a7f-b245-76ef3cfcd1e6
# ╠═5aaee16e-1121-4bf7-b4c6-60941b62414c
# ╠═4016b836-3669-4098-9725-d3708755db40
# ╠═390466c1-9b31-460b-82cd-676d85ca41da
# ╠═081adcae-8bfe-4c6f-9a17-1ddd8e9dd87e
# ╠═3a91b67e-8efa-4e84-a519-6426594ea30e
# ╠═a4d697a0-71b6-4e01-879e-6a91d9060b8d
# ╠═23ce8398-5f79-47d4-afe9-c46818011b9a
# ╠═25f512fd-e5f8-40d7-a7eb-2d04e18b884b
# ╠═6669794a-cebf-46f7-8d1b-5bb91ccb7280
# ╠═c71344c3-4faa-4849-a35a-6f511eaf4f34
# ╠═16fff5ab-601e-42d7-9f74-1670690733bb
# ╠═86179fb3-d3d0-4f5b-a8d3-14fb43dac834
# ╠═5fde98ae-bdf2-4dcb-90a0-f5017b97262a
# ╠═d90664a5-8e7d-469b-8b05-b103d829f672
# ╠═ab03ca5a-3186-493e-9e3e-13367a85e624
# ╠═cefdac37-7b63-4c8c-ad92-41ca18b78452
# ╠═140e4f5f-3c19-4c94-8829-4a46a95eb6cd
# ╠═9513c89b-dd9e-40dd-91be-c480438e3f07
# ╠═a3f63803-5803-4942-8dce-b1b90b86f6a1
# ╠═56ed77ee-fcfb-498a-ae5c-f9e77ceb971a
# ╠═e429bac0-0add-4e15-8a7a-6ab22a5c9943
# ╠═2b618c14-b478-4d5b-abfc-19bbeb073fde
# ╠═ea86f47c-ed99-4108-bb5c-c9ee225dc7ef
# ╠═7fcfeb9f-5e35-4408-8cfb-b24c5f52c6a2
# ╠═eb2f9c7c-b5bb-493d-8cc1-abdc452e646d
# ╠═ed261b2a-459a-440e-99ad-91afdeda2952
# ╠═dc85be9a-f0af-4b8b-924e-758b3ce45a7a
# ╠═bfe35e43-60fd-4ba9-9f05-f5a7f868e011
# ╠═a81dc14a-fd54-417b-bcd7-84615a62aa37
# ╠═a2dbec3c-961f-4904-961d-be2d9d538720
# ╠═7c1bce27-3c08-4b0a-9364-9e992a793a97
# ╠═bf514f88-b2f4-4a35-8c4a-08fec382e74c
# ╠═c1002495-1c24-450c-9585-1aeb7d588bb5
# ╠═823b1914-ab79-499a-98ec-348f24b0107a
# ╠═80836ca5-66c5-4f9c-a815-9c57230f6948
# ╠═66983a7b-8a6a-4b9f-91ba-bb40b5b5ebce
# ╠═6c7f7c71-ed77-46b1-97a7-0d8b9c8d758e
# ╠═29c659f1-4e5a-4c65-aa53-33dea7d730f3
# ╠═d3fc7047-56a6-4dce-96d3-61189556299e
# ╠═1d820491-8699-4324-b113-b4d7ac0a819a
# ╠═bc42854a-ecbe-46c2-b289-114b12cf6fe9
# ╠═26d0dc17-cd92-47d4-bec0-33ad54c002cd
# ╠═c2cf87ff-d11c-4bb3-b44d-6a16f4040c88
# ╠═fd549de3-8c1f-4358-8ec0-48fdd85d0bd1
# ╠═b37f9ca0-5a53-4716-9589-9b57268f72b2
# ╠═67638f93-2a44-4772-9c35-2e728cebc277
# ╠═c99dbb5d-77da-4c2e-9b7a-874f5ead6574
# ╠═409ec9ef-22c6-4192-a53c-382d059ff299
# ╠═664f5ea5-378a-45cc-a6cf-cab779b7f933
# ╠═d1c2e993-4ba4-46a5-840c-38bfbf397858
# ╠═cb81e3ba-5eea-4169-a4ff-5b7433e8b171
# ╠═07f53334-5428-4890-b02a-9e3c5de7b22c
# ╠═40bbff56-4264-4d9a-af0e-59485c923f92
# ╠═0ee4e9af-4b1e-41fb-8746-12ff51d20e27
# ╠═9503a9e1-8c7e-4dfe-aeed-015d7ec1095c
# ╠═9fbe17fb-5b2a-4cda-a70d-1598de0a6ab0
# ╠═58e326c7-e5b5-4f97-884e-7321a8084139
# ╠═e30b68c9-9aaa-44f7-86d9-c8e1ab9254dd
# ╠═3dc430dc-a218-42bf-90da-8ba8ba6dd56a
# ╠═46b1d7d3-2e9e-4d19-91d8-07150e981ef0
# ╠═809b0c5c-1538-4bbe-bb4d-293f6058b2a6
# ╠═beb418fd-d39b-4c5e-bf9d-98e936b9c819
# ╠═de7c2cc7-7931-4290-97be-39b24634d755
# ╠═088e365f-05f2-429c-9fa2-763ddf30b82b
# ╠═ea219179-f224-4aa3-bb1d-fdf494742363
# ╠═bcc0fec3-9925-4c52-8411-a59583f9171e
# ╠═b117d0ed-b44b-46a6-b9ea-204e1d5195e9
# ╠═5408e59e-53ae-4ef9-b9ae-e86be3119575
# ╠═72faa038-d71c-46f9-8110-d68b7633f2e6
# ╠═9b4cf8f2-84ae-4aba-981a-7ef0e669a453
# ╠═d7fb857c-0871-4be2-8322-e36ab95e6f19
# ╠═05ff34ad-5cfc-48a4-bfc9-4d4c7cfda4d1
# ╠═171a33c2-e532-4167-82a6-27d830402104
# ╠═657e7156-57ef-4f64-a4ba-13af4bb2784e
# ╠═89114c37-3d83-4e4c-a46f-32cd04e85bf0
# ╠═d6b27bbe-824e-492b-9aea-dd5041f4917f
# ╠═e0bfd78c-25fa-4752-a4d7-660aa4b6cd66
# ╠═b6a72db1-ec4d-4cc1-8174-c6179316dec7
# ╠═f2829ecb-cbdf-4e0d-953c-0072f03546ac
# ╠═74279cfc-ab7a-4e9a-b048-ad17a378c20f
# ╠═5b5358c8-f35c-473f-9e92-56bdb4dbc6a7
# ╠═9031f497-26d2-4b3b-ba34-3438df409fd9
# ╠═ad11f743-1eee-4298-87bf-620a0f1cd1ec
# ╠═f793c68f-731e-426f-ad6e-ef5419f16370
# ╠═cc022b63-4a04-4142-976c-02756047076d
# ╠═a6b2a25c-6fb2-4dd0-a7f1-ea16ceb63a97
# ╠═69798f9b-f9a8-41f7-a2ff-8f333a6ed728
# ╠═741a21f9-ea94-4620-83a6-a43dd5fc300d
# ╠═e3f624cc-1c45-461a-ac19-bed551b1ba4e
# ╠═5e8c4ba1-9f8d-40df-b02e-cfeec4bf4d57
# ╠═e8cee36d-6886-4e27-806e-b4df1dadff7d
# ╠═59a166ff-402f-44cd-91d4-1f58ca676e33
# ╠═25c85ec0-af9e-4f74-bb9f-df7325604100
# ╠═4cea1683-7f74-42dd-a0d5-5d66eca41253
# ╠═58a03eb1-4847-4f48-a86e-9da1233023a3
# ╠═56b54a7c-fb57-4aa0-9e84-45ade0d06b5c
# ╠═52a1db73-e957-41fd-a6dc-f483b874e9f2
# ╠═fbd6502e-efa4-4e52-81db-3bb998314b14
# ╠═bf177ece-0e74-45f2-8bb0-2b5a2ac2ddbe
# ╠═8d2aafde-d745-4a0d-bcdb-d82a7877292f
# ╠═3d10cafd-f87c-40a7-ba15-f64e924cc686
# ╠═28779e98-1a0e-481f-bb10-c6973bd98976
# ╠═40173831-c379-4e62-90fc-f7f2e55ff8c6
# ╠═95a3d0c3-09c8-490d-ae5c-993d8ec009e8
# ╠═2dbf183b-159e-44bf-87e5-ae9d6517eeb3
# ╠═4444b683-1571-4df8-9841-0e28df365dec
# ╠═757fa776-c741-4eda-9e85-c1110231fcd4
# ╠═b8437ca7-b885-4dad-a72d-dc0b9c9a6be5
# ╠═90c789e2-999c-4e6c-8d49-0926ca080e71
# ╠═f1795d88-1e1f-4aa7-b245-b2f2cef7db17
# ╠═967bcd01-03a6-4cc5-ac07-b997b65d8a8f
# ╠═322240e8-a405-4ccf-8d94-3c52aa60eb26
# ╠═e698c85d-ab27-4dda-9349-cf679b22a0f5
# ╠═4415b732-d854-4eb9-9ae6-ed872b087c3e
# ╠═8dd89e2e-78e3-4d2a-b2f8-c9d29e949752
# ╠═c817df72-596c-41be-9d55-80dcbe73c485
# ╠═76fbfe22-1223-454d-8f6c-0505cdd6c8c4
# ╠═e33af351-2a7f-45fc-be7e-b36467b3e324
# ╠═b2fabd16-f406-4b70-b5ab-343ff87b6a05
# ╠═7b353ba1-cc33-405b-962e-5b7f82d1323b
# ╠═cb6f8399-addd-4095-98bf-915f1b38ee4d
# ╠═3cbff0f7-eb89-49e5-88e4-43ec57ab87b3
# ╠═46ce6ac5-c128-46df-8d38-90776f42861b
# ╠═ca061644-6cf0-4547-bb0d-c21642f76f3f
# ╠═d636baaf-b65f-4c5e-a9ed-775a7f1c9909
# ╠═134e34a4-e467-49e6-9c39-12123401d1b7
# ╠═d62aed96-41e6-4a01-a614-8f7c2d440b38
# ╠═1bf15795-87e8-4cf7-92b1-9db7b102450a
# ╠═d9be422c-3630-4c86-859a-3df1cd725405
# ╠═10a0495d-b861-4959-a5c2-1d74b6a0f5c6
# ╠═0e7cfe6c-a700-4ec1-b8dc-a67f246782b3
# ╠═3b6c3b40-d54e-41c2-afac-6bf8234ae648
# ╠═f97950ff-8c2b-4fc4-99e2-51c29bf7b489
# ╠═b9ee5aa2-ced4-4f71-b43a-bbfd6d5afef9
# ╠═bd94a0c6-0e84-4699-b449-36e45b2c2c23
# ╠═b79d6323-c6be-4b18-ae93-1a1ce48a1db3
# ╠═ca38417c-1d5c-45c6-816c-d2d6898af7db
# ╠═692af519-1ac2-483c-86c5-f1801d8426bd
# ╠═76d3d1c4-0e01-4360-bc3d-bc274ff1c39c
# ╠═3bfb00e8-f963-4c03-825f-71c4343074fa
# ╠═1c622f87-b763-4f8b-a253-fd38a1fd2830
# ╠═45a5a8f7-a1bb-434c-b882-3b24d3efa412
# ╠═96482819-8af4-4ce7-bfda-2e20ac05a66c
# ╠═6a70233a-4754-49e8-bd81-292d364b4031
# ╠═445fe109-2361-4bea-bd03-74f3f12dca65
# ╠═f7a4a1b1-5970-4e1c-a37a-5db901c1115c
# ╠═4605495f-6df1-44fb-8a33-6e84bebe83df
# ╠═38f8fbed-20a9-46ab-92e9-edf3cb693bfa
# ╠═f84bfc2a-c037-4627-b0c6-4f47d4f45a1b
# ╠═a0798ff2-c72f-4de5-b31a-3c8a3bacf117
# ╠═007e6893-0718-41c0-8ca0-03fd367f6cc7
# ╠═49252d74-9a99-4795-9527-3e7eccf20eed
# ╠═e9788e02-f17f-4bd8-b344-caa18863e92e
# ╠═9ba8c503-f59e-4d6a-b97c-9826279d744a
# ╠═3266d3d6-64a2-4936-91ac-059a1ccf76c4
# ╠═da65e2fe-54a8-4ca6-b599-411ee2a3e9d4
# ╠═7212801b-ee05-4699-b4fb-c02377710e10
# ╠═aa699be6-0fa9-468d-847c-b7b3b7f91c72
# ╠═6a423afd-67f4-4152-9aa8-7d69c4262f63
# ╠═c02b568d-e607-4729-9cef-1d189ae22bc1
# ╠═7aa12ba0-e6e7-4c36-9620-7ed443013d40
# ╠═a11ab145-e6d2-4452-8df8-6976f5f67453
# ╠═046608c2-5c79-4651-bcfe-1b8b965e63c2
# ╠═e2bf25a2-b5c2-48be-95b3-9ee2905b7b35
# ╠═5a1dfded-7f69-4fae-ad5c-9242686026d8
# ╠═77dbbadb-2588-4db7-8e84-5b375dbec467
# ╠═ce3be394-4d52-40c6-9df2-4cad1ca55fd7
# ╠═9f184944-6093-487f-a095-98a8572af620
# ╠═ff570f16-c8fb-4e0d-8851-a0fb93b23645
# ╠═b23610f5-240b-4385-955b-216b3221d1d8
# ╠═119b9407-c1cd-4eff-a0a6-993b9b05f50b
# ╠═f208fdf1-e949-468a-a987-d798baa939ca
# ╠═434ef740-5386-4a64-9ebe-8211e40b0c2b
# ╠═f9dbfebc-9c8e-451a-a46c-03b548b630de
# ╠═43603d8b-82dd-4670-a999-51802c6173df
# ╠═d5a93508-e559-4767-a078-54f9bb8fb3bf
# ╠═4ad4579a-4b64-460b-9ec4-a80f86ee391a
# ╠═e86d7b59-449e-4e8f-8ad2-1c708ea15076
# ╠═6412e0da-abe4-44fb-9e06-7cce1835d686
# ╠═2d9caf6d-154d-408e-8e4d-aef4f1c3ff78
# ╠═1a8aac68-4308-4700-ad2d-8e6989757416
# ╠═4b5089dd-74aa-4dbf-8e26-bdaecc0b5b36
# ╠═7361c452-d2c9-4793-855e-055911702d8a
# ╠═ae549f94-6e59-48f2-aedd-5703b7301202
# ╠═0a46f8ed-2a03-4b04-98be-5d85473343ea
# ╠═3c36cc71-db84-4f59-81fe-e4dd7761f72b
# ╠═e9d02b89-5a52-424f-a371-d0ca248bf316
# ╠═913d87d1-6d0a-4db3-9e46-f8a10d247dd5
# ╠═faddf0dd-f931-4379-a13a-021bd5b23cae
# ╠═c4a7d77d-68cd-40ef-a168-6dbf9dc410b8
# ╠═d3567161-530b-476a-b339-a0b2d415da7e
# ╠═ff2cea1e-4bb4-4cd2-a1f4-a68cfeb6efde
# ╠═26299f99-3851-4130-84b6-eaac2c59fbd8
# ╠═da893d5b-1fb4-41f1-a8c9-16c09196739a
# ╠═f474351b-6b02-4402-a722-9bd03ac6b1b9
# ╠═78ac5d01-2ded-41c5-9b6b-64dfc94a2141
# ╠═69b7aa3a-7d0f-4af9-8c71-f1479d60225b
# ╠═577ba306-497a-4cfc-8ea4-bdbc8afb2076
# ╠═f621f4d2-eb36-4733-84a1-9a77f3d88cb9
# ╠═31aa6d10-3ba3-4af0-9458-2c6a60aa01d1
# ╠═552238b4-26f1-4a3c-9706-65410a8ea4d1
# ╠═966e302e-5ab9-4806-ac72-244eb4be2aae
# ╠═eb32157e-74aa-4e61-9b4c-3b57bfa1cb0a
# ╠═734aebf9-6da9-4e6e-bc79-2905bcd75e7a
# ╠═a6b67788-7967-47fc-9d47-bc4e7c202f5a
# ╠═9a545506-f35d-48b8-aeff-ca24d8ab6442
# ╠═3b9fd69d-ca22-4f05-a563-3492ddd66285
# ╠═4167bc84-b881-4eb0-b140-aa35c0acd365
# ╠═c92f80c5-6f47-47c3-92d8-8d1cc9fd6c42
# ╠═50443a47-3600-4733-8b77-1d56717b00b4
# ╠═3060c2aa-078b-43e3-8451-b17642f470eb
# ╠═a65930bd-58ca-4326-a04d-0416614b7abb
# ╠═5b380673-9c93-4e6f-a158-efe469ab8c52
# ╠═e265dcee-67b7-4238-987a-19242ffe9f03
# ╠═d3dd2610-f200-4b22-bfb3-6d28df96de13
# ╠═e6e36c0d-f522-4959-a979-20573ffda8c2
# ╠═55243de1-4f6e-4184-a1c3-40b3e24b1a1b
# ╠═ee927604-3622-4dd0-8f02-b6838e036fc1
# ╠═2254d304-dc0a-4dae-87d7-f3cd18ebe152
# ╠═4e1e1f43-de2e-4584-ac04-31bcfd7036ad
# ╠═5fbffc42-e89e-47d8-b493-e297ad66ca0b
# ╠═8201e87f-19c7-4b61-bfb8-46ed75bfedcc
# ╠═3dd015d7-3c68-4b7b-9331-d6bd51c3702d
# ╠═e9c15162-c3d4-449f-88db-0280b5874a69
# ╠═70267a1f-34da-49c6-a09f-1c214d0dfecb
# ╠═c15eca41-5fda-4869-a298-4991fa4ad72d
# ╠═a0e6818d-4d4d-4644-84ab-a4c45c9dcb26
# ╠═36477658-b2f5-4ff4-b57a-72367ed61819
# ╠═4593f1eb-c752-4043-98d9-1a00954d738e
# ╠═0f4380ab-a77f-4413-823f-c16044efd721
# ╠═91d20be7-a7a4-4120-8756-e8e5656ec804
# ╠═0e7d33ba-f19d-45b3-8606-0fe9eb1845c3
# ╠═bcd1df63-e0b6-467b-9858-de5a15a61765
# ╠═05558267-7c22-4480-99f7-847b02a273cd
# ╠═ce1658a1-67fd-4743-90f3-9dc914d857be
# ╠═a1052731-3d45-42e5-b851-9dc0a8cf8177
# ╠═749d91db-cebd-43f2-a22b-248fd9b8587d
# ╠═730a04cc-b8db-406a-94ad-f951baf08adf
# ╠═8805231b-7e9f-4dec-a999-ac8e750a49b9
# ╠═a6b2b943-c350-4d99-9124-cfaaed1ddd85
