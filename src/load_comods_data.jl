include("lib.jl")

using FASTX
using AlgebraOfGraphics
using SQLite
import MultipleTesting
import LinearAlgebra
import HypothesisTests

mod_ref = CSV.read("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/combined_mod_ref.tsv", DataFrame, delim = '\t')

gtf = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/gencode.v41.annotation.gtf", delim = '\t', comment = "#", header = [:chr, :source, :feature, :start, :end, :score, :strand, :frame, :attribute]))
gtf[!, :transcript_id] = map(a -> replace(split(split(a, "; ")[2], " ")[2], "\"" => ""), gtf.attribute)

glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k.bed"))
# glori[!, "pos"] = glori.start .+ 2
glori[:, :pos] = glori.start
glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
glori[!, "modified"] .= 1

tx_ref = Dict()
FASTAReader(open("/work/mzdravkov/gencode.v41.transcripts.fa")) do reader
	for record in reader
		tx_ref[identifier(record)] = sequence(record)
	end
end

ivt = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_nanocompore_results.tsv",
				   "GMM_chi2_qvalue",
				   shift=4,
				   genomic_collapse=false,
				   LOR_threshold = 0.8)
ivt[!, :pos] .+= 4
ivt[!, :mod_ratio] = sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))

peaks_ivt = peaks(ivt, 4)
sig_peaks_ivt = peaks_ivt[peaks_ivt.predicted .!== missing .&& peaks_ivt.predicted .>= -log10(0.01), :]

an_sig_peaks_ivt = annotate_peaks(sig_peaks_ivt, mod_ref, writer_motifs, writer_mods)

comods = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/comod_50perc_contingency.tsv")) # this run always has fisher's test
comods[!, :pos1] .+= 4
comods[!, :pos2] .+= 4
comods[!, :chi2_pvalue] = comods.pvalue
comods[!, :pvalue] = ifelse.(comods.fisher_pvalue .!== missing,
							 comods.fisher_pvalue,
							 comods.pvalue)
comods[!, :pvalue] = clamp.(comods.pvalue, 1e-300, 1)
comods[!, :relation] = map(assoc_type, comods.residual00, comods.residual01, comods.residual10, comods.residual11)
comods[!, :phi] .*= ifelse.(comods.relation .== :comod, 1, -1)
comods[:, :qvalue] = MultipleTesting.adjust(comods.pvalue, MultipleTesting.BenjaminiHochberg())
comods = comods[abs.(comods.pos1 .- comods.pos2) .> 8, :]

# add genomic positions
comods = leftjoin(comods, ivt[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
				  on=[:reference => :ref_id, :pos1 => :pos])
rename!(comods, :genomicPos => :genomicPos1)
comods = leftjoin(comods, ivt[:, [:ref_id, :pos, :genomicPos]],
				  on=[:reference => :ref_id, :pos2 => :pos])
rename!(comods, :genomicPos => :genomicPos2)


sig_comods = comods[comods.qvalue .< 0.05, :]
an_sig_comods = annotate_comods(an_sig_peaks_ivt, sig_comods)

