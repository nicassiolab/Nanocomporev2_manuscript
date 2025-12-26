include("lib.jl")
using AlgebraOfGraphics
import MultipleTesting

v1_results = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1_seed_72/outnanocompore_results.tsv", delim = '\t'))
v1_results[!, :GMM_logit_pvalue] = [p === NaN ? 1 : p for p in v1_results.GMM_logit_pvalue]
v1_results[!, :GMM_LOR] = [v == "NC" ? 0 : parse(Float64, v) for v in v1_results.Logit_LOR]
v1_results[!, :GMM_logit_pvalue] = clamp.(v1_results.GMM_logit_pvalue, 1e-300, 1)

rna002_tx = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv", delim = '\t'))
rna002_tx[!, :GMM_chi2_pvalue] = clamp.(coalesce.(rna002_tx.GMM_chi2_pvalue, 1), 1e-300, 1)
rna002_tx[!, :GMM_chi2_qvalue] = clamp.(coalesce.(rna002_tx.GMM_chi2_qvalue, 1), 1e-300, 1)
rna002_tx[!, :GMM_LOR] = coalesce.(rna002_tx.GMM_LOR, 0)


glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
glori[:, :pos] = glori.start
glori = rename!(glori, :chr => :chrom)
glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
glori[!, "modified"] .= 1

binned_glori = unique(glori[!, ["chrom", "strand", "bin", "modified"]])
binned_glori = combine(groupby(glori, [:chrom, :strand, :bin, :modified]),
					   :ratio => maximum => :ratio)
binned_glori = rename!(binned_glori, :chrom => :chr)

v1_v2_common = innerjoin(v1_results, rna002_tx,
					     on = [:ref_id, :pos],
						 renamecols = "_v1" => "_v2")


common = copy(v1_v2_common)
common[:, :bin] = map(pos -> div(pos, BIN_SIZE), common.genomicPos_v2)
# In v1, the pvalue is already corrected with Benjamini-Hohchberg in place.
# In v2 we keep the pvalue unmodified and save the corrected one in the qvalue column
# However, due to changes in the coverage filtering of the positions, v2
# includes more positions than v1, so the multiple test correction is more
# punishing. To make them more fairly comparable, we apply the multiple
# test correction for v2 only on the common sites.
common[:, :GMM_chi2_pvalue_v2] = MultipleTesting.adjust(common.GMM_chi2_pvalue_v2, MultipleTesting.BenjaminiHochberg())
# Instead of setting the pvalues of positions with low LOR to 1, which creates
# a big jump in the precision-recall curve, we divide the log-transformed p-value
# to make them insignificant while keeping their order.
common[:, :predicted_v1] = ifelse.(abs.(common.GMM_LOR_v1) .< 0.5, -log10.(common.GMM_logit_pvalue_v1) ./ 300, -log10.(common.GMM_logit_pvalue_v1))
common[:, :predicted_v2] = ifelse.(abs.(common.GMM_LOR_v2) .< 0.5, -log10.(common.GMM_chi2_pvalue_v2) ./ 300, -log10.(common.GMM_chi2_pvalue_v2))
col_selector = (pred, col) -> col[argmax(pred)]
binned_v1_v2_common = combine(groupby(common, [:chr_v2, :strand_v2, :bin]),
	:predicted_v1 => maximum => :predicted_v1,
	:predicted_v2 => maximum => :predicted_v2,
	[:predicted_v1, :GMM_LOR_v1] => col_selector => :GMM_LOR_v1,
	[:predicted_v1, :GMM_logit_pvalue_v1] => col_selector => :GMM_pvalue_v1,
	[:predicted_v2, :GMM_LOR_v2] => col_selector => :GMM_LOR_v2,
	[:predicted_v2, :GMM_chi2_pvalue_v2] => col_selector => :GMM_pvalue_v2)
binned_v1_v2_common = rename!(binned_v1_v2_common, :chr_v2 => :chr, :strand_v2 => :strand)

binned_v1_v2_common = leftjoin(binned_v1_v2_common, binned_glori, on = [:chr, :strand, :bin])
binned_v1_v2_common[:, :modified] = coalesce.(binned_v1_v2_common.modified, 0)
binned_v1_v2_common[:, :ratio] = coalesce.(binned_v1_v2_common.ratio, 0)
binned_v1_v2_common = dropmissing(binned_v1_v2_common, disallowmissing = true)

pvals1 = common.predicted_v1
pvals2 = common.predicted_v2

println("Correlation: ", cor(pvals1, pvals2))


fig = Figure(size = (1200, 600))

df = DataFrame(p1=Float64.(pvals1), p2=Float64.(pvals2))
# df = df[df.p1 .> 1 .|| df.p2 .> 1, :]
plt1 = data(df) *
		mapping(:p1, :p2) *
		visual(Scatter; markersize=5, rasterize=2)
draw!(fig[1, 1], plt1;
	 axis = (;
			 aspect = 1,
			 xlabel = "v1: -log₁₀(𝑃-value)",
			 ylabel = "v2: -log₁₀(𝑃-value)",
			 xlabelsize=18,
			 ylabelsize=18))

	

rocs1 = roc(binned_v1_v2_common.modified, binned_v1_v2_common.predicted_v1, nrow(binned_v1_v2_common))
precisions1 = vcat([0], map(precision, rocs1), [1])
recalls1 = vcat([1], map(recall, rocs1), [0])
rocs2 = roc(binned_v1_v2_common.modified, binned_v1_v2_common.predicted_v2, nrow(binned_v1_v2_common))
precisions2 = vcat([0], map(precision, rocs2), [1])
recalls2 = vcat([1], map(recall, rocs2), [0])

df1 = DataFrame(precision = precisions1, recall = recalls1)
df1[:, :version] .= "v1"
df2 = DataFrame(precision = precisions2, recall = recalls2)
df2[:, :version] .= "v2"

df = vcat(df1, df2)

p001_v1 = precision(roc(binned_v1_v2_common.modified, binned_v1_v2_common.predicted_v1 .>= -log10(0.01)))
r001_v1 = recall(roc(binned_v1_v2_common.modified, binned_v1_v2_common.predicted_v1 .>= -log10(0.01)))

p001_v2 = precision(roc(binned_v1_v2_common.modified, binned_v1_v2_common.predicted_v2 .>= -log10(0.01)))
r001_v2 = recall(roc(binned_v1_v2_common.modified, binned_v1_v2_common.predicted_v2 .>= -log10(0.01)))



plt2 = (
		data(df) *
			mapping(:recall => "Recall", :precision => "Precision", color=:version, linestyle=:version) *
			visual(Lines; rasterize=4)
		+
		data(DataFrame(precision=[p001_v1, p001_v2],
					   recall=[r001_v1, r001_v2],
					   version=["v1", "v2"])) *
			mapping(:recall => "Recall", :precision => "Precision", color=:version) *
			visual(Scatter)
	   )
grid = draw!(fig[1, 2], plt2, scales(Color = (; palette = [:black, COL_RNA002])); axis = (; aspect = 1, xlabelsize=18, ylabelsize=18))
legend!(fig[1, 2], grid;
        tellheight=false,
		tellwidth=false,
		halign=:right,
		valign=:top,
		margin=(10, 30, 10, 10),
		framevisible=true,
		labelsize=18,
		titlesize=18)


for (label, layout) in zip(["A", "B"],
						   [fig[1, 1], fig[1, 2]])
    Label(layout[1, 1, TopLeft()],
		  label,
		  fontsize=26,
		  font=:bold,
		  padding=(0, 5, 5, 0),
		  halign=:right)
end

save("supplemental_fig2_v1_vs_v2.png", fig)
save("supplemental_fig2_v1_vs_v2.pdf", fig)
