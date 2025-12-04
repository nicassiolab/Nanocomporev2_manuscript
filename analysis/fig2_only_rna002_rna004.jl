include("lib.jl")
import MultipleTesting

binned_glori = get_binned_glori("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed")

glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
# glori[!, "pos"] = glori.start .+ 2
glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.start)
glori[!, "modified"] .= 1

ref002 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv"))
ref002 = combine(groupby(ref002, ["chr", "strand", "genomicPos"]),
                 :modified => maximum => :modified)

ref004 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv"))
ref004 = combine(groupby(ref004, ["chr", "strand", "genomicPos"]),
                 :modified => maximum => :modified)

rna002 = annotate_results(
  read_results(
    "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv",
    "GMM_chi2_qvalue",
    shift=2,
    genomic_collapse=true),
  ref002)

binned_rna002 = annotate_binned_results(bin(rna002, BIN_SIZE), 
                                        binned_glori)

rna004 = annotate_results(
  read_results(
    "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_test_v2_release/out_nanocompore_results.tsv",
    "GMM_chi2_qvalue",
    shift=4,
    genomic_collapse=true),
  ref004)

binned_rna004 = annotate_binned_results(bin(rna004, BIN_SIZE),
                                        binned_glori)

rna002_100 = read_results(
  "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_depth_100_v2.0.0/out_nanocompore_results.tsv",
  "GMM_chi2_qvalue",
  shift=2,
  genomic_collapse=true)
binned_rna002_100 = annotate_binned_results(bin(rna002_100, BIN_SIZE),
                                            binned_glori)

rna002_500 = read_results(
  "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_depth_500_v2.0.0/out_nanocompore_results.tsv",
  "GMM_chi2_qvalue",
  shift=2,
  genomic_collapse=true)
binned_rna002_500 = annotate_binned_results(bin(rna002_500, BIN_SIZE),
                                            binned_glori)

rna002_1000 = read_results(
  "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_depth_1000_v2.0.0/out_nanocompore_results.tsv",
  "GMM_chi2_qvalue",
  shift=2,
  genomic_collapse=true)
binned_rna002_1000 = annotate_binned_results(bin(rna002_1000, BIN_SIZE),
                                             binned_glori)

rna004_100 = read_results(
  "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_100_v2.0.0/out_nanocompore_results.tsv",
  "GMM_chi2_qvalue",
  shift=4,
  genomic_collapse=true)
binned_rna004_100 = annotate_binned_results(bin(rna004_100, BIN_SIZE),
                                            binned_glori)

rna004_500 = read_results(
  "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_500_v2.0.0/out_nanocompore_results.tsv",
  "GMM_chi2_qvalue",
  shift=4,
  genomic_collapse=true)
binned_rna004_500 = annotate_binned_results(bin(rna004_500, BIN_SIZE),
                                            binned_glori)

rna004_1000 = read_results(
  "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_1000_v2.0.0/out_nanocompore_results.tsv",
  "GMM_chi2_qvalue",
  shift=4,
  genomic_collapse=true)
binned_rna004_1000 = annotate_binned_results(bin(rna004_1000, BIN_SIZE),
                                             binned_glori)


# Take only the common sites
binned_rna002_100 = innerjoin(binned_rna002_100, binned_rna004_100[:, [:chr, :strand, :bin]],
			      on = [:chr, :strand, :bin])
binned_rna002_500 = innerjoin(binned_rna002_500, binned_rna004_500[:, [:chr, :strand, :bin]],
			      on = [:chr, :strand, :bin])
binned_rna002_1000 = innerjoin(binned_rna002_1000, binned_rna004_1000[:, [:chr, :strand, :bin]],
		 	       on = [:chr, :strand, :bin])
binned_rna004_100 = innerjoin(binned_rna004_100, binned_rna002_100[:, [:chr, :strand, :bin]],
			      on = [:chr, :strand, :bin])
binned_rna004_500 = innerjoin(binned_rna004_500, binned_rna002_500[:, [:chr, :strand, :bin]],
			      on = [:chr, :strand, :bin])
binned_rna004_1000 = innerjoin(binned_rna004_1000, binned_rna002_1000[:, [:chr, :strand, :bin]],
		 	       on = [:chr, :strand, :bin])


roc002 = roc(binned_rna002.modified .== 1,
	  	   binned_rna002.predicted .>= -log10(0.01))
roc004 = roc(binned_rna004.modified .== 1,
		   binned_rna004.predicted .>= -log10(0.01))


f = Figure(size = (730, 800))

g1 = f[1, 1] = GridLayout()
g2 = f[2, 1] = GridLayout()
g3 = f[3, 1] = GridLayout()


# ==== A: GLORI barplot ====

ax_glori = Axis(g1[1, 1],
	        yticks = (1:2, ["RNA004", "RNA002"]),
	        xlabel = "m⁶A sites in GLORI",
	        xticks = [0, 50000, 100000],
	        xtickformat = "{:n}",
	        height = 75)
xlims!(ax_glori, (0, 177000))

glori_not_covered_rna002 = size(
      antijoin(glori, binned_rna002,
	       on = [:chr, :strand, :bin]), 1)
glori_not_covered_rna004 = size(
      antijoin(glori, binned_rna004,
	       on = [:chr, :strand, :bin]), 1)

df = DataFrame(
     #kit = ["RNA004", "RNA004", "RNA004", "RNA002", "RNA002", "RNA002"],
     #class = ["Unmapped", "Missed", "Detected", "Unmapped", "Missed", "Detected"],
     kit = [1, 1, 1, 2, 2, 2],
     class = [0, 1, 2, 0, 1, 2],
     number = [glori_not_covered_rna004,
	       false_negative(roc004),
	       true_positive(roc004),
	       glori_not_covered_rna002,
	       false_negative(roc002),
	       true_positive(roc002)])

colormap = Dict(0 => :black,
		1 => :grey,
		2 => :lightgreen)

barplot!(ax_glori, df.kit, df.number,
	 stack = df.class,
	 color = map(c -> colormap[c], df.class),
	 direction = :x,
	 bar_labels = [(@sprintf "%d" glori_not_covered_rna004),
		       (@sprintf "%d" false_negative(roc004)),
		       (@sprintf "%d" true_positive(roc004)),
		       (@sprintf "%d" glori_not_covered_rna002),
		       (@sprintf "%d" false_negative(roc002)),
		       (@sprintf "%d" true_positive(roc002))],
	 label_position = :center,
	 label_size = [12, 12, 12, 12, 12, 12],
	 label_color = [:white, :white, :black, :white, :white, :black])

Legend(g1[2, 1],
       [PolyElement(color = :black),
	PolyElement(color = :grey),
	PolyElement(color = :lightgreen)],
       ["no coverage", "false negatives", "true positives"],
       orientation = :horizontal,
       framevisible = false)


# ==== B: Precision/recall curve ====

ax_b = Axis(g2[1:2, 1:2],
	    title = "Precision/recall for v2 by chemistry",
	    aspect = 1,
	    xlabel = "Recall",
	    ylabel = "Precision")

a, p, r = auprc(binned_rna002.predicted,
		Int.(binned_rna002.modified),
		Set([1]))
lines!(ax_b, r[r .> 1e-200], p[r .> 1e-200],
       label = "RNA002\nAUC=$(round(a; digits = 2))",
       color = COL_RNA002)
p001 = precision(roc(binned_rna002.modified,
		     binned_rna002.predicted .>= -log10(0.01)))
r001 = recall(roc(binned_rna002.modified,
		  binned_rna002.predicted .>= -log10(0.01)))
scatter!(ax_b, [r001], [p001], color = COL_RNA002)
println("RNA002 precision = $p001, recall = $r001")

a, p, r = auprc(binned_rna004.predicted,
		Int.(binned_rna004.modified),
		Set([1]))
lines!(ax_b, r[r .> 1e-200], p[r .> 1e-200],
       label = "RNA004\nAUC=$(round(a; digits = 2))",
       color = COL_RNA004)
p001 = precision(roc(binned_rna004.modified,
		     binned_rna004.predicted .>= -log10(0.01)))
r001 = recall(roc(binned_rna004.modified,
		  binned_rna004.predicted .>= -log10(0.01)))
scatter!(ax_b, [r001], [p001], color = COL_RNA004)
println("RNA004 precision = $p001, recall = $r001")

axislegend(ax_b, "Chemistry", position = :rt, labelsize = 12)


# ==== C: Same coverage =====
ax_c1 = Axis(g2[1, 3:4],
	     title = "Precision at equal coverage",
	     xlabel = "Depth",
	     xticks = (1:3, ["100", "500", "1000"]),
	     ylabel = "Precision",
	     yticks = [0, 0.25, 0.5, 0.75, 1.0],
	     height=100)
ylims!(ax_c1, [0, 1.15])
ax_c2 = Axis(g2[2, 3:4],
	     title = "Recall at equal coverage",
	     xlabel = "Depth",
	     xticks = (1:3, ["100", "500", "1000"]),
	     ylabel = "Recall",
	     yticks = [0, 0.25, 0.5, 0.75, 1.0],
	     height=100)
ylims!(ax_c2, [0, 1.15])

roc_rna002_100 = roc(binned_rna002_100.modified,
                     binned_rna002_100.predicted .>= -log10(0.01))
roc_rna002_500 = roc(binned_rna002_500.modified,
                     binned_rna002_500.predicted .>= -log10(0.01))
roc_rna002_1000 = roc(binned_rna002_1000.modified,
                      binned_rna002_1000.predicted .>= -log10(0.01))
roc_rna004_100 = roc(binned_rna004_100.modified,
                     binned_rna004_100.predicted .>= -log10(0.01))
roc_rna004_500 = roc(binned_rna004_500.modified,
                     binned_rna004_500.predicted .>= -log10(0.01))
roc_rna004_1000 = roc(binned_rna004_1000.modified,
                      binned_rna004_1000.predicted .>= -log10(0.01))

df = DataFrame(
     # kit = ["RNA002", "RNA002", "RNA002", "RNA004", "RNA004", "RNA004"],
     kit = [1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2],
     # depth = [100, 500, 1000, 100, 500, 1000],
     depth = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
     metric = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
     values = [precision(roc_rna002_100),
	       precision(roc_rna002_500),
	       precision(roc_rna002_1000),
	       precision(roc_rna004_100),
	       precision(roc_rna004_500),
	       precision(roc_rna004_1000),
	       recall(roc_rna002_100),
	       recall(roc_rna002_500),
	       recall(roc_rna002_1000),
	       recall(roc_rna004_100),
	       recall(roc_rna004_500),
	       recall(roc_rna004_1000)])

colormap = Dict(1 => "#56B4E9",
		2 => "#009E73")

barplot!(ax_c1,
	 df.depth[df.metric .== 1],
	 df.values[df.metric .== 1],
	 dodge = df.kit[df.metric .== 1],
	 color = map(k -> colormap[k], df.kit[df.metric .== 1]),
	 bar_labels = map(v -> (@sprintf("%.2f", v)), df.values[df.metric .== 1]),
	 label_size = 11)

barplot!(ax_c2,
	 df.depth[df.metric .== 2],
	 df.values[df.metric .== 2],
	 dodge = df.kit[df.metric .== 2],
	 color = map(k -> colormap[k], df.kit[df.metric .== 2]),
	 bar_labels = map(v -> (@sprintf("%.2f", v)), df.values[df.metric .== 2]),
	 label_size = 11)


# ==== D: density ====

ax_stoich = Axis(g3[1, 1],
		 title = "Modification stoichiometry density",
		 xlabel = "Modification stoichiometry",
		 ylabel = "Density",)

density!(ax_stoich,
	 glori.ratio,
	 color = (COL_RNA004, 0),
	 strokecolor = :black,
	 strokewidth = 1,
	 linestyle = :solid,
	 label = "Overall")
density!(ax_stoich,
	 leftjoin(rna002[rna002.modified .== 1, :],
		  glori,
		  on = [:chr, :strand, :genomicPos => :start],
		  makeunique = true).ratio,
	 color = (COL_RNA002, 0),
	 strokecolor = COL_RNA002,
	 strokewidth = 1,
	 linestyle = :solid,
	 label = "RNA002")
density!(ax_stoich,
	 leftjoin(rna002[rna002.modified .== 1 .&& rna002.predicted .>= -log10(0.01), :],
		  glori,
		  on = [:chr, :strand, :genomicPos => :start],
		  makeunique = true).ratio,
	 color = (COL_RNA002, 0),
	 strokecolor = COL_RNA002,
	 strokewidth = 2,
	 linestyle = :dash,
	 label = "RNA002")
density!(ax_stoich,
	 leftjoin(rna004[rna004.modified .== 1, :],
		  glori,
		  on = [:chr, :strand, :genomicPos => :start],
		  makeunique = true).ratio,
	 color = (COL_RNA004, 0),
	 strokecolor = COL_RNA004,
	 strokewidth = 1,
	 linestyle = :solid,
	 label = "RNA004")
density!(ax_stoich,
	 leftjoin(rna004[rna004.modified .== 1 .&& rna004.predicted .>= -log10(0.01), :],
		  glori,
		  on = [:chr, :strand, :genomicPos => :start],
		  makeunique = true).ratio,
	 color = (COL_RNA004, 0),
	 strokecolor = COL_RNA004,
	 strokewidth = 2,
	 linestyle = :dash,
	 label = "RNA004")


# ==== E: venn ====

ax_venn = Axis(g3[1, 2],
	       title = "Overlap of true positives",
	       height = 76)

rna002_tps = binned_rna002[binned_rna002.modified .== 1 .&& binned_rna002.predicted .>= -log10(0.01), :]
rna004_tps = binned_rna004[binned_rna004.modified .== 1 .&& binned_rna004.predicted .>= -log10(0.01), :]

common = size(innerjoin(rna002_tps, rna004_tps,
			on = [:chr, :strand, :bin],
			makeunique = true), 1)
rna002_only = size(antijoin(rna002_tps, rna004_tps,
				  on = [:chr, :strand, :bin],
				  makeunique = true), 1)
rna004_only = size(antijoin(rna004_tps, rna002_tps,
				  on = [:chr, :strand, :bin],
				  makeunique = true), 1)



poly!(ax_venn, Circle(Point2f(-60, 0), sqrt(9851/pi)), color = (:white, 0), linestyle = :solid, strokewidth = 2, strokecolor = COL_RNA002)
poly!(ax_venn, Circle(Point2f(0, 0), sqrt(21627/pi)), color = (:white, 0), linestyle = :solid, strokewidth = 2, strokecolor = COL_RNA004)

text!(ax_venn, -95, 0, text = "$rna002_only", align = (:center, :center), fontsize = 12)
text!(ax_venn, -50, 0, text = "$common", align = (:center, :center), fontsize = 12)
text!(ax_venn, 35, 0, text = "$rna004_only", align = (:center, :center), fontsize = 12)


Legend(g3[2, 1],
       [LineElement(color = :black, linestyle = :solid),
        LineElement(color = :black, linestyle = :dash)],
       ["GLORI", "Nanocompore"],
       valign = :top,
       halign = :center,
       orientation = :horizontal,
       framevisible=false)

Legend(g3[2, 2],
       [PolyElement(color = COL_RNA002),
	PolyElement(color = COL_RNA004)],
       ["RNA002", "RNA004"],
       valign = :top,
       halign = :center,
       orientation = :horizontal,
       framevisible=false)

hidedecorations!(ax_venn)
hidespines!(ax_venn)


for (label, layout) in zip(["A", "B", "C", "D", "E"],
			   [g1[1, 1], g2[1, 1], g2[1, 3], g3[1, 1], g3[1, 2]])
  Label(layout[1, 1, TopLeft()],
        label,
        fontsize = 18,
        font = :bold,
        padding = (50, 30, 5, 0),
        halign = :right)
end

# rowgap!(f.layout, -10)
rowgap!(g1, 0)
colgap!(g2, -15)
rowgap!(g2, 0)
colgap!(g3, -15)
rowgap!(g3, 0)

save("fig_rna002_rna004_comparison.png", f)

