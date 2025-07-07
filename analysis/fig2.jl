include("lib.jl")

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


f = Figure(size = (820, 500))
ax_a = Axis(f[1, 1:2],
	       aspect = 1,
	       title = "Overall precision and recall",
	       xlabel = "Recall",
	       ylabel = "Precision")
xlims!(ax_a, [0, 1])
ylims!(ax_a, [0, 1])

ax_c1 = Axis(f[2, 1:5][1, 1],
             title = "Precision at equal coverage",
             xlabel = "Depth",
             xticks = (1:3, ["100", "500", "1000"]),
             ylabel = "Precision",
             yticks = [0, 0.25, 0.5, 0.75, 1.0])
ylims!(ax_c1, [0, 1.3])
ax_c2 = Axis(f[2, 1:5][1, 2],
             title = "Recall at equal coverage",
             xlabel = "Depth",
             xticks = (1:3, ["100", "500", "1000"]),
             ylabel = "Recall",
             yticks = [0, 0.25, 0.5, 0.75, 1.0])
ylims!(ax_c2, [0, 1.3])

a, p, r = auprc(binned_rna002.predicted,
                Int.(binned_rna002.modified),
                Set([1]))
lines!(ax_a, r, p,
       label="RNA002\nAUC=$(round(a; digits = 2))",
       color="#56b4e9")
p001 = precision(roc(binned_rna002.modified,
		     binned_rna002.predicted .>= -log10(0.01)))
r001 = recall(roc(binned_rna002.modified,
    		  binned_rna002.predicted .>= -log10(0.01)))
scatter!(ax_a, [r001], [p001], color = "#56b4e9")

a, p, r = auprc(binned_rna004.predicted,
                Int.(binned_rna004.modified),
                Set([1]))
lines!(ax_a, r, p,
       label="RNA004\nAUC=$(round(a; digits = 2))",
       color="#009e73")
p001 = precision(roc(binned_rna004.modified,
		     binned_rna004.predicted .>= -log10(0.01)))
r001 = recall(roc(binned_rna004.modified,
    		  binned_rna004.predicted .>= -log10(0.01)))
scatter!(ax_a, [r001], [p001], color = "#009e73")



axislegend(ax_a, "Chemistry", position = :rt, labelsize = 12)

# ==== Same coverage (c) =====

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
         bar_labels = map(v -> (@sprintf("%.2f", v)), df.values[df.metric .== 1]))

barplot!(ax_c2,
         df.depth[df.metric .== 2],
         df.values[df.metric .== 2],
         dodge = df.kit[df.metric .== 2],
         color = map(k -> colormap[k], df.kit[df.metric .== 2]),
         bar_labels = map(v -> (@sprintf("%.2f", v)), df.values[df.metric .== 2]))

Legend(f[2, 6],
       [PolyElement(color = "#56B4E9"), PolyElement(color = "#009E73")],
       ["RNA002", "RNA004"])



# ==== top right ====
ax11 = Axis(f[1, 3:6][1, 1:2],
	    yticks = (1:2, ["RNA002", "RNA004"]),
	    xlabel = "m6A sites in GLORI",
	    xticks = [0, 50000, 100000],
	    xtickformat = "{:n}",
	    height = 75)
xlims!(ax11, (0, 160000))

ax21 = Axis(f[1, 3:6][2, 1],
            title = "Overlap of true positives",
            height = 100)

ax22 = Axis(f[1, 3:6][2, 2],
            xlabel = "GLORI modification stoichiometry",
            ylabel = "Density",
            height = 100)

glori_not_covered_rna002 = size(
      antijoin(glori, binned_rna002,
	       on = [:chr, :strand, :bin]), 1)
glori_not_covered_rna004 = size(
      antijoin(glori, binned_rna004,
	       on = [:chr, :strand, :bin]), 1)
roc002 = roc(binned_rna002.modified .== 1, binned_rna002.predicted .>= -log10(0.01))
roc004 = roc(binned_rna004.modified .== 1, binned_rna004.predicted .>= -log10(0.01))

df = DataFrame(
     #kit = ["RNA002", "RNA002", "RNA002", "RNA004", "RNA004", "RNA004"],
     #class = ["Unmapped", "Missed", "Detected", "Unmapped", "Missed", "Detected"],
     kit = [1, 1, 1, 2, 2, 2],
     class = [0, 1, 2, 0, 1, 2],
     number = [glori_not_covered_rna002,
	       false_negative(roc002),
	       true_positive(roc002),
	       glori_not_covered_rna004,
	       false_negative(roc004),
	       true_positive(roc004)])

colormap = Dict(0 => :black,
                1 => :grey,
                2 => :lightgreen)

barplot!(ax11, df.kit, df.number,
         stack = df.class,
         color = map(c -> colormap[c], df.class),
         direction = :x,
         bar_labels = ["insufficient coverage",
                       "undetected",
                       (@sprintf "%d" true_positive(roc002)),
                       "insufficient coverage",
                       "undetected",
                       (@sprintf "%d" true_positive(roc004))],
         label_position = :end, #[:start, :start, :end, :start, :center, :end],
         # label_align = :lb,
         label_offset = [-160, -70, 5, -135, -80, 5],
         label_color = [:white, :white, :black, :white, :white, :black])
# bar_labels = map(n -> (@sprintf "%d" n), df.number),
# label_position = :center,
# label_color = :white)

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

axislegend(ax11, [PolyElement(color = :lightgreen)], ["detected"]; position=:rt)

rowsize!(f.layout, 1, Relative(0.75))

poly!(ax21, Circle(Point2f(-60, 0), sqrt(9851/pi)), color = (:white, 0), linestyle = :solid, strokewidth = 2, strokecolor = COL_RNA002)
poly!(ax21, Circle(Point2f(0, 0), sqrt(21627/pi)), color = (:white, 0), linestyle = :solid, strokewidth = 2, strokecolor = COL_RNA004)

# xlims!(ax21, (-1.6, 2))
# ylims!(ax21, (-1.6, 2))

text!(ax21, -95, 0, text = "$rna002_only", align = (:center, :center), fontsize = 10)
text!(ax21, -50, 0, text = "$common", align = (:center, :center), fontsize = 10)
text!(ax21, 35, 0, text = "$rna004_only", align = (:center, :center), fontsize = 10)

# poly!(ax21, Circle(Point2f(-0.9, 0), 1.3*2.777/(2*pi)),
#       color = (:white, 0),
#       linestyle = :solid,
#       strokewidth = 2,
#       strokecolor = COL_RNA002)
# poly!(ax21, Circle(Point2f(0.4, 0), 9.043/(2*pi)),
#       color = (:white, 0),
#       linestyle = :solid,
#       strokewidth = 2,
#       strokecolor = COL_RNA004)
# 
# xlims!(ax21, (-1.6, 2))
# ylims!(ax21, (-1.6, 2))
#
# text!(ax21, -1.25, 0, text = "$rna002_only", align = (:center, :center), fontsize = 10)
# text!(ax21, -0.74, 0, text = "$common", align = (:center, :center), fontsize = 10)
# text!(ax21, 0.3, 0, text = "$rna004_only", align = (:center, :center), fontsize = 10)


density!(ax22,
	 glori.ratio,
         color = (COL_RNA004, 0),
         strokecolor = :black,
         strokewidth = 1,
         # strokearound = true,
	 linestyle = :solid,
         label = "Overall")
density!(ax22,
         leftjoin(rna002[rna002.modified .== 1, :],
                  glori,
                  on = [:chr, :strand, :genomicPos => :start],
                  makeunique = true).ratio,
         color = (COL_RNA002, 0),
         strokecolor = COL_RNA002,
         strokewidth = 1,
         # strokearound = true,
	 linestyle = :solid,
         # boundary = (0, 1),
         label = "RNA002")
density!(ax22,
         leftjoin(rna002[rna002.modified .== 1 .&& rna002.predicted .>= -log10(0.01), :],
                  glori,
                  on = [:chr, :strand, :genomicPos => :start],
                  makeunique = true).ratio,
         color = (COL_RNA002, 0),
         strokecolor = COL_RNA002,
         strokewidth = 2,
         # strokearound = true,
	 linestyle = :dash,
         # boundary = (0, 1),
         label = "RNA002")
density!(ax22,
         leftjoin(rna004[rna004.modified .== 1, :],
                  glori,
                  on = [:chr, :strand, :genomicPos => :start],
                  makeunique = true).ratio,
         color = (COL_RNA004, 0),
         strokecolor = COL_RNA004,
         strokewidth = 1,
         # strokearound = true,
	 linestyle = :solid,
         label = "RNA004")
density!(ax22,
         leftjoin(rna004[rna004.modified .== 1 .&& rna004.predicted .>= -log10(0.01), :],
                  glori,
                  on = [:chr, :strand, :genomicPos => :start],
                  makeunique = true).ratio,
         color = (COL_RNA004, 0),
         strokecolor = COL_RNA004,
         strokewidth = 2,
         # strokearound = true,
	 linestyle = :dash,
         label = "RNA004")
# density!(ax22,
# glori.ratio,
# color = (:grey, 0),
# linestyle = :dash,
# strokecolor = :grey,
# strokewidth = 2,
# strokearound = true)

hidedecorations!(ax21)
hidespines!(ax21)

save("fig_rna002_rna004_comparison.png", f)

