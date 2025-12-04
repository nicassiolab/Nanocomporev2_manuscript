include("lib.jl")


println("Loading GLORI")

# binned_glori = get_binned_glori("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed")

glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
glori[:, :bin] = map(pos -> div(pos, BIN_SIZE), glori.start)
glori[:, :modified] .= 1
binned_glori = unique(glori[:, [:chr, :strand, :bin, :modified]])


println("Loading reference")
ref = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv"))
ref = combine(groupby(ref, ["chr", "strand", "genomicPos"]),
              :modified => maximum => :modified)

println("Loading eventalign results")
eventalign = annotate_results(
  read_results(
    "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_test_v2_release/out_nanocompore_results.tsv",
    "GMM_chi2_qvalue",
    shift=4,
    genomic_collapse=true),
  ref)

binned_eventalign = annotate_binned_results(bin(eventalign, BIN_SIZE),
                                            binned_glori)

println("Loading uncalled4 results")
uncalled4 = annotate_results(
  read_results(
    "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_uncalled4_v2.0.0/out_nanocompore_results.tsv",
    "GMM_chi2_qvalue",
    shift=4,
    genomic_collapse=true),
  ref)

binned_uncalled4 = annotate_binned_results(bin(uncalled4, BIN_SIZE),
                                           binned_glori)

println("Loading remora results")
remora = annotate_results(
  read_results(
     "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_remora_v2.0.0/out_nanocompore_results.tsv",
     "GMM_chi2_qvalue",
     shift=0,
     genomic_collapse=true),
  ref)

binned_remora = annotate_binned_results(bin(remora, BIN_SIZE),
                                        binned_glori)

println("Joining pairs of results")
eventalign_uncalled4 = innerjoin(binned_eventalign,
                                 binned_uncalled4,
                                 on = "bin",
                                 makeunique = true)

eventalign_remora = innerjoin(binned_eventalign,
                              binned_remora,
                              on = "bin",
                              makeunique = true)

uncalled4_remora = innerjoin(binned_uncalled4,
                             binned_remora,
                             on = "bin",
                             makeunique = true)

println("Plotting")

resq_fig = Figure(size = (550, 1500))

prc_ax = Axis(resq_fig[1, 1],
              aspect = 1,
              title = "Precision and recall of different resquigglers",
              xlabel = "Recall",
              ylabel = "Precision",
              height = 460)
right_grid = resq_fig[2, 1] = GridLayout()
xlims!(prc_ax, [0, 1])
ylims!(prc_ax, [0, 1])
a, p, r = auprc(disallowmissing(binned_eventalign[!, "predicted"]),
                binned_eventalign[!, "modified"],
                Set([1]))
lines!(prc_ax, r, p, label="Eventalign (AUC=$(round(a; digits = 2)))", color="#E69F00")
p001 = precision(roc(binned_eventalign.modified,
		     binned_eventalign.predicted .>= -log10(0.01)))
r001 = recall(roc(binned_eventalign.modified,
    		  binned_eventalign.predicted .>= -log10(0.01)))
scatter!(prc_ax, [r001], [p001], color = "#E69F00")
println("Eventalign predicted modified $(sum(binned_eventalign[!, :predicted] .>= -log10(0.01)))")

a, p, r = auprc(disallowmissing(binned_uncalled4[!, "predicted"]),
                binned_uncalled4[!, "modified"],
                Set([1]))
lines!(prc_ax, r, p, label="Uncalled4 (AUC=$(round(a; digits = 2)))", color="#D55E00")
p001 = precision(roc(binned_uncalled4.modified,
		     binned_uncalled4.predicted .>= -log10(0.01)))
r001 = recall(roc(binned_uncalled4.modified,
    		  binned_uncalled4.predicted .>= -log10(0.01)))
scatter!(prc_ax, [r001], [p001], color = "#D55E00")
println("Uncalled4 predicted modified $(sum(binned_uncalled4[!, :predicted] .>= -log10(0.01)))")

a, p, r = auprc(disallowmissing(binned_remora[!, "predicted"]),
                binned_remora[!, "modified"],
                Set([1]))
lines!(prc_ax, r, p, label="Remora (AUC=$(round(a; digits = 2)))", color="#CC79A7")
p001 = precision(roc(binned_remora.modified,
		     binned_remora.predicted .>= -log10(0.01)))
r001 = recall(roc(binned_remora.modified,
    		  binned_remora.predicted .>= -log10(0.01)))
scatter!(prc_ax, [r001], [p001], color = "#CC79A7")
println("Remora predicted modified $(sum(binned_remora[!, :predicted] .>= -log10(0.01)))")

axislegend("Resquiggler", position = :rt, labelsize = 14, framevisible = false, titlesize = 0)

sf_ax_ev = Axis(right_grid[1, 1],
                title = "Eventalign",
		height=200,
                aspect = 1)
xlims!(sf_ax_ev, (0, 5))
sharkfin(sf_ax_ev,
         discard_low_sig_positions(eventalign, 0.005), #[1:10000, :],
         "predicted_raw",
         "GMM_LOR_raw",
         markersize = 4)

sf_ax_u4 = Axis(right_grid[2, 1],
                title = "Uncalled4",
		height=200,
                aspect = 1)
xlims!(sf_ax_u4, (0, 5))
sharkfin(sf_ax_u4,
         discard_low_sig_positions(uncalled4, 0.005), #[1:10000, :], 
         "predicted_raw",
         "GMM_LOR_raw",
         markersize = 4)

sf_ax_re = Axis(right_grid[3, 1],
                title = "Remora",
		height=200,
                aspect = 1)
xlims!(sf_ax_re, (0, 5))
sharkfin(sf_ax_re,
         discard_low_sig_positions(remora, 0.005), #[1:10000, :],
         "predicted_raw",
         "GMM_LOR_raw",
         markersize = 4)

evu4_cor = round(cor(eventalign_uncalled4.predicted,
                     eventalign_uncalled4.predicted_1);
                 digits = 2)
println("Eventalign/Uncalled4 correlation r = ", evu4_cor)
cor_ax_evu4 = Axis(right_grid[1, 2],
                   title = "Eventalign/Uncalled4", # (r=$evu4_cor)",
                   aspect = 1)
corr_plot(cor_ax_evu4,
          eventalign_uncalled4, #[1:10000, :],
          "Eventalign",
          "Uncalled4";
          markersize = 4)

evre_cor = round(cor(eventalign_remora.predicted,
                     eventalign_remora.predicted_1);
                 digits = 2)
println("Eventalign/Remora correlation r = ", evre_cor)
cor_ax_evre = Axis(right_grid[2, 2],
                   title = "Eventalign/Remora", # (r=$evre_cor)",
                   aspect = 1)
corr_plot(cor_ax_evre,
          eventalign_remora, #[1:10000, :],
          "Eventalign",
          "Remora";
          markersize = 4)

u4re_cor = round(cor(uncalled4_remora.predicted,
                     uncalled4_remora.predicted_1);
                 digits = 2)
println("Uncalled4/Remora correlation r = ", u4re_cor)
cor_ax_u4re = Axis(right_grid[3, 2],
                   title = "Uncalled4/Remora", #(r=$u4re_cor)",
                   aspect = 1)
corr_plot(cor_ax_u4re,
          uncalled4_remora, #[1:10000, :],
          "Uncalled4",
          "Remora";
          markersize = 4)

linkaxes!(sf_ax_ev, sf_ax_u4, sf_ax_re)

Legend(resq_fig[3, 1],
       [PolyElement(color = COL_GLORI_NEG),
        PolyElement(color = COL_GLORI_POS)],
       ["GLORI-", "GLORI+"],
       orientation = :horizontal,
       valign = 1,
       framevisible=false)

# rowgap!(resq_fig.layout, 0)
rowgap!(resq_fig.layout, 2, -10)
# colgap!(right_grid, -10)


for (label, layout) in zip(["A", "B", "C"],
			   [resq_fig[1, 1], right_grid[1, 1], right_grid[1, 2]])
  Label(layout[1, 1, TopLeft()],
        label,
        fontsize = 16,
        font = :bold,
        padding = (25, 25, 5, 0),
        halign = :right)
end


println("Saving the figure")
save("fig_resquigglers.png", resq_fig)
println("Done")

