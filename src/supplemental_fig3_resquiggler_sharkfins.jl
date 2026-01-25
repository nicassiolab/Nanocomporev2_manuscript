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
                                 on="bin",
                                 makeunique=true)

eventalign_remora = innerjoin(binned_eventalign,
                              binned_remora,
                              on="bin",
                              makeunique=true)

uncalled4_remora = innerjoin(binned_uncalled4,
                             binned_remora,
                             on="bin",
                             makeunique=true)

println("Plotting")

fig = Figure(size=(1000, 600), fontsize=18)


sf_ax_ev_p = Axis(fig[1, 1],
				  title="Eventalign",
				  aspect=1)
sharkfin(sf_ax_ev_p,
         discard_low_sig_positions(eventalign[eventalign.modified .== 1, :], 0.005),
         "GMM_chi2_pvalue",
         "GMM_LOR_raw",
         markersize=4,
		 rasterize=2)

sf_ax_ev_n = Axis(fig[2, 1],
				  aspect=1)
sharkfin(sf_ax_ev_n,
         discard_low_sig_positions(eventalign[eventalign.modified .== 0, :], 0.005),
         "GMM_chi2_pvalue",
         "GMM_LOR_raw",
         markersize=4,
		 rasterize=2)
xlims!(sf_ax_ev_p, (0, 5))
xlims!(sf_ax_ev_n, (0, 5))

sf_ax_u4_p = Axis(fig[1, 2],
				  title="Uncalled4",
				  aspect=1)
sharkfin(sf_ax_u4_p,
         discard_low_sig_positions(uncalled4[uncalled4.modified .== 1, :], 0.005),
         "GMM_chi2_pvalue",
         "GMM_LOR_raw",
         markersize=4,
		 rasterize=2)

sf_ax_u4_n = Axis(fig[2, 2],
				  aspect=1)
sharkfin(sf_ax_u4_n,
         discard_low_sig_positions(uncalled4[uncalled4.modified .== 0, :], 0.005),
         "GMM_chi2_pvalue",
         "GMM_LOR_raw",
         markersize=4,
		 rasterize=2)
xlims!(sf_ax_u4_p, (0, 5))
xlims!(sf_ax_u4_n, (0, 5))

sf_ax_re_p = Axis(fig[1, 3],
				  title="Remora",
				  aspect=1)
sharkfin(sf_ax_re_p,
         discard_low_sig_positions(remora[remora.modified .== 1, :], 0.005),
         "GMM_chi2_pvalue",
         "GMM_LOR_raw",
         markersize=4,
		 rasterize=2)

sf_ax_re_n = Axis(fig[2, 3],
				  aspect=1)
sharkfin(sf_ax_re_n,
         discard_low_sig_positions(remora[remora.modified .== 0, :], 0.005),
         "GMM_chi2_pvalue",
         "GMM_LOR_raw",
         markersize=4,
		 rasterize=2)
xlims!(sf_ax_re_p, (0, 5))
xlims!(sf_ax_re_n, (0, 5))


linkaxes!(sf_ax_ev_p, sf_ax_u4_p, sf_ax_re_p, sf_ax_ev_n, sf_ax_u4_n, sf_ax_re_n)

Legend(fig[3, 1:3],
       [PolyElement(color=COL_GLORI_NEG),
        PolyElement(color=COL_GLORI_POS)],
       ["GLORI-", "GLORI+"],
       orientation=:horizontal,
       # valign=1,
       framevisible=false)

# rowgap!(fig.layout, 0)
# rowgap!(fig.layout, 2, -10)
# rowgap!(fig, 1, -7)
# rowgap!(fig, 2, -7)
# colgap!(fig, -5)
# colsize!(fig.layout, 1, 500)
rowsize!(fig.layout, 3, 6)


for (label, layout) in zip(["A", "B", "C"],
						   [fig[1, 1], fig[1, 2], fig[1, 3]])
  Label(layout[1, 1, TopLeft()],
        label,
        fontsize=18,
        font=:bold,
        padding=(15, 40, 0, 0),
        halign=:right)
end


println("Saving the figure")
save("supplemental_fig3_resquigglers_sharkfins.png", fig)
save("supplemental_fig3_resquigglers_sharkfins.pdf", fig)
println("Done")

