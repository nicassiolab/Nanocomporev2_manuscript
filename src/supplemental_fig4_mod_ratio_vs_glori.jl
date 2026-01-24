include("lib.jl")
using AlgebraOfGraphics


ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=false,
		LOR_threshold = 0.8)
ivt[!, :pos] .+= 4
ivt[!, :mod_ratio] = sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))

peaks_ivt = peaks(ivt, 4)
sig_peaks_ivt = peaks_ivt[peaks_ivt.predicted .!== missing .&& peaks_ivt.predicted .>= -log10(0.01), :]


glori_raw1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv"))
glori_raw2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv"))

glori_raw = innerjoin(glori_raw1, glori_raw2,
					  on=[:Chr, :Strand, :Sites],
					  makeunique=true)


f = Figure(size=(500, 500))


glori_ratios = Dict(map(r -> (r.Chr, r.Strand, r.Sites) => (r.Ratio + r.Ratio_1)/2,
			eachrow(glori_raw)))

df = copy(sig_peaks_ivt)
closest_m6A = map(
	r -> let offset_results = [(get(glori_ratios, (r.chr, r.strand, r.genomicPos + offset), missing), offset)
					   for offset in [0, 1, -1, 2, -2, 3, -3, 4, -4]]
		ratios = map(first, offset_results)
		offsets = map(last, offset_results)
		non_missing = ratios .!== missing
		# existing = collect(skipmissing(ratios))
		any(non_missing) ? (first(ratios[non_missing]), first(offsets[non_missing])) : (0.0, missing)
	end,
	eachrow(df))
df[!, :glori_mod_ratio] = map(first, closest_m6A)
df[!, :distance] = abs.(map(last, closest_m6A))

df = df[df.glori_mod_ratio .> 0, :]

c = round(cor(df.glori_mod_ratio, df.mod_ratio), digits = 2)

plt1 = data(df) *
	mapping(:glori_mod_ratio, :mod_ratio => "IVT ratio", color = :distance => "Distance to the GLORI site") *
	visual(Scatter, markersize = 5, alpha = 0.5, rasterize=2) |>
	draw(; axis=(; title="Nanocompore (Native/IVT) vs. GLORI", aspect=1, xlabel="GLORI mod. ratio", ylabel="Nanocompore mod. ratio", xticks=0:0.1:1, yticks=0:0.1:1), figure=(size=(550, 500),))


save("supplemental_fig4_mod_ratio_vs_GLORI.png", plt1)
save("supplemental_fig4_mod_ratio_vs_GLORI.pdf", plt1)
