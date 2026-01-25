include("load_comods_data.jl")


glori_raw1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv"))
glori_raw2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv"))
glori_raw = innerjoin(glori_raw1, glori_raw2,
					  on = [:Chr, :Strand, :Sites],
					  makeunique = true)


cancer_genes = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/cancerGeneList.tsv", delim = '\t', normalizenames = true))


LABELPAD = 0

f = Figure(size = (1600, 1250), fontsize=26)

# PANEL A
panel_a = f[1, 1] = GridLayout()

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
	mapping(:glori_mod_ratio,
			:mod_ratio => "IVT ratio",
			color = :distance => "Distance to the GLORI site") *
	visual(Scatter, markersize = 7, alpha = 0.5, rasterize=2)
draw!(panel_a[1, 1],
	  plt1;
	  axis=(;
			#title = "Nanocompore (Native/IVT) vs. GLORI",
			aspect = 1,
			xlabel = "GLORI stoichiometry",
			ylabel = "Nanocompore stoichiometry",
			xticks = 0:0.2:1,
			yticks = 0:0.2:1))

Colorbar(panel_a[1, 2], limits = [0, 4], colormap = :viridis,
vertical = true, label = "Distance to the GLORI site", labelpadding = LABELPAD)

# PANEL B

sig_color = "#2E2585"
insig_color = "#DDDDDD"


df = copy(comods)
df[!, :significant] = df.qvalue .< 0.05
plt = data(df) *
mapping(:phi, :pvalue => (p -> -log10(p)), color = :significant) *
visual(Scatter; markersize=7, rasterize=2)
draw!(f[2, 1],
  plt,
  scales(Color = (; palette = [insig_color, sig_color]));
  axis = (; #title = "All modification pairs",
			xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
			ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
			limits = ((-1, 1), (-2, 310))))


# PANEL C

panel_c = f[1:2, 2] = GridLayout()

df = filter(r -> abs(r.phi) > 0.1, sig_comods)
df[!, :gene] = [x[6] for x in split.(df.reference, "|")]

df_cancer = innerjoin(df, cancer_genes, on = [:gene => :Hugo_Symbol])

df[!, :log10_distance] = log10.(abs.(df.pos1 .- df.pos2))
df[!, :log10_distance_bin] = round.(df.log10_distance ./ 0.5; digits = 0)
df[!, :combined_score] = 100 .* abs.(df.phi) .+ clamp.(-log10.(df.pvalue), 0, 100)
labels_per_bin = 5
labeled = combine(groupby(df, [:log10_distance_bin]),
					# df -> vcat(sort(df[df.phi .>= 0, :], :combined_score, rev = true)[1:(df[1, :log10_distance_bin] == 2 ? 2*labels_per_bin : labels_per_bin), :],
					# 		   sort(df[df.phi .< 0, :], :combined_score, rev = true)[1:(df[1, :log10_distance_bin] == 2 ? 2*labels_per_bin : labels_per_bin), :]))
					df -> vcat(sort(df[df.phi .>= 0, :], :phi, by = abs, rev = true)[1:(df[1, :log10_distance_bin] <= 2 ? 0 : labels_per_bin), :],
							   sort(df[df.phi .< 0, :], :phi, by = abs, rev = true)[1:(df[1, :log10_distance_bin] == 2 ? labels_per_bin : labels_per_bin), :]))
high_significance = df[df.pvalue .< 1e-100 .&&
						 (.! (df.log10_distance_bin .<= 2 .&&
							 df.phi .> 0 .&&
							 df.phi .< 0.6)), :]
labeled = unique(vcat(labeled, high_significance))
cancer_geneset = Set(cancer_genes.Hugo_Symbol)
labeled[!, :cancer] = [gene in cancer_geneset for gene in labeled.gene]



ax = Axis(panel_c[1, 1],
			xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
			ylabel = "Distance between associated modifications",
			yticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
			yminorticksvisible = true,
			yminorgridvisible = true,
			yminorticksize = 2,
			yminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]))
plt = data(df) *
mapping(:phi,
		(:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))),
		color = :pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)") *
visual(Scatter; markersize = 8, rasterize=2)
draw!(ax, #panel_c[1, 1],
  plt,
 scales(Color = (; colormap = :seaborn_crest_gradient, colorrange = (0, 100), highclip = :black)))
 #  axis = (;
	# yticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
	# yminorticksvisible = true,
 #    yminorgridvisible = true,
	# yminorticksize = 2,
	# yminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]),
	# limits = ((-1, 1), (0.9, 4.2))))
Colorbar(panel_c[1, 2], limits = [0, 100], colormap = :seaborn_crest_gradient, #:viridis,
vertical = true, label = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", labelpadding = LABELPAD)
colsize!(panel_c, 1, Relative(8/9))

points = Point2f.(labeled.phi, labeled.log10_distance)
labels = labeled.gene
label_positions = place_labels_nonoverlapping(points, labels; offset = 0.1, gene_offsets = Dict("SRSF2" => 0.45), height=0.075, width_factor=0.0475)
for (p, lp, label, cancer) in zip(points, label_positions, labels, labeled.cancer)
		xalign = (lp[1] > p[1] ? :left : :right)
		text!(ax, label, position=lp, align=(xalign, :center), color=:black, fontsize = 14, font = cancer ? :bold : :italic)
		lines!(ax, [p, lp], color=:gray, linewidth=0.5)
end
xlims!(ax, -1.55, 1.55)


colsize!(f.layout, 1, Relative(8/18))
colgap!(f.layout, 0)


for (label, layout) in zip(["A", "B", "C"],
					   [f[1, 1], f[2, 1], f[1, 2]])
Label(layout[1, 1, TopLeft()], label,
	fontsize = 26,
	font = :bold,
	padding = (0, 5, 5, 0),
	halign = :right)
end

	

save("fig5_comods.png", f)
save("fig5_comods.pdf", f)
