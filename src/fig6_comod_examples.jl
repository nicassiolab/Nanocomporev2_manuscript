include("load_comods_data.jl")

f = Figure(size = (1600, 1250), fontsize=22)
trow = f[1, 1] = GridLayout()
brow = f[2, 1] = GridLayout()

LABELPAD = 0

function contingency_to_df(M)
    rows = Tuple{Int,Int}[]
    for i in 1:2, j in 1:2
        append!(rows, fill((i-1, j-1), M[i,j]))
    end
    return DataFrame(rows, [:x, :y])
end

# PANEL A

ax = Axis(trow[1, 1])


plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, an_sig_comods), gtf, "ENST00000234875.9"; minheight = 5, maxheight = 150, dotsize=15, xlimits=(1170, 1820), legend=true, colormap = :seaborn_crest_gradient, highclip=:darkblue, fig=trow[1, 1], axis = ax, fontsize=20) # RPL22
ylims!(ax, -180, 150)
text!(ax, 0, -30; text = "RPL22", align = (:left, :top))
text!(ax, 1720, 130; text = "↑ Co-occurring", align = (:left, :top))
text!(ax, 1720, -130; text = "↓ Mutually exclusive", align = (:left, :bottom))
ax.xlabelpadding = -20

inset_ax = Axis(trow[1, 1],
					  width=Relative(0.5),
					  height=Relative(0.2),
					  halign=0.05,
					  valign=0.95)
hidedecorations!(inset_ax)
tx = "ENST00000234875.9"
name = "RPL22"
gene_start = minimum(gtf[gtf.transcript_id .== tx, :start])

plot_isoforms_model!(inset_ax, "RPL22"; transcripts = Set([tx]), colors = Dict([tx => :black]), rename = Dict([tx => name]), focus=(gene_start+1170, gene_start+1820), fontsize=16, lpadding=2000)


# lines!(ax, [2200, 3200], [1000, 1000], color = :black)
Colorbar(trow[1, 2], limits = [0, 100], colormap = :seaborn_crest_gradient, #:viridis,
vertical = true, label = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", labelpadding = LABELPAD)

ax.xticks = ([1200, 1800],
	 [format_with_commas(Integer(round(gene_start+1200; digits = 0))),
	  format_with_commas(Integer(round(gene_start+1800; digits = 0)))])
hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)


# PANELS B, C, and D


ax = Axis(brow[1, 1:2],
		  xtickformat = "{:,d}")# (vals -> map(format_with_commas, vals)))
		
gpos1 = 76066564
gpos2 = 76066573
row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
tx1 = split(row1.reference, "|")[1]
tx2 = split(row2.reference, "|")[1]
name1 = split(row1.reference, "|")[5]
name2 = split(row2.reference, "|")[5]

blue = "#00aaff"

plot_isoforms_model!(ax, "MDH"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => blue, tx2 => :orange]), focus = (gpos1-100, gpos2+100), rename = Dict([tx1 => name1, tx2 => name2]), lpadding=1650)

tx1_mask = occursin.(tx1, an_sig_comods.reference)
tx2_mask = occursin.(tx2, an_sig_comods.reference)
colorrange = (0, 100)
seq = tx_ref[row1.reference][(row1.pos1-10):(row1.pos2+10)]
println(seq)
mod_symbols = Dict(["m⁵C" => :utriangle,
						  "m⁶A" => :cross,
						  "?" => :circle,
						  "Ψ" => :rect])
mod_colors = Dict(["m⁵C" => "#e69f00",
						 "m⁶A" => "#009e73",
						 "?" => "#0072b2",
						 "Ψ" => "#cc79a7"])
ax1, max_radius1 = arc_plot_genomic2(brow[2, 1],
				  renamemods(an_sig_comods[tx1_mask, :]),
				  an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
				  range=(row1.pos1-2000, row1.pos2+2000),
				  grange=(gpos1-10, gpos2+10),
				  colorrange=colorrange,
				  title=name1,
				  highlight=(row1.pos1, row1.pos2),
				  spinecolor=blue,
				  ticks = 3,
				  sequence=seq,
				  mod_symbols=mod_symbols,
				  mod_colors=mod_colors)
text!(ax1, 1, 14; text = "↑ Co-occurring", align = (:left, :top))
text!(ax1, 1, -14; text = "↓ Mutually exclusive", align = (:left, :bottom))
Legend(brow[2, 1],
   [MarkerElement(marker=mod_symbols["?"], markersize=16, color=mod_colors["?"]),
		MarkerElement(marker=mod_symbols["m⁵C"], markersize=16, color=mod_colors["m⁵C"])],
   ["Unclassified", "m⁵C"],
	   "Modification";
	   tellheight=false,
	   tellwidth=false,
	   halign=:right,
	   valign=:bottom,
	   orientation=:horizontal,
	   framevisible=false,
	   titlegap=2,
	   colgap=4,
	   patchlabelgap=2)
ax2, max_radius2 = arc_plot_genomic2(brow[3, 1],
				  renamemods(an_sig_comods[tx2_mask, :]),
				  an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
				  range=(row2.pos1-2000, row2.pos2+2000),
				  grange=(gpos1-10, gpos2+10),
				  colorrange=colorrange,
				  title=name2,
				  highlight=(row2.pos1, row2.pos2),
				  spinecolor=:orange,
				  ticks = 3,
				  sequence=seq,
				  mod_symbols=mod_symbols,
				  mod_colors=mod_colors)
text!(ax2, 1, 14; text = "↑ Co-occurring", align = (:left, :top))
text!(ax2, 1, -14; text = "↓ Mutually exclusive", align = (:left, :bottom))
Legend(brow[3, 1],
   [MarkerElement(marker=mod_symbols["?"], markersize=16, color=mod_colors["?"]),
		MarkerElement(marker=mod_symbols["m⁵C"], markersize=16, color=mod_colors["m⁵C"])],
   ["Unclassified", "m⁵C"],
	   "Modification";
	   tellheight=false,
	   tellwidth=false,
	   halign=:right,
	   valign=:bottom,
	   orientation=:horizontal,
	   framevisible=false,
	   titlegap=2,
	   colgap=4,
	   patchlabelgap=2)
max_radius = max(max_radius1, max_radius2) * 2
ylims!(ax1, -max_radius, max_radius)
ylims!(ax2, -max_radius, max_radius)
Colorbar(brow[2:3, 2], limits = colorrange, colormap = :seaborn_crest_gradient, #viridis
vertical = true, label = "-log₁₀(𝑃-value)", labelpadding = LABELPAD)


# PANEL E

ref = row1.reference
# probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
# df1 = DataFrame(Tables.table(probs[:, [row1.pos1+1-4, row1.pos2+1-4]] .> 0.5))
# dropmissing!(df1)
# observed1 = contingency(df1.Column1, df1.Column2) .+ 1

observed1 = [row1.a row1.b;
			 row1.c row1.d]
df1 = contingency_to_df(observed1)
res1 = HypothesisTests.ChisqTest(observed1)
data1 = map(r -> join(Int.(collect(r)), "-"), eachrow(df1))
counts1 = data1 |> countmap
println("$name1 observations:")
println(observed1)

ref = row2.reference
# probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
# df2 = DataFrame(Tables.table(probs[:, [row2.pos1+1-4, row2.pos2+1-4]] .> 0.5))
# dropmissing!(df2)
# observed2 = contingency(df2.Column1, df2.Column2) .+ 1

observed2 = [row2.a row2.b;
			 row2.c row2.d]
df2 = contingency_to_df(observed2)
res2 = HypothesisTests.ChisqTest(observed2)
data2 = map(r -> join(Int.(collect(r)), "-"), eachrow(df2))
counts2 = data2 |> countmap
println("$name2 observations:")
println(observed2)


# fig = Figure(size = (600, 600))

grid = brow[1:3, 3] = GridLayout()	
top = Axis(grid[1, 1],
				 xlabel = "",
				 xticks = (1:4, ["", "", "", ""]),
				 ylabel = "Pattern occurrences")
xlims!(top, 0.5, 4.5)
ylims!(top, 0, 2000)

df1 = DataFrame(tx = name1, pattern = data1)
df2 = DataFrame(tx = name2, pattern = data2)
df = vcat(df1, df2)
df = combine(groupby(df, [:tx, :pattern]), nrow => :count)
barplt = data(df) *
	mapping(:pattern, :count, dodge = :tx, color = :tx) *
	visual(BarPlot; width = 0.7, dodge_gap = 0.15)
draw!(top, barplt, scales(Color = (; palette = [blue, :orange])))

expected = DataFrame(tx = [name1, name1, name1, name1,
								 name2, name2, name2, name2],
						   pattern = ["0-0", "0-1", "1-0", "1-1",
									  "0-0", "0-1", "1-0", "1-1"],
						   count = [res1.expected[1, 1],
								    res1.expected[1, 2],
								    res1.expected[2, 1],
								    res1.expected[2, 2],
								    res2.expected[1, 1],
								    res2.expected[1, 2],
								    res2.expected[2, 1],
								    res2.expected[2, 2]])
println("Expected: \n", expected)
barplt = data(expected) *
	mapping(:pattern, :count, dodge = :tx) *
	visual(BarPlot;
		   width = 0.7,
		   dodge_gap = 0.15,
		   color = (:white, 0),
		   strokewidth = 1.5,
		   strokecolor = :black)
draw!(top, barplt)

order = ["0-0", "0-1", "1-0", "1-1"]
bottom = Axis(grid[2, 1],
			yticks = (1:2, map(format_with_commas, [gpos2, gpos1])))
xlims!(bottom, 0.5, 4.5)
ylims!(bottom, 0.5, 2.5)
scatter!(bottom,
		 [1, 1, 2, 2, 3, 3, 4, 4],
		 [1, 2, 1, 2, 1, 2, 1, 2],
		 # [Point2f(1, 1),
		 #  Point2f(2, 1),
		 #  Point2f(1, 2),
		 #  Point2f(2, 2),
		 #  Point2f(1, 3),
		 #  Point2f(2, 3),
		 #  Point2f(1, 4),
		 #  Point2f(2, 4)];
		 color = [:white, :white,
				  :white, :black,
				  :black, :white,
				  :black, :black],
		 strokecolor = :black,
		 strokewidth = 1.5,
		 markersize = 25)
# hidexdecorations!(top)
hidexdecorations!(bottom)
hidespines!(bottom)

rowsize!(grid, 1, Relative(5/6))
rowsize!(grid, 2, Relative(1/6))

Legend(
grid[1, 1],
	[PolyElement(color=blue),
	 PolyElement(color=:orange),
	 PolyElement(color=(:white, 0), strokecolor=:black, strokewidth=1.5)],
	[name1, name2, "Expected"],
# "$ha & $va",
tellheight = false,
tellwidth = false,
	halign = :left,
	valign = :top,
margin = (10, 10, 10, 10),
)


colsize!(brow, 1, Relative(12/18))
colsize!(brow, 2, Relative(0.25/18))
colsize!(brow, 3, Relative(5.75/18))

rowsize!(f.layout, 1, Relative(3/8))
rowsize!(f.layout, 2, Relative(5/8))

rowsize!(brow, 1, Relative(1/7))
rowsize!(brow, 2, Relative(3/7))
rowsize!(brow, 3, Relative(3/7))


for (label, layout) in zip(["A", "B", "C", "D", "E"],
						   [trow[1, 1], brow[1, 1], brow[2, 1], brow[3, 1], brow[1, 3]])
    Label(layout[1, 1, TopLeft()], label,
	fontsize = 26,
	font = :bold,
	padding = (0, 5, 5, 0),
	halign = :right)
end


save("fig6_comod_examples.png", f)
save("fig6_comod_examples.pdf", f)
