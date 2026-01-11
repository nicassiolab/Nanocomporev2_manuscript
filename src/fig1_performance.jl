using CairoMakie
using Printf

f = Figure(size = (640, 340))

ax_collapse = Axis(f[1, 1],
		   aspect = 1/2.6,
		   title = "Preprocessing\ntime",
		   #titlesize = 13,
		   ylabel = "Execution time (s) / transcript",
		   xticks = ([1, 2], ["v1", "v2"]),
		   yticks = 0:0.5:2.5)
ylims!(ax_collapse, (-0.1, 2.7))
ax_space = Axis(f[1, 2],
		aspect = 1/2.6,
		title = "Preprocessing\ndisk space",
		ylabel = "Disk space (GB) / 1 million reads",
		xticks = ([1, 2], ["v1", "v2"]),
		yticks = 10:10:60)
ylims!(ax_space, (-1, 67))
ax_comp = Axis(f[1, 3],
	       aspect = 1/2.6,
	       title = "Analysis\ntime",
	       ylabel = "Execution time (s) / transcript",
	       xticks = ([1, 2], ["v1", "v2"]),
	       yticks = 0:5:30)
ylims!(ax_comp, (-1, 31))
ax_space2 = Axis(f[1, 4],
		 aspect = 1/2.6,
		 title = "Analysis\ndisk space",
		 ylabel = "Disk space (MB) / 1,000 tested sites",
		 xticks = ([1, 2], ["v1", "v2"]),
		 yticks = 0:5:20)
ylims!(ax_space2, (-0.05, 22.2))

t1 = (35*60+3)/1000
t2 = (7*60)/1000
barplot!(ax_collapse,
	 #[(35*60+3)/60, (7*60)/60],
	 [t1, t2],
	 color = :black,
	 bar_labels = ["$(@sprintf("%.2f", t1))s",
		       "$(@sprintf("%.2f", t2))s"])

t1 = (59*60+57)/136
t2 = (3*60+21)/136
barplot!(ax_comp,
	 # [59u"minute" + 57u"s", 3u"minute" + 21u"s"],
	 # [(59*60+57)/60, (3*60+21)/60],
	 [t1, t2],
	 color = :black,
	 bar_labels = ["$(@sprintf("%.2f", t1))s",
		       "$(@sprintf("%.2f", t2))s"])

barplot!(ax_space,
	 [55.52, 19.43],
	 color = :black,
	 bar_labels = ["56GB", "19GB"])

barplot!(ax_space2,
	 [17, 1],
	 color = :black,
	 bar_labels = ["17MB", "1MB"])

colgap!(f.layout, 0.7)

for (label, layout) in zip(["A", "B", "C", "D"],
			   [f[1, 1], f[1, 2], f[1, 3], f[1, 4]])
  Label(layout[1, 1, TopLeft()],
        label,
        fontsize = 16,
        font = :bold,
        padding = (25, 25, 5, 0),
        halign = :right)
end

save("fig1_performance.png", f)
save("fig1_performance.pdf", f)


