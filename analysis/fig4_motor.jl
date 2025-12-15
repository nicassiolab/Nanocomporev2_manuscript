include("lib.jl")

using AlgebraOfGraphics


rna002_tx = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv",
						 "GMM_chi2_qvalue",
						 shift=2,
						 genomic_collapse=false)
rna002_p10_tx = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_p10_motor_v2.0.0/out_nanocompore_results.tsv",
							 "GMM_chi2_qvalue",
							 shift=2,
							 genomic_collapse=false) 

rna004_tx = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_fix_mod_clust_inferring/out_nanocompore_results.tsv",
						 "GMM_chi2_qvalue",
						 shift=4,
						 genomic_collapse=false) 

rna004_p12_tx = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_p12_motor_fix_gap_predictions_logdwell/out_nanocompore_results_gx.tsv",
							 "GMM_chi2_qvalue",
							 shift=4,
							 genomic_collapse=false)

function get_all_positions(res)
	all_positions = combine(groupby(res, :ref_id),
							:pos => (p -> minimum(p):maximum(p)) => :pos)
	all_positions = leftjoin(all_positions, rna004_tx, on = [:ref_id, :pos])
	sort(all_positions, [:ref_id, :pos])
end

function get_ref_crosscorrs(all_positions)
	fn = function(i, d)
		zip(-20:20,
		    crosscor(-log10.(coalesce.(i, 1)),
					 -log10.(coalesce.(d, 1)),
					 -20:20))
	end
	valid_refs = filter(r -> r.count > 40,
						combine(groupby(all_positions, :ref_id), nrow => :count)).ref_id |> Set
	valid_mask = map(r -> in(r, valid_refs), all_positions.ref_id)
	combine(groupby(all_positions[valid_mask, :], :ref_id),
			[:KS_intensity_pvalue, :KS_dwell_pvalue] => fn => [:offset, :xcorr])
end

function get_sig_corrs(res)
	sig_corrs = []
	t = copy(res)
	for i in -20:20
		t[!, "motor_pos"] = disallowmissing(t.pos) .+ i
		tj = innerjoin(t[t.KS_intensity_qvalue .!== missing .&& t.KS_intensity_qvalue .<= 0.01, :], t,
					   on=[:ref_id => :ref_id,
						   :motor_pos => :pos],
					   makeunique=true)
		tj = tj[tj.KS_dwell_pvalue_1 .!== missing, :]
		push!(sig_corrs, cor(tj.KS_intensity_pvalue, tj.KS_dwell_pvalue_1))
	end
	sig_corrs
end


function peak_annotate(df, ref, peak_radius)
	df_peaks = peaks(df, peak_radius)
	collapsed_peaks = combine(groupby(df_peaks, [:chr, :strand, :genomicPos]),
							  :predicted => maximum => :predicted)
	result = leftjoin(ref, collapsed_peaks,
					  on = [:chr, :strand, :genomicPos])
	result[!, :predicted] = coalesce.(result.predicted, 0)
	sort(result, [:chr, :strand, :genomicPos])
end

function plot_prc!(ax, gt, pred, thres; label = "", color = :black, linestyle = :dash, dot_thres = 2, linewidth = 2)
	rocs = roc(gt, pred, thres)
	precisions = map(precision, rocs) |> reverse
	recalls = map(recall, rocs) |> reverse
	push!(precisions, 0)
	push!(recalls, last(recalls))
	push!(precisions, 0)
	push!(recalls, 1)
	area = auc(recalls, precisions)
	lines!(ax, recalls, precisions, color = color, linestyle = linestyle, label = "$label (AUC=$(round(area; digits = 2)))", linewidth = linewidth)

	local evaluation = roc(gt, pred .> dot_thres)
	local p001 = precision(evaluation)
	local r001 = recall(evaluation)
	scatter!(ax, [r001], [p001], color = color)
end

function trapezoidal_rule(x, y)
	auc = 0.0
	for i in 2:length(x)
		auc += (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2
	end
	return auc
end

function auc(recall, precision)
    trapezoidal_rule(recall, precision)
end

function plot_peak_dist_hist!(fig, gridpos, res, reference, peaks, motor_peaks; title = "", radius = 15, col_palette = [:black, :red], style_palette = [:solid, :dash], titlesize = 28, linewidth = 2, maxy=nothing)
	ref = innerjoin(reference, res[:, [:chr, :strand, :genomicPos]],
						  on = [:chr, :strand, :genomicPos])
	mods = unique(ref[ref.modified .== 1, :])

	peaks2d = peaks[peaks.predicted .> 2, :]

	distances2d = map(r -> begin
					      distances = peaks2d[peaks2d.chr .== r.chr .&& peaks2d.strand .== r.strand, :genomicPos] .- r.genomicPos
						  push!(distances, 1000)
						  distances[argmin(abs.(distances))]
					  end,
					  eachrow(mods))

	peaks3d = motor_peaks[motor_peaks.predicted .> 2, :]

	distances3d = map(r -> begin
					      distances = peaks3d[peaks3d.chr .== r.chr .&& peaks3d.strand .== r.strand, :genomicPos] .- r.genomicPos
					  	  push!(distances, 1000)
					  	  distances[argmin(abs.(distances))]
					  end,
					  eachrow(mods))

	distances2d = distances2d[abs.(distances2d) .<= radius]
	distances3d = distances3d[abs.(distances3d) .<= radius]

	distances = vcat(DataFrame(distance=distances2d, type="MD−"),
					 DataFrame(distance=distances3d, type="MD+"))

	
	plt = data(distances) *
		mapping(:distance => "Distance",
			    color = :type  => presorted => "GMM test",
				linestyle = :type => presorted => "GMM test") * 
		AlgebraOfGraphics.histogram(Stairs;
									normalization=:pdf,
									bins=-radius:(radius+1)) *
		visual(; linewidth = linewidth)

	local f = draw!(gridpos,
					plt,
					scales(Color=(; palette = col_palette),
						   LineStyle=(; palette = style_palette));
					axis = (; title=title,
							  aspect=1,
							  titlesize=titlesize,
							  xticks=((-radius + 0.5):(radius+0.5),
									   [string(t) for t in -radius:radius]),
				   			  ylabel="PDF",
						      limits=(nothing, (0, maxy))))
	
	legend!(gridpos, f;
			tellwidth=false,
			halign=:right,
			valign=:top,
			margin=(10, 10, 10, 10),
			patchsize=(40, 40))
end

all_positions_002 = get_all_positions(rna002_tx)
all_positions_004 = get_all_positions(rna004_tx)

ref_crosscorrs_002 = get_ref_crosscorrs(all_positions_002)
ref_crosscorrs_004 = get_ref_crosscorrs(all_positions_004)


mean_xcorrs_002 = combine(groupby(ref_crosscorrs_002, :offset),
						  :xcorr => (vs -> mean(filter(v -> !isnan(v), vs))) => :mean_xcorr)
mean_xcorrs_004 = combine(groupby(ref_crosscorrs_004, :offset),
					  	  :xcorr => (vs -> mean(filter(v -> !isnan(v), vs))) => :mean_xcorr)

rna002_corrs_sig = get_sig_corrs(rna002_tx)
rna004_corrs_sig = get_sig_corrs(rna004_tx)

ref002 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv"))
ref002 = combine(groupby(ref002, ["chr", "strand", "genomicPos"]), :modified => maximum => :modified)
peakannot_002 = peak_annotate(rna002_tx, ref002, 4)

ref004 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv"))
ref004 = combine(groupby(ref004, ["chr", "strand", "genomicPos"]), :modified => maximum => :modified)
peakannot_004 = peak_annotate(rna004_tx, ref004, 4)

peakannot_002_binned = copy(peakannot_002)
peakannot_002_binned[!, :bin] = div.(peakannot_002.genomicPos, BIN_SIZE)
peakannot_002_binned = combine(groupby(peakannot_002_binned, [:chr, :strand, :bin]),
							   :modified => maximum => :modified,
							   :predicted => maximum =>  :predicted)
peakannot_002_p10 = peak_annotate(rna002_p10_tx, ref002, 4)
peakannot_002_p10_binned = copy(peakannot_002_p10)
peakannot_002_p10_binned[!, :bin] = div.(peakannot_002_p10.genomicPos, BIN_SIZE)
peakannot_002_p10_binned = combine(groupby(peakannot_002_p10_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)

peakannot_004_binned = copy(peakannot_004)
peakannot_004_binned[!, :bin] = div.(peakannot_004.genomicPos, BIN_SIZE)
peakannot_004_binned = combine(groupby(peakannot_004_binned, [:chr, :strand, :bin]),
							   :modified => maximum => :modified,
							   :predicted => maximum =>  :predicted)
peakannot_004_p12 = peak_annotate(rna004_p12_tx, ref004, 4)
peakannot_004_p12_binned = copy(peakannot_004_p12)
peakannot_004_p12_binned[!, :bin] = div.(peakannot_004_p12.genomicPos, BIN_SIZE)
peakannot_004_p12_binned = combine(groupby(peakannot_004_p12_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)

COL_ALL_SITES = :black # "#7f7f7f22"
COL_MOD_SITES = "#00A0A9" #"#FFD700"
COL_2D = "#004488"  # "#3E8241" # "#17becf"
COL_3D = "#D75F00" # "#6948A3" # "#bcbd22"

fig = Figure(size = (1150, 1850), fontsize = 26)

cor_maxy = 0.21
markersize = 9

ax = Axis(fig[1, 1],
		  title="RNA002",
		  titlesize=34)
hidedecorations!(ax)
hidespines!(ax)
ax = Axis(fig[1, 2],
		  title="RNA004",
		  titlesize=34)
hidedecorations!(ax)
hidespines!(ax)


ax = Axis(fig[2, 1],
		  titlesize=20,
		  aspect=1,
		  xlabel="Offset",
		  ylabel="Correlation",
		  yticks=0:0.02:cor_maxy)
ylims!(ax, (-0.01, cor_maxy))

scatter!(ax, mean_xcorrs_002.offset, mean_xcorrs_002.mean_xcorr, markersize=markersize, color=COL_ALL_SITES, label="All sites")
lines!(ax, mean_xcorrs_002.offset, mean_xcorrs_002.mean_xcorr, color=COL_ALL_SITES, linewidth=2.5)

scatter!(ax, -20:20, rna002_corrs_sig, markersize=markersize, color=COL_MOD_SITES, label="Modified sites")
lines!(ax, -20:20, rna002_corrs_sig, color=COL_MOD_SITES, linestyle=:dash, linewidth=2.5)

legend_elem_all = [LineElement(color=COL_ALL_SITES, linestyle=:solid),
				   MarkerElement(color=COL_ALL_SITES, marker=:circle, markersize=15, strokecolor=COL_ALL_SITES)]
legend_elem_mod = [LineElement(color=COL_MOD_SITES, linestyle=:dash),
				   MarkerElement(color=COL_MOD_SITES, marker=:circle, markersize=15, strokecolor=COL_MOD_SITES)]
Legend(fig[2, 1],
	   [legend_elem_all, legend_elem_mod],
	   ["All sites", "Modified sites"],
	   patchsize=(40, 40),
	   rowgap=5,
	   tellwidth=false,
	   halign=:right,
	   valign=:top,
	   margin=(10, 10, 10, 10))

ax = Axis(fig[2, 2],
		  titlesize=20,
		  aspect=1,
		  xlabel="Offset",
		  ylabel="Correlation",
		  yticks=0:0.02:cor_maxy)
ylims!(ax, (-0.01, cor_maxy))

scatter!(ax, mean_xcorrs_004.offset, mean_xcorrs_004.mean_xcorr, markersize=markersize, color=COL_ALL_SITES, label="All sites")
lines!(ax, mean_xcorrs_004.offset, mean_xcorrs_004.mean_xcorr, color=COL_ALL_SITES, linewidth=2.5)

scatter!(ax, -20:20, rna004_corrs_sig, markersize=markersize, color=COL_MOD_SITES, label="Modified sites")
lines!(ax, -20:20, rna004_corrs_sig, color=COL_MOD_SITES, linestyle=:dash, linewidth=2.5)

Legend(fig[2, 2],
	   [legend_elem_all, legend_elem_mod],
       ["All sites", "Modified sites"],
       patchsize=(40, 40),
       rowgap=5,
       tellwidth=false,
       halign=:right,
       valign=:top,
       margin=(10, 10, 10, 10))

ax = Axis(fig[3, 1],
		  titlesize=20,
		  aspect=1,
		  xlabel="Recall",
		  ylabel="Precision")
plot_prc!(ax,
		  peakannot_002.modified,
		  peakannot_002.predicted,
		  sort(filter(v -> v > 0, peakannot_002.predicted));
		  linestyle=:solid,
		  linewidth=2.5,
		  color=COL_2D,
		  label="MD−")
plot_prc!(ax,
		  peakannot_002_p10.modified,
		  peakannot_002_p10.predicted,
		  sort(filter(v -> v > 0, peakannot_002_p10.predicted));
		  label="MD+",
		  linestyle=:dash,
		  linewidth=2.5,
		  color=COL_3D)

axislegend(ax, "GMM test", patchsize=(40, 40))

ax = Axis(fig[3, 2],
		  titlesize=20,
		  aspect=1,
		  xlabel="Recall",
		  ylabel="Precision")
plot_prc!(ax,
		  peakannot_004.modified,
		  peakannot_004.predicted,
		  sort(filter(v -> v > 0, peakannot_004.predicted));
		  linestyle=:solid,
		  linewidth=2.5,
		  color=COL_2D,
		  label="MD−")
plot_prc!(ax,
		  peakannot_004_p12.modified,
		  peakannot_004_p12.predicted,
		  sort(filter(v -> v > 0, peakannot_004_p12.predicted));
		  label="MD+",
		  linestyle=:dash,
		  linewidth=2.5,
		  color=COL_3D)

axislegend(ax, "GMM test", patchsize=(40, 40))

plot_peak_dist_hist!(fig, fig[4, 1], rna002_tx, ref002, peakannot_002, peakannot_002_p10; radius=5, col_palette=[COL_2D, COL_3D], titlesize=28, linewidth=2.5, maxy=0.85)
plot_peak_dist_hist!(fig, fig[4, 2], rna004_tx, ref004, peakannot_004, peakannot_004_p12; radius=5, col_palette=[COL_2D, COL_3D], titlesize=28, linewidth=2.5, maxy=0.85)

rowsize!(fig.layout, 1, Fixed(20))

for (label, layout) in zip(["B", "C", "D"],
						   [fig[2, 1], fig[3, 1], fig[4, 1]])
	Label(layout[1, 1, TopLeft()], label,
		  fontsize=26,
		  font=:bold,
		  padding=(50, 30, 5, 0),
		  halign=:right)
end

save("fig4_motor.png", fig)
save("fig4_motor.pdf", fig)


f = Figure(size = (900, 450))


ax = Axis(f[1, 1],
		  title="RNA002",
		  aspect=1)
plot_prc!(ax,
	      peakannot_002_binned.modified,
	      peakannot_002_binned.predicted,
	      sort(filter(v -> v > 0, peakannot_002_binned.predicted));
	      linestyle=:solid,
	      label="MD−",
	      color=COL_2D)
plot_prc!(ax,
	      peakannot_002_p10_binned.modified,
	      peakannot_002_p10_binned.predicted,
	      sort(filter(v -> v > 0, peakannot_002_p10_binned.predicted));
	      label="MD+",
	      linestyle=:dash,
	      color=COL_3D)

axislegend(ax, "GMM test}")


ax = Axis(f[1, 2],
		  title="RNA004",
		  aspect=1)
plot_prc!(ax,
	      peakannot_004_binned.modified,
	      peakannot_004_binned.predicted,
	      sort(filter(v -> v > 0, peakannot_004_binned.predicted));
	      linestyle=:solid,
	      label="MD−",
	      color=COL_2D)
plot_prc!(ax,
	      peakannot_004_p12_binned.modified,
	      peakannot_004_p12_binned.predicted,
	      sort(filter(v -> v > 0, peakannot_004_p12_binned.predicted));
	      label="MD+",
	      linestyle=:dash,
	      color=COL_3D)

axislegend(ax, "GMM test")

save("supplemental_fig3_binned_motor.png", f)
save("supplemental_fig3_binned_motor.pdf", f)
