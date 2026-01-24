using Pkg

Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("MLBase")
Pkg.add("CairoMakie")
Pkg.add("StatsBase")
Pkg.add("Unitful")
Pkg.add("Peaks")


using DataFrames
using CSV
using MLBase
using CairoMakie
using Random
using StatsBase
using Unitful
using Printf
using Peaks


BIN_SIZE = 9
LOR_THRESHOLD = 0.8

COL_RNA002 = "#56B4E9"
COL_RNA004 = "#009E73"

COL_GLORI_POS = "#2E2585"
COL_GLORI_NEG = "#DDDDDD"


function get_binned_glori(path)
  glori = DataFrame(CSV.File(path))
  # glori[!, "pos"] = glori.start .+ 2
  glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.start)
  glori[!, "modified"] .= 1
  binned_glori = unique(glori[!, ["chr", "strand", "bin", "modified"]])
  binned_glori
end

function read_results(filepath, col; shift = 0, genomic_collapse = true, LOR_threshold = LOR_THRESHOLD)
  df = DataFrame(CSV.File(filepath))
  if shift != 0
  	df[!, "genomicPos"] .+= ifelse.(df[!, "strand"] .== "+", shift, -shift)
  end
  # df = df[.! ismissing.(df[!, col]), :]
  df[:, col] = coalesce.(df[:, col], 1)

  # df[df[!, col] .== 0, col] .= eps(Float64)
  df[df[!, col] .== 0, col] .= 1e-300
  df[ismissing.(df[!, "GMM_LOR"]), "GMM_LOR"] .= 0
  df[!, "predicted_raw"] = -log10.(df[!, col])
  lor_corrected_significance = ifelse.(abs.(df[!, "GMM_LOR"]) .>= LOR_threshold, df[!, col], 1.0)
  df[!, "predicted"] = -log10.(lor_corrected_significance)
  additional_cols = []
  col_selector = (predicted, col) -> col[findmax(predicted)[2]]
  if any(map(n -> occursin("_mod", n), names(df)))
  	df[!, "mean_cov"] = sum.(eachrow(df[:, filter(n -> occursin("_mod", n) || occursin("_unmod", n), names(df))]))
  	push!(additional_cols, [:predicted_raw, :mean_cov] => col_selector => :mean_cov)
  end
  qvals = filter(n -> occursin("qvalue", n) || occursin("pvalue", n), names(df))
  if genomic_collapse
    return combine(groupby(df, ["chr", "strand", "genomicPos"]),
                   :predicted => maximum => :predicted,
                   [:predicted, :GMM_LOR] => col_selector => :GMM_LOR,
                   :predicted_raw => maximum => :predicted_raw,
                   [:predicted_raw, :GMM_LOR] => col_selector => :GMM_LOR_raw,
                   [n => (vs -> minimum(skipmissing(vs), init = 1.0)) => n for n in qvals]...,
                   additional_cols...)
  end
  df
end

function annotate_results(results, ref)
	joined = leftjoin(ref,
			  results,
			  on=["chr", "strand", "genomicPos"],
			  makeunique=true)
	joined[ismissing.(joined.predicted), "predicted"] .= 0
	joined[!, "predicted"] = disallowmissing(joined.predicted)
	joined[ismissing.(joined.GMM_LOR), "GMM_LOR"] .= 0
	joined[!, "GMM_LOR"] = disallowmissing(joined.GMM_LOR)
	joined[ismissing.(joined.predicted_raw), "predicted_raw"] .= 0
	joined[!, "predicted_raw"] = disallowmissing(joined.predicted_raw)
	if "GMM_LOR_raw" in names(joined)
		joined[ismissing.(joined.GMM_LOR_raw), "GMM_LOR_raw"] .= 0
		joined[!, "GMM_LOR_raw"] = disallowmissing(joined.GMM_LOR_raw)
	end
	joined
end

function annotate_binned_results(binned_results, binned_glori)
	joined = leftjoin(binned_results,
			  binned_glori,
			  on=[:chr, :strand, :bin])
	joined[ismissing.(joined[!, "modified"]), "modified"] .= 0
	joined[!, "modified"] = disallowmissing(joined.modified)
	sort(joined, [:chr, :strand, :bin])
end

function bin(df, bin_size = BIN_SIZE)
	df[!, "bin"] = div.(df.genomicPos, bin_size)
	# columns = [:modified => maximum => :modified]
	columns = []
	for col in names(df)
	  if occursin("pvalue", col) || occursin("qvalue", col)
	    push!(columns, col => minimum => col)
	  end
	end
	if "predicted" in names(df)
		push!(columns, :predicted => maximum => :predicted)
	end
	if "predicted_raw" in names(df)
		push!(columns, :predicted_raw => maximum => :predicted_raw)
	end
	if "mean_cov" in names(df)
		push!(columns, :mean_cov => mean => :mean_cov)
	end
	combine(groupby(df, ["chr", "strand", "bin"]), columns...)
end


function sharkfin(ax, df, col, lor_col; markersize=7, rasterize=false)
	mods = df[!, "modified"] .== 1
	# nomods = .! mods
	# scatter!(ax, abs.(df[nomods, lor_col]), df[nomods, col], color="#DDDDDD", alpha=1, markersize=markersize)
	# scatter!(ax, abs.(df[mods, lor_col]), df[mods, col], color="#2E2585", alpha=0.8, markersize=markersize)
	color = ifelse.(mods, "#2E2585", "#DDDDDD")
	score = -log10.(clamp.(coalesce.(df[:, col], 1), 1e-300, 1))
	scatter!(ax, abs.(df[:, lor_col]), score, color=color, alpha=0.8, markersize=markersize, rasterize=rasterize)
	ax.xlabel = "|LOR|"
	ax.ylabel = "-log₁₀(𝑃-value)"
end


function corr_plot(ax, df, xlabel, ylabel; markersize=7, rasterize=false)
	mods = df[!, "modified"] .== 1
	nomods = .! mods
	scores1 = -log10.(clamp.(coalesce.(df.GMM_chi2_pvalue, 1), 1e-300, 1))
	scores2 = -log10.(clamp.(coalesce.(df.GMM_chi2_pvalue_1, 1), 1e-300, 1))
	scatter!(ax,
			 scores1[nomods],
			 scores2[nomods],
			 color="#DDDDDD",
			 alpha=1,
			 markersize=markersize,
			 rasterize=rasterize)
	scatter!(ax,
			 scores1[mods],
			 scores2[mods],
			 color="#2E2585",
			 alpha=0.8,
			 markersize=markersize,
			 rasterize=rasterize)
	ax.xlabel = xlabel
	ax.ylabel = ylabel
end

function discard_low_sig_positions(df, prob=0.1)
	low_sig = df.predicted .== 0 .&& abs.(df.GMM_LOR) .< LOR_THRESHOLD
	indices = 1:size(df, 1)
	selected_low_prob = randsubseq(indices[low_sig], prob)
	mask = repeat([false], size(df, 1))
	mask[selected_low_prob] .= true
	high_sig = .! low_sig
	selected = high_sig .|| mask
	df[selected, :]
end

function Base.count(labels::AbstractArray{Int64}, pos_labels::Set{Int64})
	num_pos, num_neg = 0, 0
	for label in labels
		if label in pos_labels
			num_pos += 1
		else
			num_neg += 1
		end
	end
	num_pos, num_neg
end

"""
auprc(scores::AbstractArray{Float64}, classes::AbstractArray{Symbol}, pos_labels::Set{Symbol})

Computes the area under the Precision-Recall curve using a lower
trapezoidal estimator, which is more accurate for skewed datasets.
"""
function auprc(scores::AbstractArray{Float64}, classes::AbstractArray{Int64}, pos_labels::Set{Int64})
    num_scores = length(scores) + 1
	ordering = sortperm(scores, rev=true)
	labels = classes[ordering]
	num_pos, num_neg = count(labels, pos_labels)

	tn, fn, tp, fp = 0, 0, num_pos, num_neg

	p = Array{Float64}(undef, num_scores)
	r = Array{Float64}(undef, num_scores)
	p[num_scores] = tp/(tp+fp)
	r[num_scores] = tp/(tp+fn)
	auprc, prev_r = 0.0, r[num_scores]
	pmin, pmax = p[num_scores], p[num_scores]

	# traverse scores from lowest to highest
	for i in num_scores-1:-1:1
		dtn = labels[i] in pos_labels ? 0 : 1
		tn += dtn
		fn += 1-dtn
		tp = num_pos - fn
		fp = num_neg - tn
		p[i] = (tp+fp) == 0 ? 1-dtn : tp/(tp+fp)
		r[i] = tp/(tp+fn)

		# update max precision observed for current recall value
		if r[i] == prev_r
			pmax = p[i]
		else
			pmin = p[i] # min precision is always at recall switch
			auprc += (pmin + pmax)/2*(prev_r - r[i])
			prev_r = r[i]
			pmax = p[i]
		end
	end
	auprc, p, r
end

function auprc(scores::AbstractArray{Float64}, classes::AbstractArray{Bool}, pos_labels::Set{Int64})
	auprc(scores, Int.(classes), pos_labels)
end

function peaks(df, radius)
	all_positions = DataFrames.combine(groupby(df, :ref_id),
									   :pos => (p -> 1:maximum(p)) => :pos)
	all_positions = leftjoin(all_positions, df,
			 				 on = [:ref_id, :pos])
	all_positions = sort(all_positions, [:ref_id, :pos])

	peaks = DataFrames.combine(groupby(all_positions, :ref_id),
				   	 	       :predicted_raw => (p -> argmaxima(p, radius, strict = false)) => :pos)
	sort(leftjoin(peaks, df, on = [:ref_id, :pos]),
		 [:ref_id, :pos])
end


function assoc_type(r00, r01, r10, r11)
	if abs(r00 + r11) > abs(r01 + r10)
		ifelse(r00 + r11 > 0, :comod, :excl)
	else
		ifelse(r01 + r10 > 0, :excl, :comod)
	end
end

function place_labels_nonoverlapping(points, labels; offset=0.05, padding=0.01, width_factor = 0.05, height = 0.05, gene_offsets = Dict())
    n = length(points)
    label_positions = Vector{Point2f}(undef, n)

    # Initial placement: offset from each point
    for i in 1:n
		final_offset = sign(points[i][1]) * get(gene_offsets, labels[i], offset)
        label_positions[i] = points[i] .+ (final_offset, 0)# (rand()-0.5)/10)
    end

    # Function to compute label extents (simplified)
    function compute_extents(pos, label)
        # Approximate text width/height (adjust as needed)
        width = width_factor * length(label)
        Rect2f(pos[1] - width/2, pos[2] - height/2, width, height)
    end

	function intersects(a::Rect2f, b::Rect2f)
	    a_min_x = a.origin[1]
	    a_max_x = a.origin[1] + a.widths[1]
	    a_min_y = a.origin[2]
	    a_max_y = a.origin[2] + a.widths[2]
	
	    b_min_x = b.origin[1]
	    b_max_x = b.origin[1] + b.widths[1]
	    b_min_y = b.origin[2]
	    b_max_y = b.origin[2] + b.widths[2]
	
	    return !(a_max_x < b_min_x || a_min_x > b_max_x ||
	             a_max_y < b_min_y || a_min_y > b_max_y)
	end

    # Iteratively adjust positions to avoid overlap
    for _ in 1:30  # Max iterations
        for i in 1:n
            ext_i = compute_extents(label_positions[i], labels[i])
            for j in 1:n
                if i == j
                    continue
                end
                ext_j = compute_extents(label_positions[j], labels[j])
                if intersects(ext_i, ext_j)
                    # Move label i away from label j
                    dir = LinearAlgebra.normalize(label_positions[i] - label_positions[j])
                    label_positions[i] += dir * 0.02
                end
            end
        end
    end

    return label_positions
end

function plot_gene_comods(sig_comods, gtf, gene; axis = nothing, colormap = :seaborn_crest_gradient, colorrange = (0, 100), highclip = :black, min_w = 2, max_w = 8)
	pairs = sig_comods[occursin.(gene, sig_comods.reference), :]
	exons = gtf[gtf.feature .== "exon" .&& occursin.(gene, gtf.attribute), :]


	shift = minimum(exons.start)
	# minx = minimum(exons.start)
	minx = 0
	maxx = maximum(exons.end) - shift

	if isnothing(axis)
		fig = Figure(size = (1200, 600))
		ax = Axis(fig[1, 1])
	else
		ax = axis
	end

	ticks = round.(range(minx, maxx; length = 4); digits = 0)
	ax.xticks = (ticks, [format_with_commas(Int(t)) for t in ticks .+ shift])

	chr = pairs[1, :chr]
	# text!(ax, -(maxx-minx)/30, 0; text = chr, align = (:right, :center))
	ax.xlabel = chr
	ax.xlabelsize = LABELSIZE

	lines!(ax, [minx, maxx], [0, 0], color = :black, linewidth = 1)
	scatter!(ax, range(minx, maxx; length = 40), zeros(40);
			 marker = pairs[1, :strand] == "+" ? :rtriangle : :ltriangle,
			 color = :black)

	for exon in eachrow(exons)
		lines!(ax, [exon.start - shift, exon.end - shift], [0, 0], color = :black, linewidth = 15)
	end

	pairs[!, :neglogp] = -log10.(pairs.fisher_pvalue)
	# tx_data[‘neglogp’] = -np.log10(tx_data[‘p_value’])
    # min_w, max_w = 0.5, 6
    # w_norm = (pairs.neglogp .- minimum(pairs.neglogp)) ./
                 # (maximum(pairs.neglogp) - minimum(pairs.neglogp) + 1e-9)
	w_norm = pairs.phi
	pairs[!, :linewidth] = min_w .+ (max_w - min_w) .* w_norm

	for pair in eachrow(pairs)
		alpha = 0.25 + (1 - 0.25) * abs(pair.phi)
		left = min(pair.genomicPos1, pair.genomicPos2) - shift
		right = max(pair.genomicPos1, pair.genomicPos2) - shift
		distance = right - left
		radius = distance/2
		center = left + radius
		orientation = sign(pair.phi)
		arc!(ax, Point2f(center, 0), radius, 0, orientation*π,
			 color = pair.neglogp, #pair.phi,
			 linewidth = pair.linewidth,
			 alpha = alpha,
			 colormap = colormap,
			 colorrange = colorrange,#(minimum(pairs.phi), maximum(pairs.phi)),
			 highclip = highclip,
			 resolution = 10000)
		println(pair.phi)
	end

	hideydecorations!(ax)

	if isnothing(axis)
		fig
	else
		axis
	end
end

function plot_gene_comods_elliptical_arcs(sig_comods, gtf, gene; axis = nothing, colormap = :seaborn_crest_gradient, colorrange = (0, 100), highclip = :darkblue, min_w = 2, max_w = 8, xlimits = nothing, minheight=missing, maxheight=missing, dotsize=10, fig=nothing, legend=true, fontsize=14)
	pairs = copy(sig_comods[occursin.(gene, sig_comods.reference), :])
	exons = gtf[gtf.feature .== "exon" .&& occursin.(gene, gtf.attribute), :]
	cds = gtf[gtf.feature .== "CDS" .&& occursin.(gene, gtf.attribute), :]


	shift = minimum(exons.start)
	# minx = minimum(exons.start)
	minx = 0
	maxx = maximum(exons.end) - shift

	if isnothing(axis)
		fig = Figure(size = (1200, 600))
		ax = Axis(fig[1, 1])
	else
		ax = axis
	end

	ticks = round.(range(minx, maxx; length = 4); digits = 0)
	ax.xticks = (ticks, [format_with_commas(Int(t)) for t in ticks .+ shift])

	chr = pairs[1, :chr]
	# text!(ax, -(maxx-minx)/30, 0; text = chr, align = (:right, :center))
	ax.xlabel = chr
	# ax.xlabelsize = LABELSIZE

	lines!(ax, [minx, maxx], [0, 0], color = :black, linewidth = 1)
	scatter!(ax, range(minx, maxx; length = 40), zeros(40);
			 marker = pairs[1, :strand] == "+" ? :rtriangle : :ltriangle,
			 color = :black)

	for exon in eachrow(exons)
		lines!(ax, [exon.start - shift, exon.end - shift], [0, 0], color = :black, linewidth = 10)
	end
	for region in eachrow(cds)
		lines!(ax, [region.start - shift, region.end - shift], [0, 0], color = :black, linewidth = 20)
	end

	pairs[!, :neglogp] = -log10.(pairs.fisher_pvalue)
	# tx_data[‘neglogp’] = -np.log10(tx_data[‘p_value’])
    # min_w, max_w = 0.5, 6
    # w_norm = (pairs.neglogp .- minimum(pairs.neglogp)) ./
                 # (maximum(pairs.neglogp) - minimum(pairs.neglogp) + 1e-9)
	w_norm = pairs.phi
	pairs[!, :linewidth] = min_w .+ (max_w - min_w) .* w_norm

	pairs = renamemods(pairs)

	mods1 = pairs[:, [:genomicPos1, :mod1]]
	rename!(mods1, [:pos, :mod])
	mods2 = pairs[:, [:genomicPos2, :mod2]]
	rename!(mods2, [:pos, :mod])
	mods = vcat(mods1, mods2)
	mods[!, :y] .= 0

	# scatter!(ax, mods.pos, zeros(nrow(mods)), color=coalesce.(mods.mod, "missing"))

	previous_pos = chr == "+" ? 0 : Inf
	close_num = 0
	for (i, pair) in enumerate(eachrow(pairs))
		alpha = 0.25 + (1 - 0.25) * abs(pair.phi)
		left = min(pair.genomicPos1, pair.genomicPos2) - shift
		right = max(pair.genomicPos1, pair.genomicPos2) - shift
		distance = right - left
		radius = distance/2

		if !isnothing(xlimits) && (left < xlimits[1] || right > xlimits[2])
			continue
		end
		
		center = left + radius
		orientation = sign(pair.phi)

		h = ismissing(minheight) ? radius : max(minheight, radius)
		h = ismissing(maxheight) ? h : min(h, maxheight)
		current_pos = chr == "+" ? left : right
		if abs(previous_pos - current_pos) > (maxx - minx)/10
			previous_pos = current_pos
			close_num = 0
		else
			close_num += 1
		end
		h += close_num*coalesce(maxheight/10, 100)

		t = range(0, orientation*π, 100)
		x = radius * cos.(t)
		y = minimum(h) * sin.(t) .+ orientation .* 5
		lines!(ax, center.+x, y,
			   color = pair.neglogp, #pair.phi,
			   # linewidth = pair.linewidth,
			   linewidth = 2,
			   # alpha = alpha,
			   colormap = colormap,
			   colorrange = colorrange,#(minimum(pairs.phi), maximum(pairs.phi)),
			   highclip = highclip)
		
	end

	local dots = data(mods) *
		mapping(:pos => (p -> p - shift), :y, color=:mod => "Modification", marker=:mod => "Modification") *
		visual(Scatter; markersize=dotsize)
	local grid = draw!(ax, dots)
	if legend && !isnothing(fig)
		legend!(fig[1, 1], grid; tellheight=false, tellwidth=false, halign=:center, valign=:bottom, orientation=:horizontal, framevisible=false, labelsize=fontsize, titlesize=fontsize)
	end

	hideydecorations!(ax)

	if !isnothing(xlimits)
		xlims!(ax, xlimits)
	end

	if isnothing(axis)
		fig
	else
		axis
	end
end

function plot_isoforms_model!(ax, gene; transcripts = nothing, colors = Dict(), focus = nothing, rename = nothing, fontsize=12, lpadding=0)
	exons = gtf[occursin.(gene, gtf.attribute) .&& occursin.("transcript_type \"protein_coding\"", gtf.attribute) .&& gtf.feature .== "exon", :]
	cdses = gtf[occursin.(gene, gtf.attribute) .&& occursin.("transcript_type \"protein_coding\"", gtf.attribute) .&& gtf.feature .== "CDS", :]
	refs = exons.transcript_id |> unique |> collect |> sort
	if !isnothing(transcripts)
		refs = transcripts |> unique |> collect |> sort
	end
	exons = exons[in.(exons.transcript_id, Ref(refs)), :]
	leftmost = minimum(exons.start)
	rightmost = maximum(exons.end)
	distance = rightmost - leftmost
	# refs = filter(r -> occursin(gene, r), keys(gtf.attribute)) |> unique |> collect |> sort)
	N = length(refs)
	i = length(refs) - 1
	max_label_len = 0
	for ref in refs
		ref_exons = exons[exons.transcript_id .== ref, :]
		ref_cdses = cdses[cdses.transcript_id .== ref, :]
		color = get(colors, ref, :black)
		for exon in eachrow(ref_exons)
			lines!(ax, [exon.start, exon.end], [i, i], linewidth = 12, color = color)
		end
		for cds in eachrow(ref_cdses)
			lines!(ax, [cds.start, cds.end], [i, i], linewidth = 20, color = color)
		end
		lines!(ax, [minimum(ref_exons.start), maximum(ref_exons.end)], [i, i], linewidth = 1, color = color)
		label = isnothing(rename) ? ref : get(rename, ref, ref)
		text!(ax, leftmost - distance * 0.02, i; text = label, align = (:right, :center), fontsize=fontsize)
		max_label_len = max(max_label_len, length(label))

		scatter!(ax, range(leftmost, rightmost; length = 40), repeat([i], 40);
			 marker = exons[1, :strand] == "+" ? :rtriangle : :ltriangle,
			 color = color)
		
		i -= 1
	end

	if !isnothing(focus)
		(left, right) = focus
		poly!(Point2f[(left, N-1+0.45), (right, N-1+0.45), (right, -0.45), (left, -0.45)], strokecolor = :red, strokewidth = 2, color=("#F0F0F0", 0.25))
	end
	
	ylims!(ax, -0.5, length(refs) - 0.5)
	xlims!(ax, low=leftmost - distance * 0.02 - lpadding, high=rightmost)
	ax.xlabel = exons[1, :chr]
	# ax.xlabelsize = LABELSIZE
	hidespines!(ax)
	# hidexdecorations!(ax)
	hideydecorations!(ax)
end

function arc_plot_genomic2(gridpos, df, allmods; range=nothing, grange=nothing, colorrange=nothing, title=nothing, highlight=nothing, spinecolor = :black, ticks = 5, sequence=missing, mod_symbols=Dict(), mod_colors=Dict())
	df = copy(df)
	allmods = copy(allmods)
	minpos = min(minimum(df.genomicPos1), minimum(df.genomicPos2))
	maxpos = max(maximum(df.genomicPos1), maximum(df.genomicPos2))
	if grange !== nothing
		minpos, maxpos = grange
	end
	df[!, :genomicPos1] .-= minpos
	df[!, :genomicPos2] .-= minpos
	allmods[!, :genomicPos] .-= minpos
	positions = unique(union(df.genomicPos1, df.genomicPos2))

	tstart = Int64(round((maxpos - minpos)*0.1; digits=0))
	tend = Int64(round((maxpos - minpos)*0.9; digits=0))
	println(tstart, " ", tend, " ", minpos)
	xticks = (Base.range(tstart, tend, length=ticks),
			  map(t -> format_with_commas(Integer(round(t+minpos; digits = 0))), Base.range(tstart, tend, length=ticks)))

	
	# f = Figure(size = (1400, 600))
	ax = Axis(gridpos,
		  	  title=title,
			  xlabel="Chromosome position",
			  # xlabelsize = 18,
			  xticks=xticks,
			  yticksvisible = false,
			  yticklabelsvisible = false,
			  spinewidth = 2,
			  topspinecolor = spinecolor,
			  bottomspinecolor = spinecolor,
			  leftspinecolor = spinecolor,
			  rightspinecolor = spinecolor)
	

	if grange !== nothing
		xlims!(ax, (0, maxpos-minpos))
	end

	local mods = Dict()
	for row in eachrow(df)
		if row.genomicPos1 ∉ keys(mods) || mods[row.genomicPos1] == "?"
			mods[row.genomicPos1] = row.mod1 === missing ? "?" : row.mod1
		end
		if row.genomicPos2 ∉ keys(mods) || mods[row.genomicPos2] === missing
			mods[row.genomicPos2] = row.mod2 === missing ? "?" : row.mod2
		end
	end
	if range !== nothing
		allmods = allmods[between.(allmods.pos, range[1], range[2]), :]
		df = df[between.(df.pos1, range[1], range[2]) .&& between.(df.pos2, range[1], range[2]), :]
	end

	max_radius = maximum(abs.(df.genomicPos2 .- df.genomicPos1))/2
	arc_offset = log10(max_radius)

	if max_radius < 16
		ylims!(ax, (-16, 16))
	end

	hlines!(ax, [0, 0], [minimum(keys(mods)), maximum(keys(mods))], color = :lightgrey, alpha = 0.5, linewidth = 16)

	if colorrange == nothing
		colorrange = (-log10(maximum(df.pvalue)), -log10(minimum(df.pvalue)))
	end

	highlighted_gpos = []

	max_vis_radius = 0
	for (i, row) in enumerate(eachrow(df))
		radius = abs(row.genomicPos2 - row.genomicPos1)/2
		if between(row.genomicPos1, 0, maxpos-minpos) && between(row.genomicPos2, 0, maxpos-minpos)
			max_vis_radius = max(radius, max_vis_radius)
		end
		center = min(row.genomicPos1, row.genomicPos2) + radius
		angle = row.relation == :comod ? π : -π
		y = row.relation == :comod ? arc_offset : -arc_offset
		# Vector.(eachrow(df[1:i, [:pos1, :pos2]]))
		prev = df[1:i, :]
		if y > 0
			y = 2 # max(y, max_radius/16)
		end
		if y < 0
			y = -2 # min(y, -max_radius/16)
		end

		linewidth = 2
		if highlight != nothing && (row.pos1, row.pos2) == highlight
			linewidth = 5
			highlighted_gpos = [row.genomicPos1, row.genomicPos2] 
		end

		style = :solid
		arc!(ax, Point2f(center, y), radius, 0, angle,
			 color = -log10(row.pvalue),
			 colorrange = colorrange,
			 colormap = :seaborn_crest_gradient,
			 linewidth = linewidth)		
	end

	highlight_mask = [x in highlighted_gpos for x in keys(mods)]
	mod_labels = ifelse.(values(mods) .=== missing, "?", values(mods))

	if sequence !== missing
		for (pos, letter) in zip(grange[1]:grange[2], sequence)
			mod = get(mods, pos - minpos, nothing)
			if !isnothing(mod)
				println(mod_symbols, " ", mod)
				println(mod_colors, " ", mod)
				marker = get(mod_symbols, mod, :circle)
				color = get(mod_colors, mod, :black)
				scatter!(ax, [pos - minpos], [0], marker=marker, color=color, markersize=16)
			else
				text!(ax, [Point2f(pos - minpos, 0)]; text=string(letter), align=(:center, :center), fontsize=12)
			end
			
		end
	else
		text!(ax,
			  [Point2f(x, 0) for x in collect(keys(mods))[.! highlight_mask]];
			  text = mod_labels[.! highlight_mask],
			  align = (:center, :center),
			  fontsize=12)
		
		text!(ax,
			  [Point2f(x, 0) for x in collect(keys(mods))[highlight_mask]];
			  text = collect(values(mods))[highlight_mask],
			  font = :bold,
			  align = (:center, :center),
			  fontsize=14)
	end
	
	
	hideydecorations!(ax)
	hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
	
	ax, max_vis_radius
end

function renamemods(df)
	t = copy(df)
	local mapping = Dict(
		missing => "Unclassified",
		# "m6A" => "ₘ⁶A",
		# "m1A" => "ₘ¹A",
		# "m7G" => "ₘ⁷G",
		# "m5C" => "ₘ⁵C",
		"m6A" => "m⁶A",
		"m1A" => "m¹A",
		"m7G" => "m⁷G",
		"m5C" => "m⁵C",
		"Y" => "Ψ",
		"Nm" => "Nm"
	)
	t[!, :mod1] = map(m -> get!(mapping, m, m), t.mod1)
	t[!, :mod2] = map(m -> get!(mapping, m, m), t.mod2)
	t
end



function contingency(a, b)::Matrix{Integer}
    [sum(a .== 0 .&& b .== 0) sum(a .== 0 .&& b .== 1)
     sum(a .== 1 .&& b .== 0) sum(a .== 1 .&& b .== 1)]
end

function format_with_commas(x)
    s = string(Int(round(x; digits=2)))
	return replace(s, r"(\d)(?=(\d{3})+(?!\d))" => s"\1,")
end

function infer_mods_by_motifs(peak_df)
	matching_mods = []
	for row in eachrow(peak_df)
		seq = tx_ref[row.ref_id]
		mod = missing
		for (writer, motif) in writer_motifs
			area = seq[max(row.pos + 1 - 5, 1):min(row.pos + 1 + 5, length(seq))]
			if occursin(motif, area)
				mod = writer_mods[writer]
				break
			end
		end
		push!(matching_mods, mod)
	end
	matching_mods
end

function infer_mods_by_motifs(refs, positions)
	matching_mods = []
	for (ref, pos) in zip(refs, positions)
		seq = tx_ref[ref]
		mod = missing
		for (writer, motif) in writer_motifs
			area = seq[max(pos + 1 - 5, 1):min(pos + 1 + 5, length(seq))]
			if occursin(motif, area)
				mod = writer_mods[writer]
				break
			end
		end
		push!(matching_mods, mod)
	end
	matching_mods
end

function annotate_peaks(peaks, mod_ref, writer_motifs, writer_mods)
	local df = copy(peaks)
	df[:, :peak_id] = 1:nrow(df)
	df[:, :start] = max.(df.genomicPos .- 4, 0)
	# end is non-inclusive in BED files
	df[:, :end] = df.genomicPos .+ 5
	local ranges = df[:, [:peak_id, :chr, :strand]]
	ranges[:, :pos] = range.(df.start, df.end)
	local mod_overlaps = innerjoin(flatten(ranges, :pos), mod_ref,
			  					   on = [:chr, :strand, :pos])
	rename!(mod_overlaps, :pos => :mod_pos)
	annot_peaks = leftjoin(df, mod_overlaps,
			 			   on = [:peak_id, :chr, :strand])
	# We may have multiple reference mods overlapping the peak,
	# we select the one that is from the closest cell type
	function select_closest_mod(sources)
		argmax([
			if missing === s
				0
			elseif occursin("HEK293T", s)
				3
			elseif occursin("HEK293", s)
				2
			else
				1
			end for s in sources])
	end
	annot_peaks = combine(groupby(annot_peaks, :peak_id),
						  (df -> df[select_closest_mod(df.source), :]))


	motif_mods = infer_mods_by_motifs(annot_peaks)
	mod_category = map((annot_mod, motif_mod) -> if annot_mod !== missing && annot_mod === motif_mod
			annot_mod, :motif_and_annot
		elseif annot_mod !== missing
			annot_mod, :annot
		elseif motif_mod !== missing
			motif_mod, :motif
		else
			missing, missing
		end,
	   annot_peaks.mod, motif_mods)
	annot_peaks[!, :mod] = map(first, mod_category)
	annot_peaks[!, :category] = map(last, mod_category)
	
	annot_peaks
end

function high_conf_not_m6A(comods)
	not_DRACH1 = infer_mods_by_motifs(comods.reference, comods.pos1) .!== "m6A"
	not_DRACH2 = infer_mods_by_motifs(comods.reference, comods.pos2) .!== "m6A"

	glori_m6As = Set(zip(glori.chrom, glori.strand, glori.pos))

	not_annotated1 = map(r -> !any([(r.chr, r.strand, r.genomicPos1 + offset) in glori_m6As for offset in -5:5]),
						 eachrow(comods))
	not_annotated2 = map(r -> !any([(r.chr, r.strand, r.genomicPos2 + offset) in glori_m6As for offset in -5:5]),
						 eachrow(comods))

	tuple.(not_DRACH1 .&& not_annotated1, not_DRACH2 .&& not_annotated2)
end

function annotate_comods(annotated_res, comods)
	mods = copy(unique(annotated_res[:, [:ref_id, :pos, :mod]]))
	mods = Dict(zip(mods.ref_id, mods.pos) .=> mods.mod)
	df = copy(comods)
	df[!, :mod1] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos1 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))
	df[!, :mod2] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos2 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))


	# Annotate if we're quite certain that the mods are not m6A
	# (not in GLORI and not a DRACH)
	local not_m6A = high_conf_not_m6A(comods)
	df[!, :pos1_not_m6A] = first.(not_m6A)
	df[!, :pos2_not_m6A] = last.(not_m6A)
	
	df
end

function get_read_mods(dbfile, ref)
	local db = SQLite.DB(dbfile)

	local reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame

	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	local probs = reduce(hcat, Array.(reads.probs))
	probs = ifelse.(probs .== -1, missing, probs)
	transpose(probs ./= 100)
end

function between(x, s, e)
	s <= x < e
end


writer_motifs = Dict(
	# m6A
	# Note m6A motif is first, so we match it with highest priority.
	"METTL3" => r"[GAT][AG]AC[ATC]",

	# m1A
	"TRMT6-61A" => r"GTTC[AG]A",
	
	# pseudouridine
	"PUS4" => r"[AG][AG]TTC[AGCT]A", # https://doi.org/10.1038/s41587-021-00915-6
	"PUS7" => r"T[ACTG]TA[AG]",
	
	# m5C
	"NSUN6" => r"CTCC[AT]",
	
	# m3C
	"METTL8_1" => r"GCCTCC[AG]",
	"METTL8_2" => r"AGGC[AT]GG",
	
	# m7G
	"METTL1" => r"[AG]AGGT",
	
	# Nm
	"FTSJ3" => r"[ACTG]AGATC[AG][AG]", #[AG][AG][AG],
	
	# ac4C
	# "NAT10" =>  r"C[ACTG][ACTG]C"
)


writer_mods = Dict(
	"TRMT6-61A" => "m1A",
	"PUS4" => "Y",
	"PUS7" => "Y",
	"NSUN6" => "m5C",
	"METTL8_1" => "m3C",
	"METTL8_2" => "m3C",
	"METTL1" => "m7G",
	"FTSJ3" => "Nm",
	"NAT10" => "ac4C",
	"METTL3" => "m6A"
)

