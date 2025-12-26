include("lib.jl")
using FASTX
using AlgebraOfGraphics
using SQLite
import MultipleTesting
import LinearAlgebra
import HypothesisTests

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
	ax.xlabelsize = LABELSIZE

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
	ax.xlabelsize = LABELSIZE
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
			  xlabelsize = 18,
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
	# for row in eachrow(allmods)
	# 	if row.genomicPos ∉ keys(mods) || mods[row.genomicPos] === missing
	# 		mods[row.genomicPos] = "?"
	# 	end
	# end

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
	
	
	# Colorbar(f[1:2, 2], limits = colorrange, colormap = :viridis,
    # vertical = true, label = "-log₁₀(q-value)")

	# Legend(f[2, 1],
	# 	   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
	# 		LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	#        ["Positive", "Negative"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)

	hideydecorations!(ax)
	hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
	
	ax, max_vis_radius
end

function renamemods(df)
	t = copy(df)
	local mapping = Dict(
		missing => "Unknown",
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

mod_ref = CSV.read("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/combined_mod_ref.tsv", DataFrame, delim = '\t')

gtf = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/gencode.v41.annotation.gtf", delim = '\t', comment = "#", header = [:chr, :source, :feature, :start, :end, :score, :strand, :frame, :attribute]))
gtf[!, :transcript_id] = map(a -> replace(split(split(a, "; ")[2], " ")[2], "\"" => ""), gtf.attribute)

glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k.bed"))
# glori[!, "pos"] = glori.start .+ 2
glori[:, :pos] = glori.start
glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
glori[!, "modified"] .= 1

tx_ref = Dict()
FASTAReader(open("/work/mzdravkov/gencode.v41.transcripts.fa")) do reader
	for record in reader
		tx_ref[identifier(record)] = sequence(record)
	end
end

ivt = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_nanocompore_results.tsv",
				   "GMM_chi2_qvalue",
				   shift=4,
				   genomic_collapse=false,
				   LOR_threshold = 0.8)
ivt[!, :pos] .+= 4
ivt[!, :mod_ratio] = sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))

peaks_ivt = peaks(ivt, 4)
sig_peaks_ivt = peaks_ivt[peaks_ivt.predicted .!== missing .&& peaks_ivt.predicted .>= -log10(0.01), :]

an_sig_peaks_ivt = annotate_peaks(sig_peaks_ivt, mod_ref, writer_motifs, writer_mods)

comods = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/comod_50perc_contingency.tsv")) # this run always has fisher's test
comods[!, :pos1] .+= 4
comods[!, :pos2] .+= 4
comods[!, :chi2_pvalue] = comods.pvalue
comods[!, :pvalue] = ifelse.(comods.fisher_pvalue .!== missing,
							 comods.fisher_pvalue,
							 comods.pvalue)
comods[!, :pvalue] = clamp.(comods.pvalue, 1e-300, 1)
comods[!, :relation] = map(assoc_type, comods.residual00, comods.residual01, comods.residual10, comods.residual11)
comods[!, :phi] .*= ifelse.(comods.relation .== :comod, 1, -1)
comods[:, :qvalue] = MultipleTesting.adjust(comods.pvalue, MultipleTesting.BenjaminiHochberg())
comods = comods[abs.(comods.pos1 .- comods.pos2) .> 8, :]

# add genomic positions
comods = leftjoin(comods, ivt[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
				  on=[:reference => :ref_id, :pos1 => :pos])
rename!(comods, :genomicPos => :genomicPos1)
comods = leftjoin(comods, ivt[:, [:ref_id, :pos, :genomicPos]],
				  on=[:reference => :ref_id, :pos2 => :pos])
rename!(comods, :genomicPos => :genomicPos2)


sig_comods = comods[comods.qvalue .< 0.05, :]
an_sig_comods = annotate_comods(an_sig_peaks_ivt, sig_comods)


cancer_genes = DataFrame(CSV.File("/home/mzdravkov/cancerGeneList.tsv", delim = '\t', normalizenames = true))


LABELSIZE = 18
LABELPAD = 0


f = Figure(size=(1600, 1100))
trow = f[1, 1] = GridLayout()
brow = f[2, 1] = GridLayout()

sig_color = "#2E2585"
insig_color = "#DDDDDD"


df = copy(comods)#[1:100, :]
df[!, :significant] = df.qvalue .< 0.05
plt = data(df) *
	mapping(:phi, :pvalue => (p -> -log10(p)), color=:significant) *
	visual(Scatter; markersize=5, rasterize=2)
draw!(trow[1, 1],
	  plt,
	  scales(Color=(; palette=[insig_color, sig_color]));
	  axis=(; #title="All modification pairs",
				xlabel="← Mutually exclusive    𝛗    Co-occurring →        ",
				ylabel="-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
			    xlabelsize=LABELSIZE,
				ylabelsize=LABELSIZE,
				limits=((-1, 1), (-2, 310))))


# PANEL B

panel_b = brow[1:3, 1] = GridLayout()

df = copy(sig_comods)

df[!, :gene] = [x[6] for x in split.(df.reference, "|")]

df_cancer = innerjoin(df, cancer_genes, on=[:gene => :Hugo_Symbol])
println(size(df_cancer))

df[!, :log10_distance] = log10.(abs.(df.pos1 .- df.pos2))
df[!, :log10_distance_bin] = round.(df.log10_distance ./ 0.5; digits=0)
labels_per_bin = 5
labeled = combine(groupby(df, [:log10_distance_bin]),
				  df -> vcat(sort(df[df.phi .>= 0, :], :phi, by=abs, rev=true)[1:(df[1, :log10_distance_bin] == 2 ? 0 : labels_per_bin), :],
							 sort(df[df.phi .< 0, :], :phi, by=abs, rev=true)[1:(df[1, :log10_distance_bin] == 2 ? labels_per_bin : labels_per_bin), :]))
# srsf2 = df[df.gene .== "SRSF2", :]
# labeled = vcat(labeled, srsf2[argmax(srsf2.phi):argmax(srsf2.phi), :])

high_significance = df[df.pvalue .< 1e-100 .&&
					   (.! (df.log10_distance_bin .<= 2 .&&
						    df.phi .> 0 .&&
						    df.phi .< 0.6)), :]
labeled = unique(vcat(labeled, high_significance))

cancer_geneset = Set(cancer_genes.Hugo_Symbol)
labeled[!, :cancer] = [gene in cancer_geneset for gene in labeled.gene]



ax = Axis(panel_b[1, 1],
		  xlabel="← Mutually exclusive    𝛗    Co-occurring →                     ",
		  ylabel="Distance between associated modifications",
		  xlabelsize=LABELSIZE,
		  ylabelsize=LABELSIZE,
		  yticks=(1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
		  yminorticksvisible=true,
		  yminorgridvisible=true,
		  yminorticksize=2,
		  yminorticks=map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]))
plt = data(df) *
	mapping(:phi,
			(:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))),
			color=:pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)") *
	visual(Scatter; markersize=5, rasterize=2)
draw!(ax,
	  plt,
	  scales(Color=(; colormap=:seaborn_crest_gradient, colorrange=(0, 100), highclip=:black)))
Colorbar(panel_b[1, 2], limits=[0, 100], colormap=:seaborn_crest_gradient, #:viridis,
		 vertical=true, label="-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", labelsize=LABELSIZE, labelpadding=LABELPAD)
colsize!(panel_b, 1, Relative(8/9))

points = Point2f.(labeled.phi, labeled.log10_distance)
labels = labeled.gene
println(labeled)
label_positions = place_labels_nonoverlapping(points, labels; offset=0.1, gene_offsets=Dict("SRSF2" => 0.45))
for (p, lp, label, cancer) in zip(points, label_positions, labels, labeled.cancer)
	xalign = (lp[1] > p[1] ? :left : :right)
    text!(ax, label, position=lp, align=(xalign, :center), color=:black, fontsize=11, font=cancer ? :bold : :italic)
    lines!(ax, [p, lp], color=:gray, linewidth=0.5)
end
xlims!(ax, -1.2, 1.55)


# PANEL C

ax = Axis(trow[1, 2])

# Plot SRSF2
# plot_gene_comods(sig_comods, gtf, "ENST00000359995.10"; axis=ax, colormap=:seaborn_crest_gradient, highclip=:black)
# text!(ax, 1600, -60; text="SRSF2", align=(:center, :top))
# Colorbar(trow[1, 3], limits=[0, 100], colormap=:seaborn_crest_gradient, #:viridis,
# vertical=true, label="-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", labelsize=LABELSIZE, labelpadding=LABELPAD)

plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, an_sig_comods), gtf, "ENST00000234875.9"; minheight = 5, maxheight = 150, dotsize=15, xlimits=(1170, 1820), legend=true, colormap = :seaborn_crest_gradient, highclip=:darkblue, fig=trow[1, 2], axis = ax, fontsize=20) # RPL22
ylims!(ax, -150, 150)
text!(ax, 0, -30; text = "RPL22", align = (:left, :top))
text!(ax, 1720, 130; text = "↑ Co-occurring", align = (:left, :top))
text!(ax, 1720, -130; text = "↓ Mutually exclusive", align = (:left, :bottom))

inset_ax = Axis(trow[1, 2],
			    width=Relative(0.5),
				height=Relative(0.2),
				halign=0.05,
				valign=0.95)
hidedecorations!(inset_ax)
tx = "ENST00000234875.9"
name = "RPL22"
gene_start = minimum(gtf[gtf.transcript_id .== tx, :start])

plot_isoforms_model!(inset_ax, "RPL22"; transcripts = Set([tx]), colors = Dict([tx => :black]), rename = Dict([tx => name]), focus=(gene_start+1170, gene_start+1820), fontsize=16, lpadding=2000)

Colorbar(trow[1, 3], limits = [0, 100], colormap = :seaborn_crest_gradient, vertical = true, label = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", labelsize = LABELSIZE, labelpadding = LABELPAD)



# PLOT D


ax = Axis(brow[1, 2:4], xtickformat="{:,d}")
		
gpos1 = 76066564
gpos2 = 76066573
row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
tx1 = split(row1.reference, "|")[1]
tx2 = split(row2.reference, "|")[1]
name1 = split(row1.reference, "|")[5]
name2 = split(row2.reference, "|")[5]

blue = "#00aaff"

plot_isoforms_model!(ax, "MDH"; transcripts=Set([tx1, tx2]), colors=Dict([tx1 => blue, tx2 => :orange]), focus=(gpos1-100, gpos2+100), rename=Dict([tx1 => name1, tx2 => name2]), lpadding=1650)

tx1_mask = occursin.(tx1, an_sig_comods.reference)
tx2_mask = occursin.(tx2, an_sig_comods.reference)
colorrange = (0, 100)
seq = tx_ref[row1.reference][(row1.pos1-10):(row1.pos2+10)]
mod_symbols = Dict(["m⁵C" => :utriangle,
					"m⁶A" => :cross,
					"Unknown" => :circle,
					"Ψ" => :rect])
mod_colors = Dict(["m⁵C" => "#e69f00",
				   "m⁶A" => "#009e73",
				   "Unknown" => "#0072b2",
				   "Ψ" => "#cc79a7"])
ax1, max_radius1 = arc_plot_genomic2(brow[2, 2],
				  					 renamemods(an_sig_comods[tx1_mask, :]),
				  					 an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
				  					 range=(row1.pos1-2000, row1.pos2+2000),
				  					 grange=(gpos1-10, gpos2+10),
				  					 colorrange=colorrange,
				  					 title=name1,
				  					 highlight=(row1.pos1, row1.pos2),
				  					 spinecolor=blue,
				  					 ticks=3,
									 sequence=seq,
									 mod_symbols=mod_symbols,
									 mod_colors=mod_colors)
text!(ax1, 1, 13; text="↑ Co-occurring", align=(:left, :top))
text!(ax1, 1, -13; text="↓ Mutually exclusive", align=(:left, :bottom))
Legend(brow[2, 2],
    	   [MarkerElement(marker=mod_symbols["Unknown"], markersize=16, color=mod_colors["Unknown"]),
	 		MarkerElement(marker=mod_symbols["m⁵C"], markersize=16, color=mod_colors["m⁵C"])],
    	   ["Unknown", "m⁵C"],
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
ax2, max_radius2 = arc_plot_genomic2(brow[3, 2],
				  				     renamemods(an_sig_comods[tx2_mask, :]),
				  				     an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
				  				     range=(row2.pos1-2000, row2.pos2+2000),
				  				     grange=(gpos1-10, gpos2+10),
				  				     colorrange=colorrange,
				  				     title=name2,
				  				     highlight=(row2.pos1, row2.pos2),
				  				     spinecolor=:orange,
				  				     ticks=3,
									 sequence=seq,
									 mod_symbols=mod_symbols,
									 mod_colors=mod_colors)
text!(ax2, 1, 13; text="↑ Co-occurring", align=(:left, :top))
text!(ax2, 1, -13; text="↓ Mutually exclusive", align=(:left, :bottom))
Legend(brow[3, 2],
    	   [MarkerElement(marker=mod_symbols["Unknown"], markersize=16, color=mod_colors["Unknown"]),
	 		MarkerElement(marker=mod_symbols["m⁵C"], markersize=16, color=mod_colors["m⁵C"])],
    	   ["Unknown", "m⁵C"],
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
Colorbar(brow[2:3, 3], limits=colorrange, colormap=:seaborn_crest_gradient, #viridis
		 vertical=true, label="-log₁₀(𝑃-value)", labelsize=LABELSIZE, labelpadding=LABELPAD)


ref = row1.reference
probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
df1 = DataFrame(Tables.table(probs[:, [row1.pos1+1-4, row1.pos2+1-4]] .> 0.5))
dropmissing!(df1)
observed1 = contingency(df1.Column1, df1.Column2) .+ 1
res1 = HypothesisTests.ChisqTest(observed1)
data1 = map(r -> join(Int.(collect(r)), "-"), eachrow(df1))
counts1 = data1 |> countmap

ref = row2.reference
probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
df2 = DataFrame(Tables.table(probs[:, [row2.pos1+1-4, row2.pos2+1-4]] .> 0.5))
dropmissing!(df2)
observed2 = contingency(df2.Column1, df2.Column2) .+ 1
res2 = HypothesisTests.ChisqTest(observed2)
data2 = map(r -> join(Int.(collect(r)), "-"), eachrow(df2))
counts2 = data2 |> countmap



grid = brow[2:3, 4] = GridLayout()	
top = Axis(grid[1, 1],
		   xlabel="",
		   xticks=(1:4, ["", "", "", ""]),
		   ylabel="Pattern occurrences",
		   ylabelsize=LABELSIZE)
xlims!(top, 0.5, 4.5)
ylims!(top, 0, 2350)

df1 = DataFrame(tx=name1, pattern=data1)
df2 = DataFrame(tx=name2, pattern=data2)
df = vcat(df1, df2)
df = combine(groupby(df, [:tx, :pattern]), nrow => :count)
barplt = data(df) *
	mapping(:pattern, :count, dodge=:tx, color=:tx) *
	visual(BarPlot; width=0.7, dodge_gap=0.15)
draw!(top, barplt, scales(Color=(; palette=[blue, :orange])))

expected = DataFrame(tx=[name1, name1, name1, name1, name2, name2, name2, name2],
				  	 pattern=["0-0", "0-1", "1-0", "1-1",
							  "0-0", "0-1", "1-0", "1-1"],
				     count=[res1.expected[1, 1],
							res1.expected[1, 2],
							res1.expected[2, 1],
							res1.expected[2, 2],
							res2.expected[1, 1],
							res2.expected[1, 2],
							res2.expected[2, 1],
							res2.expected[2, 2]])
barplt = data(expected) *
	mapping(:pattern, :count, dodge=:tx) *
	visual(BarPlot;
		   width=0.7,
		   dodge_gap=0.15,
		   color=(:white, 0),
		   strokewidth=1.5,
		   strokecolor=:black)
draw!(top, barplt)

order = ["0-0", "0-1", "1-0", "1-1"]
bottom = Axis(grid[2, 1],
			  yticks=(1:2, map(format_with_commas, [gpos2, gpos1])))
xlims!(bottom, 0.5, 4.5)
ylims!(bottom, 0.5, 2.5)
scatter!(bottom,
		 [1, 1, 2, 2, 3, 3, 4, 4],
		 [1, 2, 1, 2, 1, 2, 1, 2],
		 color=[:white, :white,
				:white, :black,
				:black, :white,
				:black, :black],
		 strokecolor=:black,
		 strokewidth=1.5,
		 markersize=25)
hidexdecorations!(bottom)
hidespines!(bottom)

rowsize!(grid, 1, Relative(5/6))
rowsize!(grid, 2, Relative(1/6))

Legend(grid[1, 1],
	   [PolyElement(color=blue),
		PolyElement(color=:orange),
		PolyElement(color=(:white, 0), strokecolor=:black, strokewidth=1.5)],
	   [name1, name2, "Expected"],
	   tellheight=false,
	   tellwidth=false,
	   halign=:left,
	   valign=:top,
	   margin=(10, 10, 10, 10))


colsize!(trow, 1, Relative(6/18))

colsize!(brow, 1, Relative(9/18))
colsize!(brow, 2, Relative(5/18))
colsize!(brow, 3, Relative(0.25/18))
colsize!(brow, 4, Relative(3.75/18))

rowsize!(f.layout, 1, Relative(3/8))
rowsize!(f.layout, 2, Relative(5/8))

rowsize!(brow, 1, Relative(1/7))
rowsize!(brow, 2, Relative(3/7))
rowsize!(brow, 3, Relative(3/7))




for (label, layout) in zip(["A", "B", "C", "D", "E"],
						   [trow[1, 1], brow[1, 1], trow[1, 2], brow[1, 2], brow[2, 4]])
    Label(layout[1, 1, TopLeft()],
		  label,
		  fontsize=26,
		  font=:bold,
		  padding=(0, 5, 5, 0),
		  halign=:right)
end

save("fig5_comods.png", f)
save("fig5_comods.pdf", f)
