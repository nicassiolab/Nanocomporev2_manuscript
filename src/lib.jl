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
