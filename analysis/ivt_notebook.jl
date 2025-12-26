### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 41818e42-3cb8-4d48-93d9-b53e0eea7873
using Pkg

# ╔═╡ 6e8e2f5d-d998-4a4c-a14f-dc63e247b8a7
Pkg.add("Associations")

# ╔═╡ 82dba028-3fbe-44e2-9e8f-377f738363a7
Pkg.add("Peaks")

# ╔═╡ ab03ca5a-3186-493e-9e3e-13367a85e624
Pkg.add("FASTX")

# ╔═╡ 9fbe17fb-5b2a-4cda-a70d-1598de0a6ab0
Pkg.add("Clustering")

# ╔═╡ 692af519-1ac2-483c-86c5-f1801d8426bd
Pkg.add("AlgebraOfGraphics")

# ╔═╡ 9ba8c503-f59e-4d6a-b97c-9826279d744a
Pkg.add("GenomicAnnotations")

# ╔═╡ bcd1df63-e0b6-467b-9858-de5a15a61765
Pkg.add("StatsPlots")

# ╔═╡ f4100e9f-9456-4663-9328-485adc82baf2
Pkg.add("CategoricalArrays")

# ╔═╡ 705da010-749c-4841-9259-c323ce711464
Pkg.add("MultivariateStats")

# ╔═╡ a731b6b7-bc7e-4e92-a54d-d0ccc8ff678a
Pkg.add("UMAP")

# ╔═╡ 93db1a5d-4de2-482b-a4e2-540db2059689
begin
	using DataFrames
	using CSV
	using MLBase
	using CairoMakie
	using Random
	using StatsBase
	using Printf
	# using MultipleTesting
end

# ╔═╡ c88509a9-5a1d-4895-9ae7-5174cc199772
using Peaks

# ╔═╡ cefdac37-7b63-4c8c-ad92-41ca18b78452
using FASTX

# ╔═╡ dd730961-b1f9-4b2b-af69-84b2ac751117
using SQLite

# ╔═╡ 76d3d1c4-0e01-4360-bc3d-bc274ff1c39c
using AlgebraOfGraphics

# ╔═╡ 164f8adc-3bba-11f0-3c64-19ee3bf9097e
include("lib.jl")

# ╔═╡ e6362147-b298-4b52-a2bd-b38dee879b49
import Associations

# ╔═╡ bc77ac69-225a-4f81-b219-15683bb5eb4a
# ╠═╡ disabled = true
#=╠═╡
function simple_peaks(df, radius)
	N = size(df, 1)
	keep = repeat([true], N)
	for i in 1:N
		for j in max(1, i - radius):min(N, i + radius)
			if df[i, :ref_id] == df[j, :ref_id] &&
			   abs(df[i, :pos] - df[j, :pos]) <= radius &&
			   df[i, :predicted] > df[j, :predicted]
				keep[j] = false
			end
		end
	end
	return df[keep, :]
end
  ╠═╡ =#

# ╔═╡ dedfd9b6-9b9e-4149-834a-547cbac99873
# ╠═╡ disabled = true
#=╠═╡
# function peaks(df, radius)
# 	all_positions = DataFrames.combine(groupby(df, :ref_id),
# 									   :pos => (p -> minimum(p):maximum(p)) => :pos)
# 	all_positions = leftjoin(all_positions, df,
# 			 				 on = [:ref_id, :pos])
# 	all_positions = sort(all_positions, [:ref_id, :pos])
# 	# all_positions[:, :predicted] = coalesce.(all_positions.predicted, 0)
# 	peak_indices = argmaxima(all_positions.predicted, radius, strict = false)
# 	all_positions[peak_indices, :]
# end
  ╠═╡ =#

# ╔═╡ 1795bc42-bdc8-4685-bb52-39271581bc06
function peaks(df, radius)
	all_positions = DataFrames.combine(groupby(df, :ref_id),
									   :pos => (p -> 1:maximum(p)) => :pos)
	all_positions = leftjoin(all_positions, df,
			 				 on = [:ref_id, :pos])
	all_positions = sort(all_positions, [:ref_id, :pos])

	peaks = DataFrames.combine(groupby(all_positions, :ref_id),
				   	 	       :predicted_raw => (p -> argmaxima(p, radius, strict = false)) => :pos)
	sort(leftjoin(peaks, df,
			 	  on = [:ref_id, :pos]),
		 [:ref_id, :pos])
end

# ╔═╡ 5a2ad3e8-9f14-49ad-a78b-6fdd4be05220
# ╠═╡ disabled = true
#=╠═╡
function expand_pos(df, pos_col, left, right)
	N = left + right + 1
	expanded_df = repeat(df, N)
	pos_offsets = sort(repeat(-left:right, size(df, 1)))
	expanded_df[:, pos_col] += pos_offsets
	expanded_df
end
  ╠═╡ =#

# ╔═╡ 8445b562-cb03-4a7e-911b-5dc44778ee5c
# ╠═╡ disabled = true
#=╠═╡
begin
	glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/high_confidence_GLORI_sites.bed"))
	glori[!, "pos"] = glori.start .+ 2
	glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
	glori[!, "modified"] .= 1
	glori
end
  ╠═╡ =#

# ╔═╡ 8b930a4c-cef4-49e6-8a76-d59d4df8a349
	begin
	glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k.bed"))
	# glori[!, "pos"] = glori.start .+ 2
	glori[:, :pos] = glori.start
	glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
	glori[!, "modified"] .= 1
	glori
end

# ╔═╡ a6d02ee4-2e77-4520-b140-36fb2150b1b0
binned_glori = unique(glori[!, ["chrom", "strand", "bin", "modified"]])

# ╔═╡ b89d56f0-5104-4b0d-a1ae-6090d262db25
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
	"FTSJ3" => r"[ACTG]AGATC[AG][AG]", #[AG][AG][AG]",
	
	# ac4C
	# "NAT10" =>  r"C[ACTG][ACTG]C"
)

# ╔═╡ 566437e3-1abe-49bf-a262-e9a4c702432e
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

# ╔═╡ b5d8eafa-9ea4-4c59-bef8-3070126bb095
rmbase_psu = DataFrame(CSV.File("/home/mzdravkov/rmbase_mods/human.hg38.Pseudo.result.col29.bed",
				   header = ["chr", "start", "end", "mod_id", "score", "strand", "mod_type", "support_num", "support_list", "support_list_sub", "pub_list", "cell_list", "seq_type_list", "gene_id", "transcript_id", "gene_name", "gene_type", "region", "seq", "motif_score", "conserved_sites_list", "snoRNA_detail_info", "snoRNA_guide_site", "snoRNA_name_list", "snoRNA_database", "writer_loc", "writer_id", "writer_name_list", "source"]))

# ╔═╡ c9fabd63-2131-41fc-a50a-bda27319bc6d
begin
	atoi = DataFrame(CSV.File(
		"/home/mzdravkov/rmbase_mods/human.hg38.RNA-editing.result.col29.bed",
		delim='\t',
		header = vcat(["chr", "start", "end", "mod_id", "score", "strand", "mod_type"], repeat(["."], 22))))[:, 1:7]
	for col in ["support_num", "support_list", "support_list_sub", "pub_list", "cell_list"]
		atoi[:, col] .= "na"
	end
	atoi
end

# ╔═╡ 43a5d6b6-652d-403c-b4fd-945ba2b80ead
begin
	rmbase_all = copy(rmbase_psu[:, 1:12])
	for file in ["human.hg38.m1A.result.col29.bed",
				 "human.hg38.m5C.result.col29.bed",
				 "human.hg38.m7G.result.col29.bed",
				 "human.hg38.Nm.result.col29.bed",
				 "human.hg38.otherMod.result.col29.bed"]
				 # "human.hg38.RNA-editing.result.col29.bed"]
		local mods = DataFrame(CSV.read("/home/mzdravkov/rmbase_mods/$file",
										DataFrame,
										select = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
									    header = ["chr", "start", "end", "mod_id", "score", "strand", "mod_type", "support_num", "support_list", "support_list_sub", "pub_list", "cell_list"],
									    ntasks = 1))[:, 1:12]
		rmbase_all = vcat(rmbase_all, mods)
	end
	# rmbase_all = rmbase_all[occursin.("HEK293T", rmbase_all.cell_list), :]
	rmbase_all
end

# ╔═╡ 5abecd8c-e06c-477c-8938-545f64591385
hela_ac4C = DataFrame(CSV.File("/home/mzdravkov/GSE162043_RedaCTseq_ac4C_sites_GSE162043_RedaCTseq_ac4C_sites_hg38_lift.bed", header = [:chr, :start, :end, :name, :score, :strand]))

# ╔═╡ b0a1b932-7f0b-40ac-a2d6-e4867d2b54da
begin
	m7GDB_ngs = sort(DataFrame(CSV.File("/home/mzdravkov/base_human_filtered.csv")),
					 [:seqnames, :start])
	m7GDB_ngs.start .-= 1
end

# ╔═╡ 7ca021f1-532f-4a3c-98b6-0ff08c27d01e
m5C_atlas = DataFrame(CSV.File("/home/mzdravkov/hglift_m5C_to_grch38.bed",
							   header = [:chr, :start, :end, :name, :score, :strand]))

# ╔═╡ 33a21a05-dac1-40e0-97a4-8dae37fa8bfc
begin
	local glori_sel = glori[:, [:chrom, :strand, :pos]]
	rename!(glori_sel, [:chr, :strand, :pos])
	glori_sel[:, :mod] .= "m6A"
	glori_sel[:, :source] .= "HEK293T"
	
	local rmbase_sel = rmbase_all[:, [:chr, :strand, :start, :mod_type, :cell_list]]
	rename!(rmbase_sel, :mod_type => :mod, :start => :pos, :cell_list => :source)

	local atoi_sel = atoi[:, [:chr, :strand, :start, :mod_type, :cell_list]]
	rename!(atoi_sel, :start => :pos, :mod_type => :mod, :cell_list => :source)

	local m7g_sel = m7GDB_ngs[:, [:seqnames, :strand, :start, :source]]
	m7g_sel[:, :mod] .= "m7G"
	rename!(m7g_sel, :seqnames => :chr, :start => :pos)

	local m5c_sel = m5C_atlas[:, [:chr, :strand, :start, :name]]
	m5c_sel[:, :mod] .= "m5C"
	rename!(m5c_sel, :start => :pos, :name => :source)

	local ac4c_sel = hela_ac4C[:, [:chr, :strand, :start]]
	ac4c_sel[:, :mod] .= "ac4C"
	ac4c_sel[:, :source] .= "HeLa"
	rename!(ac4c_sel, :start => :pos)

	mod_ref = vcat(glori_sel, rmbase_sel, atoi_sel, m7g_sel, m5c_sel, ac4c_sel)

	mod_ref[:, :mod] = ifelse.(map(m -> m in Set(["Um", "Am", "Cm", "Gm"]),  mod_ref.mod), "Nm", mod_ref.mod)
	mod_ref
end

# ╔═╡ 9232e724-d8df-4970-931c-f1c11f0ce859
mod_ref |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/combined_mod_ref.tsv", delim = '\t')

# ╔═╡ 82a9954d-2c23-46bc-a0bb-c2fc078eabb1
begin
	local counts = DataFrame([(k, v) for (k, v) in countmap(mod_ref.mod)],
							 [:mod, :count])
	counts = counts[counts.count .> 5, :]
	sort!(counts, :count, rev = true)
	
	local f = Figure()
	local ax = Axis(f[1, 1],
				    xlabel = "Modification type",
				    ylabel = "Log count",
				    xticks = (1:nrow(counts), counts.mod),
				    xticklabelrotation = 45,
				    yscale = log10)
	
	barplot!(ax, 1:nrow(counts), counts.count)
	
	f
end

# ╔═╡ 57a4e1e9-59cb-4180-b735-9da13feb64a6
begin
	local df = copy(mod_ref)
	rename!(df, :pos => :start)
	df[:, :end] = df.start .+ 1
	df[:, :score] .= 0
	local mod_ref_bed = df[:, [:chr, :start, :end, :mod, :score, :strand, :source]]
	CSV.write("/home/mzdravkov/mod_ref.bed", mod_ref_bed, delim = '\t')
end

# ╔═╡ fbd57196-5c5b-49a9-a5f5-a39d25a33c87
begin
	ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
		LOR_threshold = 0.8)
	ivt[!, :pos] .+= 4
	ivt[!, :mod_ratio] = sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(ivt[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	ivt
end

# ╔═╡ d41911f6-ab79-4a7f-b245-76ef3cfcd1e6
begin
	peaks_ivt = peaks(ivt, 4)
	sig_peaks_ivt = peaks_ivt[peaks_ivt.predicted .!== missing .&&
						 	  peaks_ivt.predicted .>= -log10(0.01), :]
end

# ╔═╡ 5aaee16e-1121-4bf7-b4c6-60941b62414c
begin
	local df = copy(sig_peaks_ivt)
	df[:, :start] = max.(df.genomicPos .- 4, 0)
	# end is non-inclusive in BED files
	df[:, :end] = df.genomicPos .+ 5
	local peaks_bed = sort(df[:, [:chr, :start, :end, :ref_id, :predicted, :strand]],
		 				   [:chr, :start])
	CSV.write("/home/mzdravkov/ivt_sig_peaks.bed", peaks_bed, delim = '\t')
end

# ╔═╡ f3dba5a1-0352-4334-ab79-89b73172a289
sig_peaks_ivt.ref_id |> unique |> length

# ╔═╡ 9f605ad2-6d05-49d5-a743-2dbf960b1a51


# ╔═╡ fe4155bb-d187-4113-9e85-dfc3fedee4ad
mod_ref

# ╔═╡ a4d697a0-71b6-4e01-879e-6a91d9060b8d
function format_power(base, pow)
	superscript_digits = Dict(
	    '0' => '⁰',
	    '1' => '¹',
	    '2' => '²',
	    '3' => '³',
	    '4' => '⁴',
	    '5' => '⁵',
	    '6' => '⁶',
	    '7' => '⁷',
	    '8' => '⁸',
	    '9' => '⁹'
	)
	
	superscript_number = join([superscript_digits[digit]
							   for digit in string(pow)])
	"$(base)$(superscript_number)"
end

# ╔═╡ 16109b94-06e5-4285-9333-5b8a98761ef7
begin
	local df = copy(mod_ref)
	df[:, :category] = map(r -> if r.source !== missing &&
								occursin("HEK293T", r.source)
									1
								elseif r.source !== missing && occursin("HEK293", r.source)
									2
								else
							        3
								end,
								eachrow(df))
	local counts = DataFrames.combine(groupby(df, [:mod, :category]), nrow => :count)
	local totals = DataFrames.combine(groupby(counts, :mod), :count => sum => :count)
	local common_mods = totals[totals.count .>= 10, :mod] |> Set
	counts = counts[map(m -> m in common_mods, counts.mod), :]
	totals = Dict(zip(totals.mod, totals.count))
	counts[:, :total] = map(m -> totals[m], counts.mod)
	sort!(counts, :total, rev = true)
	local mods = unique(counts.mod)
	local mods_ids = Dict(collect(zip(mods, 1:length(mods))))
	local colormap = Dict(1 => :green,
						  2 => :orange,
						  3 => :lightgrey,
						  4 => :lightyellow)
	local f = Figure()
	local max_pow = Int(ceil(log10.(max(values(totals)...))))
	local ax = Axis(f[1, 1],
					title = "Modifications in the reference",
				    xlabel = "Modification type",
				    ylabel = "log₁₀(count)",
				    xticks = (1:length(unique(counts.mod)), unique(counts.mod)),
				    xticklabelrotation = 45,
				    # yscale = log10)
				    yticks = (0:max_pow, [format_power(10, i) for i in 0:max_pow]))
	barplot!(ax,
			 map(m -> mods_ids[m], counts.mod),
			 map(m -> log10(totals[m]), counts.mod) .* counts.count ./ map(m -> totals[m], counts.mod),
	         stack = counts.category,
	         color = map(c -> colormap[c], counts.category),
			 strokecolor = :black,
			 strokewidth = 0.25)

	Legend(f[2, 1],
	       [PolyElement(color = :green),
			PolyElement(color = :orange),
		    PolyElement(color = :lightgrey)],
	       ["HEK293T", "HEK293", "Other/Unspecified"],
		   orientation = :horizontal)
	
	f
end

# ╔═╡ 23ce8398-5f79-47d4-afe9-c46818011b9a
begin
	local df = innerjoin(mod_ref, ivt,
		  				 on = [:chr, :strand, :pos => :genomicPos],
		  				 makeunique = true)
	println(nrow(df))
	df[:, :category] = map(r -> if r.source !== missing &&
								occursin("HEK293T", r.source)
									1
								elseif r.source !== missing && occursin("HEK293", r.source)
									2
								else
							        3
								end,
								eachrow(df))
	local counts = DataFrames.combine(groupby(df, [:mod, :category]), nrow => :count)
	local totals = DataFrames.combine(groupby(counts, :mod), :count => sum => :count)
	local common_mods = totals[totals.count .>= 10, :mod] |> Set
	counts = counts[map(m -> m in common_mods, counts.mod), :]
	totals = Dict(zip(totals.mod, totals.count))
	counts[:, :total] = map(m -> totals[m], counts.mod)
	sort!(counts, :total, rev = true)
	local mods = unique(counts.mod)
	local mods_ids = Dict(collect(zip(mods, 1:length(mods))))
	local colormap = Dict(1 => :green,
						  2 => :orange,
						  3 => :lightgrey,
						  4 => :lightyellow)
	local f = Figure()
	local max_pow = Int(ceil(log10.(max(values(totals)...))))
	local ax = Axis(f[1, 1],
					title = "Reference modifications that are covered by Nanocompore",
				    xlabel = "Modification type",
				    ylabel = "log₁₀(count)",
				    xticks = (1:length(unique(counts.mod)), unique(counts.mod)),
				    xticklabelrotation = 45,
				    # yscale = log10)
				    yticks = (0:max_pow, [format_power(10, i) for i in 0:max_pow]))
	ylims!(ax, (0, max_pow + 0.5))
	local max_categories = DataFrames.combine(groupby(counts, :mod),
											  :category => maximum => :max_cat)
	local max_categories = Dict(zip(max_categories.mod, max_categories.max_cat))
	local bar_labels = ifelse.(counts.category .== map(m -> max_categories[m], counts.mod),
							   map(t -> "$(Int(t))", counts.total),
							   "")
	barplot!(ax,
			 map(m -> mods_ids[m], counts.mod),
			 map(m -> log10(totals[m]), counts.mod) .* counts.count ./ map(m -> totals[m], counts.mod),
	         stack = counts.category,
	         color = map(c -> colormap[c], counts.category),
			 bar_labels = bar_labels,
			 strokecolor = :black,
			 strokewidth = 0.25)

	Legend(f[2, 1],
	       [PolyElement(color = :green),
			PolyElement(color = :orange),
		    PolyElement(color = :lightgrey)],
	       ["HEK293T", "HEK293", "Other/Unspecified"],
		   orientation = :horizontal)
	
	f
end

# ╔═╡ 6669794a-cebf-46f7-8d1b-5bb91ccb7280
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT/WT MA plot - peaks",
				    xlabel = "Total coverage in both conditions",
					yticks = -9:1:9,
				    ylabel = "LOR",
				    xscale = log10)
	xlims!(ax, (10, 12000))
	ylims!(ax, (-9, 9))
	local df = an_peaks_ivt
	scatter!(ax, df.mean_cov, df.GMM_LOR,
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ c71344c3-4faa-4849-a35a-6f511eaf4f34
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT/WT qq-plot",
				    xlabel = "Expected p-value (-log₁₀-scale)",
				    ylabel = "Observed p-value (-log₁₀-scale)")
	local df = peaks_ivt[peaks_ivt.GMM_chi2_pvalue .!== missing, :]
	local expected = -log10.((1:nrow(df))./nrow(df))
	local observed = -log10.(sort(df.GMM_chi2_pvalue .+ 1e-310))
	
	qqplot!(ax,
			expected,
			observed,
		    qqline = :identity,
		    markercolor = :blue,
			markersize = 5,
		    color = :red)
	
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	
	local median_observed_chisq = median(quantile.(Distributions.Chisq(n), df.GMM_chi2_pvalue))
	local expected_chisq = quantile(Distributions.Chisq(n), 0.5)
	local λ = median_observed_chisq/expected_chisq
	println("Lambda gc: ", median_observed_chisq/expected_chisq)

	text!(ax, 0.1, 280; text = "λgc = $(round(λ; digits = 4))")
	
	f
end
  ╠═╡ =#

# ╔═╡ 86179fb3-d3d0-4f5b-a8d3-14fb43dac834
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
				    xlabel = "Coverage",
				    ylabel = "# positions",
				    # xscale = log10)
				    xticks = 0:2000:10000)
	hist!(ax, ivt.mean_cov, bins = 100)
	vlines!(ax, [200], color = :red, linestyle = :dash, linewidth = 1)
	f
end

# ╔═╡ 5fde98ae-bdf2-4dcb-90a0-f5017b97262a


# ╔═╡ d90664a5-8e7d-469b-8b05-b103d829f672
# ╠═╡ disabled = true
#=╠═╡
tx_ref["ENST00000001008.6|ENSG00000004478.8|OTTHUMG00000090429.4|OTTHUMT00000206861.3|FKBP4-201|FKBP4|3715|protein_coding|"][(1763 + 1 + 4):(1763 + 1 + 8 + 4)]
  ╠═╡ =#

# ╔═╡ 140e4f5f-3c19-4c94-8829-4a46a95eb6cd
begin
	tx_ref = Dict()
	FASTAReader(open("/work/mzdravkov/gencode.v41.transcripts.fa")) do reader
	   for record in reader
		   tx_ref[identifier(record)] = sequence(record)
	   end
	end
end

# ╔═╡ 74f8fb7a-e501-490a-a530-224f51d6312c
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

# ╔═╡ 7a783418-ca3a-4ac7-a8a1-8968a5bc9096
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

# ╔═╡ 88dac2fb-7062-4e75-aaaf-c2596d90b299
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
	
	# for (writer, motif) in writer_motifs
	# 	# matched = annot_peaks.mod .=== missing .&&
	# 	# 		  occursin.(motif, annot_peaks.ref_kmer)
	# 	# annot_peaks[matched, :mod] .= writer_mods[writer]
	# 	# annot_peaks[matched, :source] .= "motif"
	# 	matched = occursin.(motif, annot_peaks.ref_kmer)
		
	# 	mod_category = map((match, mod) -> if match && mod === missing
	# 		writer_mods[writer], :motif
	# 	elseif match && mod === writer_mods[writer]
	# 		mod, :motif_and_annot
	# 	elseif mod !== missing
	# 		mod, :annot
	# 	else
	# 		missing, missing
	# 	end,
	# 	   matched,
	# 	   annot_peaks.mod)
	# 	annot_peaks[!, :mod] = map(first, mod_category)
	# 	annot_peaks[!, :category] = map(last, mod_category)
	# end
	annot_peaks
end

# ╔═╡ 4016b836-3669-4098-9725-d3708755db40
an_peaks_ivt = annotate_peaks(peaks_ivt, mod_ref, writer_motifs, writer_mods)

# ╔═╡ 25f512fd-e5f8-40d7-a7eb-2d04e18b884b
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT/WT peaks",
				    xlabel = "|LOR| > 0.5",
				    ylabel = "-log₁₀(q-value)",
				    xticks = 0:0.5:9,
				    yticks = 0:25:300)
	xlims!(ax, (0, 9))
	local df = an_peaks_ivt[an_peaks_ivt.GMM_chi2_pvalue .!== missing .&& an_peaks_ivt.GMM_chi2_qvalue .!= 1, :]
	scatter!(ax,
			 abs.(df.GMM_LOR),
			 -log10.(df.GMM_chi2_qvalue),
			 color = map(r -> if r.mod === missing
				 				:grey
							  elseif occursin("motif", r.source)
								:orange
							  else
								:red
							  end, eachrow(df)),
			 markersize = 5)
	vlines!(ax, [0.8], linestyle = :dash, color = :black, linewidth = 0.8)
	hlines!(ax, [2], linestyle = :dash, color = :black, linewidth = 0.8)
	Legend(f[2, 1],
	       [PolyElement(color = :grey),
			PolyElement(color = :orange),
			PolyElement(color = :red)],
	       ["No reference overlaps", "Motif match", "Overlaps reference mods"],
		   orientation = :horizontal)
	f
end

# ╔═╡ 390466c1-9b31-460b-82cd-676d85ca41da
begin
	an_sig_peaks_ivt = annotate_peaks(sig_peaks_ivt, mod_ref, writer_motifs, writer_mods)
	
	# local modkit_mods = Dict(zip(wt1mk.ref_id, wt1mk.pos) .=> wt1mk.mod)
	# local mods = map(r -> if r.mod !== missing
	# 		r.mod
	# 	else
	# 		m = get!(modkit_mods, (r.ref_id, r.pos), "")
	# 		Dict("a" => "m6A", "m" => "m5C", "17802" => "Y", "17596" => "A-I", "" => missing)[m]
	# 	end,
	# 	eachrow(an_sig_peaks_ivt))
	# # an_sig_peaks_ivt[!, :mod] = ifelse.(mods .== "", missing, mods)
	# an_sig_peaks_ivt[!, :mod] = mods
	# an_sig_peaks_ivt
end

# ╔═╡ 9539bce6-7042-4396-b4f0-1badf268b068
(an_sig_peaks_ivt.mod .!== missing) |> countmap

# ╔═╡ 081adcae-8bfe-4c6f-9a17-1ddd8e9dd87e
begin
	local bed = an_sig_peaks_ivt[:, [:chr, :genomicPos, :mod, :predicted, :strand]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, :mod] = coalesce.(bed.mod, ".")
	bed[:, [:chr, :genomicPos, :end, :mod, :predicted, :strand]] |> CSV.write("/home/mzdravkov/ivt_mods.bed", delim = '\t', header = false)
end

# ╔═╡ 3a91b67e-8efa-4e84-a519-6426594ea30e
begin
	local bed = an_sig_peaks_ivt[:, [:chr, :genomicPos, :predicted]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, [:chr, :genomicPos, :end, :predicted]] |> CSV.write("/home/mzdravkov/ivt_mods.bedGraph", delim = '\t', header = false)
end

# ╔═╡ 16fff5ab-601e-42d7-9f74-1670690733bb
begin
	local g = DataFrames.combine(
		groupby(an_sig_peaks_ivt, [:peak_id, :mod]),
		:source => (
		  source -> min(map(s -> if occursin("HEK293T", s)
				  					1
								 elseif occursin("HEK293", s)
									2
		  						 elseif occursin("motif", s)
									4
								 else
									3
								 end,
							ifelse.(source .=== missing, "", source))...)
		  # s -> any(occursin.("HEK293T",
				#  			 ifelse.(s .=== missing, "", s)))
		) => :category)
	local counts = DataFrames.combine(groupby(g, [:mod, :category]),
									  nrow => :count)
	counts[:, :mod] = ifelse.(counts.mod .=== missing,
							  "Not matched",
							  counts.mod)
	local totals = DataFrames.combine(groupby(counts, :mod), :count => sum => :count)
	local totals = Dict(zip(totals.mod, totals.count))
	counts[:, :total] = map(m -> totals[m], counts.mod)
	sort!(counts, :total, rev = true)
	local mods = unique(counts.mod)
	local mods_ids = Dict(collect(zip(mods, 1:length(mods))))
	local f = Figure()
	local max_pow = Int(ceil(log10.(max(values(totals)...))))
	local colormap = Dict(1 => :green,
						  2 => :orange,
						  3 => :lightgrey,
						  4 => :lightyellow,
						  5 => :black)
	local ax = Axis(f[1, 1],
					title = "Putative modifications in significant peaks",
				    xlabel = "Modification type",
				    ylabel = "log₁₀(count)",
				    xticks = (1:length(unique(counts.mod)), unique(counts.mod)),
				    xticklabelrotation = 45,
				    # yscale = log10)
				    yticks = (0:max_pow, [format_power(10, i) for i in 0:max_pow]))
	ylims!(ax, (0, max_pow + 1))
	local max_categories = DataFrames.combine(groupby(counts, :mod),
											  :category => maximum => :max_cat)
	local max_categories = Dict(zip(max_categories.mod, max_categories.max_cat))
	local bar_labels = ifelse.(counts.category .== map(m -> max_categories[m], counts.mod),
							   map(t -> "$(Int(t))", counts.total),
							   "")
	counts[counts.mod .== "Not matched", :category] .= 5
	barplot!(ax,
			 map(m -> mods_ids[m], counts.mod),
			 map(m -> log10(totals[m]), counts.mod) .* counts.count ./ map(m -> totals[m], counts.mod),
	         stack = counts.category,
	         color = map(c -> colormap[c], counts.category),
			 bar_labels = bar_labels,
			 strokecolor = :black,
			 strokewidth = 0.25)

	Legend(f[2, 1],
	       [PolyElement(color = :green),
			PolyElement(color = :orange),
		    PolyElement(color = :lightgrey),
		    PolyElement(color = :lightyellow)],
	       ["HEK293T", "HEK293", "Other/Unspecified", "Motif"],
		   "Attribution source",
		   orientation = :horizontal)
	# local f = Figure()
	# local ax = Axis(f[1, 1],
	# 			    title = "Detected modifications (WT / IVT)",
	# 			    xticks = (1:length(mods), unique(counts.mod)),
	# 			    yticks = (0:500:12000, [if mod(n, 1000) == 0 "$n" else "" end
	# 										for n in 0:500:12000]))
	# 			    # ytickformat = "{:d}")
	# local colormap = Dict(true => :green,
	# 					  # 2 => :orange,
	# 					  false => :lightgrey,
	# 					  4 => :lightyellow)
	# barplot!(ax,
	# 		 map(m -> mods_ids[m], counts.mod),
	# 		 counts.count,
	#          stack = counts.category,
	#          color = map(c -> colormap[c], counts.category),
	# 		 strokecolor = :black,
	# 		 strokewidth = 0.25)
	# Legend(f[2, 1],
	# 	   [PolyElement(color = :green),
	# 		PolyElement(color = :orange),
	# 		PolyElement(color = :lightgrey),
	# 		PolyElement(color = :lightyellow)],
	#        ["HEK293T", "HEK293", "Unspecified/Other cell line", "Motif"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	f
end

# ╔═╡ 14c823fa-6412-4f38-88c2-66dfa4bc4842
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

# ╔═╡ 4cb861eb-f76b-4871-a783-554c85642e42
begin
	local df = copy(an_sig_peaks_ivt)

	local covered_ref = innerjoin(ivt[:, [:chr, :strand, :genomicPos]], mod_ref,
		  						  on = [:chr, :strand, :genomicPos => :pos])

	local modifications = []
	local total_mods = []
	local perc_by_annot = []
	local perc_by_motif = []
	local perc_by_both = []
	local total_annot = []
	local total_motifs = []
	for mod in unique(df.mod)
		push!(modifications, mod)
		subset = df[df.mod .=== mod, :]
		total = nrow(subset)
		push!(total_mods, total)
		category_counts = countmap(subset.category)
		by_annot = get(category_counts, :annot, 0)/ total
		push!(perc_by_annot, by_annot)
		by_motif = get(category_counts, :motif, 0) / total
		push!(perc_by_motif, by_motif)
		by_both = get(category_counts, :motif_and_annot, 0) / total
		push!(perc_by_both, by_both)
		push!(total_annot, sum(covered_ref.mod .== mod))
		motif_occurences = 0
		for writer in keys(filter(p -> p[2] === mod, writer_mods))
			if !haskey(writer_motifs, writer)
				continue
			end
			motif = writer_motifs[writer]
			# motif_occurences += sum(occursin.(motif, ivt.ref_kmer))
			refs = combine(groupby(ivt, :ref_id),
						   :pos => minimum => :minpos,
					   	   :pos => maximum => :maxpos)
			for row in eachrow(refs)
				seq = tx_ref[row.ref_id]
				area = seq[row.minpos+1:row.maxpos+1]
				motif_occurences += length(collect(eachmatch(motif, area)))
			end
		end
		push!(total_motifs, motif_occurences)
	end
	
	sort(DataFrame(modification = modifications, total = total_mods, perc_annot = perc_by_annot, perc_motif = perc_by_motif, perc_both = perc_by_both, annotated = total_annot, motif_occurences = total_motifs), :total, rev = true)
end

# ╔═╡ a3f63803-5803-4942-8dce-b1b90b86f6a1
function between(x, s, e)
	s <= x < e
end

# ╔═╡ 9513c89b-dd9e-40dd-91be-c480438e3f07
function overlaps(mod_ref, row)
	matches = mod_ref[mod_ref.chr .== row.chr .&&
		    		  mod_ref.strand .== row.strand .&&
		    		  between.(mod_ref.pos, row.start, row.end), :]
	Dict(zip(matches.mod, matches.pos))
end

# ╔═╡ 56ed77ee-fcfb-498a-ae5c-f9e77ceb971a
begin
	local df = copy(mod_ref)
	df = df[map(m -> m in Set(["m6A", "m5C", "Nm", "m7G", "m1A", "ac4C", "A-I", "Y", "m3C"]), df.mod), :]
	df[:, :source] = map(r -> if occursin("HEK293T", r.source)
							   		 "HEK293T"
								   elseif occursin("HEK293", r.source)
									 "HEK293"
								   else
								   	 "Other"
								   end, eachrow(df))
		
	local table = combine(
			groupby(df, :mod),
			nrow => :count,
		    :source => (s -> 100*sum(s .== "HEK293T")/length(s)) => :HEK293T,
		    :source => (s -> 100*sum(s .== "HEK293")/length(s)) => :HEK293,
		    :source => (s -> 100*sum(s .== "Other")/length(s)) => :Other)

	local covered_ref_mods = combine(
		groupby(innerjoin(mod_ref, ivt,
						  on = [:chr, :strand, :pos => :genomicPos],
						  makeunique = true),
				:mod),
		nrow => :ivt_covered)

	table = leftjoin(table, covered_ref_mods, on = :mod)

	local mod_peaks = combine(groupby(an_sig_peaks_ivt, :mod),
							  nrow => :peaks)

	table = outerjoin(table, mod_peaks,
					  on = :mod,
			 		  matchmissing = :equal)
	sort!(table, :count)
	ivt_mod_table = table
end

# ╔═╡ e429bac0-0add-4e15-8a7a-6ab22a5c9943
begin
	local ref = copy(mod_ref[:, [:mod, :source]])
	ref[:, :source] = map(
		s -> if occursin("HEK293T", s)
				"HEK293T"
			 elseif occursin("HEK293", s)
				"HEK293"
			 else
				"Other"
			 end,
		mod_ref.source)
	local ref_counts = combine(groupby(ref, [:mod, :source]), nrow => :count)


	local covered_ref = innerjoin(mod_ref, ivt,
						  on = [:chr, :strand, :pos => :genomicPos],
						  makeunique = true)[:, [:mod, :source]]
	covered_ref[:, :source] = map(
		s -> if occursin("HEK293T", s)
				"HEK293T"
			 elseif occursin("HEK293", s)
				"HEK293"
			 else
				"Other"
			 end,
		covered_ref.source)
	local covered_ref_counts = combine(groupby(covered_ref, [:mod, :source]),
									   nrow => :count)

	local an_sig_peaks_ivt_cp = copy(an_sig_peaks_ivt[:, [:mod, :source]])
	an_sig_peaks_ivt_cp[:, :source] = map(
		s -> if s === missing
			    "Unclassified"
			 elseif s === "motif"
				"Motif"
			 elseif occursin("HEK293T", s)
				"HEK293T"
			 elseif occursin("HEK293", s)
				"HEK293"
			 else
				"Other"
			 end,
		an_sig_peaks_ivt_cp.source)
	local detected_counts = combine(groupby(an_sig_peaks_ivt_cp[an_sig_peaks_ivt_cp.mod .!== missing, :], [:mod, :source]),
								    nrow => :count)

	local colors = Dict(
		"HEK293T" => :green,
		"HEK293" => :orange,
		"Motif" => :lightyellow,
		"Other" => :lightgrey,
		"Unclassified" => :black
	)
	
	local unclassified = ivt_mod_table[ivt_mod_table.mod .=== missing, :peaks][1]
	local table = ivt_mod_table[ivt_mod_table.mod .!== missing, :]
	
	local f = Figure(size = (800, 600))
	local ax = Axis(f[1, 1])
				    # aspect = 1)
	# xlims!(ax, (6, 44))
	local rows = 0:10:(10+10*nrow(table))
	for (h, x) in zip(["Modification", "Reference", "IVT covered", "In peaks"], [0, 25, 55, 85])
		if h == "Modification"
			text!(ax, x, (10*nrow(table));
				  text = h,
				  font = :bold,
				  align = (:left, :center),
				  fontsize = 16)
		else
			text!(ax, x, (10*nrow(table));
				  text = h,
				  font = :bold,
				  align = (:left, :center),
				  fontsize = 16)
		end
	end
	for (mod, y) in zip(table.mod, rows)
		text!(ax, 0, y;
			  text = coalesce(mod, "Unclassified"),
			  align = (:left, :center),
			  fontsize = 15)
		hlines!(ax, [y - 5], color = :black, linewidth = 0.5)
		hlines!(ax, [y + 5], color = :black, linewidth = 0.5)
	end
	for (mod, c, y) in zip(table.mod, table.count, rows)
		text!(ax, 25, y;
			  text = ifelse(c === missing, "", string(c)),
			  align = (:right, :center),
			  fontsize = 15)
		rc = ref_counts[ref_counts.mod .=== mod, :]
		pie!(ax, 33, y, rc.count,
			 color = map(s -> colors[s], rc.source),
		 	 radius = 4)
	end
	for (mod, c, y) in zip(table.mod, table.ivt_covered, rows)
		text!(ax, 55, y;
			  text = ifelse(c === missing, "", string(c)),
			  align = (:right, :center),
			  fontsize = 15)
		rc = covered_ref_counts[covered_ref_counts.mod .=== mod, :]
		pie!(ax, 63, y, rc.count,
			 color = map(s -> colors[s], rc.source),
		 	 radius = 4)
	end
	for (mod, c, y) in zip(table.mod, table.peaks, rows)
		text!(ax, 85, y;
			  text = ifelse(c === missing, "", string(c)),
			  align = (:right, :center),
			  fontsize = 15)
		rc = detected_counts[detected_counts.mod .=== mod, :]
		pie!(ax, 93, y, rc.count,
			 color = map(s -> colors[s], rc.source),
		 	 radius = 4)
	end

	text!(ax, 0, -12,
		  text = "Total",
		  fontsize = 15,
		  align = (:left, :center))
	text!(ax, 25, -12,
		  text = string(sum(ref_counts.count)),
		  fontsize = 15,
		  align = (:right, :center))
	text!(ax, 55, -12,
		  text = string(sum(covered_ref_counts.count)),
		  fontsize = 15,
		  align = (:right, :center))
	text!(ax, 85, -12,
		  text = string(sum(detected_counts.count)),
		  fontsize = 15,
		  align = (:right, :center))
	
	text!(ax, 0, -22,
		  text = "Unclassified peaks",
		  fontsize = 15,
		  align = (:left, :center))
	text!(ax, 85, -22,
		  text = string(unclassified),
		  fontsize = 15,
		  align = (:right, :center))
	
	Legend(f[2, 1],
		   [PolyElement(color = :green),
			PolyElement(color = :orange),
			PolyElement(color = :lightgrey),
			PolyElement(color = :lightyellow)],
	       ["HEK293T", "HEK293", "Unspecified/Other cell line", "Motif"],
		   "Attribution source",
		   orientation = :horizontal,
		   framevisible = false)
	
	hidedecorations!(ax)
	hidespines!(ax)
	f
end

# ╔═╡ 2b618c14-b478-4d5b-abfc-19bbeb073fde
md"""
## STORM/WT
"""

# ╔═╡ ea86f47c-ed99-4108-bb5c-c9ee225dc7ef
begin
	storm = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_fix_mod_clust_inferring/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
	    LOR_threshold = 0.8)
	storm[!, :pos] .+= 4
end

# ╔═╡ 7fcfeb9f-5e35-4408-8cfb-b24c5f52c6a2
begin
	peaks_storm = peaks(storm, 4)
	sig_peaks_storm = peaks_storm[peaks_storm.predicted .!== missing .&&
						  	      peaks_storm.predicted .>= -log10(0.01), :]
end

# ╔═╡ eb2f9c7c-b5bb-493d-8cc1-abdc452e646d
an_peaks_storm = annotate_peaks(peaks_storm, mod_ref, writer_motifs, writer_mods)

# ╔═╡ ed261b2a-459a-440e-99ad-91afdeda2952
an_sig_peaks_storm = annotate_peaks(sig_peaks_storm, mod_ref, writer_motifs, writer_mods)

# ╔═╡ dc85be9a-f0af-4b8b-924e-758b3ce45a7a
begin
	local bed = an_sig_peaks_storm[:, [:chr, :genomicPos, :mod, :predicted, :strand]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, :mod] = coalesce.(bed.mod, ".")
	bed[:, [:chr, :genomicPos, :end, :mod, :predicted, :strand]] |> CSV.write("/home/mzdravkov/storm_mods.bed", delim = '\t', header = false)
end

# ╔═╡ bfe35e43-60fd-4ba9-9f05-f5a7f868e011
begin
	local bed = an_sig_peaks_storm[:, [:chr, :genomicPos, :predicted]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, [:chr, :genomicPos, :end, :predicted]] |> CSV.write("/home/mzdravkov/storm_mods.bedGraph", delim = '\t', header = false)
end

# ╔═╡ a81dc14a-fd54-417b-bcd7-84615a62aa37
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "STORM/WT peaks",
				    xlabel = "|LOR| > 0.5",
				    ylabel = "-log₁₀(q-value)",
					xticks = 0:0.5:9,
				    yticks = 0:25:300)
	xlims!(ax, (0, 9))
	local df = an_peaks_storm[an_peaks_storm.GMM_chi2_pvalue .!== missing .&& an_peaks_storm.GMM_chi2_qvalue .!= 1, :]
	scatter!(ax,
			 abs.(df.GMM_LOR),
			 -log10.(df.GMM_chi2_qvalue),
			 color = map(r -> if r.mod === missing
				 				:grey
							  elseif occursin("motif", r.source)
								:orange
							  else
								:red
							  end, eachrow(df)),
			 markersize = 5)
	vlines!(ax, [0.8], linestyle = :dash, color = :black, linewidth = 0.8)
	hlines!(ax, [2], linestyle = :dash, color = :black, linewidth = 0.8)
	Legend(f[2, 1],
	       [PolyElement(color = :grey),
			PolyElement(color = :orange),
			PolyElement(color = :red)],
	       ["No reference overlaps", "Motif match", "Overlaps reference mods"],
		   orientation = :horizontal)
	f
end
  ╠═╡ =#

# ╔═╡ a2dbec3c-961f-4904-961d-be2d9d538720
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "STORM/WT MA plot - peaks",
				    xlabel = "Total coverage in both conditions",
					yticks = -9:1:9,
				    ylabel = "LOR",
				    xscale = log10)
	xlims!(ax, (10, 12000))
	ylims!(ax, (-9, 9))
	local df = an_peaks_storm
	scatter!(ax, df.mean_cov, df.GMM_LOR,
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ 7c1bce27-3c08-4b0a-9364-9e992a793a97
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "STORM/WT qq-plot",
				    xlabel = "Expected p-value (-log₁₀-scale)",
				    ylabel = "Observed p-value (-log₁₀-scale)")
	local df = peaks_storm[peaks_storm.GMM_chi2_pvalue .!== missing, :]
	local expected = -log10.((1:nrow(df))./nrow(df))
	local observed = -log10.(sort(df.GMM_chi2_pvalue .+ 1e-310))
	
	qqplot!(ax,
			expected,
			observed,
		    qqline = :identity,
		    markercolor = :blue,
			markersize = 5,
		    color = :red)
	
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	
	local median_observed_chisq = median(quantile.(Distributions.Chisq(n), df.GMM_chi2_pvalue))
	local expected_chisq = quantile(Distributions.Chisq(n), 0.5)
	local λ = median_observed_chisq/expected_chisq
	println("Lambda gc: ", median_observed_chisq/expected_chisq)

	text!(ax, 0.1, 280; text = "λgc = $(round(λ; digits = 4))")
	
	f
end
  ╠═╡ =#

# ╔═╡ bf514f88-b2f4-4a35-8c4a-08fec382e74c
md"""
## Single molecule IVT/WT
"""

# ╔═╡ c1002495-1c24-450c-9585-1aeb7d588bb5
import MultipleTesting

# ╔═╡ 823b1914-ab79-499a-98ec-348f24b0107a
begin
	# comod = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/single_molecule_comodification2.tsv"))
	comod = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/single_molecule_comodifications.tsv"))
	comod[!, :pos1] .+= 4
	comod[!, :pos2] .+= 4
	comod = comod[.! isnan.(comod.pvalue) .&& comod.expected .!= 0 .&& comod.observed .>= 5, :]
	comod[:, :qvalue] = MultipleTesting.adjust(comod.pvalue, MultipleTesting.BenjaminiHochberg())
	comod = comod[abs.(comod.pos1 .- comod.pos2) .> 8, :]
	comod
end

# ╔═╡ 80836ca5-66c5-4f9c-a815-9c57230f6948
sig_comod = comod[comod.qvalue .<= 0.01, :]

# ╔═╡ 66983a7b-8a6a-4b9f-91ba-bb40b5b5ebce
combine(groupby(sig_comod, [:reference, :pos1, :pos2]), nrow => :count)

# ╔═╡ 237e9e39-95c8-40a1-b8d3-306a7b98ce0f
function relation_type(r)
	s = (r.state1, r.state2)
	if (s in [(1, 1), (0, 0)] && r.observed > r.expected) || (s in [(0, 1), (1, 0)] && r.observed < r.expected)
		:positive
	elseif (s in [(1, 1), (0, 0)] && r.observed < r.expected) || (s in [(0, 1), (1, 0)] && r.observed > r.expected)
		:negative
	else
		:weird
	end
end

# ╔═╡ cf152b80-2eea-4ed4-8d7f-d0d30a1aec2f
begin
	local df = copy(sig_comod)
	df[!, :relation] = map(relation_type, eachrow(df))
	combine(groupby(df, [:reference, :pos1, :pos2]),
		    :relation => length∘Set => :relations).relations |> countmap
end

# ╔═╡ 6c7f7c71-ed77-46b1-97a7-0d8b9c8d758e
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure(size = (1200, 800))
	local ax1 = Axis(f[1, 1],
					 title = "Co-unmodified (0-0)",
					 xlabel = "log₂(Obs/Exp)",
					 ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)")
	local ax2 = Axis(f[1, 2],
					 title = "Mutually exclusive (0-1/1-0)",
					 xlabel = "log₂(Obs/Exp)",
					 ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)")
	local ax3 = Axis(f[1, 3],
					 title = "Co-modified (1-1)",
					 xlabel = "log₂(Obs/Exp)",
					 ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)")

	local s00 = comod.state1 .== 0 .&& comod.state2 .== 0
	local s01 = comod.state1 .== 0 .&& comod.state2 .== 1
	local s10 = comod.state1 .== 1 .&& comod.state2 .== 0
	local s11 = comod.state1 .== 1 .&& comod.state2 .== 1

	local effsize = log2.(comod.observed ./ comod.expected)
	
	scatter!(ax1,
			 effsize[s00],
			 -log10.(comod[s00, :pvalue]),
			 color = ifelse.(comod[s00, :qvalue] .<= 0.01, :red, :grey),
			 markersize = 5)
	
	local selected = s00 .&& -log10.(comod.qvalue) .> 60
	text!(ax1,
		  [Point2f(x, y) for (x, y) in zip(effsize[selected], -log10.(comod.pvalue[selected]))];
		  text = map(r -> split(r, "|")[5], comod.reference[selected]),
		  align = (:left, :bottom),
		  fontsize = 9)

	println("Significant 0-0: ", nrow(comod[s00 .&& comod.qvalue .<= 0.01, :]))

	scatter!(ax2,
			 effsize[s01 .|| s10],
			 -log10.(comod[s01 .|| s10, :pvalue]),
			 color = ifelse.(comod[s01 .|| s10, :qvalue] .<= 0.01, :red, :grey),
			 markersize = 5)
	local selected = (s01 .|| s10) .&& -log10.(comod.qvalue) .> 110
	local offsets = Dict("SERBP1-203" => (0.3, 3),
						 "YBX3-201" => (0, 2),
						 "TUBA1B-201" => (0.01, -4),
						 "PGAM1-201" => (0.01, -4),
						 "HSPD1-202" => (-0.04, 3),
						 "EEF2-201" => (0.08, -2))
	local offset = m -> m in keys(offsets) ? offsets[m] : (0, 0)
	text!(ax2,
		  [Point2f(x + 0.01 + offset(m)[1],
				   y + 1 + offset(m)[2])
		   for (x, y, m) in zip(effsize[selected],
							    -log10.(comod.pvalue[selected]),
							    map(r -> String(split(r, "|")[5]), comod.reference[selected]))];
		  text = map(r -> split(r, "|")[5], comod.reference[selected]),
		  align = (:left, :bottom),
		  fontsize = 9)
	println("Significant 0-1: ", nrow(comod[s01 .&& comod.qvalue .<= 0.01, :]))
	println("Significant 1-0: ", nrow(comod[s10 .&& comod.qvalue .<= 0.01, :]))

	scatter!(ax3,
			 effsize[s11],
			 -log10.(comod[s11, :pvalue]),
			 color = ifelse.(comod[s11, :qvalue] .<= 0.01, :red, :grey),
			 markersize = 5)
	local selected = s11 .&& -log10.(comod.qvalue) .> 75
	text!(ax3,
		  [Point2f(x, y) for (x, y) in zip(effsize[selected], -log10.(comod.pvalue[selected]))];
		  text = map(r -> split(r, "|")[5], comod.reference[selected]),
		  align = (:left, :bottom),
		  fontsize = 9)

	println("Significant 1-1: ", nrow(comod[s11 .&& comod.qvalue .<= 0.01, :]))

	linkaxes!(ax1, ax2, ax3)

	Legend(f[2, 1:3],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ 29c659f1-4e5a-4c65-aa53-33dea7d730f3
sig_comod[-log10.(sig_comod.qvalue) .> 120, :]

# ╔═╡ d3fc7047-56a6-4dce-96d3-61189556299e
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Mod. associations qq-plot",
				    xlabel = "Expected p-value (-log₁₀-scale)",
				    ylabel = "Observed p-value (-log₁₀-scale)")
	local df = comod
	local expected = -log10.((1:nrow(df))./nrow(df))
	local observed = -log10.(sort(df.pvalue))
	
	qqplot!(ax,
			expected,
			observed,
		    qqline = :identity,
		    markercolor = :blue,
			markersize = 5,
		    color = :red)
	
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	
	local median_observed_chisq = median(quantile.(Distributions.Chisq(n), df.pvalue))
	local expected_chisq = quantile(Distributions.Chisq(n), 0.5)
	local λ = median_observed_chisq/expected_chisq
	println("Lambda gc: ", median_observed_chisq/expected_chisq)

	text!(ax, 0.1, 280; text = "λgc = $(round(λ; digits = 4))")
	
	f
end
  ╠═╡ =#

# ╔═╡ 1d820491-8699-4324-b113-b4d7ac0a819a
begin
	local f = Figure(size = (1200, 800))
	local ax = Axis(f[1, 1],
					title = "Distance between significant associated sites",
				    xlabel = "Distance (nt)",
				    ylabel = "Count",
				    xticks = 0:50:1400)

	df = sig_comod[:, [:pos1, :pos2]] |> unique

	hist!(ax, abs.(df.pos2 .- df.pos1), bins=nrow(df), color = :black)
	f
end

# ╔═╡ da1c3021-eccf-40bd-8988-4a0d068eb0ba
# ╠═╡ disabled = true
#=╠═╡
an_ivt = annotate_peaks(ivt, mod_ref, writer_motifs, writer_mods)
  ╠═╡ =#

# ╔═╡ bc42854a-ecbe-46c2-b289-114b12cf6fe9
#=╠═╡
begin
	# local p1 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	# rename!(p1, :mod => :mod1)
	# local p2 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	# rename!(p2, :mod => :mod2)

	# an_sig_comod = leftjoin(leftjoin(sig_comod, p1,
	# 	 	 		        		 on = [:reference => :ref_id, :pos1 => :pos]),
	# 		     		    p2,
	# 		 				on = [:reference => :ref_id, :pos2 => :pos])

	local mods = copy(unique(an_ivt[:, [:ref_id, :pos, :mod]]))
	mods = Dict(zip(mods.ref_id, mods.pos) .=> mods.mod)

	local df = copy(sig_comod)
	df[!, :mod1] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos1 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))
	df[!, :mod2] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos2 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))
	an_sig_comod = df
end
  ╠═╡ =#

# ╔═╡ 3de4cb79-b2c6-4c12-9e8c-8a78a3af4457
#=╠═╡
begin
	local df = copy(an_sig_comod)
	df[!, :relation] = map(r -> let s = (r.state1, r.state2)
									if (s in [(1, 1), (0, 0)] && r.observed > r.expected) || (s in [(0, 1), (1, 0)] && r.observed < r.expected)
										:positive
									elseif (s in [(1, 1), (0, 0)] && r.observed < r.expected) || (s in [(0, 1), (1, 0)] && r.observed > r.expected)
										:negative
									else
										:weird
									end
								end,
						   eachrow(df))
	df[!, "mods"] = map(r -> join(sort([coalesce(r.mod1, ""), coalesce(r.mod2, "")]), "-"),
						 eachrow(df))
	local selected = map(r -> r.mods in ["m6A-m7G", "m5C-m6A", "Y-m6A", "Y-m5C"], eachrow(df))
	sort(combine(groupby(df[selected, :], [:mods, :relation]), nrow => :count))
end
  ╠═╡ =#

# ╔═╡ be74ab35-49fe-401e-9829-6dbc8d90f37f
#=╠═╡
begin
	local df = copy(an_sig_comod)
	df[!, :relation] = map(r -> let s = (r.state1, r.state2)
									if (s in [(1, 1), (0, 0)] && r.observed > r.expected) || (s in [(0, 1), (1, 0)] && r.observed < r.expected)
										:positive
									elseif (s in [(1, 1), (0, 0)] && r.observed < r.expected) || (s in [(0, 1), (1, 0)] && r.observed > r.expected)
										:negative
									else
										:weird
									end
								end,
						   eachrow(df))
	df[!, "mods"] = map(r -> join(sort([coalesce(r.mod1, ""), coalesce(r.mod2, "")]), "-"),
						 eachrow(df))
	local selected = map(r -> r.mods in ["m6A-m7G", "m5C-m6A", "Y-m6A", "Y-m5C"], eachrow(df))

	println(sort(combine(groupby(df, :mods), nrow => :count)))

	
	# sort(combine(groupby(df[selected, :], [:mods, :relation]), nrow => :count))

	data(df[selected, :]) *
		mapping(:mods, color=:relation, dodge=:relation) *
		frequency() |> draw
end
  ╠═╡ =#

# ╔═╡ 97711651-3a07-4c6b-b749-572dddb4a145
#=╠═╡
begin
	local df = copy(an_sig_comod)
	df[!, :relation] = map(r -> let s = (r.state1, r.state2)
									if (s in [(1, 1), (0, 0)] && r.observed > r.expected) || (s in [(0, 1), (1, 0)] && r.observed < r.expected)
										:positive
									elseif (s in [(1, 1), (0, 0)] && r.observed < r.expected) || (s in [(0, 1), (1, 0)] && r.observed > r.expected)
										:negative
									else
										:weird
									end
								end,
						   eachrow(df))
	df[!, "mods"] = map(r -> join(sort([coalesce(r.mod1, ""), coalesce(r.mod2, "")]), "-"),
						 eachrow(df))
	local selected = map(r -> r.mods in ["m6A-m7G", "m5C-m6A", "Y-m6A", "Y-m5C"], eachrow(df))
	df[!, :distance] = abs.(df.pos1 .- df.pos2)
	
	# sort(combine(groupby(df[selected, :], [:mods, :relation]), nrow => :count))

	data(df[selected, :]) *
		mapping(:mods => "Modifications pair",
				:distance => "Distance",
				color=:relation => "Association",
				dodge=:relation) *
		visual(BoxPlot) |> draw
end
  ╠═╡ =#

# ╔═╡ d8feca5f-d4e8-4e78-8d87-4f0f7a4755c4
#=╠═╡
sort(an_sig_comod, :qvalue)
  ╠═╡ =#

# ╔═╡ ca2f63ea-677c-42e7-a736-b12a280cad22
#=╠═╡
an_sig_comod_gx[occursin.("TUBA1B", an_sig_comod.reference), :]
  ╠═╡ =#

# ╔═╡ f28e2414-efe5-4a73-91e3-79ed437d34e5
#=╠═╡
begin
	local df = copy(an_sig_comod_gx)
	df[!, :relation] = map(relation_type, eachrow(df))
	df = combine(groupby(df, [:reference, :pos1, :pos2, :relation, :chr, :strand, :genomicPos1, :genomicPos2, :mod1, :mod2]),
			:pvalue => minimum => :pvalue,
			:qvalue => minimum => :qvalue)
	# sort(df[occursin.("TUBA1B", df.reference), :])
	simple_comods = sort(df)
end
  ╠═╡ =#

# ╔═╡ cf8c009c-a109-4121-8e39-b3e96cbd4092
#=╠═╡
begin
	local df = leftjoin(an_sig_comod, ivt[:, [:ref_id, :pos, :chr, :strand, :genomicPos]], on=[:reference => :ref_id, :pos1 => :pos])
	rename!(df, :genomicPos => :genomicPos1)
	df = leftjoin(df, ivt[:, [:ref_id, :pos, :genomicPos]], on=[:reference => :ref_id, :pos2 => :pos])
	rename!(df, :genomicPos => :genomicPos2)
	an_sig_comod_gx = df
	# sort(df[occursin.("TUBA1B", df.reference), :])
end
  ╠═╡ =#

# ╔═╡ a7fe57e3-cfcc-4973-810e-98edacaf6fdb
nrow(sig_peaks_ivt)

# ╔═╡ 39ad0f2b-299e-4b35-9ecf-8cc83b55cd3f
an_sig_peaks_ivt.mod |> countmap

# ╔═╡ e113b1dc-a755-438a-8c08-6747cde5afde
sig_peaks_ivt.ref_id |> unique |> length

# ╔═╡ 414f6a39-e6c8-4dbe-9f35-41f83a065040
filter(r -> r.count > 1, eachrow(combine(groupby(sig_peaks_ivt, :ref_id), nrow => :count))) |> nrow

# ╔═╡ 98ecad14-b670-4ba6-bb55-8adbef018a95
#=╠═╡
simple_comods.reference |> unique |> length
  ╠═╡ =#

# ╔═╡ 435b9d3a-1eb7-4755-8778-05ce6f7def8c


# ╔═╡ 1aa2a4e9-fb40-4573-982b-cb6a0a8c2b59
#=╠═╡
vcat(collect(zip(simple_comods.reference, simple_comods.pos1)),
	 collect(zip(simple_comods.reference, simple_comods.pos2))) |> unique |> length
  ╠═╡ =#

# ╔═╡ f203a977-1521-460d-969d-796e5e4fecd2
#=╠═╡
begin
	local df = copy(simple_comods[simple_comods.mod1 .!== missing .&& simple_comods.mod2 .!== missing, :])
	df[!, :same_mod] = df.mod1 .== df.mod2
	
	sort(combine(groupby(df, [:same_mod, :relation]), nrow => :count))
end
  ╠═╡ =#

# ╔═╡ 2825eb3b-12a2-44ec-9146-06d95381dcbe
#=╠═╡
simple_comods[occursin.("TUBA1B", simple_comods.reference), :]
  ╠═╡ =#

# ╔═╡ 87144656-e38e-4e91-aa08-33a9c7b790e8
#=╠═╡
arc_plot_genomic(simple_comods[occursin.("ENST00000332858.10", simple_comods.reference), :], an_sig_peaks_ivt[occursin.("ENST00000332858.10", an_sig_peaks_ivt.ref_id), :]; range=(1000, 3200), grange=(49127780, 49128020))
  ╠═╡ =#

# ╔═╡ f064b8c4-d0ad-4e00-8296-4bbed1316256
#=╠═╡
arc_plot_genomic(simple_comods[occursin.("ENST00000336023.9", simple_comods.reference), :], an_sig_peaks_ivt[occursin.("ENST00000336023.9", an_sig_peaks_ivt.ref_id), :]; range=(1300, 1700), grange=(49127780, 49128020))
  ╠═╡ =#

# ╔═╡ 134194d3-c29f-4eb3-878d-a3c24892dd6c
function save_bedpe(comods, reference, filename)
	local df = copy(comods[comods.reference .== reference, :])
	# df[!, :chr2] = df.chr
	# df[!, :strand] = ifelse.(df.strand .== "+", 0, 1)
	# df[!, :strand2] = df.strand
	# df[!, :frag1] .= 0
	# df[!, :frag2] .= 1
	# df = rename!(df, :strand => :strand1, :chr => :chr1, :genomicPos1 => :pos1, :genomicPos2 => :pos2)
	# df[:, [:chr1, :pos1,  :frag1, :strand2, :chr2, :pos2, :frag2]] |> CSV.write("/home/mzdravkov/$filename", delim='\t', header=["# chrom1", :start1, :stop1, :chrom2, :start2, :stop2, :name, :qual, :strand1, :strand2, :filters, :info])
	# nothing
	header = "# chrom1 start1 stop1 chrom2 start2 stop2 name qual strand1 strand2 filters info"
	lines = map(r -> join([r.chr, r.genomicPos1, r.genomicPos1, r.chr, r.genomicPos2, r.genomicPos2, split(r.reference, "|")[5], 0, r.strand, r.strand, ".", r.relation], '\t'), eachrow(df))
	open(filename, "w") do f
		write(f, join([header, lines...], "\n"))
	end
	nothing
end

# ╔═╡ d2c86575-c50f-4069-a48b-cced8084a84b
ivt[ivt.ref_id .== "ENST00000332858.10|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409006.1|TUBA1B-201|TUBA1B|3198|retained_intron|", :]

# ╔═╡ 70c7949f-915c-4ded-96fd-bf158d9a8109
#=╠═╡
begin
	save_bedpe(simple_comods, "ENST00000332858.10|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409006.1|TUBA1B-201|TUBA1B|3198|retained_intron|", "/home/mzdravkov/ENST00000332858.10.bedpe")
	save_bedpe(simple_comods, "ENST00000336023.9|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409005.1|TUBA1B-202|TUBA1B|1627|protein_coding|", "/home/mzdravkov/ENST00000336023.9.bedpe")
end
  ╠═╡ =#

# ╔═╡ 02599f81-c0c5-453a-a278-4dfa7104a001
function get_read_mods(dbfile, ref)
	local db = SQLite.DB(dbfile)

	local reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame

	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	local probs = reduce(hcat, Array.(reads.probs))
	probs = ifelse.(probs .== -1, missing, probs)
	transpose(probs ./= 100)
end

# ╔═╡ 5542defd-747d-462e-abf2-6894ddb54081
function get_read_mods_with_ids(dbfile, ref)
	local db = SQLite.DB(dbfile)

	local reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame

	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	local probs = reduce(hcat, Array.(reads.probs))
	probs = ifelse.(probs .== -1, missing, probs)
	reads.read, transpose(probs ./= 100)
end

# ╔═╡ ac85532f-1b1d-4571-a187-e4304882ba27
begin
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000332858.10|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409006.1|TUBA1B-201|TUBA1B|3198|retained_intron|")
	local f = Figure(size = (600, 600))
	local ax = Axis(f[1, 1],
				    aspect = 0.3,
				    xticks=(1:2, ["Y", "m5C"]),
				    ylabel="Reads")

	local df = DataFrame(Tables.table(probs[:, [3155, 3172]] .> 0.75))
	dropmissing!(df)
	heatmap!(ax, transpose(Array(sort(df))), colormap=[:white, :black])
	f
end

# ╔═╡ 47ff4422-b91d-4748-9937-b0ddc8e9ecd5
begin
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000336023.9|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409005.1|TUBA1B-202|TUBA1B|1627|protein_coding|")
	local f = Figure(size = (600, 600))
	local ax = Axis(f[1, 1],
				    aspect = 0.3,
					xticks=(1:2, ["Y", "m5C"]),
				    ylabel="Reads")

	local df = DataFrame(Tables.table(probs[:, [1580, 1592]] .> 0.75))
	dropmissing!(df)
	heatmap!(ax, transpose(Array(sort(df))), colormap=[:white, :black])
	f
end

# ╔═╡ 9fc1f07c-48e1-420b-9335-4bfa031e86cf
ivt[ivt.genomicPos .== 49127828, :]

# ╔═╡ fc0494df-e2fa-4a36-a0e0-dc6de09350cb
ivt[ivt.genomicPos .== 49127807, :]

# ╔═╡ fd8b9aeb-d76b-431d-b75f-ef86bfa4c7de
ivt[ivt.ref_id .== "ENST00000336023.9|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409005.1|TUBA1B-202|TUBA1B|1627|protein_coding|" .&& ivt.pos .> 1500, :]

# ╔═╡ 89aa206f-14bd-485d-ae1c-d7b6cabdfee6
#=╠═╡
an_sig_comod[((an_sig_comod.mod1 .=== "m6A" .&& an_sig_comod.mod2 .=== "m5C") .|| (an_sig_comod.mod1 .=== "m5C" .&& an_sig_comod.mod2 .=== "m6A")) .&& an_sig_comod.observed .> an_sig_comod.expected, :]
  ╠═╡ =#

# ╔═╡ 5c70bdac-40a0-4378-9d00-1bc8966da523
[1 1 1 1 1 1 1 1 1;
 1 2 2 2 2 2 2 2 1;
 1 2 3 3 3 3 3 2 1;
 1 2 3 4 4 4 3 2 1;
 1 2 3 4 5 4 3 2 1;
 1 2 3 4 4 4 3 2 1;
 1 2 3 3 3 3 3 2 1;
 1 2 2 2 2 2 2 2 1;
 1 1 1 1 1 1 1 1 1]

# ╔═╡ 694914ad-e7d0-4cbf-a9db-55fb49e393ac
#=╠═╡
begin
	local df = copy(simple_comods)
	df[!, :gene] = [r[2] for r in split.(df.reference, "|")]
	local area = [1 1 1 1 1 1 1 1 1;
				  1 2 2 2 2 2 2 2 1;
				  1 2 3 3 3 3 3 2 1;
				  1 2 3 4 4 4 3 2 1;
				  1 2 3 4 5 4 3 2 1;
				  1 2 3 4 4 4 3 2 1;
				  1 2 3 3 3 3 3 2 1;
				  1 2 2 2 2 2 2 2 1;
				  1 1 1 1 1 1 1 1 1]
	local t = nothing
	# for gene in filter(g -> length(unique(g.reference)) > 1, groupby(df, :gene))
		
	# 	t = groupby(gene, :reference)
	# 	break
	# end
	# t
	println(length(filter(g -> length(unique(g.reference)) == 5, groupby(df, :gene))))
	local t = filter(g -> length(unique(g.reference)) > 1, groupby(df, :gene))[9]
	local minpos = min(minimum(t.genomicPos1), minimum(t.genomicPos2))
	local maxpos = max(maximum(t.genomicPos1), maximum(t.genomicPos2))
	local N = maxpos - minpos + 1
	local matrices = []
	for ref in unique(t.reference)
		m = zeros((N, N))
		for row in eachrow(t[t.reference .== ref, :])
			value = row.relation == :positive ? 1 : -1
			centerr = row.genomicPos1 - minpos + 1
			centerc = row.genomicPos2 - minpos + 1
			left = max(centerc - 4, 1)
			right = min(centerc + 4, N)
			up = max(centerr - 4, 1)
			down = min(centerr + 4, N)
			println("$centerr $centerc $left $right $up $down")
			m[up:down, left:right] = value * area[(5 - (centerr - up)):(5 + (down - centerr)), (5 - (centerc - left)):(5 + (right - centerc))]
		end
		push!(matrices, m)
	end
	show(t)
	t
	heatmap(matrices[2])
	heatmap(std(cat(matrices..., dims=3), dims=3)[:, :, 1])
end
  ╠═╡ =#

# ╔═╡ 968ac33d-b86d-46cc-820f-c951f61a8e24
#=╠═╡
unique(map(r -> String(r[5]), split.(simple_comods.reference, "|"))) |> length
  ╠═╡ =#

# ╔═╡ 3b5d190e-5b56-4b7c-b389-68bbbb81f836
#=╠═╡
begin
	local df = copy(dropmissing(simple_comods))
	df[!, :gene] = [r[2] for r in split.(df.reference, "|")]
	# for gene in filter(g -> length(unique(g.reference)) > 1, groupby(df, :gene))
		
	# 	t = groupby(gene, :reference)
	# 	break
	# end
	local genes = []
	local scores = []
	local mod_counts = []
	local isoforms = []
	# local gene_scores = Dict()
	# local gene_mod_counts = Dict()
	for t in filter(g -> length(unique(g.reference)) > 1, groupby(df, :gene))
		gene = String(split(t[1, :reference], "|")[6])
		mods = sort(unique(vcat(t.genomicPos1, t.genomicPos2)))
		for i in 1:(length(mods)-1)
			if abs(mods[i] - mods[i + 1]) <= 4
				mods[i] = -1
			end
		end
		mods = mods[mods .!= -1]
		N = length(mods)
		T = length(unique(t.reference))
		M = zeros((N, N, T))
		for (i, ref) in enumerate(unique(t.reference))
			for row in eachrow(t[t.reference .== ref, :])
				m1 = [abs(m - row.genomicPos1) <= 4 for m in mods]
				m2 = [abs(m - row.genomicPos2) <= 4 for m in mods]
				M[m1, m2, i] .= row.relation == :positive ? 1 : -1
			end
		end
		score = sum(std(M, dims=3) .> 0) / (N*N)
		push!(scores, score)
		push!(genes, gene)
		push!(mod_counts, N)
		push!(isoforms, length(unique(t.reference)))
		# gene_scores[gene] = score
		# gene_mod_counts[gene] = N
	end
	# genes[sortperm(scores)]
	sort(DataFrame(gene=genes, score=scores, mods=mod_counts, isoforms=isoforms), :score)
end
  ╠═╡ =#

# ╔═╡ a617f69a-669d-499d-a40c-ef8f8b747e49
#=╠═╡
begin
	local df = copy(dropmissing(simple_comods))
	df[!, :gene] = [r[6] for r in split.(df.reference, "|")]

	# combine(groupby(df, :gene), :reference => length ∘ unique => :refs)
	local j = innerjoin(df, df,
						on=[:gene, :genomicPos1, :genomicPos2], 
						makeunique=true)
	sort(j[j.relation .!= j.relation_1, :], :gene)
end
  ╠═╡ =#

# ╔═╡ dd97de04-b552-4d11-a9a4-f46e0f5d670c
#=╠═╡
arc_plot_genomic(simple_comods[occursin.("ENST00000571732.5", simple_comods.reference), :], an_sig_peaks_ivt[occursin.("ENST00000571732.5", an_sig_peaks_ivt.ref_id), :]; range=(500, 1800), grange=(1344520, 1344980))
  ╠═╡ =#

# ╔═╡ 8ad73f0f-7e4f-4848-ae37-4e5780a1b0da
renamemods = df -> begin
	t = copy(df)
	local mapping = Dict(
		missing => missing,
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

# ╔═╡ 31ef15f2-8cf6-4711-82b8-59f34d2e043e
renamemods2 = df -> begin
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

# ╔═╡ 82a98275-f6d6-4ed1-9aa7-dd8ab59d9370
#=╠═╡
arc_plot_genomic(renamemods(simple_comods[occursin.("ENST00000264335.13", simple_comods.reference), :]), an_sig_peaks_ivt[occursin.("ENST00000264335.13", an_sig_peaks_ivt.ref_id), :]; range=(500, 1800), grange=(1344520, 1344980))
  ╠═╡ =#

# ╔═╡ ec42c786-4e6d-41e0-a40b-8694f4f99590
begin
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000571732.5|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000437356.1|YWHAE-208|YWHAE|1818|protein_coding|")
	local f = Figure(size = (600, 600))
	local ax = Axis(f[1, 1],
				    aspect = 0.3,
				    xticks=(1:2, ["m6A", "m5C"]),
				    ylabel="Reads")

	local df = DataFrame(Tables.table(probs[:, [1416, 1615]] .> 0.75))
	dropmissing!(df)
	println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	heatmap!(ax, transpose(Array(sort(df))), colormap=[:white, :black])
	f
end

# ╔═╡ 61e39473-08e6-45fb-bfd2-545001a2ecef
ivt[ivt.genomicPos .== 1344937, :]

# ╔═╡ dd6d235f-b102-46b0-9ea7-71857a4db0ba
begin
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000264335.13|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000259354.4|YWHAE-201|YWHAE|2052|protein_coding|")
	local f = Figure(size = (600, 600))
	local ax = Axis(f[1, 1],
				    aspect = 0.3,
				    xticks=(1:2, ["m6A", "m5C"]),
				    ylabel="Reads")

	local df = DataFrame(Tables.table(probs[:, [1385, 1584]] .> 0.75))
	dropmissing!(df)
	println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	println(nrow(df))
	println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	heatmap!(ax, transpose(Array(sort(df))), colormap=["#F0F0F0", :black])
	f
end

# ╔═╡ 1513ad3e-b7c9-4a14-8374-27d79eb6630e
#=╠═╡
begin
	local f = Figure(size = (1200, 1450))
	local ga = f[1, 1] = GridLayout()
	local gb = f[2, 1] = GridLayout()
	local gc = f[3, 1] = GridLayout()


	# local df = renamemods(an_sig_comod)
	# df[!, :relation] = map(r -> let s = (r.state1, r.state2)
	# 								if (s in [(1, 1), (0, 0)] && r.observed > r.expected) || (s in [(0, 1), (1, 0)] && r.observed < r.expected)
	# 									:positive
	# 								elseif (s in [(1, 1), (0, 0)] && r.observed < r.expected) || (s in [(0, 1), (1, 0)] && r.observed > r.expected)
	# 									:negative
	# 								else
	# 									:weird
	# 								end
	# 							end,
	# 					   eachrow(df))
	local df = renamemods(simple_comods)
	df[!, "mods"] = map(r -> join(sort([coalesce(r.mod1, ""), coalesce(r.mod2, "")]), "-"),
						 eachrow(df))
	println(sort(combine(groupby(df, :mods), nrow => :count), :count, rev=true))
	df[!, :distance] = abs.(df.pos1 .- df.pos2)
	# local selected = map(r -> r.mods in ["m⁶A-m⁷G", "m⁵C-m⁶A", "m⁶A-Ψ", "m⁵C-Ψ"], eachrow(df))
	local common = filter(r -> r.count > 40, combine(groupby(df, :mods), nrow => :count)).mods |> Set
	println(common)
	local selected = [non_trivial_simple(r) && r.mod1 !== missing && r.mod2 !== missing && r.mods in common && r.mod1 != r.mod2 for r in eachrow(df)]
	
	# sort(combine(groupby(df[selected, :], [:mods, :relation]), nrow => :count))

	local plt = data(df[selected, :]) *
		mapping(:mods => "Modification pair", color=:relation, dodge=:relation) *
		frequency()
	draw!(ga[1, 1], plt)
	
	local plt = data(df[selected, :]) *
		mapping(:mods => "Modifications pair",
				:distance => "Distance",
				color=:relation => "Association",
				dodge=:relation) *
		visual(BoxPlot)
	draw!(ga[1, 2], plt)

	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.("ENST00000264335.13", simple_comods.reference)
	local tx2_mask = occursin.("ENST00000571732.5", simple_comods.reference)
	local colorrange = (-log10(maximum(simple_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(simple_comods[tx1_mask .|| tx2_mask, :pvalue])))
	arc_plot_genomic2(gb[1, 1],
					  renamemods(simple_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.("ENST00000264335.13", an_sig_peaks_ivt.ref_id), :];
					  range=(500, 1800),
					  grange=(1344545, 1344950),
					  colorrange=colorrange,
					  title="YWHAE-201",
					  highlight=(1388, 1587))
	arc_plot_genomic2(gb[2, 1],
					  renamemods(simple_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.("ENST00000571732.5", an_sig_peaks_ivt.ref_id), :];
					  range=(500, 1800),
					  grange=(1344545, 1344950),
					  colorrange=colorrange,
					  title="YWHAE-208",
					  highlight=(1419, 1618))
	Colorbar(gb[1:2, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000264335.13|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000259354.4|YWHAE-201|YWHAE|2052|protein_coding|")
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [1385, 1584]] .> 0.75))
	dropmissing!(df)
	local n2 = nrow(df)
	local ax1 = Axis(gc[1, 1],
					 title="YWHAE-201",
				     # aspect = 0.3,
				     yticks=(1:2, ["m⁶A", "m⁵C"]),
					 xticks=(0:(n2/10):n2, ["$x%" for x in 0:10:100]),
				     xlabel="")
	println("n2: $n2")
	println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	heatmap!(ax1, Array(sort(df)), colormap=["#dbdbdb", :black])

	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000571732.5|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000437356.1|YWHAE-208|YWHAE|1818|protein_coding|")
	# local f = Figure(size = (600, 600))

	local df = DataFrame(Tables.table(probs[:, [1416, 1615]] .> 0.75))
	dropmissing!(df)
	local n1 = nrow(df)
	local ax2 = Axis(gc[2, 1],
				 title="YWHAE-208",
				 # aspect = 0.3,
				 yticks=(1:2, ["m⁶A", "m⁵C"]),
				 xticks=(0:(n1/10):n1, ["$x%" for x in 0:10:100]),
				 xlabel="Reads")
	println("n1: $n1")
	println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	heatmap!(ax2, Array(sort(df)), colormap=["#dbdbdb", :black])
	
	# xlims!(ax1, (0, max(n1, n2)))
	# xlims!(ax2, (0, max(n1, n2)))
	

	# colsize!(f.layout, 2, Relative(7/10))
	# colsize!(f.layout, 3, Relative(1/10))
	rowsize!(f.layout, 1, Relative(1.8/7))
	rowsize!(f.layout, 2, Relative(3.9/7))
	rowsize!(f.layout, 3, Relative(1.3/7))

	rowsize!(gb, 1, Relative(1/2))
	rowsize!(gb, 2, Relative(1/2))

	Legend(ga[2, 1:2],
		   [PolyElement(color=Makie.wong_colors()[1]),
		    PolyElement(color=Makie.wong_colors()[2])],
	       ["Negative association", "Positive association"],
		   orientation = :horizontal,
		   framevisible = false)
	Legend(gb[3, 1],
		   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π])],
	       ["Negative association", "Positive association"],
		   orientation = :horizontal,
		   framevisible = false)
	Legend(gc[3, 1],
		   [PolyElement(color=:black),
			PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	       ["Modified", "Not modified"],
		   orientation = :horizontal,
		   framevisible = false)


	for (label, layout) in zip(["A", "B", "C", "D"],
							   [ga[1, 1], ga[1, 2], gb[1, 1], gc[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end
  ╠═╡ =#

# ╔═╡ 40eaa862-f1ee-4fc1-ad89-5897f2938f0d
Makie.wong_colors()[1]

# ╔═╡ 3e577c30-ae50-4deb-a361-83912f0904fc
an_sig_peaks_ivt.mod |> countmap

# ╔═╡ 619bda6e-9209-4634-88bb-a3fbff14f74b
(an_sig_peaks_ivt.mod .!== missing) |> countmap

# ╔═╡ 301d30fe-4485-47a7-97e5-d093d2b8feff
#=╠═╡
an_sig_comod.reference |> unique |> length
  ╠═╡ =#

# ╔═╡ 73ddf0c4-5961-4c81-b5d6-833cd0c06f0d
#=╠═╡
filter(non_trivial, an_sig_comod).reference |> unique |> length
  ╠═╡ =#

# ╔═╡ e3665e85-ba56-4929-bd10-20c949e951e7
#=╠═╡
sort(combine(groupby(an_sig_comod, [:mod1, :mod2]), nrow => :count), :count, rev=true)
  ╠═╡ =#

# ╔═╡ 79736277-c65c-4ec5-94b0-2268b38bc095
#=╠═╡
sort(combine(groupby(filter(non_trivial_simple, simple_comods), [:mod1, :mod2]), nrow => :count), :count, rev=true)
  ╠═╡ =#

# ╔═╡ f58d900d-a5f3-40eb-8533-95542d92de30
#=╠═╡
begin
	save_bedpe(simple_comods, "ENST00000571732.5|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000437356.1|YWHAE-208|YWHAE|1818|protein_coding|", "/home/mzdravkov/ENST00000571732.5.bedpe")
	save_bedpe(simple_comods, "ENST00000264335.13|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000259354.4|YWHAE-201|YWHAE|2052|protein_coding|", "/home/mzdravkov/ENST00000264335.13.bedpe")
end
  ╠═╡ =#

# ╔═╡ a5a081d7-2ea0-4f7f-afd8-4e2f70ff3345
#=╠═╡
simple_comods[occursin.("TMEM241", simple_comods.reference), :]
  ╠═╡ =#

# ╔═╡ c90a46c3-e3ca-4088-83af-1a1368559f20
function get_transcript_signals(files, ref)
	samples = []
	reads = []
	intensities = []
	dwells = []
	for (sample, file) in files
		db = SQLite.DB(file)
		tx_id = DBInterface.execute(db, "SELECT * FROM transcripts WHERE name = \"$ref\"") |> DataFrame
		tx_id = tx_id[1, :id]
		
		t = DBInterface.execute(db, "SELECT read_id, intensity, dwell FROM signal_data s JOIN transcripts t ON t.id = s.transcript_id WHERE s.transcript_id = $tx_id") |> DataFrame
		intensity = reinterpret.(Float32, t.intensity)
		dwell = reinterpret.(Float32, t.dwell)

		append!(samples, repeat([sample], length(t.read_id)))
		append!(reads, t.read_id)
		append!(intensities, intensity)
		append!(dwells, dwell)
	end
	DataFrame(sample = samples, read = reads, intensity = intensities, dwell = dwells)
end

# ╔═╡ 38cc2f30-0ece-4dd7-ae64-c5b3a28b99fd
dbs = Dict(
	:WT_1 => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/WT_1_eventalign_collapsed.sqlite",
	:WT_2 => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/WT_2_eventalign_collapsed.sqlite",
	:WT_SPK => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/WT_SPK_eventalign_collapsed.sqlite",
	:STORM_1 => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/STORM_1_eventalign_collapsed.sqlite",
	:STORM_2 => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/STORM_2_eventalign_collapsed.sqlite",
	:STORM_SPK => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/STORM_SPK_eventalign_collapsed.sqlite",
	:IVT_1 => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/IVT_tailed_1_eventalign_collapsed.sqlite",
	:IVT_2 => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/IVT_tailed_2_eventalign_collapsed.sqlite",
	:IVT_3 => "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/IVT_notail_eventalign_collapsed.sqlite"
)
	

# ╔═╡ 8c658c04-0948-4d33-83fe-5d9756f43552
function get_pos_data(tx_data, pos)
	pos_data = copy(tx_data)
	# pos_data[:, :sample] = map(s -> sample_ids[s], pos_data.sample)
	pos_data[:, :intensity_p] = map(r -> Float64(r[pos]), pos_data.intensity)
	pos_data[:, :dwell_p] = map(r -> log10(Float64(r[pos])), pos_data.dwell)
	pos_data = pos_data[.! isnan.(pos_data.intensity_p), :]
	pos_data = pos_data[:, [:sample, :read, :intensity_p, :dwell_p]]
	rename!(pos_data, :intensity_p => :intensity, :dwell_p => :dwell)
	pos_data[:, :condition] = map(string, first.(split.(string.(pos_data.sample), "_")))
	pos_data
end

# ╔═╡ 5ed770d3-da06-4f86-b8fc-3425ddc7bdfc
tuba1b201_data = get_transcript_signals(dbs, "ENST00000332858.10|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409006.1|TUBA1B-201|TUBA1B|3198|retained_intron|")

# ╔═╡ 5be0dc53-9039-4324-b1e2-e24c344f069c
tuba1b202_data = get_transcript_signals(dbs, "ENST00000336023.9|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409005.1|TUBA1B-202|TUBA1B|1627|protein_coding|")

# ╔═╡ 46068dad-f73f-4660-9c09-2fbc85be0a09
(
	data(get_pos_data(tuba1b201_data[occursin.("WT", String.(tuba1b201_data.sample)), :], 3155)) *
		mapping(:dwell => "log10(dwell)", :intensity, color = :condition) *
		(visual(Scatter, markersize=3.5, alpha=0.45)) +
	data(get_pos_data(tuba1b201_data[occursin.("IVT", String.(tuba1b201_data.sample)), :], 3155)) *
		mapping(:dwell => "log10(dwell)", :intensity, color = :condition) *
		(visual(Scatter, markersize=3.5, alpha=0.45))) |> draw(scales(Color = (; palette = [:green, :purple]));
																																												   axis = (; title = "TUBA1B-201, position 3155 (Y)", limits = ((-3, 0.5), (80, 110))))

# ╔═╡ 6ba634cb-e770-44ba-a044-578029cfd8ec
begin
	local f = Figure()
	local ga = f[1, 1] = GridLayout()

	# (data(get_pos_data(tuba1b201_data[occursin.("WT", String.(tuba1b201_data.sample)), :], 3155)) *
	# 	mapping(:dwell => "log10(dwell)", :intensity, color = :condition) *
	# 	(visual(Scatter, markersize=3.5, alpha=0.45)) +
	# data(get_pos_data(tuba1b201_data[occursin.("IVT", String.(tuba1b201_data.sample)), :], 3155)) *
	# 	mapping(:dwell => "log10(dwell)", :intensity, color = :condition) *
	# 	(visual(Scatter, markersize=3.5, alpha=0.45))) |> draw!(f, scales(Color = (; palette = [:green, :purple]));
	# 																																											   axis = (; title = "TUBA1B-201, position 3155 (Y)", limits = ((-3, 0.5), (80, 110))))


	local t201p3155 = get_pos_data(tuba1b201_data, 3155)
	
	local ax = Axis(ga[1, 1])
	scatter!(ax, t201p3155.dwell, t201p3155.intensity, color=ifelse.(t201p3155.condition .== "WT", :purple, :green), markersize=3.5, alpha=0.4)
	xlims!(ax, (-3, 0.5))
	ylims!(ax, (75, 115))

	local t201p3172 = get_pos_data(tuba1b201_data, 3167)
	
	local ax = Axis(ga[1, 2])
	scatter!(ax, t201p3172.dwell, t201p3172.intensity, color=ifelse.(t201p3172.condition .== "WT", :purple, :green), markersize=3.5, alpha=0.4)
	xlims!(ax, (-3, 0.5))
	ylims!(ax, (50, 100))

	
	local t202p1580 = get_pos_data(tuba1b202_data, 1580)

	local ax = Axis(ga[2, 1])
	scatter!(ax, t202p1580.dwell, t202p1580.intensity, color=ifelse.(t202p1580.condition .== "WT", :purple, :green), markersize=3.5, alpha=0.4)
	xlims!(ax, (-3, 0.5))
	ylims!(ax, (75, 115))


	local t202p1592 = get_pos_data(tuba1b202_data, 1592)

	local ax = Axis(ga[2, 2])
	scatter!(ax, t202p1592.dwell, t202p1592.intensity, color=ifelse.(t202p1592.condition .== "WT", :purple, :green), markersize=3.5, alpha=0.4)
	xlims!(ax, (-3, 0.5))
	ylims!(ax, (50, 100))
	
	f
	
end

# ╔═╡ d4443059-dade-4a9a-aeef-75408913ce16


# ╔═╡ 10c9acd5-47d0-4f2a-a619-44fd903c4b36


# ╔═╡ 26d0dc17-cd92-47d4-bec0-33ad54c002cd
#=╠═╡
begin
	local f = Figure()


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)

	local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
	# mods = ifelse.(mods .=== missing, "Unknown", mods)
	local M = length(mods)

	local ax = Axis(f[1, 1],
					title = "Number of associated modification pair sites by modification types",
				    xticks = (1:M, mods),
				    yticks = (1:M, mods))

	df = unique(df[:, [:pos1, :pos2, :mod1, :mod2]])

	local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)

	local t = 0
	local counts = zeros((length(mods), length(mods)))
	for (i, m1) in enumerate(mods)
		for (j, m2) in enumerate(mods)
			if j < i
				continue
			end
			c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
						(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
			counts[i, j] = c
			t += c
		end
	end
	println(t)

	heatmap!(ax, 1:M, 1:M, counts, colormap = :heat)
	text!(ax,
		  [string(Int(counts[i, j])) for i in 1:M for j in i:M],
		  position = [Point2f(x, y) for x in 1:M for y in x:M],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[i, j] > 2000 ? :white : :black for i in 1:M for j in i:M])

	f
end
  ╠═╡ =#

# ╔═╡ c2cf87ff-d11c-4bb3-b44d-6a16f4040c88
#=╠═╡
begin
	local mods = ("m5C", "Y")
	local df = an_sig_comod[(an_sig_comod.mod1 .=== mods[1] .&& an_sig_comod.mod2 .=== mods[2]) .||
		(an_sig_comod.mod1 .=== mods[2] .&& an_sig_comod.mod2 .=== mods[1]), :]
	df[:, :dir] = df.observed .> df.expected
	combine(groupby(df, [:state1, :state2, :dir]), nrow => :count)
end
  ╠═╡ =#

# ╔═╡ fd549de3-8c1f-4358-8ec0-48fdd85d0bd1
function non_trivial(row)
	if row.mod1 !== row.mod2
		true
	elseif (row.state1, row.state2) == (0, 0) && row.observed > row.expected
		false
	elseif (row.state1, row.state2) == (0, 1) && row.observed < row.expected
		false
	elseif (row.state1, row.state2) == (1, 0) && row.observed < row.expected
		false
	elseif (row.state1, row.state2) == (1, 1) && row.observed > row.expected
		false
	else
		true
	end
end

# ╔═╡ cbe49401-5a36-4f40-9ece-4cccfc89f877
function non_trivial_simple(row)
	row.mod1 !== row.mod2 || row.relation == :negative
end

# ╔═╡ b37f9ca0-5a53-4716-9589-9b57268f72b2
#=╠═╡
begin
	local f = Figure()


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)

	df = filter(non_trivial, df)

	local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
	# mods = ifelse.(mods .=== missing, "Unknown", mods)
	local M = length(mods)

	local ax = Axis(f[1, 1],
					title = "Number of non-trivial associated modification pair sites by modification types",
				    xticks = (1:M, mods),
				    yticks = (1:M, mods))

	df = unique(df[:, [:pos1, :pos2, :mod1, :mod2]])
	
	local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)

	local t = 0
	local counts = zeros((length(mods), length(mods)))
	for (i, m1) in enumerate(mods)
		for (j, m2) in enumerate(mods)
			if j < i
				continue
			end
			c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
						(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
			counts[i, j] = c
			t += c
		end
	end
	println(t)

	heatmap!(ax, 1:M, 1:M, counts, colormap = :heat)
	text!(ax,
		  [string(Int(counts[i, j])) for i in 1:M for j in i:M],
		  position = [Point2f(x, y) for x in 1:M for y in x:M],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[i, j] > 1500 ? :white : :black for i in 1:M for j in i:M])

	f
end
  ╠═╡ =#

# ╔═╡ 67638f93-2a44-4772-9c35-2e728cebc277


# ╔═╡ c99dbb5d-77da-4c2e-9b7a-874f5ead6574
#=╠═╡
map(println, sort(unique(map(r -> r[6], split.(filter(r -> non_trivial(r) && r.mod1 !== missing && r.mod2 !== missing, an_sig_comod).reference, "|")))))
  ╠═╡ =#

# ╔═╡ 409ec9ef-22c6-4192-a53c-382d059ff299
#=╠═╡
begin
	local f = Figure(size = (1200, 800))
	local ax = Axis(f[1, 1],
					title = "Distance between significant non-trivially associated sites",
				    xlabel = "Distance (nt)",
				    ylabel = "Count",
				    xticks = 0:50:1400)
	
	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)
	df = filter(non_trivial, df)

	hist!(ax, abs.(df.pos2 .- df.pos1), bins=1400)
	f
end
  ╠═╡ =#

# ╔═╡ 664f5ea5-378a-45cc-a6cf-cab779b7f933
#=╠═╡
begin
	# local df = copy(an_sig_comod)
	
	# sort(combine(groupby(df, [:pattern, :direction]), nrow => :count), [:pattern, :direction])
	
	local f = Figure()


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)
	df[:, :pattern] = map(r -> Dict((0, 0) => "00", (0, 1) => "01", (1, 0) => "01", (1, 1) => "11")[(r.state1, r.state2)], eachrow(df))
	df[:, :direction] = sign.(log2.(df.observed ./ df.expected))
	df = df[df.pattern .== "01" .&& df.direction .== 1, :]

	local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
	# mods = ifelse.(mods .=== missing, "Unknown", mods)
	local M = length(mods)

	local ax = Axis(f[1, 1],
					title = "Number of associated modification pair sites by modification types",
				    xticks = (1:M, mods),
				    yticks = (1:M, mods))

	local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)

	local t = 0
	local counts = zeros((length(mods), length(mods)))
	for (i, m1) in enumerate(mods)
		for (j, m2) in enumerate(mods)
			if j < i
				continue
			end
			c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
						(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
			counts[i, j] = c
			t += c
		end
	end
	println(t)

	heatmap!(ax, 1:M, 1:M, counts, colormap = :heat)
	text!(ax,
		  [string(Int(counts[i, j])) for i in 1:M for j in i:M],
		  position = [Point2f(x, y) for x in 1:M for y in x:M],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[i, j] > 4000 ? :white : :black for i in 1:M for j in i:M])

	f
end
  ╠═╡ =#

# ╔═╡ d1c2e993-4ba4-46a5-840c-38bfbf397858
#=╠═╡
begin
	# local df = copy(an_sig_comod)
	
	# sort(combine(groupby(df, [:pattern, :direction]), nrow => :count), [:pattern, :direction])
	
	local f = Figure(size = (1600, 950))

	for (c, pattern) in enumerate(["00", "01", "11"]), (r, direction) in enumerate([1, -1])

		local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)
		df[:, :pattern] = map(r -> Dict((0, 0) => "00", (0, 1) => "01", (1, 0) => "01", (1, 1) => "11")[(r.state1, r.state2)], eachrow(df))
		df[:, :direction] = sign.(log2.(df.observed ./ df.expected))
		
	
		local mods = Vector(union(unique(df.mod1), unique(df.mod2)))
		# mods = ifelse.(mods .=== missing, "Unknown", mods)
		local M = length(mods)

		dir = direction == 1 ? "common" : "rare"
		local ax = Axis(f[r, c],
						title = "# mod-associated pairs | pattern: $pattern direction: $dir",
					    xticks = (1:M, mods),
					    yticks = (1:M, mods))
	
		local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)
	
		local t = 0
		local counts = zeros((length(mods), length(mods)))
		for (i, m1) in enumerate(mods)
			for (j, m2) in enumerate(mods)
				if j < i
					continue
				end
				c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
							(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
				counts[i, j] = c
				t += c
			end
		end
		println(t)
	
		df = df[df.pattern .== pattern .&& df.direction .== direction, :]
		local grouped = combine(groupby(df, [:mod1, :mod2]), nrow => :count)
	
		local t = 0
		local subcounts = zeros((length(mods), length(mods)))
		for (i, m1) in enumerate(mods)
			for (j, m2) in enumerate(mods)
				if j < i
					continue
				end
				c = grouped[(grouped.mod1 .== m1 .&& grouped.mod2 .== m2) .||
							(grouped.mod1 .== m2 .&& grouped.mod2 .== m1), :count] |> sum
				subcounts[i, j] = c
				t += c
			end
		end
	
		heatmap!(ax, 1:M, 1:M, subcounts./(counts.+1), colormap = :heat, colorrange = [0, 1])
		text!(ax,
			  [string(Int(subcounts[i, j])) for i in 1:M for j in i:M],
			  position = [Point2f(x, y) for x in 1:M for y in x:M],
			  align = (:center, :center),
			  fontsize = 14,
			  color = [subcounts[i, j]/(counts[i, j] + 1) > 0.7 ? :white : :black for i in 1:M for j in i:M])
	end

	Colorbar(f[1:2, 4], colorrange = (0, 1), colormap = :heat)

	f
end
  ╠═╡ =#

# ╔═╡ cb81e3ba-5eea-4169-a4ff-5b7433e8b171
#=╠═╡
groupby(an_sig_comod[map(v -> v == ["m6A", "Y"], Vector.(eachrow(coalesce.(an_sig_comod[:, [:mod1, :mod2]], "Unknown")))), :], [:reference, :pos1, :pos2]) |> println
  ╠═╡ =#

# ╔═╡ 07f53334-5428-4890-b02a-9e3c5de7b22c
#=╠═╡
begin
	local f = Figure(size = (500, 1000))


	local df = ifelse.(an_sig_comod .=== missing, "Unknown", an_sig_comod)

	axes = []

	for (i, (mod, mod2)) in enumerate([("m6A", "m5C"),
									   ("m6A", "Y"),
									   ("m6A", "m1A"),
									   ("m5C", "Y"),
									   ("m6A", "Nm")])
		# "m5C", "m1A", "Y", "Nm"])
		ax = Axis(f[i, 1],
				  title = "Distance between associated $mod-$mod2 pairs",
				  xticks = 0:100:1200)
		instances = df[(df.mod1 .== mod2 .&& df.mod2 .== mod) .||
		   			   (df.mod1 .== mod .&& df.mod2 .== mod2), [:pos1, :pos2]]

		instances = unique(instances)


		hist!(ax, abs.(instances.pos2 .- instances.pos1), bins = nrow(instances))
		push!(axes, ax)
	end

	linkxaxes!(axes...)
	
	f
end
  ╠═╡ =#

# ╔═╡ 40bbff56-4264-4d9a-af0e-59485c923f92
begin
	local f = Figure(size = (400, 1200))

	local df = sig_comod

	# local axes = []
	

	local colors = [:orange, :green, :purple, :grey]
	for (i, pattern) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)])
		ax = Axis(f[i, 1],
		  title = "Distance between associated pairs with pattern $pattern")

		instances = df[df.state1 .== pattern[1] .&&
		   			   df.state2 .== pattern[2], :]

		instances = instances[abs.(instances.pos2 .- instances.pos1) .> 10, :]

		if sum(instances.observed .> instances.expected) > nrow(instances)/2
			hist!(ax, abs.(instances[instances.observed .> instances.expected, :].pos2 .- instances[instances.observed .> instances.expected, :].pos1),
					 label = "obs>exp",
					 color = (:orange, 0.5),
				   	 bins = nrow(instances[instances.observed .> instances.expected, :]))
			hist!(ax, abs.(instances[instances.observed .< instances.expected, :].pos2 .- instances[instances.observed .< instances.expected, :].pos1),
					 label = "obs<exp",
					 color = (:blue, 0.5),
				     bins = nrow(instances[instances.observed .< instances.expected, :]))
		else
			hist!(ax, abs.(instances[instances.observed .< instances.expected, :].pos2 .- instances[instances.observed .< instances.expected, :].pos1),
				 label = "obs<exp",
				 color = (:blue, 0.5),
			     bins = nrow(instances[instances.observed .< instances.expected, :]))
			hist!(ax, abs.(instances[instances.observed .> instances.expected, :].pos2 .- instances[instances.observed .> instances.expected, :].pos1),
				 label = "obs>exp",
				 color = (:orange, 0.5),
			   	 bins = nrow(instances[instances.observed .> instances.expected, :]))
		end
		# hist!(ax, abs.(instances.pos2 .- instances.pos1),
		# 	  	 # label = ifelse.(instances.observed .> instances.expected,
		# 					# 	 "obs>exp",
		# 					# 	 "obs<exp"),
		# 		 color = ifelse.(instances.observed .> instances.expected,
		# 						 :orange,
		# 						 :blue),

		# 	     bins = nrow(instances))

		axislegend(ax)
		# xlims!(ax, (10, 500))

		# push!(axes, ax)
	end

	# linkxaxes!(axes...)
	# linkyaxes!(axes...)

	f
end

# ╔═╡ 0ee4e9af-4b1e-41fb-8746-12ff51d20e27
#=╠═╡
an_sig_comod[occursin.("lnc", an_sig_comod.reference), :]
  ╠═╡ =#

# ╔═╡ 9503a9e1-8c7e-4dfe-aeed-015d7ec1095c
# ╠═╡ disabled = true
#=╠═╡
import SQLite
  ╠═╡ =#

# ╔═╡ 58e326c7-e5b5-4f97-884e-7321a8084139
import Clustering

# ╔═╡ e30b68c9-9aaa-44f7-86d9-c8e1ab9254dd
map(r -> println(split(r, "|")[6]), unique(sig_comod.reference))

# ╔═╡ 3dc430dc-a218-42bf-90da-8ba8ba6dd56a
begin
	local db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))


	local reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	local probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	probs = ifelse.(probs .== -1, 0, probs)
	probs ./= 100

	# local c = Clustering.kmeans(probs, 4)
	# local idx =  sortperm(assignments(c))
	# # M = M[:,idx]
	# println(idx)

	local N = size(probs, 2)
	local D = size(probs, 1)
	
	dicer_probs = probs
	dicer_distances = zeros((N, N))
	for i in 1:N
		for j in 1:N
			dicer_distances[i, j] = sum((probs[:, i] .- probs[:, j]) .^ 2)/D
		end
	end
end

# ╔═╡ 46b1d7d3-2e9e-4d19-91d8-07150e981ef0
begin
	local f = Figure()

	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))

	local ax = Axis(f[1, 1],
					title = "Modification probability for significant DANCR sites",
					xlabel = "Significant positions",
					ylabel = "Reads",
				    xticks = (1:length(positions), map(string, positions)),
				    xticklabelrotation = 45)

	local clust = Clustering.hclust(dicer_distances, branchorder = :barjoseph)
	heatmap!(ax, dicer_probs[:, clust.order])

	Colorbar(f[1, 2], colorrange = [0, 1])

	f
end

# ╔═╡ 809b0c5c-1538-4bbe-bb4d-293f6058b2a6
#=╠═╡
function arc_plot2(ref; range = nothing)
	name = split(ref, "|")[5]

	df = an_sig_comod[an_sig_comod.reference .== ref, :]
	if range !== nothing
		df = df[between.(df.pos1, range[1], range[2]) .&& between.(df.pos2, range[1], range[2]), :]
	end
	positions = unique(union(df.pos1, df.pos2))
	
	f = Figure(size = (1000, 700))
	ax = Axis(f[1, 1],
			  title = "$name modification site associations",
			  xlabel = "Reference position",
			  # xticks = 700:20:960,
			  yticksvisible = false,
			  yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	# local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end
	allmods = comod[comod.reference .== ref, :]
	if range !== nothing
		allmods = allmods[between.(allmods.pos1, range[1], range[2]) .&& between.(allmods.pos2, range[1], range[2]), :]
	end
	for row in eachrow(allmods)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for row in eachrow(df)
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			colormap = :reds
			colorrange = [-200, 10]
		elseif row.state1 == 0 && row.state2 == 1
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [0, 10]
		elseif row.state1 == 1 && row.state2 == 0
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [10, 0]
		else
			colormap = :greens
			colorrange = [-200, 10]
		end
		direction = row.observed > row.expected ? 1 : -1
		y = row.observed > row.expected ? 2.5 : -2.5
		style = row.observed > row.expected ? :solid : :dash
		arcsize = π/10
		for i in 0:9
			arcstart = direction*i*arcsize
			arcend = direction*(i+1)*arcsize
			arc!(ax, Point2f(center, y), radius, arcstart, arcend,
				 color = i,
				 colorrange = colorrange,
				 colormap = colormap,
				 linestyle = :solid,
				 linewidth = 2)
				 # linewidth = -log10(row.qvalue))
		end
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :red),
			PolyElement(color = :green),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["Non-modified", "Modified", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ ab71927a-b85c-4834-a20a-cebe1f734b7b
#=╠═╡
an_sig_comod_gx[occursin.("ENST00000332858.10", an_sig_comod_gx.reference), :]
  ╠═╡ =#

# ╔═╡ 3a9274bf-3af7-46d0-880e-1ea0e8b421de
#=╠═╡
arc_plot("ENST00000332858.10|ENSG00000123416.15|OTTHUMG00000170410.2|OTTHUMT00000409006.1|TUBA1B-201|TUBA1B|3198|retained_intron|"; range=(3000, 4000))
  ╠═╡ =#

# ╔═╡ b5423ce9-cad4-48ec-88ec-33d12319d478
function arc_plot_genomic(df, allmods; range=nothing, grange=nothing)
	#name = split(ref, "|")[5]
	# df = an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
	#df = an_sig_comod[an_sig_comod.reference .== ref, :]
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
	if range !== nothing
		df = df[between.(df.pos1, range[1], range[2]) .&& between.(df.pos2, range[1], range[2]), :]
	end
	positions = unique(union(df.genomicPos1, df.genomicPos2))

	tstart = (maxpos - minpos)*0.1
	tend = (maxpos - minpos)*0.9
	xticks = (Base.range(tstart, tend, length=5),
			  map(t -> string(Integer(t+minpos)), Base.range(tstart, tend, length=5)))

	
	f = Figure(size = (1400, 600))
	ax = Axis(f[1, 1],
		  	  title = "Modification site associations",
			  xlabel = "Chromosome position",
			  xticks = xticks,
			  yticksvisible = false,
			  yticklabelsvisible = false)

	if grange !== nothing
		xlims!(ax, (0, maxpos-minpos))
	end

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	# local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	local mods = Dict()
	for row in eachrow(df)
		if row.genomicPos1 ∉ keys(mods) || mods[row.genomicPos1] === missing
			mods[row.genomicPos1] = row.mod1
		end
		if row.genomicPos2 ∉ keys(mods) || mods[row.genomicPos2] === missing
			mods[row.genomicPos2] = row.mod2
		end
	end
	# allmods = comod[comod.reference .== ref, :]
	if range !== nothing
		allmods = allmods[between.(allmods.pos, range[1], range[2]), :]
	end
	for row in eachrow(allmods)
		if row.genomicPos ∉ keys(mods) || mods[row.genomicPos] === missing
			mods[row.genomicPos] = "?"
		end
		# if row.genomicPos2 ∉ keys(mods) || mods[row.genomicPos2] === missing
		# 	mods[row.genomicPos2] = "?"
		# end
	end

	max_radius = maximum(abs.(df.genomicPos2 .- df.genomicPos1))/2
	arc_offset = log10(max_radius)

	if max_radius < 16
		ylims!(ax, (-16, 16))
	end

	hlines!(ax, [0, 0], [minimum(keys(mods)), maximum(keys(mods))], color = :lightgrey, alpha = 0.5, linewidth = 16)

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	colorrange = (-log10(maximum(df.qvalue)), -log10(minimum(df.qvalue)))
	println(colorrange)
	
	for (i, row) in enumerate(eachrow(df))
		radius = abs(row.genomicPos2 - row.genomicPos1)/2
		center = min(row.genomicPos1, row.genomicPos2) + radius
		# if row.relation == :positive
		# 	color = :green
		# else
		# 	color = :orange
		# end
		# if row.state1 == 0 && row.state2 == 0
		# 	color = :blue
		# elseif row.state1 == 0 && row.state2 == 1
		# 	color = :orange
		# elseif row.state1 == 1 && row.state2 == 0
		# 	color = :orange
		# else
		# 	color = :purple
		# end
		# angle = row.observed > row.expected ? π : -π
		angle = row.relation == :positive ? π : -π
		# y = row.observed > row.expected ? arc_offset : -arc_offset
		y = row.relation == :positive ? arc_offset : -arc_offset
		Vector.(eachrow(df[1:i, [:pos1, :pos2]]))
		prev = df[1:i, :]
		# prevcount = sum(prev.genomicPos1 .== row.genomicPos1 .&& prev.genomicPos2 .== row.genomicPos2 .&& (prev.observed .> prev.expected) .== (row.observed > row.expected))
		# prevcount = 1
		
		# y *= prevcount * 1.5*arc_offset
		if y > 0
			y = max(y, max_radius/16)
		end
		if y < 0
			y = min(y, -max_radius/16)
		end

		# style = row.observed > row.expected ? :solid : :dash
		style = :solid
		arc!(ax, Point2f(center, y), radius, 0, angle,
			 color = -log10(row.qvalue),
			 colorrange = colorrange,
			 linewidth = 2)
			 # linewidth = -log10(row.qvalue))
		
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center),
		  fontsize=12)
	
	Colorbar(f[1:2, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(q-value)")

	Legend(f[2, 1],
		   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["Positive", "Negative"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 5229ac7e-5a3c-4339-8c06-ce4e07f379c3
function format_with_commas(x)
    s = string(Int(round(x; digits=2)))
    return replace(s, r"(\d)(?=(\d{3})+(?!\d))" => s"\1,")
end

# ╔═╡ 25e2f5f3-0612-48ea-83fc-3c4e85eb793d
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

# ╔═╡ beb418fd-d39b-4c5e-bf9d-98e936b9c819
#=╠═╡
arc_plot2("ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|")
  ╠═╡ =#

# ╔═╡ de7c2cc7-7931-4290-97be-39b24634d755
#=╠═╡
arc_plot2("ENST00000370990.5|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000095785.1|SERBP1-202|SERBP1|1621|protein_coding|"; range = (1100, 1450))
  ╠═╡ =#

# ╔═╡ 088e365f-05f2-429c-9fa2-763ddf30b82b
#=╠═╡
arc_plot2("ENST00000605788.6|ENSG00000166197.17|OTTHUMG00000018944.5|OTTHUMT00000050012.2|NOLC1-209|NOLC1|3723|protein_coding|"; range = (3650, Inf))
  ╠═╡ =#

# ╔═╡ ea219179-f224-4aa3-bb1d-fdf494742363
#=╠═╡
begin
	local df = an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
	local positions = unique(union(df.pos1, df.pos2))
	
	local f = Figure(size = (1000, 700))
	local ax = Axis(f[1, 1],
				    title = "DANCR (lncRNA) modification site associations",
				    xlabel = "Reference position",
					xticks = 700:20:960,
				    yticksvisible = false,
				    yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	local mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end
	for row in eachrow(comod[comod.reference .== ref, :])
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	hlines!(ax, [0, 0], [minimum(keys(mods)), maximum(keys(mods))], color = :lightgrey, alpha = 0.5, linewidth = 16)

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for row in eachrow(df)
		println(row.pos1, " ", row.pos2, " ", row.state1, " ", row.state2, " ", row.expected, " ", row.observed)
		# angle = row.state1 == row.state2 == 0 ? -π : π
		# y = row.state1 == row.state2 == 0 ? -2.5 : 2.5
		# color = row.observed > row.expected ? :green : :red
		# color, cmap = colormap[(row.state1, row.state2)]
		# style = row.state1 == row.state2 == 0 ? :dash : :solid
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			colormap = :reds
			colorrange = [-200, 10]
		elseif row.state1 == 0 && row.state2 == 1
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [0, 10]
		elseif row.state1 == 1 && row.state2 == 0
			colormap = :diverging_gkr_60_10_c40_n256
			colorrange = [10, 0]
		else
			colormap = :greens
			colorrange = [-200, 10]
		end
		direction = row.observed > row.expected ? 1 : -1
		y = row.observed > row.expected ? 2.5 : -2.5
		style = row.observed > row.expected ? :solid : :dash
		arcsize = π/10
		for i in 0:9
			arcstart = direction*i*arcsize
			arcend = direction*(i+1)*arcsize
			arc!(ax, Point2f(center, y), radius, arcstart, arcend,
				 color = i,
				 colorrange = colorrange,
				 colormap = colormap,
				 linestyle = :solid,
				 linewidth = -log10(row.qvalue))
		end
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :red),
			PolyElement(color = :green),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["Non-modified", "Modified", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ bcc0fec3-9925-4c52-8411-a59583f9171e
#=╠═╡
begin
	local df = an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
	local positions = unique(union(df.pos1, df.pos2))
	
	local f = Figure(size = (1000, 700))
	local ax = Axis(f[1, 1],
				    title = "DANCR (lncRNA) modification site associations",
				    xlabel = "Reference position",
					xticks = 700:20:960,
				    yticksvisible = false,
				    yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	local mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end

	for row in eachrow(comod[comod.reference .== ref, :])
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for row in eachrow(df)
		println(-log10(row.qvalue))
		println(row.pos1, " ", row.pos2, " ", row.state1, " ", row.state2, " ", row.expected, " ", row.observed)
		# angle = row.state1 == row.state2 == 0 ? -π : π
		# y = row.state1 == row.state2 == 0 ? -2.5 : 2.5
		# color = row.observed > row.expected ? :green : :red
		# color, cmap = colormap[(row.state1, row.state2)]
		# style = row.state1 == row.state2 == 0 ? :dash : :solid
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			color = :blue
		elseif row.state1 == 0 && row.state2 == 1
			color = :orange
		elseif row.state1 == 1 && row.state2 == 0
			color = :orange
		else
			color = :purple
		end
		angle = row.observed > row.expected ? π : -π
		y = row.observed > row.expected ? 2.5 : -2.5
		style = row.observed > row.expected ? :solid : :dash
		arc!(ax, Point2f(center, y), radius, 0, angle,
			 color = color,
			 # linewidth = 2,
			 linewidth = -log10(row.qvalue))
		
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :blue),
			PolyElement(color = :orange),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["State 0-0", "State 0-1/1-0", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ b117d0ed-b44b-46a6-b9ea-204e1d5195e9
#=╠═╡
function arc_plot(ref; range = nothing)
	name = split(ref, "|")[5]
	# df = an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
	df = an_sig_comod[an_sig_comod.reference .== ref, :]
	if range !== nothing
		df = df[between.(df.pos1, range[1], range[2]) .&& between.(df.pos2, range[1], range[2]), :]
	end
	positions = unique(union(df.pos1, df.pos2))
	
	f = Figure(size = (1000, 700))
	ax = Axis(f[1, 1],
		  	  title = "$name modification site associations",
			  xlabel = "Reference position",
			  # xticks = 700:20:960,
			  yticksvisible = false,
			  yticklabelsvisible = false)

	# local colormap = Dict(
	# 	(0, 0) => ([0, 0], :grays),
	# 	(0, 1) => (0:1:1, :redgreensplit),
	# 	(1, 0) => (1:-1:0, :redgreensplit),
	# 	(1, 1) => (0:1, :green)
	# )

	# local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"

	local mods = Dict()
	for row in eachrow(df)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = row.mod1
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = row.mod2
		end
	end
	allmods = comod[comod.reference .== ref, :]
	if range !== nothing
		allmods = allmods[between.(allmods.pos1, range[1], range[2]) .&& between.(allmods.pos2, range[1], range[2]), :]
	end
	for row in eachrow(allmods)
		if row.pos1 ∉ keys(mods) || mods[row.pos1] === missing
			mods[row.pos1] = "?"
		end
		if row.pos2 ∉ keys(mods) || mods[row.pos2] === missing
			mods[row.pos2] = "?"
		end
	end

	max_radius = maximum(df.pos2 .- df.pos1)/2
	arc_offset = log10(max_radius)

	if max_radius < 16
		ylims!(ax, (-16, 16))
	end

	hlines!(ax, [0, 0], [minimum(keys(mods)), maximum(keys(mods))], color = :lightgrey, alpha = 0.5, linewidth = 16)

	# local mods = sort(an_peaks_ivt[an_peaks_ivt.ref_id .== ref .&& an_peaks_ivt.pos .> 760 .&& an_peaks_ivt.pos .< 960, [:pos, :mod]]) |> unique
	
	for (i, row) in enumerate(eachrow(df))
		radius = (row.pos2 - row.pos1)/2
		center = row.pos1 + radius
		if row.state1 == 0 && row.state2 == 0
			color = :blue
		elseif row.state1 == 0 && row.state2 == 1
			color = :orange
		elseif row.state1 == 1 && row.state2 == 0
			color = :orange
		else
			color = :purple
		end
		angle = row.observed > row.expected ? π : -π
		y = row.observed > row.expected ? arc_offset : -arc_offset
		Vector.(eachrow(df[1:i, [:pos1, :pos2]]))
		prev = df[1:i, :]
		prevcount = sum(prev.pos1 .== row.pos1 .&& prev.pos2 .== row.pos2 .&& (prev.observed .> prev.expected) .== (row.observed > row.expected))
		y *= prevcount * 1.5*arc_offset
		if y > 0
			y = max(y, max_radius/16)
		end
		if y < 0
			y = min(y, -max_radius/16)
		end

		style = row.observed > row.expected ? :solid : :dash
		arc!(ax, Point2f(center, y), radius, 0, angle,
			 color = color,
			 linewidth = 2)
			 # linewidth = -log10(row.qvalue))
		
	end

	text!(ax,
		  [Point2f(x, 0) for x in keys(mods)];
		  text = ifelse.(values(mods) .=== missing, "?", values(mods)),
		  align = (:center, :center))

	Legend(f[2, 1],
		   [PolyElement(color = :blue),
			PolyElement(color = :orange),
			PolyElement(color = :purple),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π])],
	       ["Co-unmodified (0-0)", "Mutually exclusive (0-1/1-0)", "Co-modified (1-1)", "Observed > Expected", "Observed < Expected"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ 5408e59e-53ae-4ef9-b9ae-e86be3119575
#=╠═╡
arc_plot("ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|")
  ╠═╡ =#

# ╔═╡ 72faa038-d71c-46f9-8110-d68b7633f2e6
#=╠═╡
arc_plot("ENST00000370990.5|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000095785.1|SERBP1-202|SERBP1|1621|protein_coding|"; range = (1100, 1500))
  ╠═╡ =#

# ╔═╡ 9b4cf8f2-84ae-4aba-981a-7ef0e669a453
#=╠═╡
an_sig_comod[occursin.("SERBP1-202", an_sig_comod.reference), :]
  ╠═╡ =#

# ╔═╡ d7fb857c-0871-4be2-8322-e36ab95e6f19
#=╠═╡
arc_plot("ENST00000370994.8|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000025983.2|SERBP1-203|SERBP1|6676|protein_coding|"; range = (1120, 1650))
  ╠═╡ =#

# ╔═╡ 05ff34ad-5cfc-48a4-bfc9-4d4c7cfda4d1
#=╠═╡
arc_plot("ENST00000552461.5|ENSG00000089157.16|OTTHUMG00000169317.10|OTTHUMT00000403466.1|RPLP0-226|RPLP0|2591|retained_intron|"; range = (2500, Inf))
  ╠═╡ =#

# ╔═╡ 171a33c2-e532-4167-82a6-27d830402104
#=╠═╡
arc_plot("ENST00000419869.7|ENSG00000023734.11|OTTHUMG00000168791.2|OTTHUMT00000401114.2|STRAP-202|STRAP|1874|protein_coding|", range = (1800, Inf))
  ╠═╡ =#

# ╔═╡ 657e7156-57ef-4f64-a4ba-13af4bb2784e
#=╠═╡
arc_plot("ENST00000264720.8|ENSG00000115207.15|OTTHUMG00000097785.9|-|GTF3C2-201|GTF3C2|3607|protein_coding|", range = (2000, Inf))
  ╠═╡ =#

# ╔═╡ 89114c37-3d83-4e4c-a46f-32cd04e85bf0
#=╠═╡
sort(an_sig_comod[an_sig_comod.mod1 .!== an_sig_comod.mod2 .&& an_sig_comod.mod1 .!== missing .&& an_sig_comod.mod2 .!== missing, :], :qvalue)
  ╠═╡ =#

# ╔═╡ d6b27bbe-824e-492b-9aea-dd5041f4917f
#=╠═╡
sort(innerjoin(innerjoin(an_sig_comod[an_sig_comod.mod1 .!== an_sig_comod.mod2 .&& an_sig_comod.mod1 .!== missing .&& an_sig_comod.mod2 .!== missing .&& (.! occursin.("GTF3C2", an_sig_comod.reference))
, :],
		  	an_sig_peaks_ivt[an_sig_peaks_ivt.source .!== "motif", [:ref_id, :pos]],
		  	on = [:reference => :ref_id, :pos1 => :pos]),
		  an_sig_peaks_ivt[an_sig_peaks_ivt.source .!== "motif", [:ref_id, :pos]],
		  on = [:reference => :ref_id, :pos2 => :pos]), :qvalue)
  ╠═╡ =#

# ╔═╡ e0bfd78c-25fa-4752-a4d7-660aa4b6cd66
begin
	local db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	local ref = "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]


	local reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	local orig_probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	local probs = ifelse.(orig_probs .== -1, 0, orig_probs)
	probs ./= 100

	# local c = Clustering.kmeans(probs, 4)
	# local idx =  sortperm(assignments(c))
	# # M = M[:,idx]
	# println(idx)

	local N = size(probs, 2)
	local D = size(probs, 1)
	
	slc5a6_probs = probs
	slc5a6_distances = zeros((N, N))
	for i in 1:N
		for j in 1:N
			# slc5a6_distances[i, j] = sqrt(sum((probs[:, i] .- probs[:, j]) .^ 2))/D
			slc5a6_distances[i, j] = sum(abs.(probs[:, i] .- probs[:, j]))/D
		end
	end
	slc5a6_probs = ifelse.(orig_probs .== -1, missing, probs)
end

# ╔═╡ b6a72db1-ec4d-4cc1-8174-c6179316dec7
function get_probs(ref)
	
	db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]


	reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	orig_probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	probs = ifelse.(orig_probs .== -1, 0, orig_probs)
	probs ./= 100

	# local c = Clustering.kmeans(probs, 4)
	# local idx =  sortperm(assignments(c))
	# # M = M[:,idx]
	# println(idx)

	N = size(probs, 2)
	D = size(probs, 1)
	
	distances = zeros((N, N))
	for i in 1:N
		for j in 1:N
			# slc5a6_distances[i, j] = sqrt(sum((probs[:, i] .- probs[:, j]) .^ 2))/D
			distances[i, j] = sum(abs.(probs[:, i] .- probs[:, j]))/D
		end
	end
	probs = ifelse.(orig_probs .== -1, missing, probs)
end

# ╔═╡ f2829ecb-cbdf-4e0d-953c-0072f03546ac
# ╠═╡ disabled = true
#=╠═╡
begin
	local ref = "ENST00000411630.7|ENSG00000226950.8|OTTHUMG00000154670.4|OTTHUMT00000336530.2|DANCR-201|DANCR|987|lncRNA|"
	local dancr_probs = get_probs(ref)
	DataFrame(transpose(dancr_probs), map(string, unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2])))) |> CSV.write("/home/mzdravkov/dancr_probs.tsv", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ 74279cfc-ab7a-4e9a-b048-ad17a378c20f
# ╠═╡ disabled = true
#=╠═╡
begin
	local ref = "ENST00000253024.10|ENSG00000130726.12|OTTHUMG00000183546.3|OTTHUMT00000467074.2|TRIM28-201|TRIM28|3364|protein_coding|"
	local trim28_probs = get_probs(ref)
	DataFrame(transpose(trim28_probs), map(string, unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2])))) |> CSV.write("/home/mzdravkov/trim28_probs.tsv", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ 5b5358c8-f35c-473f-9e92-56bdb4dbc6a7
arc_plot("ENST00000370990.5|ENSG00000142864.15|OTTHUMG00000009372.3|OTTHUMT00000095785.1|SERBP1-202|SERBP1|1621|protein_coding|")#
, range = (3200, Inf))

# ╔═╡ 9031f497-26d2-4b3b-ba34-3438df409fd9
#=╠═╡
an_sig_comod[an_sig_comod.reference .== "ENST00000316448.10|ENSG00000179218.15|OTTHUMG00000180574.3|OTTHUMT00000451952.2|CALR-201|CALR|1901|protein_coding|"


, :]
  ╠═╡ =#

# ╔═╡ ad11f743-1eee-4298-87bf-620a0f1cd1ec
an_sig_peaks_ivt[an_sig_peaks_ivt.ref_id .== "ENST00000264720.8|ENSG00000115207.15|OTTHUMG00000097785.9|-|GTF3C2-201|GTF3C2|3607|protein_coding|"


, :]

# ╔═╡ f793c68f-731e-426f-ad6e-ef5419f16370
#=╠═╡
an_sig_comod[occursin.("DANCR", an_sig_comod.reference), :]
  ╠═╡ =#

# ╔═╡ cc022b63-4a04-4142-976c-02756047076d
begin
	local f = Figure()

	local ref = "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|"
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	local positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]

	local ax = Axis(f[1, 1],
					title = "Modification probability for significant SLC5A6 sites",
					xlabel = "Significant positions",
					ylabel = "Reads",
				    xticks = (1:length(positions), map(string, positions)),
				    xticklabelrotation = 45)

	local clust = Clustering.hclust(slc5a6_distances, branchorder = :barjoseph)
	heatmap!(ax, slc5a6_probs[:, clust.order])

	Colorbar(f[1, 2], colorrange = [0, 1])

	f
	# clust
end

# ╔═╡ a6b2a25c-6fb2-4dd0-a7f1-ea16ceb63a97
DataFrame(transpose(slc5a6_probs), map(string, unique(union(comod[comod.reference .== "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|", :pos1], comod[comod.reference .== "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|", :pos2])))) |> CSV.write("/home/mzdravkov/slc5a6_probs.tsv", delim = '\t')

# ╔═╡ 69798f9b-f9a8-41f7-a2ff-8f333a6ed728
begin
	local df = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/datasets/miR-eCLIP/sample_comparisons/HEK293T_DMSO_vs_HEK293T_STM.miRNA_target_sample_comparisons.final.tsv"))
	rename!(df, [:chr, :start, :end, :stran, :gene, :gene_id, :feature, :mirna, :avg_ctrl, :avg_storm, :lfc, :pscore, :qscore])
	mireclip = df
end

# ╔═╡ 741a21f9-ea94-4620-83a6-a43dd5fc300d
#=╠═╡
begin
	local df = copy(an_sig_comod)
	df[:, :gene_id] = map(rx -> rx[2], split.(df.reference, "|"))
	df = innerjoin(
			innerjoin(df, rename(ivt[:, [:ref_id, :pos, :genomicPos]], :genomicPos => :genomicPos1),
				   	  on = [:reference => :ref_id, :pos1 => :pos]),
			rename(ivt[:, [:ref_id, :pos, :genomicPos]], :genomicPos => :genomicPos2),
			on =  [:reference => :ref_id, :pos2 => :pos])
	df = innerjoin(df, mireclip,
			  	   on = [:gene_id])

	# df[between.(df.genomicPos1, df.start, df.end) .|| between.(df.genomicPos2, df.start, df.end), :]

	df[df.qscore .>= 2, :]
end
  ╠═╡ =#

# ╔═╡ e3f624cc-1c45-461a-ac19-bed551b1ba4e
md"""
## STORM resistant mods
"""

# ╔═╡ 5e8c4ba1-9f8d-40df-b02e-cfeec4bf4d57
# ╠═╡ disabled = true
#=╠═╡
begin
	# local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A" .&& an_sig_peaks_ivt.source .=== "motif", :]
	local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A", :]
	println(length(unique(ivt_m6a.peak_id)))
	local refjoin = innerjoin(ivt_m6a, storm,
			  				  on = [:ref_id],
			  				  makeunique = true)
	refjoin = refjoin[between.(refjoin.genomicPos_1, refjoin.start, refjoin.end), :]

	local agg_ivt_peaks = combine(groupby(refjoin, :peak_id),
							   nrow => :overlaps,
							   :GMM_chi2_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_pval)
	storm_res_peaks = agg_ivt_peaks[agg_ivt_peaks.overlaps .== 9 .&& agg_ivt_peaks.min_storm_pval .> 0.2, :peak_id] |> unique
end
  ╠═╡ =#

# ╔═╡ e8cee36d-6886-4e27-806e-b4df1dadff7d
begin
	# local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A" .&& an_sig_peaks_ivt.source .=== "motif", :]
	local ivt_m6a = an_sig_peaks_ivt[an_sig_peaks_ivt.mod .=== "m6A", :]
	println(length(unique(ivt_m6a.peak_id)))
	local refjoin = innerjoin(ivt_m6a, storm,
			  				  on = [:ref_id],
			  				  makeunique = true)
	refjoin = refjoin[between.(refjoin.genomicPos_1, refjoin.start, refjoin.end), :]

	local agg_ivt_peaks = combine(groupby(refjoin, :peak_id),
							   nrow => :overlaps,
							   :GMM_chi2_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_pval,
							   :KS_intensity_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_ksi_pval,
							   :KS_dwell_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_ksd_pval)
	storm_res_peaks = agg_ivt_peaks[agg_ivt_peaks.overlaps .== 9 .&& agg_ivt_peaks.min_storm_pval .> 0.2 .&& agg_ivt_peaks.min_storm_ksi_pval .> 0.1 .&& agg_ivt_peaks.min_storm_ksd_pval .> 0.1, :peak_id] |> unique
end

# ╔═╡ 59a166ff-402f-44cd-91d4-1f58ca676e33
storm[storm.GMM_chi2_pvalue .=== missing, :]

# ╔═╡ 25c85ec0-af9e-4f74-bb9f-df7325604100
an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :source] |> countmap

# ╔═╡ 4cea1683-7f74-42dd-a0d5-5d66eca41253
map(r -> split(r, "|")[6], an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&& an_sig_peaks_ivt.source .== "HEK293T", :ref_id]) |> unique |> length

# ╔═╡ 58a03eb1-4847-4f48-a86e-9da1233023a3
map(r -> println(r[6]), split.(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&& an_sig_peaks_ivt.source .== "HEK293T", :ref_id] |> unique, "|")) 

# ╔═╡ 56b54a7c-fb57-4aa0-9e84-45ade0d06b5c
storm[storm.ref_id .== "ENST00000372586.4|ENSG00000175283.8|OTTHUMG00000020765.2|OTTHUMT00000054515.2|DOLK-201|DOLK|2074|protein_coding|" .&& storm.pos .== 1709, :]

# ╔═╡ 52a1db73-e957-41fd-a6dc-f483b874e9f2
begin
	local df = an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :]

	

	# for (i, sample) in enumerate(["WT_1", "IVT_1", "WT_2", "IVT_2", "WT_SPK", "IVT_3"])
	# 	violin!(ax, i, df[:, "$(sample)_mod"], label = sample)
	# end
	
	# f

	# df = stack(df, [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod])
	df = stack(df, [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod])

	df[:, :cluster] = ifelse.(occursin.("_mod", df.variable), "modified", "unmodified")
	df[:, :sample] = replace.(df.variable, "_mod" => "", "_unmod" => "")

	local samples = df.sample |> unique
	local sample_ids = Dict(zip(samples, 1:length(samples)))

	local f = Figure()
	local ax = Axis(f[1, 1],
				    xticks = (1:length(samples), samples))

	ylims!(ax, (0, 300))
	
	violin!(ax,
		    map(s -> sample_ids[s], df.sample),
		    df.value,
		    side = ifelse.(df.cluster .== "unmodified", :left, :right),
		    color = ifelse.(df.cluster .== "unmodified", :blue, :orange))
	
	f
end

# ╔═╡ fbd6502e-efa4-4e52-81db-3bb998314b14
an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&& an_sig_peaks_ivt.ref_id .== "ENST00000372586.4|ENSG00000175283.8|OTTHUMG00000020765.2|OTTHUMT00000054515.2|DOLK-201|DOLK|2074|protein_coding|" .&& an_sig_peaks_ivt.pos .== 1709, :]

# ╔═╡ bf177ece-0e74-45f2-8bb0-2b5a2ac2ddbe
begin
	local df = an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :]
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, [:WT_1_mod, :WT_1_unmod]] ./= sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, [:WT_2_mod, :WT_2_unmod]] ./= sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, [:WT_SPK_mod, :WT_SPK_unmod]] ./= sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	counts[:, [:IVT_1_mod, :IVT_1_unmod]] ./= sum.(eachrow(counts[:, [:IVT_1_mod, :IVT_1_unmod]]))
	counts[:, [:IVT_2_mod, :IVT_2_unmod]] ./= sum.(eachrow(counts[:, [:IVT_2_mod, :IVT_2_unmod]]))
	counts[:, [:IVT_3_mod, :IVT_3_unmod]] ./= sum.(eachrow(counts[:, [:IVT_3_mod, :IVT_3_unmod]]))

	# local sums = sum.(eachrow(counts))

	# counts ./ sums

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Per sample proportion of reads")

	f
	# transpose(Array(counts))
end

# ╔═╡ 8d2aafde-d745-4a0d-bcdb-d82a7877292f
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :])

	df[:, :mean_wt_mod] = sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	df[:, :glori] = df.source .== "HEK293T"
	df = sort(df, [:glori, :mean_wt_mod])

	local non_glori_count = sum(.! df.glori)
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, :WT_1] = counts.WT_1_mod ./ sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, :WT_2] = counts.WT_2_mod ./ sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, :WT_SPK] = counts.WT_SPK_mod ./ sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	counts[:, :IVT_1] = counts.IVT_1_mod ./ sum.(eachrow(counts[:, [:IVT_1_mod, :IVT_1_unmod]]))
	counts[:, :IVT_2] = counts.IVT_2_mod ./ sum.(eachrow(counts[:, [:IVT_2_mod, :IVT_2_unmod]]))
	counts[:, :IVT_3] = counts.IVT_3_mod ./ sum.(eachrow(counts[:, [:IVT_3_mod, :IVT_3_unmod]]))

	cols = [:WT_1, :WT_2, :WT_SPK, :IVT_1, :IVT_2, :IVT_3]
	counts = counts[:, cols]
	# local counts = df[:, cols]

	# local sums = sum.(eachrow(counts))


	# counts ./ sums

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
					yticks = ([900, non_glori_count, 2200, nrow(counts)],
							  ["DRACH", "$non_glori_count", "GLORI+", "$(nrow(counts))"]),
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	lines!(ax, [0.5, 6.5], [non_glori_count, non_glori_count], color = :red)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Proportion of reads modified")

	f
	# transpose(Array(counts))
end

# ╔═╡ 3d10cafd-f87c-40a7-ba15-f64e924cc686
md"""
## STORM/IVT
"""

# ╔═╡ 28779e98-1a0e-481f-bb10-c6973bd98976
# ╠═╡ disabled = true
#=╠═╡
storm_ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/STORM_IVT_eventalign_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
		LOR_threshold = 0.8)
  ╠═╡ =#

# ╔═╡ 40173831-c379-4e62-90fc-f7f2e55ff8c6
#=╠═╡
begin
	peaks_storm_ivt = peaks(storm_ivt, 4)
	sig_peaks_storm_ivt = peaks_storm_ivt[peaks_storm_ivt.predicted .!== missing .&&
						 	  peaks_storm_ivt.predicted .>= -log10(0.01), :]
end
  ╠═╡ =#

# ╔═╡ 95a3d0c3-09c8-490d-ae5c-993d8ec009e8
#=╠═╡
an_peaks_storm_ivt = annotate_peaks(peaks_storm_ivt, mod_ref, writer_motifs, writer_mods)
  ╠═╡ =#

# ╔═╡ 7bede161-f145-44f4-b6de-9fa00a72be52
#=╠═╡
an_sig_peaks_storm_ivt = annotate_peaks(sig_peaks_storm_ivt, mod_ref, writer_motifs, writer_mods)
  ╠═╡ =#

# ╔═╡ a927ecd1-b616-4479-88b3-a070e1cc95b8
#=╠═╡
begin
	local bed = an_sig_peaks_storm_ivt[:, [:chr, :genomicPos, :mod, :predicted, :strand]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, :mod] = coalesce.(bed.mod, ".")
	bed[:, [:chr, :genomicPos, :end, :mod, :predicted, :strand]] |> CSV.write("/home/mzdravkov/storm_ivt_mods.bed", delim = '\t', header = false)
end
  ╠═╡ =#

# ╔═╡ 25c51988-445b-4253-9763-22e1dd8ee210
#=╠═╡
begin
	local bed = an_sig_peaks_storm_ivt[:, [:chr, :genomicPos, :predicted]]
	bed[:, :end] = bed.genomicPos .+ 1
	bed[:, [:chr, :genomicPos, :end, :predicted]] |> CSV.write("/home/mzdravkov/storm_ivt_mods.bedGraph", delim = '\t', header = false)
end
  ╠═╡ =#

# ╔═╡ 4444b683-1571-4df8-9841-0e28df365dec
#=╠═╡
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :])

	local r = copy(an_peaks_storm_ivt)
	r[:, :range] = range.(r.pos .- 4, r.pos .+ 4)
	r = flatten(r, :range)

	local j = leftjoin(df, r,
			 		   on = [:ref_id, :pos => :range],
			 		   makeunique = true)

	(j[j.source .=== "HEK293T", :].GMM_chi2_qvalue_1 .<= 0.01) |> countmap
end
  ╠═╡ =#

# ╔═╡ 757fa776-c741-4eda-9e85-c1110231fcd4
#=╠═╡
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :])
	
	local refjoin = innerjoin(df, storm_ivt,
			  				  on = [:ref_id],
			  				  makeunique = true)
	refjoin = refjoin[between.(refjoin.genomicPos_1, refjoin.start, refjoin.end), :]

	local df_with_storm_ivt = combine(groupby(refjoin, :peak_id),
							   nrow => :overlaps,
							   :GMM_chi2_pvalue_1 => (p -> minimum(coalesce.(p, 1))) => :min_storm_ivt_pval)

	local j = leftjoin(df, df_with_storm_ivt,
		     		   on = [:peak_id])
	storm_resistant_storm_ivt = j
	storm_resistant_storm_ivt_sig = j[j.min_storm_ivt_pval .< 0.01, :peak_id] |> unique

	j[:, :glori] = j.source .=== "HEK293T"
	j[:, :storm_ivt_sig] = j.min_storm_ivt_pval .< 0.01
	combine(groupby(j, [:glori, :storm_ivt_sig]), nrow => :count)
	# (leftjoin(df[df.source .=== "HEK293T", :], df_with_storm_ivt,
	# 	      on = [:peak_id]).min_storm_ivt_pval .< 0.01) |> countmap
	# df_with_storm_ivt[df_with_storm_ivt.overlaps .== 9 .&& df_with_storm_ivt.min_storm_pval .> 0.2, :peak_id] |> unique |> length
end
  ╠═╡ =#

# ╔═╡ b8437ca7-b885-4dad-a72d-dc0b9c9a6be5


# ╔═╡ 90c789e2-999c-4e6c-8d49-0926ca080e71
#=╠═╡
begin
	local df = copy(an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id) .&&
		map(p -> p in storm_resistant_storm_ivt_sig, an_sig_peaks_ivt.peak_id), :])

	df[:, :mean_wt_mod] = sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(df[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	df[:, :glori] = df.source .== "HEK293T"
	df = sort(df, [:glori, :mean_wt_mod])

	local non_glori_count = sum(.! df.glori)
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod, :IVT_1_mod, :IVT_1_unmod, :IVT_2_mod, :IVT_2_unmod, :IVT_3_mod, :IVT_3_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, :WT_1] = counts.WT_1_mod ./ sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, :WT_2] = counts.WT_2_mod ./ sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, :WT_SPK] = counts.WT_SPK_mod ./ sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	counts[:, :IVT_1] = counts.IVT_1_mod ./ sum.(eachrow(counts[:, [:IVT_1_mod, :IVT_1_unmod]]))
	counts[:, :IVT_2] = counts.IVT_2_mod ./ sum.(eachrow(counts[:, [:IVT_2_mod, :IVT_2_unmod]]))
	counts[:, :IVT_3] = counts.IVT_3_mod ./ sum.(eachrow(counts[:, [:IVT_3_mod, :IVT_3_unmod]]))

	cols = [:WT_1, :WT_2, :WT_SPK, :IVT_1, :IVT_2, :IVT_3]
	counts = counts[:, cols]
	# local counts = df[:, cols]

	# local sums = sum.(eachrow(counts))


	# counts ./ sums

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
					yticks = ([900, non_glori_count, 2200, nrow(counts)],
							  ["DRACH", "$non_glori_count", "GLORI+", "$(nrow(counts))"]),
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	lines!(ax, [0.5, 6.5], [non_glori_count, non_glori_count], color = :red)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Proportion of reads modified")

	f
	# transpose(Array(counts))
end
  ╠═╡ =#

# ╔═╡ f1795d88-1e1f-4aa7-b245-b2f2cef7db17
begin
	local df = an_sig_peaks_ivt[map(p -> p in storm_res_peaks, an_sig_peaks_ivt.peak_id), :]
end

# ╔═╡ 967bcd01-03a6-4cc5-ac07-b997b65d8a8f
#=╠═╡
storm_resistant_storm_ivt
  ╠═╡ =#

# ╔═╡ 322240e8-a405-4ccf-8d94-3c52aa60eb26
glori_stm1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432635_GLORI-293T-STM20-1_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ e698c85d-ab27-4dda-9349-cf679b22a0f5
glori_stm2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432636_GLORI-293T-STM20-2_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ 4415b732-d854-4eb9-9ae6-ed872b087c3e
glori_stm = innerjoin(glori_stm1, glori_stm2,
		  		      on = [:Chr, :Strand, :Sites],
		  			  makeunique = true)

# ╔═╡ 8dd89e2e-78e3-4d2a-b2f8-c9d29e949752
glori_raw1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ c817df72-596c-41be-9d55-80dcbe73c485
glori_raw2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv"))

# ╔═╡ 76fbfe22-1223-454d-8f6c-0505cdd6c8c4
glori_raw = innerjoin(glori_raw1, glori_raw2,
					  on = [:Chr, :Strand, :Sites],
					  makeunique = true)

# ╔═╡ e33af351-2a7f-45fc-be7e-b36467b3e324
# ╠═╡ disabled = true
#=╠═╡
begin
	local t = innerjoin(storm_resistant_storm_ivt, glori_stm,
	          			on = [:chr => :Chr, :strand => :Strand])

	resistant_multi_glori = innerjoin(t[between.(t.Sites, t.start, t.end), :],
			  glori_raw,
			  on = [:chr => :Chr, :strand => :Strand, :Sites],
			  renamecols = "" => "_ctrl")
	resistant_multi_glori[:, :stm_mean_ratio] = (resistant_multi_glori.Ratio .+ resistant_multi_glori.Ratio_1)./2
	resistant_multi_glori[:, :ctrl_mean_ratio] = (resistant_multi_glori.Ratio_ctrl .+ resistant_multi_glori.Ratio_1_ctrl)./2
	resistant_multi_glori
end
  ╠═╡ =#

# ╔═╡ b2fabd16-f406-4b70-b5ab-343ff87b6a05
#=╠═╡
resistant_multi_glori[log2.(resistant_multi_glori.stm_mean_ratio ./ resistant_multi_glori.ctrl_mean_ratio) .> -0.2, :]
  ╠═╡ =#

# ╔═╡ 7b353ba1-cc33-405b-962e-5b7f82d1323b
#=╠═╡
storm_resistant_storm_ivt.glori |> countmap
  ╠═╡ =#

# ╔═╡ cb6f8399-addd-4095-98bf-915f1b38ee4d
#=╠═╡
(data(resistant_multi_glori) * mapping(:ctrl_mean_ratio, :stm_mean_ratio) * visual(Scatter) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; xlabel = "CTRL mean mod ratio", ylabel = "STM mean mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
  ╠═╡ =#

# ╔═╡ 3cbff0f7-eb89-49e5-88e4-43ec57ab87b3
begin
	local t = innerjoin(glori_raw, glori_stm,
		  				on = [:Chr, :Strand, :Sites],
		  				renamecols = "" => "_stm")
	t[:, :ctrl_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :stm_mean_ratio] = (t.Ratio_stm .+ t.Ratio_1_stm)./2
	t[:, :resistant] = log2.(t.stm_mean_ratio ./ t.ctrl_mean_ratio) .> -0.1
	glori_storm_resistant = t[t.resistant, :]
	println(sum(t.resistant))
	data(t) * mapping(:ctrl_mean_ratio, :stm_mean_ratio, color = :resistant) * visual(Scatter; markersize = 5) |> draw(; axis = (; xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 46ce6ac5-c128-46df-8d38-90776f42861b
innerjoin(glori_storm_resistant, ivt,
		  on = [:Chr => :chr, :Strand => :strand, :Sites => :genomicPos])

# ╔═╡ ca061644-6cf0-4547-bb0d-c21642f76f3f
begin
	local t = innerjoin(storm[storm.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))

	c = round(cor(t.glori_mean_ratio, t.pred_ratio), digits = 2)
	println(c)
	
	(data(t) * mapping(:glori_mean_ratio, :pred_ratio) * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Nanocompore (hard assignment) and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ fa354baf-9fcc-4022-828a-3458e4cdae0c
begin
	local t = innerjoin(storm[storm.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))

	t[t.glori_mean_ratio .> 0.8 .&& t.pred_ratio .< 0.3, :]
end

# ╔═╡ 68870781-fd02-4df4-b531-1ae27f53c0fd
function is_ratio_inverted(row)
	storm_counts = row.STORM_1_mod + row.STORM_2_mod + row.STORM_SPK_mod
	wt_counts = row.WT_1_mod + row.WT_2_mod + row.WT_SPK_mod

	storm_counts / (storm_counts + wt_counts) > (wt_counts) / (storm_counts + wt_counts)
end

# ╔═╡ be8d5828-880e-41dd-9195-95908bd3e829
storm_fix_clust = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_test_fix_mod_cluster_inferring/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
	    LOR_threshold = 0.8)

# ╔═╡ d888fa77-b67a-4654-8556-2c3dc2cd9f44
begin
	local t = innerjoin(storm_fix_clust[storm_fix_clust.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	# t[:, :pred_ratio] = map(r -> is_ratio_inverted(r) ? 1 - r.pred_ratio : r.pred_ratio, eachrow(t))

	local c = round(cor(t.glori_mean_ratio, t.pred_ratio), digits = 2)
	println(c)
	println(nrow(t))
	
	(data(t) * mapping(:glori_mean_ratio, :pred_ratio) * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Nanocompore (fixed) and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 8def7cc6-83ec-406c-b3d8-24c138d78199
# ╠═╡ disabled = true
#=╠═╡
storm_fix_clust_soft = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_test_fix_mod_cluster_inferring_soft_counts/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
	    LOR_threshold = 0.8)
  ╠═╡ =#

# ╔═╡ 167ea40a-7ab7-42a3-af60-8f2b9367bf5f
Makie.density(storm_fix_clust[storm_fix_clust.GMM_chi2_qvalue .< 0.01, :].GMM_LOR)

# ╔═╡ b747c8a5-6cb8-43ae-8c77-1c41260dbe0f
#=╠═╡
begin
	local t = innerjoin(storm_fix_clust_soft[storm_fix_clust_soft.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	# t[:, :pred_ratio] = map(r -> is_ratio_inverted(r) ? 1 - r.pred_ratio : r.pred_ratio, eachrow(t))

	local c = round(cor(t.glori_mean_ratio, t.pred_ratio), digits = 2)
	println(c)

	println(nrow(t))
	
	(data(t) * mapping(:glori_mean_ratio, :pred_ratio) * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Nanocompore (fixed, soft count) and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end
  ╠═╡ =#

# ╔═╡ 451907e7-eb22-4d8f-aa4e-44d23708861c
(storm_fix_clust.GMM_chi2_pvalue .=== missing) |> countmap

# ╔═╡ dd689831-a388-44df-9eed-97add91d37bf
#=╠═╡
(storm_fix_clust_soft.GMM_chi2_pvalue .=== missing) |> countmap
  ╠═╡ =#

# ╔═╡ d636baaf-b65f-4c5e-a9ed-775a7f1c9909
storm_soft = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_soft_clustering/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false,
	    LOR_threshold = 0.8)

# ╔═╡ 134e34a4-e467-49e6-9c39-12123401d1b7
begin
	local t = innerjoin(storm_soft[storm_soft.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	t[:, :pred_ratio] = map(r -> is_ratio_inverted(r) ? 1 - r.pred_ratio : r.pred_ratio, eachrow(t))

	local c = round(cor(t.glori_mean_ratio, t.pred_ratio), digits = 2)
	println(c)
	
	(data(t) * mapping(:glori_mean_ratio, :pred_ratio) * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Nanocompore (soft assignment) and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 19e82456-11a8-441c-b753-208aba90af10


# ╔═╡ 1bf15795-87e8-4cf7-92b1-9db7b102450a
#=╠═╡
begin
	local t = copy(storm_resistant_storm_ivt)
	# t[:, :wt_mod_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	
	t = innerjoin(t, storm_ivt,
		  	      on = [:ref_id, :pos],
			      renamecols = "" => "_stm")
	t[:, :wt_mod_ratio] = calc_mod_ratio(t, ["WT_1", "WT_2", "WT_SPK"])
	t[:, :stm_mod_ratio] = calc_mod_ratio(t, ["STORM_1", "STORM_2", "STORM_SPK"])
	t[(t.wt_mod_ratio .+ t.stm_mod_ratio)./2 .> 0.7 .&& log2.(t.stm_mod_ratio ./ t.wt_mod_ratio) .> -0.2, :]
end
  ╠═╡ =#

# ╔═╡ d9be422c-3630-4c86-859a-3df1cd725405
# function calc_mod_ratio_old(df, samples)
# 	sample_mod_cols = filter(n -> any(occursin.(samples, n)) && occursin("_mod", n), names(df))
# 	mod_cols = filter(n -> any(occursin("_mod", n)), names(df))
# 	sample_cols = filter(n -> any(occursin.(samples, n)), names(df))
	
# 	sample_mod_counts = sum.(eachrow(df[:, sample_mod_cols]))
# 	mod_counts = sum.(eachrow(df[:, mod_cols]))
# 	sample_counts = sum.(eachrow(df[:, sample_cols]))
	
# 	ratios = sample_mod_counts ./ sample_counts
# 	inverted = sample_mod_counts ./ mod_counts .< 0.5
# 	ifelse.(inverted, 1 .- ratios, ratios)
	
# 	# mod_cols = filter(n -> any(occursin.(samples, n)) && occursin("_mod", n), names(df))
# 	# all_cols = filter(n -> any(occursin.(samples, n)) && (occursin("_mod", n) || occursin("_unmod", n)), names(df)) 
# 	# sum.(eachrow(df[:, mod_cols])) ./ sum.(eachrow(df[:, all_cols]))
# end


function calc_mod_ratio(df, samples)
	sample_mod_cols = filter(n -> any(occursin.(samples, n)) && occursin("_mod", n),
							 names(df))
	sample_unmod_cols = filter(n -> any(occursin.(samples, n)) && occursin("_unmod", n),
							   names(df))
	dep_mod_cols = filter(n -> !any(occursin.(samples, n)) && occursin("_mod", n),
						  names(df))
	dep_unmod_cols = filter(n -> !any(occursin.(samples, n)) && occursin("_unmod", n),
							names(df))
	
	sample_mod = sum.(eachrow(df[:, sample_mod_cols]))
	sample_unmod = sum.(eachrow(df[:, sample_unmod_cols]))
	dep_mod = sum.(eachrow(df[:, dep_mod_cols]))
	dep_unmod = sum.(eachrow(df[:, dep_unmod_cols]))

	ifelse.(sample_mod .* dep_unmod .> dep_mod .* sample_unmod,
		    sample_mod,
		    sample_unmod) ./ (sample_mod .+ sample_unmod) 
end

# ╔═╡ fbea4e5f-4c3a-4492-97cb-d2eeeb4d1c63
begin
	local t = innerjoin(storm[storm.predicted .> 2, :], glori_raw,
						on = [:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	t[:, :glori_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	# t[:, :pred_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	# t[:, :pred_ratio] = map(r -> is_ratio_inverted(r) ? 1 - r.pred_ratio : r.pred_ratio, eachrow(t))
	t[:, :pred_ratio] = calc_mod_ratio(t, ["WT_1", "WT_2", "WT_SPK"])

	local c = round(cor(t.glori_mean_ratio, t.pred_ratio), digits = 2)
	println(c)
	
	(data(t) * mapping(:glori_mean_ratio, :pred_ratio) * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio corr. between Nanocompore (hard assignment, fixed) and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

end

# ╔═╡ 33cc218c-f02f-4542-ba7a-b6075d7e4357
# ╠═╡ disabled = true
#=╠═╡
begin
	# local wtstm = innerjoin(storm, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
	# 					    on = [:ref_id, :pos])
	local df = copy(storm[storm.predicted .> 2, :])
	df[:, :wt_mod_ratio] = calc_mod_ratio(df, ["WT_1", "WT_2", "WT_SPK"])
	df[:, :wt_mod_ratio2] = calc_mod_ratio2(df, ["WT_1", "WT_2", "WT_SPK"])
	
	data(df) * mapping(:wt_mod_ratio => "WT mod. ratio (WT/STORM)", :wt_mod_ratio2 => "WT mod. ratio 2 (WT/STORM)") * visual(Scatter, markersize = 5) |> draw
end
  ╠═╡ =#

# ╔═╡ 10a0495d-b861-4959-a5c2-1d74b6a0f5c6
# ╠═╡ disabled = true
#=╠═╡
begin
	local t = copy(an_sig_peaks_ivt)
	
	
	t = innerjoin(t, storm_ivt,
			      on = [:ref_id, :pos],
			      renamecols = "" => "_stm")

	t[:, :mod_ratio] = sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod]])) ./ sum.(eachrow(t[:, [:WT_1_mod, :WT_2_mod, :WT_SPK_mod, :WT_1_unmod, :WT_2_unmod, :WT_SPK_unmod]]))
	t[:, :mod_ratio_stm] = sum.(eachrow(t[:, [:STORM_1_mod_stm, :STORM_2_mod_stm, :STORM_SPK_mod_stm]])) ./ sum.(eachrow(t[:, [:STORM_1_mod_stm, :STORM_2_mod_stm, :STORM_SPK_mod_stm, :STORM_1_unmod_stm, :STORM_2_unmod_stm, :STORM_SPK_unmod_stm]]))

	t
	# t = innerjoin(t, glori_raw,
	#           	  on = [:chr => :Chr, :strand => :Strand])

	# t = t[between.(t.Sites, t.start, t.end), :]
end
  ╠═╡ =#

# ╔═╡ 0e7cfe6c-a700-4ec1-b8dc-a67f246782b3
begin
	local wtstm = copy(storm[storm.predicted .> 2, :])
	wtstm[:, :wt_mod_ratio] = calc_mod_ratio(wtstm, ["WT_1", "WT_2", "WT_SPK"])
	local wtivt = copy(ivt[ivt.predicted .> 2, :])
	wtivt[:, :wt_mod_ratio] = calc_mod_ratio(wtivt, ["WT_1", "WT_2", "WT_SPK"])
	local t = innerjoin(wtstm, wtivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")
	
	# t = innerjoin(t, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
				  # on = [:ref_id, :pos])
	
	# data(t) * mapping(:wt_mod_ratio => "WT mod. ratio (WT/STORM)", :wt_mod_ratio_ivt => "WT mod. ratio (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
	data(t) * mapping(:GMM_LOR => abs => "|LOR| (WT/STORM)", :GMM_LOR_ivt => abs => "|LOR| (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
end

# ╔═╡ 3b6c3b40-d54e-41c2-afac-6bf8234ae648
begin
	# local wtstm = innerjoin(storm, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
	# 					    on = [:ref_id, :pos])
	local wtstm = copy(storm[storm.predicted .> 2, :])
	wtstm[:, :wt_mod_ratio] = calc_mod_ratio(wtstm, ["WT_1", "WT_2", "WT_SPK"])
	local wtivt = copy(ivt[ivt.predicted .>= 2, :])
	wtivt[:, :wt_mod_ratio] = calc_mod_ratio(wtivt, ["WT_1", "WT_2", "WT_SPK"])
	local t = innerjoin(wtstm, wtivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")

	data(t) * mapping(:wt_mod_ratio => "WT mod. ratio (WT/STORM)", :wt_mod_ratio_ivt => "WT mod. ratio (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
end

# ╔═╡ f97950ff-8c2b-4fc4-99e2-51c29bf7b489
# ╠═╡ disabled = true
#=╠═╡
begin
	local wtstm = innerjoin(storm, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
						    on = [:ref_id, :pos])
	wtstm[:, :wt_mod_ratio] = calc_mod_ratio(wtstm, ["WT_1", "WT_2", "WT_SPK"])
	local wtivt = copy(ivt[ivt.predicted .>= 0, :])
	wtivt[:, :wt_mod_ratio] = calc_mod_ratio(wtivt, ["WT_1", "WT_2", "WT_SPK"])
	local t = innerjoin(wtstm, wtivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")

	data(t) * mapping(:GMM_LOR => abs => "|LOR| (WT/STORM)", :GMM_LOR_ivt => abs => "|LOR| (WT/IVT)") * visual(Scatter, markersize = 5) |> draw
end
  ╠═╡ =#

# ╔═╡ b9ee5aa2-ced4-4f71-b43a-bbfd6d5afef9
# ╠═╡ disabled = true
#=╠═╡
begin
	local wtstm = innerjoin(storm, storm_resistant_storm_ivt[:, [:ref_id, :pos]],
						    on = [:ref_id, :pos])
	wtstm[:, :storm_mod_ratio] = calc_mod_ratio(wtstm, ["STORM_1", "STORM_2", "STORM_SPK"])
	local stmivt = copy(storm_ivt[storm_ivt.predicted .>= 0, :])
	stmivt[:, :storm_mod_ratio] = calc_mod_ratio(stmivt, ["STORM_1", "STORM_2", "STORM_SPK"])
	local t = innerjoin(wtstm, stmivt,
		      	        on = [:ref_id, :pos],
		 	  			renamecols = "" => "_ivt")

	data(t) * mapping(:storm_mod_ratio, :storm_mod_ratio_ivt) * visual(Scatter, markersize = 5) |> draw
	
end
  ╠═╡ =#

# ╔═╡ bd94a0c6-0e84-4699-b449-36e45b2c2c23
#=╠═╡
storm_resistant_storm_ivt
  ╠═╡ =#

# ╔═╡ b79d6323-c6be-4b18-ae93-1a1ce48a1db3
function get_read_probs(ref)
	
	db = SQLite.DB("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_read_level_test_all/out_sampComp_sql.db")
	
	# local df = sig_comod[sig_comod.reference .== ref, :]
	# local positions = peaks_ivt[peaks_ivt.ref_id .== ref, :pos]
	positions = unique(union(comod[comod.reference .== ref, :pos1], comod[comod.reference .== ref, :pos2]))
	# positions = [2471, 2488]


	reads = SQLite.DBInterface.execute(db, "SELECT * FROM read_results r INNER JOIN transcripts t ON t.id = r.transcript_id WHERE t.name = \"$ref\"") |> DataFrame
	
	# for i in 1:nrow(reads)
	# 	reinterpret.(Int8, reads.mod_probs[i])
	# end
	reads[:, :probs] = map(probs -> Float32.(reinterpret.(Int8, probs)), reads.mod_probs)

	orig_probs = reduce(hcat, Array.(map(probs -> probs[positions], reads.probs)))
	probs = ifelse.(orig_probs .== -1, missing, orig_probs)
	probs ./= 100
end

# ╔═╡ ca38417c-1d5c-45c6-816c-d2d6898af7db
glori_storm_resistant

# ╔═╡ 3bfb00e8-f963-4c03-825f-71c4343074fa
begin
	eclip_files = ["GSM4970401_meCLIP_HEK293_rep1.m6aList_highConfidence_3054.bed",
				   "GSM4970403_meCLIP_HEK293_rep2.m6aList_highConfidence_1775.bed",
				   "GSM4970405_meCLIP_HEK293_rep3.m6aList_highConfidence_2911.bed"]

	eclip_reps = []
	for fn in eclip_files
		df = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data//GSE147440_RAW/" * fn, header = ["chr", "start", "end", "name", "data", "strand"]))
		df = DataFrames.transform(df, :data => ByRow(x -> map(v -> parse(Float64, v), split(x, "|"))) => [:mod_ratio, :nmod, :nreads])
		df = select!(df, Not(:data))
		push!(eclip_reps, df)
	end
	eclip_reps
end

# ╔═╡ 1c622f87-b763-4f8b-a253-fd38a1fd2830
begin
	local c = copy(sig_peaks_storm)
	c[:, :start] = c.genomicPos .- 4
	c[:, :end] = c.genomicPos .+ 4
	local t = innerjoin(c, eclip_reps[1],
		  				on = [:chr, :strand],
		  				makeunique = true)
	t[between.(t.start_1, t.start, t.end), :]
end

# ╔═╡ 45a5a8f7-a1bb-434c-b882-3b24d3efa412
#=╠═╡
begin
	local t = innerjoin(storm_resistant_storm_ivt, eclip_reps[1],
		  				on = [:chr, :strand],
		  				makeunique = true)
	t[between.(t.start_1, t.start, t.end), :]
end
  ╠═╡ =#

# ╔═╡ 96482819-8af4-4ce7-bfda-2e20ac05a66c
begin
	local encore_experiments = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/ENCORE_m6A/experiment_report_2025_7_1_10h_33m.tsv", header = 2))
	encore_experiments = innerjoin(encore_experiments, DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/ENCORE_m6A/ENCORE_m6a_associated_RBPs_filtered.tsv")),
								   on = [:Accession => :accession])
	encore_experiments[:, :bed] = [string(last(split(u, "/"))) 
								   for u in encore_experiments.url]
	local encore_m6a_rbps = []
	for experiment in eachrow(encore_experiments)
		bed = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/ENCORE_m6A/$(experiment.bed)", header = [:chr, :start, :end, :name, :x, :strand, :a, :b, :c, :d], select = [:chr, :start, :end, :name, :strand]))
		bed[:, :rbp] .= experiment.rbp
		bed[:, :cell] .= experiment["Biosample term name"]
		push!(encore_m6a_rbps, bed)
	end
	encore_m6a_rbp_peaks = reduce(vcat, encore_m6a_rbps)
end

# ╔═╡ 6a70233a-4754-49e8-bd81-292d364b4031
#=╠═╡
begin
	local t = innerjoin(storm_resistant_storm_ivt, encore_m6a_rbp_peaks,
		  				on = [:chr, :strand],
		  				makeunique = true)
	combine(groupby(sort(t[between.(t.genomicPos, t.start_1, t.end_1), :], [:ref_id, :pos]), [:chr, :strand, :genomicPos, :glori]), nrow => :count)
end
  ╠═╡ =#

# ╔═╡ 445fe109-2361-4bea-bd03-74f3f12dca65
#=╠═╡
begin
	local t = innerjoin(storm_resistant_storm_ivt,
			  encore_m6a_rbp_peaks,
		  	  on = [:chr, :strand],
		  	  makeunique = true)
	eclip_overlaps = t[between.(t.genomicPos, t.start_1, t.end_1), :peak_id] |> unique

	conf_storm_resistant = copy(storm_resistant_storm_ivt)
	conf_storm_resistant[:, :overlaps_eclip] = map(p -> p in eclip_overlaps, storm_resistant_storm_ivt.peak_id)

	local conf = r -> r.storm_ivt_sig && (r.glori || r.overlaps_eclip)
	conf_storm_resistant = filter(conf, conf_storm_resistant)

	conf_storm_resistant = innerjoin(conf_storm_resistant, storm_ivt,
		 				 			 on = [:ref_id, :pos],
		  						     renamecols = "" => "_stmivt")
	conf_storm_resistant[:, :wt_mod_ratio] = calc_mod_ratio(conf_storm_resistant, ["WT_1", "WT_2", "WT_SPK"])
	conf_storm_resistant[:, :stm_mod_ratio] = calc_mod_ratio(conf_storm_resistant, ["STORM_1", "STORM_2", "STORM_SPK"])
	conf_storm_resistant
end
  ╠═╡ =#

# ╔═╡ f7a4a1b1-5970-4e1c-a37a-5db901c1115c
#=╠═╡
combine(groupby(conf_storm_resistant, [:glori, :overlaps_eclip]), nrow => :count)
  ╠═╡ =#

# ╔═╡ 4605495f-6df1-44fb-8a33-6e84bebe83df
#=╠═╡
data(conf_storm_resistant) * mapping(:wt_mod_ratio => "WT mod. ratio", :stm_mod_ratio => "STORM mod. ratio") * visual(Scatter) |> draw(; axis = (; xticks = 0:0.1:1, yticks = 0:0.1:1, limits = ((0, 1.1), (0, 1.1))))
  ╠═╡ =#

# ╔═╡ 90818cac-2f00-4fc8-9225-321eeac107e7


# ╔═╡ 38f8fbed-20a9-46ab-92e9-edf3cb693bfa
#=╠═╡
conf_storm_resistant
  ╠═╡ =#

# ╔═╡ f84bfc2a-c037-4627-b0c6-4f47d4f45a1b
#=╠═╡
(
	data(DataFrame(x = [0, 5000], y = [0, 0])) * mapping(:x, :y) * visual(Lines; color = :red) +
	data(conf_storm_resistant) * mapping((:mean_cov, :mean_cov_stmivt) => ((a, b) -> (a+b)/4), (:stm_mod_ratio, :wt_mod_ratio) => ((s, w) -> 100*(s - w))) * visual(Scatter)
	
) |> draw(; axis = (xticks = (0:500:5000, [t % 1000 == 0 ? "$t" : "" for t in 0:500:5000]), yticks = -80:10:80, ytickformat = "{}%", xlabel = "Mean coverage (per condition)", ylabel = "STORM stoichiometry - WT stoichiometry"))
  ╠═╡ =#

# ╔═╡ 71261e27-8c1b-47e0-8fdd-68fc1cecee39
#=╠═╡
(
	data(DataFrame(x = [0, 5000], y = [0, 0])) * mapping(:x, :y) * visual(Lines; color = :red) +
	data(conf_storm_resistant) * mapping((:mean_cov, :mean_cov_stmivt) => ((a, b) -> (a+b)/4), (:stm_mod_ratio, :wt_mod_ratio) => ((s, w) -> log2(s/w))) * visual(Scatter)
	
) |> draw(; axis = (xticks = (0:500:5000, [t % 1000 == 0 ? "$t" : "" for t in 0:500:5000]), yticks = -4:1:4,xlabel = "Mean coverage (per condition)", ylabel = "log₂(STORM mod. ratio / WT mod. ratio)", limits = (nothing, (-4, 4))))
  ╠═╡ =#

# ╔═╡ d836b201-4413-412a-8ec2-be15d23504c1
#=╠═╡
(
	data(DataFrame(x = [0, 5], y = [0, 0])) * mapping(:x, :y) * visual(Lines; color = :red) +
	data(conf_storm_resistant) * mapping((:mean_cov, :mean_cov_stmivt) => ((a, b) -> log10((a+b)/4)), (:stm_mod_ratio, :wt_mod_ratio) => ((s, w) -> log2(s/w))) * visual(Scatter)
	
) |> draw(; axis = (xticks = 0:1:5, yticks = -4:1:4, xlabel = "log₁₀(Mean coverage per condition)", ylabel = "log₂(STORM mod. ratio / WT mod. ratio)", limits = (nothing, (-4, 4))))
  ╠═╡ =#

# ╔═╡ 10a15475-c736-48f7-a4a4-cf972134a0f0
#=╠═╡
data(conf_storm_resistant) * mapping(:wt_mod_ratio, (:stm_mod_ratio, :wt_mod_ratio) => ((s, w) -> log2(s/w)) => "log₂(STORM mod. ratio / WT mod. ratio)") * visual(Scatter) |> draw
  ╠═╡ =#

# ╔═╡ a0798ff2-c72f-4de5-b31a-3c8a3bacf117
#=╠═╡
begin
	local breast_cancer_genes = Set(split("MDC1;SETD2;HSP90AB1;CITED2;EHMT2;FHL1;OLA1;HDLBP;HMGB2;CASC3;RPL8;PHB2;PTPRF;LAPTM4B;CHAF1A;SDF4;DHX16;GLUL;RBM5;PELP1;SRRM2;BRD2;ACTR2;HSP90AA1;SEMA6A;CSNK2A1;HLA-C;RPL13A;DDX54;LSS;SREBF2;WDR77;NME1;RHOB;CDC25B;ILF3;TCP1;RAB35;MYH9;MCM5;PTMA;GAPDH;CREB5;NOTCH2;KDM3B;HM13;GMPR2;MAZ;RPL36A;SRI;HOXB13;HNRNPDL;DVL1;TP53BP2;ABL1;SSR1;FLNA;PRAME;CSK;MAP2K7;SF3B1;RPL19;PAK4;JUN;SLC35A4;CDKN2C;XRCC5;MRPL28;NR2F2;RRP1B;DEK;FKBP1A;SQLE;HNRNPM;SLC6A8;HNRNPK;KRT18;SCD;FASN;HNRNPD;ZYX;GNAS;MAP3K10;TPT1;EIF3A", ";"))

	local t = copy(conf_storm_resistant)
	t[:, :gene] = map(r -> string(split(r, "|")[6]), t.ref_id)
	t[:, :breast_cancer] = map(g -> g in breast_cancer_genes, t.gene)
	t

	(
	data(DataFrame(x = [0, 5000], y = [0, 0])) * mapping(:x, :y) * visual(Lines; color = :red) +
	data(t) * mapping((:mean_cov, :mean_cov_stmivt) => ((a, b) -> (a+b)/4), (:stm_mod_ratio, :wt_mod_ratio) => ((s, w) -> 100*(s - w)), color = :breast_cancer) * visual(Scatter)
	
) |> draw(; axis = (xticks = (0:500:5000, [t % 1000 == 0 ? "$t" : "" for t in 0:500:5000]), yticks = -80:10:80, ytickformat = "{}%", xlabel = "Mean coverage (per condition)", ylabel = "STORM stoichiometry - WT stoichiometry"))
end
  ╠═╡ =#

# ╔═╡ 007e6893-0718-41c0-8ca0-03fd367f6cc7
#=╠═╡
begin
	local t = copy(conf_storm_resistant)
	t[:, :len] = map(r -> parse(Int64, split(r, "|")[7]), t.ref_id)
	local transcripts = map(r -> string(split(r, "|")[1]), t.ref_id)
	local all_lens = map(t -> calculate_utr_cds_lengths(t, conf_storm_resistant_gencode),
					 transcripts)
	# t[:, :utr5_len] = first.(lens)
	# t[:, :cds_len] = second.(lens)
	# t[:, :utr3_len] = last.(lens)

	t[:, :meta_pos] = [
		if pos <= utr5_len
			(pos / utr5_len)/4
		elseif pos <= utr5_len + cds_len
			1/4 + ((pos - utr5_len) / cds_len)/2
		else
			3/4 + ((pos - utr5_len - cds_len) / utr3_len)/4
		end
		for (pos, (utr5_len, cds_len, utr3_len)) in zip(t.pos, all_lens)
	]

	t = t[.! isinf.(t.meta_pos), :]
	# t

	# all_lens[5], t[5, :]

	(
		data(DataFrame(x = [1/4, 1/4], y = [0, 12])) * mapping(:x, :y) * visual(Lines, color = :red, linestyle = :dash) +
		data(DataFrame(x = [3/4, 3/4], y = [0, 12])) * mapping(:x, :y) * visual(Lines, color = :red, linestyle = :dash) +
		data(t) * mapping(:meta_pos) * AlgebraOfGraphics.histogram(bins = 50)
	)|> draw(; axis = (; xticks = ([0, 1/4, 3/4, 1], ["TSS", "Start codon", "Stop codon", "TTS"]), ylabel = "Count", yticks = 1:12))
end
  ╠═╡ =#

# ╔═╡ 49252d74-9a99-4795-9527-3e7eccf20eed
#=╠═╡
begin
	local refs = Set(unique(conf_storm_resistant.ref_id))

	local resist = combine(groupby(conf_storm_resistant, :ref_id), nrow => :res_mods)
	resist = Dict(zip(resist.ref_id, resist.res_mods))

	local t = copy(sig_peaks_storm[map(r -> r in refs, sig_peaks_storm.ref_id), :])

	local sensitive = combine(groupby(t, :ref_id), nrow => :sensitive_mods)
	sensitive = Dict(zip(sensitive.ref_id, sensitive.sensitive_mods))

	refs = [r for r in refs]
	local counts = DataFrame(ref = refs,
		 		  			 nresist = [get(resist, r, 0) for r in refs],
			  	   			 nsensitive = [get(sensitive, r, 0) for r in refs])

	data(counts) * mapping(:nsensitive => "Number of STORM-sensitive sites") * frequency() |> draw(; axis = (; title = "Transcripts with STM2457-resistant m6A sites by number of sensitive sites", ylabel = "Number of transcripts"))
	
	# combine(groupby(counts, [:nresist, :nsensitive]), nrow => :count)
	# sort(counts, [:nresist, :nsensitive])
	
end
  ╠═╡ =#

# ╔═╡ e9788e02-f17f-4bd8-b344-caa18863e92e
#=╠═╡
conf_storm_resistant_20perc = conf_storm_resistant[abs.(conf_storm_resistant.stm_mod_ratio .- conf_storm_resistant.wt_mod_ratio) .< 0.2, :]
  ╠═╡ =#

# ╔═╡ 3266d3d6-64a2-4936-91ac-059a1ccf76c4
import GenomicAnnotations

# ╔═╡ da65e2fe-54a8-4ca6-b599-411ee2a3e9d4
conf_storm_resistant_gencode = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/gencode41_annotation_for_conf_storm_resistant.gtf", header = [:chr, :source, :feature, :start, :end, :score, :strand, :frame, :attr]))

# ╔═╡ aa699be6-0fa9-468d-847c-b7b3b7f91c72
function calculate_utr_cds_lengths(transcript_id::String, gtf::DataFrame)
	transcript_rows = gtf[occursin.(transcript_id, gtf.attr), :]

	rev = transcript_rows[1, :strand] == '-'

	sort!(transcript_rows, [:start], rev = rev)

    five_prime_utr_len = 0
    cds_len = 0
    three_prime_utr_len = 0

	utr_5p = true

    for row in eachrow(transcript_rows)
        feature_type = row.feature
        start = row.start
        end_ = row.end

        if feature_type == "UTR" && utr_5p
            five_prime_utr_len += end_ - start + 1
			utr_5p = false
        elseif feature_type == "CDS"
            cds_len += end_ - start + 1
        elseif feature_type == "UTR"
            three_prime_utr_len += end_ - start + 1
        end
    end

    return (five_prime_utr_len, cds_len, three_prime_utr_len)
end

# ╔═╡ 7212801b-ee05-4699-b4fb-c02377710e10
calculate_utr_cds_lengths("ENST00000470886.1", conf_storm_resistant_gencode)

# ╔═╡ 6a423afd-67f4-4152-9aa8-7d69c4262f63


# ╔═╡ c02b568d-e607-4729-9cef-1d189ae22bc1
#=╠═╡
begin
	local filepath = "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_sites.tsv"
	local df = copy(storm_resistant_storm_ivt)
	rename!(df,
			:mean_cov => :coverage,
		    :start => :peak_start,
		    :end => :peak_end,
		    :overlaps => :storm_ivt_overlaps,
		    :min_storm_ivt_pval => :storm_ivt_min_pvalue)
	sort!(df, [:source, :ref_id, :pos])
	

	local f = open(filepath, "w")
	# local data = read(f);
	# seekstart(f);
	local comment = """
	# A table with putative STM2457-resistant sites.
	# These are peaks that were significant when comparing the IVT and WT with Nanocompore, have been inferred as m6A (because they've been detected with GLORI or due to overlapping a DRACH motif), but were not significant when comparing WT to STORM.
	# The first columns are the results from Nanocompore for the IVT/WT comparison, a few columns have been added at the end. Specifically, whether the peak was also significant when comparing the STORM and IVT conditions.
	# More technically:
	# 1. We took the IVT/WT Nonocompore results.
	# 2. We peak called with min distance between peaks equal to 5nt.
	# 3. We selected only the significant peaks at q-value < 0.01 and |GMM_LOR| > 0.8.
	# 4. We selected only the sites where the whole peak (the 9-mer centered at the peak's position) is covered by results in the Nanocompore's comparison of STORM/WT and none of them is significant even at non-strict thresholds of the non-multiple-test-corrected pvalues (GMM_chi2_pvalue > 0.2 && KS_intensity_pvalue > 0.1 && KS_dwell_pvalue 0.1).
	# 5. We combined the resulting peaks with information from Nanocompore's results from the STORM/IVT comparison. We marked the peak as storm_ivt_sig if there was at least one position overlapping the peak that was present in the STORM/IVT results with a GMM_chi2_pvalue < 0.01.
	"""
	write(f, comment)
	close(f)
	df |> CSV.write(filepath, delim = '\t', append = true, writeheader = true)
end
  ╠═╡ =#

# ╔═╡ 7aa12ba0-e6e7-4c36-9620-7ed443013d40
begin
	rna_stability = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/pelizzola_2016_rna_stability_hek293t.csv"))
	rna_stability[:, :halflife] = (log(2) ./ rna_stability.Deg) ./ 60
	rna_stability
end

# ╔═╡ a11ab145-e6d2-4452-8df8-6976f5f67453
#=╠═╡
begin
	local t = copy(conf_storm_resistant)
	local s = copy(rna_stability)
	s[:, :gene] = map(r -> string(split(r, ".")[1]), s.Gene)
	t[:, :gene] = map(r -> string(split(split(r, "|")[2], ".")[1]), t.ref_id)
	local j = innerjoin(unique(t[:, [:gene, :mod]]), s[:, [:gene, :halflife]],
			  			on = [:gene])
	println(quantile(j.halflife), mean(j.halflife))
	data(j) * mapping(:halflife => round => "Half-life (h)") * frequency() |> draw(; axis = (; ylabel = "Transcripts"))
end
  ╠═╡ =#

# ╔═╡ 046608c2-5c79-4651-bcfe-1b8b965e63c2
#=╠═╡
begin
	local t = copy(conf_storm_resistant)
	local s = copy(rna_stability)
	s[:, :gene] = map(r -> string(split(r, ".")[1]), s.Gene)
	t[:, :gene] = map(r -> string(split(split(r, "|")[2], ".")[1]), t.ref_id)
	local j = unique(leftjoin(s, t[:, [:gene, :mod]], on = :gene))
	j[:, :mod] = ifelse.(j.mod .=== missing, "No STORM-insensitive m6As", "STORM-insensitive m6As")
	data(j) * mapping(:mod => "", :halflife => "Half-life (h)") * visual(Violin) |> draw(; axis = (; yticks = 0:5:30))
end
  ╠═╡ =#

# ╔═╡ e2bf25a2-b5c2-48be-95b3-9ee2905b7b35
#=╠═╡
begin
	local t = copy(conf_storm_resistant)
	local s = copy(rna_stability)
	s[:, :gene] = map(r -> string(split(r, ".")[1]), s.Gene)
	t[:, :gene] = map(r -> string(split(split(r, "|")[2], ".")[1]), t.ref_id)
	local j = innerjoin(unique(t[:, [:gene]]), s[:, [:gene, :Deg]],
			  			on = [:gene])
	local concentration = (k, t) -> exp(-1*k*t)
	local ngenes = nrow(j)
	local concentrations = combine(
			groupby(j, :gene),
			:gene => (_ -> 0:24) => :t_h,
			:Deg => (k -> concentration.(k, 0:60:24*60)) => :concentration)

	data(combine(groupby(concentrations, :t_h), :concentration => mean => :mean_concentration)) * mapping(:t_h => "Time (h)", :mean_concentration => "Mean concentration") * visual(Lines) |> draw(; axis = (; xticks = 0:4:24))
end
  ╠═╡ =#

# ╔═╡ 5a1dfded-7f69-4fae-ad5c-9242686026d8
#=╠═╡
begin
	local t = innerjoin(conf_storm_resistant, encore_m6a_rbp_peaks,
						on = [:chr, :strand],
						renamecols = "" => "_rbp")
	t = t[between.(t.genomicPos, t.start_rbp, t.end_rbp), :]
	data(combine(groupby(t, [:peak_id, :rbp_rbp]), nrow => :count)) * mapping(:rbp_rbp => "RBP") * frequency() |> draw(; axis = (; xticklabelsize = 10))
end
  ╠═╡ =#

# ╔═╡ 77dbbadb-2588-4db7-8e84-5b375dbec467
#=╠═╡
conf_storm_resistant
  ╠═╡ =#

# ╔═╡ ce3be394-4d52-40c6-9df2-4cad1ca55fd7
#=╠═╡
begin
	local df = copy(conf_storm_resistant)

	df[:, :mean_stm_mod] = sum.(eachrow(df[:, [:STORM_1_mod_stmivt, :STORM_2_mod_stmivt, :STORM_SPK_mod_stmivt]])) ./ sum.(eachrow(df[:, [:STORM_1_mod_stmivt, :STORM_2_mod_stmivt, :STORM_SPK_mod_stmivt, :STORM_1_unmod_stmivt, :STORM_2_unmod_stmivt, :STORM_SPK_unmod_stmivt]]))
	df[:, :glori] = df.source .== "HEK293T"
	df = sort(df, [:glori, :mean_stm_mod])

	local non_glori_count = sum(.! df.glori)
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:STORM_1_mod_stmivt, :STORM_1_unmod_stmivt, :STORM_2_mod_stmivt, :STORM_2_unmod_stmivt, :STORM_SPK_mod_stmivt, :STORM_SPK_unmod_stmivt, :WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, :STORM_1] = counts.STORM_1_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_1_mod_stmivt, :STORM_1_unmod_stmivt]]))
	counts[:, :STORM_2] = counts.STORM_2_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_2_mod_stmivt, :STORM_2_unmod_stmivt]]))
	counts[:, :STORM_3] = counts.STORM_SPK_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_SPK_mod_stmivt, :STORM_SPK_unmod_stmivt]]))

	counts[:, :WT_1] = counts.WT_1_mod ./ sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, :WT_2] = counts.WT_2_mod ./ sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, :WT_3] = counts.WT_SPK_mod ./ sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	cols = [:STORM_1, :STORM_2, :STORM_3, :WT_1, :WT_2, :WT_3]
	counts = counts[:, cols]
	# local counts = df[:, cols]

	# local sums = sum.(eachrow(counts))


	# counts ./ sums

	println(cor(counts.STORM_1, counts.STORM_2))
	println(cor(counts.STORM_1, counts.STORM_3))
	println(cor(counts.STORM_2, counts.STORM_3))

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
					yticks = ([900, non_glori_count, 2200, nrow(counts)],
							  ["DRACH", "$non_glori_count", "GLORI+", "$(nrow(counts))"]),
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	lines!(ax, [0.5, 6.5], [non_glori_count, non_glori_count], color = :red)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Proportion of reads modified")

	f
	# transpose(Array(counts))
end
  ╠═╡ =#

# ╔═╡ 90cffe20-514b-4bac-8e6f-04e30bc80f67
#=╠═╡
begin
	local df = copy(conf_storm_resistant)

	# df[:, :mean_stm_mod] = sum.(eachrow(df[:, [:STORM_1_mod_stmivt, :STORM_2_mod_stmivt, :STORM_SPK_mod_stmivt]])) ./ sum.(eachrow(df[:, [:STORM_1_mod_stmivt, :STORM_2_mod_stmivt, :STORM_SPK_mod_stmivt, :STORM_1_unmod_stmivt, :STORM_2_unmod_stmivt, :STORM_SPK_unmod_stmivt]]))
	df[:, :glori] = df.source .== "HEK293T"
	df = sort(df, [:glori, :stm_mod_ratio])

	local non_glori_count = sum(.! df.glori)
	# local cols = [:WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local cols = [:STORM_1_mod_stmivt, :STORM_1_unmod_stmivt, :STORM_2_mod_stmivt, :STORM_2_unmod_stmivt, :STORM_SPK_mod_stmivt, :STORM_SPK_unmod_stmivt, :WT_1_mod, :WT_1_unmod, :WT_2_mod, :WT_2_unmod, :WT_SPK_mod, :WT_SPK_unmod]
	local counts = df[:, cols] .+ 0.0

	counts[:, :STORM_1] = calc_mod_ratio(df, ["STORM_1"])
		# counts.STORM_1_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_1_mod_stmivt, :STORM_1_unmod_stmivt]]))
	counts[:, :STORM_2] = calc_mod_ratio(df, ["STORM_2"])
		#counts.STORM_2_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_2_mod_stmivt, :STORM_2_unmod_stmivt]]))
	counts[:, :STORM_3] = calc_mod_ratio(df, ["STORM_SPK"])
	#counts.STORM_SPK_mod_stmivt ./ sum.(eachrow(counts[:, [:STORM_SPK_mod_stmivt, :STORM_SPK_unmod_stmivt]]))

	counts[:, :WT_1] = calc_mod_ratio(df, ["WT_1"])
		# counts.WT_1_mod ./ sum.(eachrow(counts[:, [:WT_1_mod, :WT_1_unmod]]))
	counts[:, :WT_2] = calc_mod_ratio(df, ["WT_2"])
		# counts.WT_2_mod ./ sum.(eachrow(counts[:, [:WT_2_mod, :WT_2_unmod]]))
	counts[:, :WT_3] = calc_mod_ratio(df, ["WT_SPK"])
		# counts.WT_SPK_mod ./ sum.(eachrow(counts[:, [:WT_SPK_mod, :WT_SPK_unmod]]))

	cols = [:STORM_1, :STORM_2, :STORM_3, :WT_1, :WT_2, :WT_3]
	counts = counts[:, cols]
	# local counts = df[:, cols]

	# local sums = sum.(eachrow(counts))


	# counts ./ sums

	println(cor(counts.STORM_1, counts.STORM_2))
	println(cor(counts.STORM_1, counts.STORM_3))
	println(cor(counts.STORM_2, counts.STORM_3))

	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "Sample",
					ylabel = "Sites",
					yticks = ([900, non_glori_count, 2200, nrow(counts)],
							  ["DRACH", "$non_glori_count", "GLORI+", "$(nrow(counts))"]),
				    xticks = (1:length(cols), map(string, cols)),
				    xticklabelrotation = 45)

	heatmap!(ax, transpose(Array(counts)), colormap = :blues)

	lines!(ax, [0.5, 6.5], [non_glori_count, non_glori_count], color = :red)

	Colorbar(f[1, 2],
			 colorrange = (0, 1),
			 colormap = :blues,
			 label = "Proportion of reads modified")

	f
	# transpose(Array(counts))
end
  ╠═╡ =#

# ╔═╡ 9f184944-6093-487f-a095-98a8572af620
#=╠═╡
conf_storm_resistant
  ╠═╡ =#

# ╔═╡ ff570f16-c8fb-4e0d-8851-a0fb93b23645
#=╠═╡
begin
	local mettl16_kd_peaks = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSE90914_m6a_peaks_mettl16_kd.csv"))
	local t = innerjoin(conf_storm_resistant, mettl16_kd_peaks,
					    on = [:chr, :strand],
					    renamecols = "" => "_m16")
	t = t[between.(t.genomicPos, t.start_m16, t.stop_m16), :]
end
  ╠═╡ =#

# ╔═╡ b23610f5-240b-4385-955b-216b3221d1d8
#=╠═╡
conf_storm_resistant
  ╠═╡ =#

# ╔═╡ 119b9407-c1cd-4eff-a0a6-993b9b05f50b
#=╠═╡
map(r -> println(split(split(r, "|")[2], ".")[1]), unique(conf_storm_resistant_20perc.ref_id)) |> length
  ╠═╡ =#

# ╔═╡ f208fdf1-e949-468a-a987-d798baa939ca
#=╠═╡
map(r -> println(split(r, "|")[6]), unique(conf_storm_resistant.ref_id)) |> length
  ╠═╡ =#

# ╔═╡ 434ef740-5386-4a64-9ebe-8211e40b0c2b
#=╠═╡
map(r -> println(">$(split(r.ref_id, "|")[5])-$(r.pos)\n$(r.ref_kmer)"), eachrow(conf_storm_resistant)) |> length
  ╠═╡ =#

# ╔═╡ f9dbfebc-9c8e-451a-a46c-03b548b630de
md"""
 I took the storm-resistant glori+ sites (43) and ran motif enrichment using https://meme-suite.org/meme/. The DRACH motif was found in 28 sites
"""

# ╔═╡ 46aee8ff-83fc-4e88-95d7-20a1c076d3b2
#=╠═╡
data(DataFrame(kmer = map(m -> m.match,
						  filter(m -> !isnothing(m),
								 match.(r"[GAT][AG]AC[ACT]",
										conf_storm_resistant.ref_kmer))))
) * mapping(:kmer) * frequency() |> draw(; axis = (; title = "STORM-insensitive sites DRACH sequences", xticklabelrotation = 45))
  ╠═╡ =#

# ╔═╡ 6a44e37f-9925-44a5-afcd-6f3379796e86
data(DataFrame(kmer = map(m -> m.match,
						  filter(m -> !isnothing(m),
								 match.(r"[GAT][AG]AC[ACT]",
										sig_peaks_storm.ref_kmer))))
) * mapping(:kmer) * frequency() |> draw(; axis = (; title = "STORM-sensitive sites DRACH sequences", xticklabelrotation = 45))

# ╔═╡ 3fd1b51b-3a5c-4d36-9cf2-1281c0e4483c
#=╠═╡
begin
	local storm_sensitive = map(m -> m.match,
						  		filter(m -> !isnothing(m),
								 	   match.(r"[GAT][AG]AC[ACT]", sig_peaks_storm.ref_kmer)))
	local storm_insensitive = map(m -> m.match,
						  		  filter(m -> !isnothing(m),
								 	     match.(r"[GAT][AG]AC[ACT]", conf_storm_resistant.ref_kmer)))

	local df = vcat(DataFrame(kmer = storm_sensitive, type = "STORM Sensitive"),
					DataFrame(kmer = storm_insensitive, type = "STORM Insensitive"))

	# data(df) * mapping(:kmer, row = :type) * frequency() |> draw(; facet = (; linkyaxes = false), axis = (; xticklabelrotation = 45, xlabel= "k-mer", ylabel = "Number of modifications"))
	data(df) * mapping(:kmer, row = :type) * visual(Hist, normalization = :pdf, gap = 0.1, bins = 0.5:1:18.5) |> draw(; axis = (; xticklabelrotation = 45, xlabel= "k-mer", ylabel = "pdf"))
					
end
  ╠═╡ =#

# ╔═╡ 2357c8ef-6b46-4f17-97da-cdb99668ca85
#=╠═╡
begin
	local storm_sensitive = map(m -> m.match,
						  		filter(m -> !isnothing(m),
								 	   match.(r"[GAT][AG]AC[ACT]", an_sig_peaks_storm[an_sig_peaks_storm.mod .=== "m6A", :ref_kmer])))
	local storm_insensitive = map(m -> m.match,
						  		  filter(m -> !isnothing(m),
								 	     match.(r"[GAT][AG]AC[ACT]", conf_storm_resistant[conf_storm_resistant.glori, :].ref_kmer)))

	local df = vcat(DataFrame(kmer = storm_sensitive, type = "STORM Sensitive"),
					DataFrame(kmer = storm_insensitive, type = "STORM Insensitive"))

	# data(df) * mapping(:kmer, row = :type) * frequency() |> draw(; facet = (; linkyaxes = false), axis = (; xticklabelrotation = 45, xlabel= "k-mer", ylabel = "Number of modifications"))
	data(df) * mapping(:kmer, row = :type) * visual(Hist, normalization = :pdf, gap = 0.1, bins = 0.5:1:18.5) |> draw(; axis = (; xticklabelrotation = 45, xlabel= "k-mer", ylabel = "pdf", title = "Only GLORI+ DRACH sites"))
					
end
  ╠═╡ =#

# ╔═╡ 29283d8b-c564-4f0a-a6a2-6d198d320162
#=╠═╡
begin
	local storm_sensitive = map(m -> m.match,
						  		filter(m -> !isnothing(m),
								 	   match.(r"[GAT][AG]AC[ACT]", sig_peaks_storm.ref_kmer)))
	local storm_insensitive = map(m -> m.match,
						  		  filter(m -> !isnothing(m),
								 	     match.(r"[GAT][AG]AC[ACT]", conf_storm_resistant[conf_storm_resistant.glori, :].ref_kmer)))
	storm_sensitive = countmap(storm_sensitive)
	storm_insensitive = countmap(storm_insensitive)
	local df = DataFrame(kmer = [k for k in union(keys(storm_sensitive), keys(storm_insensitive))])
	df[:, :sensitive] = map(k -> Float64(get!(storm_sensitive, k, 0)), df.kmer)
	df[:, :insensitive] = map(k -> Float64(get!(storm_insensitive, k, 0)), df.kmer)
	df[:, :sensitive] ./= sum(df.sensitive)
	df[:, :insensitive] ./= sum(df.insensitive)


	(data(df) * mapping(:sensitive => "Sensitive", (:sensitive, :insensitive) => ((s, i) -> i - s) => "change") * visual(Scatter)
		# +
	 # data(df) * mapping(
		#  :sensitive => (x -> x + 0.1) => "Sensitive",
		#  (:sensitive, :insensitive) => ((s, i) -> i - s + 0.1) => "change") *  visual(Scatter)
	)	|> draw

	# data(df) * mapping(:sensitive, :insensitive) * visual(Scatter) |> draw

	# data(stack(df, 2:3)) * mapping(:kmer, :value, row = :variable) * visual(BarPlot) |> draw(; axis = (; xticklabelrotation = 45))

	# data(df) * mapping([:sensitive, :insensitive] .=> "Type") * mapping(color=dims(1) => renamer(["Sensitive", "Insensitive"])) * mapping(:kmer, color = "Type") * visual(BarPlot) |> draw
	
end
  ╠═╡ =#

# ╔═╡ 452d1d7b-afae-4bf7-b0a4-88474eaf947c
#=╠═╡
begin
	local t = innerjoin(conf_storm_resistant, glori_stm,
	          			on = [:chr => :Chr, :strand => :Strand])

	local df = innerjoin(t[between.(t.Sites, t.start, t.end), :],
			  glori_raw,
			  on = [:chr => :Chr, :strand => :Strand, :Sites],
			  renamecols = "" => "_ctrl")
	df[:, :stm_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	df[:, :ctrl_mean_ratio] = (df.Ratio_ctrl .+ df.Ratio_1_ctrl)./2

	(data(df) * mapping(:ctrl_mean_ratio => "GLORI WT stoichiometry", :stm_mean_ratio => "GLORI STORM stoichiometry") * visual(Scatter)
	+
	data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines)) |> draw(; axis = (; aspect = 1, xticks = 0:0.25:1, yticks = 0:0.25:1))
end
  ╠═╡ =#

# ╔═╡ 95b209bc-4fbb-4713-b1b0-8f6b4b1a14b4
#=╠═╡
begin
	local insens_wt = conf_storm_resistant.wt_mod_ratio
	local insens_stm = conf_storm_resistant.stm_mod_ratio
	local sens_wt = calc_mod_ratio(sig_peaks_storm, ["WT_1", "WT_2", "WT_SPK"])
	local df = vcat(DataFrame(ratio = insens_wt, type = "a. STORM insensitive (WT)"),
				    DataFrame(ratio = insens_stm, type = "b. STORM insensitive (STORM)"),
				    DataFrame(ratio = sens_wt, type = "c. STORM sensitive (WT)"))
	data(df) * mapping(:ratio => "Stoichiometry", layout = :type) * visual(Density) |> draw(scales(Layout = (; palette = wrapped(cols = 2, by_col = true))))
end
  ╠═╡ =#

# ╔═╡ f652d459-817e-47c2-8456-3de6b5f89f9b
begin
	local t = innerjoin(glori_raw, glori_stm,
		  				on = [:Chr, :Strand, :Sites],
		  				renamecols = "" => "_stm")
	t[:, :ctrl_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :stm_mean_ratio] = (t.Ratio_stm .+ t.Ratio_1_stm)./2
	t[:, :resistant] = log2.(t.stm_mean_ratio ./ t.ctrl_mean_ratio) .> -0.1

	t = innerjoin(t, storm[:, [:chr, :strand, :genomicPos, :ref_kmer]],
			  	  on = [:Chr => :chr, :Strand => :strand, :Sites => :genomicPos])
	
	local storm_sensitive = map(m -> m.match,
						  		filter(m -> !isnothing(m),
								 	   match.(r"[GAT][AG]AC[ACT]", t[.! t.resistant, :ref_kmer])))
	local storm_insensitive = map(m -> m.match,
						  		  filter(m -> !isnothing(m),
								 	     match.(r"[GAT][AG]AC[ACT]", t[t.resistant, :ref_kmer])))

	local df = vcat(DataFrame(kmer = storm_sensitive, type = "STORM Sensitive"),
					DataFrame(kmer = storm_insensitive, type = "STORM Insensitive"))

	df = combine(groupby(df, [:kmer, :type]), nrow => :count)
	df[:, :prob] .= 0.0
	df[df.type .== "STORM Sensitive", :prob] = df[df.type .== "STORM Sensitive", :count] ./ sum(df[df.type .== "STORM Sensitive", :count])
	df[df.type .== "STORM Insensitive", :prob] = df[df.type .== "STORM Insensitive", :count] ./ sum(df[df.type .== "STORM Insensitive", :count])
	# df

	
	# data(df) * mapping(:kmer, row = :type) * visual(Hist, normalization = :pdf, gap = 0.1, bins = 0.5:1:18.5) |> draw(; axis = (; xticklabelrotation = 45, xlabel= "k-mer", ylabel = "pdf", title = "Using only GLORI"))

	data(df) * mapping(:kmer, :prob, color = :type, dodge = :type) * visual(BarPlot) |> draw(; axis = (; xticklabelrotation = 45, xlabel= "k-mer", ylabel = "Frequency", title = "Using only GLORI"), legend = (; position = :bottom))
					
end

# ╔═╡ 17b3496f-d5bd-4782-9ec9-65723d9ced58
#=╠═╡
begin
	local t = innerjoin(glori_raw, glori_stm,
		  				on = [:Chr, :Strand, :Sites],
		  				renamecols = "" => "_stm")
	t[:, :ctrl_mean_ratio] = (t.Ratio .+ t.Ratio_1)./2
	t[:, :stm_mean_ratio] = (t.Ratio_stm .+ t.Ratio_1_stm)./2
	t[:, :resistant] = log2.(t.stm_mean_ratio ./ t.ctrl_mean_ratio) .> -0.1

	t = innerjoin(t, storm[:, [:chr, :strand, :genomicPos, :ref_kmer]],
			  	  on = [:Chr => :chr, :Strand => :strand, :Sites => :genomicPos])
	
	local storm_sensitive = map(m -> m.match,
						  		filter(m -> !isnothing(m),
								 	   match.(r"[GAT][AG]AC[ACT]", t[.! t.resistant, :ref_kmer])))
	local storm_insensitive = map(m -> m.match,
						  		  filter(m -> !isnothing(m),
								 	     match.(r"[GAT][AG]AC[ACT]", t[t.resistant, :ref_kmer])))

	local df = vcat(DataFrame(kmer = storm_sensitive, type = :sensitive),
					DataFrame(kmer = storm_insensitive, type = :insensitive))

	df = combine(groupby(df, [:kmer, :type]), nrow => :count)
	df[:, :prob] .= 0.0
	df[df.type .== :sensitive, :prob] = df[df.type .== :sensitive, :count] ./ sum(df[df.type .== :sensitive, :count])
	df[df.type .== :insensitive, :prob] = df[df.type .== :insensitive, :count] ./ sum(df[df.type .== :insensitive, :count])

	local glori = df

	local storm_sensitive = map(m -> m.match,
						  		filter(m -> !isnothing(m),
								 	   match.(r"[GAT][AG]AC[ACT]", an_sig_peaks_storm[an_sig_peaks_storm.mod .=== "m6A", :ref_kmer])))
	local storm_insensitive = map(m -> m.match,
						  		  filter(m -> !isnothing(m),
								 	     match.(r"[GAT][AG]AC[ACT]", conf_storm_resistant[conf_storm_resistant.glori, :].ref_kmer)))

	local df = vcat(DataFrame(kmer = storm_sensitive, type = :sensitive),
					DataFrame(kmer = storm_insensitive, type = :insensitive))

	df = combine(groupby(df, [:kmer, :type]), nrow => :count)
	df[:, :prob] .= 0.0
	df[df.type .== :sensitive, :prob] = df[df.type .== :sensitive, :count] ./ sum(df[df.type .== :sensitive, :count])
	df[df.type .== :insensitive, :prob] = df[df.type .== :insensitive, :count] ./ sum(df[df.type .== :insensitive, :count])

	local nanocompore = df

	local kmer_ind = Dict(zip(sort(unique(df.kmer), rev = true), 1:length(unique(df.kmer))))
	local ind_kmer = Dict(map(reverse, collect(kmer_ind)))
	local type_ind = Dict(:sensitive => 1, :insensitive => 2)

	local yticks = (values(kmer_ind), [string(k) for k in keys(kmer_ind)])
	println(yticks)
	println(ind_kmer)

	(data(glori) * mapping(:prob, :kmer => (k -> 2*kmer_ind[k] - 0.3) => "k-mer", marker = :type) * visual(Scatter, color = :blue, label = "GLORI")
		+
		data(glori) * mapping(:prob, :kmer => (k -> 2*kmer_ind[k] - 0.3) => "k-mer", group = :kmer) * visual(Lines, linestyle = :dash, color = :blue)

	+

		data(nanocompore) * mapping(:prob, :kmer => (k -> 2*kmer_ind[k] + 0.3) => "k-mer", marker = :type) * visual(Scatter, color = :orange, label = "Nanocompore")
		+
		data(nanocompore) * mapping(:prob, :kmer => (k -> 2*kmer_ind[k] + 0.3) => "k-mer", group = :kmer) * visual(Lines, linestyle = :dash, color = :orange)
	) |> draw(scales(
        Marker = (; palette = [:xcross, :circle])
    ); axis = (; xlabel = "Frequency", yticks = (2:2:36, [string(ind_kmer[i]) for i in 1:18])))
end
  ╠═╡ =#

# ╔═╡ a8ba8f4e-96ae-49d6-91c4-52d7f839f58a
#=╠═╡
conf_storm_resistant[occursin.("YTHDF1", conf_storm_resistant.ref_id), :]
  ╠═╡ =#

# ╔═╡ dbef3c64-4989-4bdd-916c-045bbeb03263
ivt[ivt.genomicPos .== 130885771, :]

# ╔═╡ c2ca3807-116d-4251-b864-7bbf7e63d41f
storm[storm.genomicPos .== 130885771, :]

# ╔═╡ f3b640b8-ab69-4ce0-b597-d667170f9cc2
storm[occursin.("SLC5A6", storm.ref_id), :]

# ╔═╡ 43603d8b-82dd-4670-a999-51802c6173df
stm_resist_bp = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_genes_all_GO_biological_process_enirchments.txt", delim = '\t'))

# ╔═╡ d5a93508-e559-4767-a078-54f9bb8fb3bf
stm_resist_cc = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_genes_all_GO_cellular_compartment_enirchments.txt", delim = '\t'))

# ╔═╡ 4ad4579a-4b64-460b-9ec4-a80f86ee391a
stm_resist_mf = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/storm_resistant_genes_all_GO_molecular_function_enirchments.txt", delim = '\t'))

# ╔═╡ 446f9868-3f80-4200-97a3-a351de532be2
#=╠═╡
hist(map(r -> parse(Int64, split(r, "|")[7]), unique(conf_storm_resistant.ref_id)), bins = length(unique(conf_storm_resistant.ref_id)))
  ╠═╡ =#

# ╔═╡ 93e58616-9b19-42f4-84ab-da73c10561e4
#=╠═╡
begin
	local refs = unique(storm.ref_id)
	local storm_res_refs = Set(conf_storm_resistant.ref_id)
	
	local lens = map(r -> parse(Int64, split(r, "|")[7]), refs)
	local is_storm_res = map(r -> r in storm_res_refs, refs)

	local df = DataFrame(ref = refs, length = log10.(lens), storm_res = is_storm_res)
	

	# local f = data(df) * mapping(:storm_res => "Has STORM-insensitive m6As", :length => "log₁₀(Transcript length)") * visual(Violin) |> draw#(; axis =(; yticks = [0, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000], limits = (nothing, (-500, 31000))))

	local plt = data(df) * mapping(:length => "log₁₀(Transcript length)", color = :storm_res => "Has STORM-insensitive m6A") * AlgebraOfGraphics.density()

	local f = Figure()
	local gridpos = f[1, 1]

	local ax = draw!(gridpos, plt)
	legend!(gridpos, ax; tellwidth = false, halign = :left, valign = :top, margin = (10, 10, 10, 10))
	# |> draw(; legend = (; position = :bottom))

	# legend!(f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))
	f
end
  ╠═╡ =#

# ╔═╡ 5948af6a-d5df-46a5-a832-ad33d3feb8c5
#=╠═╡
begin
	local lens = map(r -> parse(Int64, split(r, "|")[7]), conf_storm_resistant.ref_id)
	local stoich_change = log2.(conf_storm_resistant.stm_mod_ratio ./ conf_storm_resistant.wt_mod_ratio)

	local df = DataFrame(length = lens, change = stoich_change)

	data(df) * mapping(:length => "Transcript length", :change => "log₂(STORM mod. ratio / WT mod. ratio)") * visual(Scatter, markersize = 5) |> draw(; axis = (; limits = (nothing, (-4, 4))))
end
  ╠═╡ =#

# ╔═╡ 820fb1ae-e3fe-4aa8-8163-e72c4748f9be
#=╠═╡
begin
	local lens = map(r -> parse(Int64, split(r, "|")[7]), conf_storm_resistant.ref_id)
	local stoich_change = log2.(conf_storm_resistant.stm_mod_ratio ./ conf_storm_resistant.wt_mod_ratio)

	local df = DataFrame(length = lens, change = stoich_change)

	data(df) * mapping(:length => log10 => "log₁₀(Transcript length)", :change => "log₂(STORM mod. ratio / WT mod. ratio)") * visual(Scatter, markersize = 5) |> draw(; axis = (; limits = (nothing, (-4, 4))))
end
  ╠═╡ =#

# ╔═╡ 5f66681b-fda7-48fe-8771-d0d44564df45
stm_expression = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/HEK_STMvsDMSOtreated.tsv", delim = '\t'))

# ╔═╡ 75c4e402-c255-49d5-a1c9-df2ffde1a578
#=╠═╡
begin
	local refs = unique(storm.ref_id)
	local storm_res_refs = Set(conf_storm_resistant.ref_id)
	
	local genes = map(r -> string(split(split(r, "|")[2], ".")[1]), refs)
	local is_storm_res = map(r -> r in storm_res_refs, refs)

	local gene_exp = Dict(zip(stm_expression.geneName, stm_expression.baseMean))

	local df = DataFrame(ref = refs,
						 exp = map(g -> get!(gene_exp, g, NaN), genes),
						 storm_res = is_storm_res)
	df = df[.! isnan.(df.exp), :]
	
	df[:, :exp] = log10.(df.exp .+ 1)
	# data(df) * mapping(:storm_res => "Has STORM-insensitive m6As", :exp => "log₁₀(count)") * visual(Violin) |> draw
	local plt = data(df) * mapping(:exp => "log₁₀(count)", color = :storm_res => "Has STORM-insensitive m6A") * AlgebraOfGraphics.density()

	local f = Figure()
	local gridpos = f[1, 1]

	local ax = draw!(gridpos, plt)
	legend!(gridpos, ax; tellwidth = false, halign = :left, valign = :top, margin = (10, 10, 10, 10))
	
	f
	
	# data(df) * mapping(:storm_res => "Has STORM-insensitive m6As", :exp => "Expression") * visual(Violin) |> draw(; axis = (; yticks = [0, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000], limits = (nothing, (-500, 31000))))
end
  ╠═╡ =#

# ╔═╡ eed7c85f-5ec5-41ef-9bd9-250dd07137db
#=╠═╡
begin
	local genes = map(r -> string(split(split(r, "|")[2], ".")[1]), conf_storm_resistant.ref_id)
	local gene_exp = Dict(zip(stm_expression.geneName, stm_expression.baseMean))
	
	local stoich_change = (conf_storm_resistant.stm_mod_ratio .- conf_storm_resistant.wt_mod_ratio) .* 100

	local df = DataFrame(count = map(g -> get!(gene_exp, g, NaN), genes),
						 change = stoich_change)

	data(df) * mapping(:count => "Expression count", :change => "Stoichiometry change") * visual(Scatter, markersize = 5) |> draw
end
  ╠═╡ =#

# ╔═╡ 6628a74b-4f8c-45eb-85f5-16681b2c849e
#=╠═╡
begin
	local t = innerjoin(conf_storm_resistant, encore_m6a_rbp_peaks,
						on = [:chr, :strand],
						renamecols = "" => "_rbp")
	t = t[between.(t.genomicPos, t.start_rbp, t.end_rbp), :]
	t = t[occursin.(r"[GAT][AG]AC[ACT]", t.ref_kmer), :]
	t[:, :motif] = [m.match for m in match.(r"[GAT][AG]AC[ACT]", t.ref_kmer)]

	data(t) * mapping(:rbp_rbp, :motif => (s -> s in Set(["AGACT", "GAACT", "GGACA", "GGACC", "GGACT"]) ? s * "*" : s * " ") => "k-mer") * frequency() |> draw(; axis = (; xticklabelrotation = 45))
end
  ╠═╡ =#

# ╔═╡ 3ef8cb41-d3de-4f38-a93e-2d0211b15263
#=╠═╡
begin
	local refs = unique(storm.ref_id)
	local storm_res_refs = Set(conf_storm_resistant.ref_id)
	
	local genes = map(r -> string(split(split(r, "|")[2], ".")[1]), refs)
	local is_storm_res = map(r -> r in storm_res_refs, refs)

	local gene_lfc = Dict(zip(stm_expression.geneName, map(v -> occursin("NA", v) ? NaN : parse(Float64, v), stm_expression.log2FoldChange)))

	local df = DataFrame(ref = refs,
						 lfc = map(g -> get!(gene_lfc, g, NaN), genes),
						 storm_res = is_storm_res)
	# df = df[.! isnan.(df.exp), :]

	# local f = data(df) * mapping(:storm_res => "Has STORM-insensitive m6As", :exp => "Expression") * visual(Violin) |> draw(; axis =(; yticks = [0, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000], limits = (nothing, (-500, 31000))))
end
  ╠═╡ =#

# ╔═╡ e86d7b59-449e-4e8f-8ad2-1c708ea15076
begin
	local cc = stm_resist_cc[stm_resist_cc.Benjamini .< 0.01, :]
	local mf = stm_resist_mf[stm_resist_mf.Benjamini .< 0.01, :]

	local Ncc = nrow(cc)
	local Nmf = nrow(mf)

	local counts = zeros((Ncc, Nmf))

	for (c, rcc) in enumerate(eachrow(cc)), (m, rmf) in enumerate(eachrow(mf))
		set_cc = Set(split(rcc.Genes, ", "))
		set_mf = Set(split(rmf.Genes, ", "))

		counts[c, m] = length(intersect(set_cc, set_mf))
	end
	local f = Figure(size = (800, 800))
	local ax = Axis(f[1, 1],
				    xticks = (1:Nmf, map(last, split.(mf.Term, "~"))),
					yticks = (1:Ncc, map(last, split.(cc.Term, "~"))),
				    xticklabelrotation = 0.7)

	heatmap!(ax, transpose(counts), colormap = :heat)

	text!(ax,
		  [string(Int(counts[j, i])) for i in 1:Nmf for j in 1:Ncc],
		  position = [Point2f(x, y) for x in 1:Nmf for y in 1:Ncc],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[j, i] > 280 ? :white : :black for i in 1:Nmf for j in 1:Ncc])
	
	
	f
end

# ╔═╡ 6412e0da-abe4-44fb-9e06-7cce1835d686
begin
	local cc = stm_resist_cc[stm_resist_cc.Benjamini .< 0.01, :]
	local bp = stm_resist_bp[stm_resist_bp.Benjamini .< 0.01, :]

	local Ncc = nrow(cc)
	local Nbp = nrow(bp)

	local counts = zeros((Ncc, Nbp))

	for (c, rcc) in enumerate(eachrow(cc)), (m, rbp) in enumerate(eachrow(bp))
		set_cc = Set(split(rcc.Genes, ", "))
		set_bp = Set(split(rbp.Genes, ", "))

		counts[c, m] = length(intersect(set_cc, set_bp))
	end
	local f = Figure(size = (800, 800))
	local ax = Axis(f[1, 1],
				    xticks = (1:Nbp, map(last, split.(bp.Term, "~"))),
					yticks = (1:Ncc, map(last, split.(cc.Term, "~"))),
				    xticklabelrotation = 0.7)

	heatmap!(ax, transpose(counts), colormap = :heat)

	text!(ax,
		  [string(Int(counts[j, i])) for i in 1:Nbp for j in 1:Ncc],
		  position = [Point2f(x, y) for x in 1:Nbp for y in 1:Ncc],
		  align = (:center, :center),
		  fontsize = 14,
		  color = [counts[j, i] > 280 ? :white : :black for i in 1:Nbp for j in 1:Ncc])
	
	
	f
end

# ╔═╡ 2d9caf6d-154d-408e-8e4d-aef4f1c3ff78
md"""
## SCRATCHPAD
"""

# ╔═╡ 1a8aac68-4308-4700-ad2d-8e6989757416
storm[storm.mean_cov .< 30, :]

# ╔═╡ 4b5089dd-74aa-4dbf-8e26-bdaecc0b5b36
begin
	local non_m6a_storm = an_sig_peaks_storm[an_sig_peaks_storm.mod .!== missing .&&
										 an_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	# local ivt_covered = with_ivt.predicted_1 .!== missing
	# with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	println(sum(with_ivt.predicted_1 .< 2))
	println(sum(with_ivt.predicted_1 .>= 2))

	with_ivt[:, :color] = map(
		s -> if occursin("HEK293T", s)
				:green
			 elseif occursin("HEK293", s)
				:orange
			 elseif s == "motif"
				:lightyellow
			 else
				:lightgrey
			 end,
		with_ivt.source)

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Non-m6A (based on GLORI) mods significant in STORM/WT",
				    xlabel = "STORM/WT -log10(q-value)",
				    ylabel = "IVT/WT -log10(q-value)")

	scatter!(ax,
			 with_ivt.predicted,
			 with_ivt.predicted_1,
			 color = with_ivt.color,
			 strokewidth = 1,
			 strokecolor = :black,
			 markersize = 7)
	text!(ax, with_ivt.predicted .+ rand(nrow(with_ivt)) .* 2, with_ivt.predicted_1 .+ rand(nrow(with_ivt)) .* 2, text = with_ivt.mod, fontsize = ifelse.(with_ivt.predicted .> 60 .|| with_ivt.predicted_1 .> 15, 8, 0))

	Legend(f[2, 1],
		   [PolyElement(color = :green),
			PolyElement(color = :orange),
			PolyElement(color = :lightgrey),
			PolyElement(color = :lightyellow)],
	       ["HEK293T", "HEK293", "Unspecified/Other cell line", "Motif"],
		   "Attribution source",
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 7361c452-d2c9-4793-855e-055911702d8a
begin
	local non_m6a_storm = an_sig_peaks_storm[an_sig_peaks_storm.mod .!== missing .&&
									  	     an_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	# local ivt_covered = with_ivt.predicted_1 .!== missing
	# with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	println(sum(with_ivt.predicted_1 .< 2))
	println(sum(with_ivt.predicted_1 .>= 2))


	local inhib = with_ivt[with_ivt.predicted_1 .< 2, :]

	inhib = combine(groupby(inhib, [:peak_id, :mod]),
					names(inhib) .=> first .=> names(inhib))
	
	# with_ivt[:, :color] = map(
	# 	s -> if occursin("HEK293T", s)
	# 			:green
	# 		 elseif occursin("HEK293", s)
	# 			:orange
	# 		 elseif s == "motif"
	# 			:lightyellow
	# 		 else
	# 			:lightgrey
	# 		 end,
	# 	with_ivt.source)
	local colormap = Dict(zip(unique(with_ivt.mod), 1:length(unique(with_ivt.mod))))

	local f = Figure(size = (1400, 1100))
	local ax = Axis(f[1, 1],
					title = "Putative non-m6A inhibited by STORM",
				    xlabel = "|LOR|",
				    ylabel = "WT/STORM -log10(q-value)")
	xlims!(ax, (0.7, 4.5))
	
	local non_drachs = .! occursin.(writer_motifs["METTL3"], inhib.ref_kmer)
	# println(inhib[non_drachs, :])

	scatter!(ax,
			 abs.(inhib.GMM_LOR),
			 inhib.predicted_raw,
			 color = map(m -> colormap[m], inhib.mod),
			 colormap = :tab10,
			 colorrange = (1, length(colormap)),
			 # strokewidth = 1,
			 # strokecolor = :black,
			 marker = ifelse.(non_drachs,
						 :circle,
						 :star5),
			 markersize = 15)

	
	local xoffsets = abs.(randn(sum(non_drachs))) ./ 3
	local yoffsets = randn(sum(non_drachs)) .* 3
	
	text!(ax,
		  abs.(inhib[non_drachs, :GMM_LOR]) .+ xoffsets,
		  disallowmissing(inhib[non_drachs, :predicted_raw]) .+ yoffsets,
		  text = inhib[non_drachs, :ref_kmer],
		  align = (:left, :center))
		  # fontsize = ifelse.(abs.(inhib.GMM_LOR) .> 1 .|| inhib.predicted_raw .> 15, 12, 0))

	for (x, y, x2, y2) in zip(abs.(inhib[non_drachs, :GMM_LOR]),
							  inhib[non_drachs, :predicted_raw],
							  abs.(inhib[non_drachs, :GMM_LOR]) .+ xoffsets,
							  inhib[non_drachs, :predicted_raw] .+ yoffsets)
		lines!(ax, [x, x2], [y, y2])
	end

	
	Legend(f[2, 1],
		   [PolyElement(color = c,
						colormap=:tab10,
						colorrange = (1, length(colormap)))
			for c in values(colormap)],
		   [m for m in keys(colormap)],
		   "Inferred modification identity",
		   orientation = :horizontal,
		   framevisible = false)
	
	f
	# inhib[non_drachs, :]
end

# ╔═╡ ae549f94-6e59-48f2-aedd-5703b7301202
groupby(storm[storm.chr .== "chr2" .&& storm.genomicPos .> 27200310 .&& storm.genomicPos .< 27200333, [:ref_id, :genomicPos, :pos, :GMM_chi2_pvalue, :GMM_LOR, :GMM_chi2_qvalue]], :ref_id)

# ╔═╡ 0a46f8ed-2a03-4b04-98be-5d85473343ea
groupby(ivt[ivt.chr .== "chr2" .&& ivt.genomicPos .> 27200305 .&& ivt.genomicPos .< 27200333, [:ref_id, :genomicPos, :pos, :GMM_chi2_qvalue, :predicted, :mean_cov, :predicted_raw, :GMM_LOR]], :ref_id)

# ╔═╡ 3c36cc71-db84-4f59-81fe-e4dd7761f72b
annot_sig_peaks[occursin.("UBAC2", annot_sig_peaks.ref_id), :]

# ╔═╡ e9d02b89-5a52-424f-a371-d0ca248bf316
annot_peaks[occursin.("UBAC2", annot_peaks.ref_id), :]

# ╔═╡ 913d87d1-6d0a-4db3-9e46-f8a10d247dd5
peaks_ivt[occursin.("UBAC2", peaks_ivt.ref_id), :]

# ╔═╡ faddf0dd-f931-4379-a13a-021bd5b23cae
peaks(ivt[occursin.("ENST00000403766.8", ivt.ref_id), :], 4)

# ╔═╡ c4a7d77d-68cd-40ef-a168-6dbf9dc410b8
ivt[occursin.("ENST00000403766.8", ivt.ref_id), :][1420:1530, [:pos, :GMM_chi2_pvalue, :GMM_chi2_qvalue, :GMM_LOR, :predicted]]

# ╔═╡ d3567161-530b-476a-b339-a0b2d415da7e
argmaxima(ivt[occursin.("ENST00000403766.8", ivt.ref_id), :][1420:1530, :predicted])

# ╔═╡ ff2cea1e-4bb4-4cd2-a1f4-a68cfeb6efde
begin
	local non_m6a_storm = an_sig_peaks_storm[an_sig_peaks_storm.mod .!== missing .&&
										 an_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	# local ivt_covered = with_ivt.predicted_1 .!== missing
	# with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	println(sum(with_ivt.predicted_1 .< 2))
	println(sum(with_ivt.predicted_1 .>= 2))


	local promote = with_ivt[with_ivt.predicted_1 .>= 2, :]

	promote = combine(groupby(promote, [:peak_id, :mod]),
					  names(promote) .=> first .=> names(promote))
	
	# with_ivt[:, :color] = map(
	# 	s -> if occursin("HEK293T", s)
	# 			:green
	# 		 elseif occursin("HEK293", s)
	# 			:orange
	# 		 elseif s == "motif"
	# 			:lightyellow
	# 		 else
	# 			:lightgrey
	# 		 end,
	# 	with_ivt.source)
	local colormap = Dict(zip(unique(with_ivt.mod), 1:length(unique(with_ivt.mod))))

	local f = Figure(size = (1400, 1100))
	local ax = Axis(f[1, 1],
					title = "Putative non-m6A promoted by STORM",
				    xlabel = "|LOR|",
				    ylabel = "WT/STORM -log10(q-value)")
	xlims!(ax, (0.7, 4.5))
	
	local non_drachs = .! occursin.(writer_motifs["METTL3"], promote.ref_kmer)

	scatter!(ax,
			 abs.(promote.GMM_LOR),
			 promote.predicted_raw,
			 color = map(m -> colormap[m], promote.mod),
			 colormap = :tab10,
			 colorrange = (1, length(colormap)),
			 # strokewidth = 1,
			 # strokecolor = :black,
			 marker = ifelse.(non_drachs,
						 :circle,
						 :star5),
			 markersize = 15)

	
	local xoffsets = abs.(randn(sum(non_drachs))) ./ 3
	local yoffsets = randn(sum(non_drachs)) .* 3
	
	text!(ax,
		  abs.(promote[non_drachs, :GMM_LOR]) .+ xoffsets,
		  disallowmissing(promote[non_drachs, :predicted_raw]) .+ yoffsets,
		  text = promote[non_drachs, :ref_kmer],
		  align = (:left, :center))

	for (x, y, x2, y2) in zip(abs.(promote[non_drachs, :GMM_LOR]),
							  promote[non_drachs, :predicted_raw],
							  abs.(promote[non_drachs, :GMM_LOR]) .+ xoffsets,
							  promote[non_drachs, :predicted_raw] .+ yoffsets)
		lines!(ax, [x, x2], [y, y2])
	end

	
	Legend(f[2, 1],
		   [PolyElement(color = c,
						colormap=:tab10,
						colorrange = (1, length(colormap)))
			for c in values(colormap)],
		   [m for m in keys(colormap)],
		   "Inferred modification identity",
		   orientation = :horizontal,
		   framevisible = false)
	
	f
	# promote[non_drachs, :]
end

# ╔═╡ 26299f99-3851-4130-84b6-eaac2c59fbd8
Set.(skipmissing.(zip([ifelse.(occursin.(motif, annot_sig_peaks.ref_kmer), writer_mods[writer], missing) for (writer, motif) in writer_motifs]...)))

# ╔═╡ da893d5b-1fb4-41f1-a8c9-16c09196739a
begin
	local motif_supported = zip(
	[ifelse.(occursin.(motif, annot_sig_peaks.ref_kmer),
			 writer_mods[writer],
			 missing)
	 for (writer, motif) in writer_motifs]...) |> 
	d -> skipmissing.(d) |>
	d -> Set.(d)
	# map(s -> annot_sig_peaks.mod in s, motif_supported)
	map(((m, s),) -> m in s, zip(annot_sig_peaks.mod, motif_supported)) |> countmap
end

# ╔═╡ f474351b-6b02-4402-a722-9bd03ac6b1b9
combine(groupby(annot_sig_peaks, [:peak_id, :mod]),
		names(annot_sig_peaks) .=> first)

# ╔═╡ 78ac5d01-2ded-41c5-9b6b-64dfc94a2141
begin
	local ambig_peaks = filter(r -> r.unique_mods > 1, eachrow(combine(groupby(annot_sig_peaks, :peak_id),
		:mod => length∘unique => :unique_mods))).peak_id |> Set
	annot_sig_peaks[map(p -> p in ambig_peaks, annot_sig_peaks.peak_id), :]
end

# ╔═╡ 69b7aa3a-7d0f-4af9-8c71-f1479d60225b
first(annot_sig_peaks)

# ╔═╡ 577ba306-497a-4cfc-8ea4-bdbc8afb2076
annot_sig_peaks.mod |> unique

# ╔═╡ f621f4d2-eb36-4733-84a1-9a77f3d88cb9
begin
	local non_m6a_storm = annot_sig_peaks_storm[annot_sig_peaks_storm.mod .!== missing .&&
										 annot_sig_peaks_storm.mod .!== "m6A", :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)

	with_ivt[with_ivt.predicted_1 .> 75, :]
end

# ╔═╡ 31aa6d10-3ba3-4af0-9458-2c6a60aa01d1
mod_storm = annotate_mods(storm)

# ╔═╡ 552238b4-26f1-4a3c-9706-65410a8ea4d1
# mod_peak_storm = simple_peaks(mod_storm, 4)
# mod_peak_storm = peaks(mod_storm, 4)

# ╔═╡ 966e302e-5ab9-4806-ac72-244eb4be2aae
mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&& mod_peak_storm.mod_type .!== missing, :]

# ╔═╡ eb32157e-74aa-4e61-9b4c-3b57bfa1cb0a
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = leftjoin(non_m6a_storm, annot_peaks,
			 				  on = [:ref_id, :pos],
							  makeunique = true)
	local ivt_covered = with_ivt.predicted_1 .!== missing
	with_ivt[:, :predicted_1] = coalesce.(with_ivt.predicted_1, 0)

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Non-m6A mods significant in STORM/WT",
				    xlabel = "STORM/WT -log10(q-value)",
				    ylabel = "IVT/WT -log10(q-value)")

	scatter!(ax, with_ivt.predicted, with_ivt.predicted_1, markersize = 5, color = ivt_covered)
	# text!(ax, with_ivt.predicted .+ rand(nrow(with_ivt)) .* 2, with_ivt.predicted_1 .+ rand(nrow(with_ivt)) .* 2, text = with_ivt.mod_type, fontsize = 2 .*log2.(with_ivt.predicted))
	f
end

# ╔═╡ 734aebf9-6da9-4e6e-bc79-2905bcd75e7a
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = innerjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
end

# ╔═╡ a6b67788-7967-47fc-9d47-bc4e7c202f5a
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = innerjoin(storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	local refs = Set(with_ivt.ref_id)
	mod_peak_storm[map(ref -> ref in refs, mod_peak_storm.ref_id) .&& mod_peak_storm.mean_m6A_level .!== missing .&& occursin.("ENST00000621592.8", mod_peak_storm.ref_id), :]
end

# ╔═╡ 9a545506-f35d-48b8-aeff-ca24d8ab6442
begin
	local non_m6a_storm = mod_peak_storm[mod_peak_storm.mean_m6A_level .=== missing .&&
										 mod_peak_storm.mod_type .!== missing, :]
	local with_ivt = leftjoin(non_m6a_storm, ivt,
			 				   on = [:ref_id, :pos],
							   makeunique = true)
	local refs = Set(with_ivt.ref_id)
	map(r -> println(split(r, "|")[6]), [r for r in refs])
	nothing
end

# ╔═╡ 3b9fd69d-ca22-4f05-a563-3492ddd66285
# ╠═╡ disabled = true
#=╠═╡
begin
	local a = glori[:, [:chrom, :start, :end, :strand]]
	rename!(a, [:chr, :start, :end, :strand])
	a[:, :name] .= "m6A"
	local b = rmbase_all[:, [:chr, :start, :end, :strand, :mod_type]]
	rename!(b, [:chr, :start, :end, :strand, :name])
	local c = hela_ac4C[:, [:chr, :start, :end, :strand]]
	c[:, :name] .= "ac4C"
	local d = atoi[:, [:chr, :start, :end, :strand]]
	d[:, :name] .= "A-I"
	local combined = vcat(a, b, c, d)
	combined[:, :score] .= 0
	combined[:, [:chr, :start, :end, :name, :score, :strand]] |> CSV.write("/home/mzdravkov/mods.bed", delim = '\t', header = false)
end
  ╠═╡ =#

# ╔═╡ 4167bc84-b881-4eb0-b140-aa35c0acd365
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(ivt_m6a_motif))
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT significant DRACH not in GLORI",
				    xlabel = "IVT/WT -log10(P-value)",
				    ylabel = "STORM/WT -log10(P-value)")
	println(sum([x >= 2 for x in t.predicted_1 if x !== missing]))
	scatter!(ax, t.predicted, t.predicted_1, markersize = 5)
	f
end

# ╔═╡ c92f80c5-6f47-47c3-92d8-8d1cc9fd6c42
# ╠═╡ disabled = true
#=╠═╡
storm_ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/STORM_IVT_eventalign_v2_0_0_mapq_gt_0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false)
  ╠═╡ =#

# ╔═╡ 50443a47-3600-4733-8b77-1d56717b00b4
#=╠═╡
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(ivt_m6a_motif), " ", nrow(t))
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT significant DRACH not in GLORI",
				    xlabel = "IVT/WT -log10(P-value)",
				    ylabel = "IVT/STORM -log10(P-value)")
	println(sum([x >= 2 for x in t.predicted_1 if x !== missing]))
	println(sum([x <= 2 for x in t.predicted_1 if x !== missing]))
	scatter!(ax, t.predicted, t.predicted_1, markersize = 5)
	f
end

  ╠═╡ =#

# ╔═╡ 3060c2aa-078b-43e3-8451-b17642f470eb
#=╠═╡
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	local res = DataFrames.combine(groupby(t, [:wt_storm, :ivt_storm]), nrow => :count)
	res[:, :ivt_wt] .= true
	res[:, [:ivt_wt, :wt_storm, :ivt_storm, :count]]
end

  ╠═╡ =#

# ╔═╡ a65930bd-58ca-4326-a04d-0416614b7abb
#=╠═╡
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, storm,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	t[t.wt_storm .=== true .&& t.ivt_storm .=== true, :]
end

  ╠═╡ =#

# ╔═╡ 5b380673-9c93-4e6f-a158-efe469ab8c52
#=╠═╡
mod_storm_ivt = annotate_mods(storm_ivt)
  ╠═╡ =#

# ╔═╡ e265dcee-67b7-4238-987a-19242ffe9f03
#=╠═╡
mod_peak_storm_ivt = peaks(mod_storm_ivt, 4)
  ╠═╡ =#

# ╔═╡ d3dd2610-f200-4b22-bfb3-6d28df96de13
#=╠═╡
begin
	local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
		  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(ivt_m6a_motif, sort(expand_labels(mod_peak_storm, :pos, 4, 4), [:ref_id, :pos]),
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, sort(expand_labels(mod_peak_storm_ivt, :pos, 4, 4), [:ref_id, :pos]),
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	DataFrames.combine(groupby(t, [:wt_storm, :ivt_storm]), nrow => :count)
end

  ╠═╡ =#

# ╔═╡ e6e36c0d-f522-4959-a979-20573ffda8c2
#=╠═╡
begin
	# local ivt_m6a_motif = mod_peak_ivt[mod_peak_ivt.motif_mod_type .=== "m6A" .&&
	# 	  							   mod_peak_ivt.mean_m6A_level .=== missing, :]
	local t = leftjoin(storm, ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	local t = leftjoin(t, storm_ivt,
			 		   on = [:ref_id, :pos],
			 		   makeunique = true)
	println(nrow(t))
	t[:, :wt_storm] = t.predicted .>= 2
	t[:, :ivt_wt] = t.predicted_1 .>= 2
	t[:, :ivt_storm] = t.predicted_2 .>= 2
	local res = DataFrames.combine(groupby(t, [:ivt_wt, :wt_storm, :ivt_storm]), nrow => :count)
	# res[:, :ivt_wt] .= true
	res[:, [:ivt_wt, :wt_storm, :ivt_storm, :count]]
end

  ╠═╡ =#

# ╔═╡ 55243de1-4f6e-4184-a1c3-40b3e24b1a1b
function genomic_peaks(df, radius)
	all_positions = DataFrames.combine(groupby(df, [:chr, :strand]), :genomicPos => (p -> minimum(p):maximum(p)) => :genomicPos)
	all_positions = leftjoin(all_positions, df,
			 				 on = [:chr, :strand, :genomicPos])
	all_positions = sort(all_positions, [:chr, :strand, :genomicPos])
	# all_positions[:, :predicted] = coalesce.(all_positions.predicted, 0)
	peak_indices = argmaxima(all_positions.predicted, radius, strict = false)
	all_positions[peak_indices, :]
end

# ╔═╡ ee927604-3622-4dd0-8f02-b6838e036fc1
mod_all_ivt = annotate_mods(ivt, only_significant = false)

# ╔═╡ 2254d304-dc0a-4dae-87d7-f3cd18ebe152
begin
	local hek293T_mods = mod_all_ivt[(mod_all_ivt.cell_list .!== missing .&& occursin.("HEK293T", mod_all_ivt.cell_list)) .|| (mod_all_ivt.mod_type .!== missing .&& mod_all_ivt.mod_type .== "m6A"), :]
	# sum(hek293T_mods.predicted .> -log10(0.01)) / nrow(hek293T_mods)
	for mod in unique(hek293T_mods.mod_type)
		t = hek293T_mods[hek293T_mods.mod_type .!== missing .&& hek293T_mods.mod_type .== mod, :]
		println(mod, " ", sum(t.predicted .> -log10(0.01)) / nrow(t))
	end
end

# ╔═╡ 4e1e1f43-de2e-4584-ac04-31bcfd7036ad
mod_all_ivt_p12 = annotate_mods(ivt_p12, only_significant = false)

# ╔═╡ 5fbffc42-e89e-47d8-b493-e297ad66ca0b
begin
	local hek293T_mods = mod_all_ivt_p12[(mod_all_ivt_p12.cell_list .!== missing .&& occursin.("HEK293T", mod_all_ivt_p12.cell_list)) .|| (mod_all_ivt_p12.mod_type .!== missing .&& mod_all_ivt_p12.mod_type .== "m6A"), :]
	for mod in unique(hek293T_mods.mod_type)
		t = hek293T_mods[hek293T_mods.mod_type .!== missing .&& hek293T_mods.mod_type .== mod, :]
		println(mod, " ", sum(t.predicted .> -log10(0.01)) / nrow(t))
	end
end

# ╔═╡ 8201e87f-19c7-4b61-bfb8-46ed75bfedcc
# genomic_peaks(ivt, 4)

# ╔═╡ 3dd015d7-3c68-4b7b-9331-d6bd51c3702d
rmbase_all.mod_type |> countmap

# ╔═╡ e9c15162-c3d4-449f-88db-0280b5874a69
md"""
## Indirect recall analysis via GLORI
"""

# ╔═╡ 70267a1f-34da-49c6-a09f-1c214d0dfecb
# ╠═╡ disabled = true
#=╠═╡
begin
	rna004_500 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_500_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=true)
	binned_rna004_500 = annotate_binned_results(bin(rna004_500, BIN_SIZE),
										        binned_glori)
	:ok
end
  ╠═╡ =#

# ╔═╡ c15eca41-5fda-4869-a298-4991fa4ad72d
# ╠═╡ disabled = true
#=╠═╡
begin
	ivt_500 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_depth_500_v2_0_0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=true)
	binned_ivt_500 = annotate_binned_results(bin(ivt_500, BIN_SIZE),
										     binned_glori)
	:ok
end
  ╠═╡ =#

# ╔═╡ a0e6818d-4d4d-4644-84ab-a4c45c9dcb26
#=╠═╡
binned_storm_ivt_common = innerjoin(binned_rna004_500, binned_ivt_500,
							  		on = [:chr, :strand, :bin],
							  		makeunique = true)

  ╠═╡ =#

# ╔═╡ 36477658-b2f5-4ff4-b57a-72367ed61819
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1], aspect = 1)
	xlims!(ax, (0, 200))
	ylims!(ax, (0, 200))
 	local t = binned_storm_ivt_common[abs.(binned_storm_ivt_common.mean_cov .- binned_storm_ivt_common.mean_cov_1) .< 50 .&& binned_storm_ivt_common.modified .== 1, :]
 	scatter!(ax, t.predicted, t.predicted_1, markersize = 5)
	f
end
  ╠═╡ =#

# ╔═╡ 4593f1eb-c752-4043-98d9-1a00954d738e
#=╠═╡
scatter(binned_storm_ivt_common.mean_cov, binned_storm_ivt_common.mean_cov_1, markersize = 5)
  ╠═╡ =#

# ╔═╡ 0f4380ab-a77f-4413-823f-c16044efd721
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT MA plot - positions",
				    xlabel = "Total coverage in both conditions",
					yticks = -9:1:9,
				    ylabel = "LOR",
				    xscale = log10)
	local df = ivt[abs.(ivt.GMM_LOR) .> 0.01, :]
	# local m = df.GMM_LOR
	# local a = log10.(df.mean_cov)
	scatter!(ax, df.mean_cov, df.GMM_LOR,
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end
  ╠═╡ =#

# ╔═╡ 91d20be7-a7a4-4120-8756-e8e5656ec804
annot_peaks

# ╔═╡ 0e7d33ba-f19d-45b3-8606-0fe9eb1845c3
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1])
	local df = annot_peaks[annot_peaks.GMM_chi2_pvalue .!== missing, :]
	local expected = -log10.((1:nrow(df))./nrow(df))
	scatter!(ax,
			 expected,
			 -log10.(sort(df.GMM_chi2_pvalue .+ 1e-310)),
			 markersize = 5)
	lines!(ax, [0, 4], [0, 4], color = :red)
	local n = nrow(df)
	local lower = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1-0.95)/2))
	local upper = -log10.(
		Distributions.quantile.(Distributions.Beta.(1:n, n .- (1:n) .+ 1), (1+0.95)/2))
	lines!(ax, expected, lower, color = :orange, alpha = 0.5)
	lines!(ax, expected, upper, color = :orange, alpha = 0.5)
	f
end
  ╠═╡ =#

# ╔═╡ 05558267-7c22-4480-99f7-847b02a273cd
import StatsPlots

# ╔═╡ a1052731-3d45-42e5-b851-9dc0a8cf8177
import Distributions

# ╔═╡ ce1658a1-67fd-4743-90f3-9dc914d857be
StatsPlots.qqplot(Distributions.Uniform(), -log10.(annot_peaks[annot_peaks.GMM_chi2_pvalue .!== missing, :GMM_chi2_pvalue]), qqline = :R)

# ╔═╡ 749d91db-cebd-43f2-a22b-248fd9b8587d
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT q-value / coverage",
				    xlabel = "Total coverage in both conditions",
					# yticks = -9:1:9,
				    ylabel = "-log10(q-value)",
				    xscale = log10)
	local df = annot_peaks
	# local m = df.GMM_LOR
	# local a = log10.(df.mean_cov)
	scatter!(ax, df.mean_cov, -log10.(df.GMM_chi2_qvalue),
			 color = ifelse.(-log10.(df.GMM_chi2_qvalue) .>= 2, :red, :lightgrey),
			 markersize = 5,
			 alpha = 0.65)
	
	Legend(f[2, 1],
		   [PolyElement(color = :lightgrey),
			PolyElement(color = :red)],
	       ["q-value > 0.01", "q-value ≤ 0.01"],
		   orientation = :horizontal,
		   framevisible = false)
	
	f
end

# ╔═╡ 730a04cc-b8db-406a-94ad-f951baf08adf
# ╠═╡ disabled = true
#=╠═╡
begin
	local glori1 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv"))
	local glori2 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GSM6432591_293T-mRNA-2_35bp_m2.totalm6A.FDR.csv"))

	local common = innerjoin(glori1, glori2,
			  				 on = [:Chr, :Strand, :Sites],
			  				 renamecols = "_1" => "_2")
	# Convert to 0-based indexing to conform to the BED file format standard
	common[:, :Sites] .-= 1
	common[:, :ratio] = (common.Ratio_1 .+ common.Ratio_2)./2
	common[:, :end] = common.Sites .+ 1

	rename!(common,
		    :Chr => :chr,
		    :Sites => :start,
		    :Strand => :strand,
		    :Gene_1 => :gene)
	
	common = common[:, [:chr, :start, :end, :gene, :ratio, :strand]]
	common |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ 8805231b-7e9f-4dec-a999-ac8e750a49b9
# ╠═╡ disabled = true
#=╠═╡
begin
	local glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
	glori[:, :modified] .= 1
	glori = glori[:, [:chr, :strand, :start, :modified, :ratio]]
	local ref002 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx.tsv"))
	local annot_ref = leftjoin(ref002, glori,
							   on = [:chr, :strand, :genomicPos => :start])
	annot_ref[:, :modified] = ifelse.(annot_ref.modified .=== missing, 0, 1)
	annot_ref[:, :ratio] = ifelse.(annot_ref.ratio .=== missing, 0, annot_ref.ratio)
	sort(annot_ref, [:ref_id, :pos]) |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ a6b2b943-c350-4d99-9124-cfaaed1ddd85
# ╠═╡ disabled = true
#=╠═╡
begin
	local glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
	glori[:, :modified] .= 1
	glori = glori[:, [:chr, :strand, :start, :modified, :ratio]]
	local ref004 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx.tsv"))
	local annot_ref = leftjoin(ref004, glori,
							   on = [:chr, :strand, :genomicPos => :start])
	annot_ref[:, :modified] = ifelse.(annot_ref.modified .=== missing, 0, 1)
	annot_ref[:, :ratio] = ifelse.(annot_ref.ratio .=== missing, 0, annot_ref.ratio)
	sort(annot_ref, [:ref_id, :pos]) |> CSV.write("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv", delim = '\t')
end
  ╠═╡ =#

# ╔═╡ 9cb352a6-ae66-49a5-8104-9dfb5efeaf33
md"""
## Plot
"""

# ╔═╡ e75d091d-a8e2-47c6-b6a6-1564bb0e196c
an_sig_peaks_ivt.mod |> countmap

# ╔═╡ 215df8c8-90a0-4edb-8a19-158c4d2ecafe
data(an_sig_peaks_ivt) * mapping(:mod_ratio, layout = :mod) * AlgebraOfGraphics.density() |> draw

# ╔═╡ e8d917c2-269f-4fa8-8aa5-98761db2168f
md"""
## Modkit
"""

# ╔═╡ 99f3cb7a-6da6-46d0-a21a-81142b26d81f
wt1mk = CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/modkit_mods/WT_1_mods.bed", delim='\t', header=[:ref_id, :pos, :end, :mod, :score, :strand, :start, :end2, :color, :nvalid, :ratio, :nmod, :ncanonical, :nother, :ndelete, :nfail, :ndiff, :nnocall], skipto=2) |> DataFrame

# ╔═╡ 5e1553f2-370d-4529-b9af-16764a2018d5
begin
	local modkit_mods = Dict(zip(wt1mk.ref_id, wt1mk.pos) .=> wt1mk.mod)
	(map(r -> if r.mod !== missing
			r.mod
		else
			get!(modkit_mods, (r.ref_id, r.pos), "")
		end,
		eachrow(an_sig_peaks_ivt)) .== "") |> countmap
end

# ╔═╡ e66bccc8-3650-45cc-93a3-b131d28aacb4
wt2mk = CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/modkit_mods/WT_2_mods.bed", delim='\t', header=[:ref_id, :pos, :end, :mod, :score, :strand, :start, :end2, :color, :nvalid, :ratio, :nmod, :ncanonical, :nother, :ndelete, :nfail, :ndiff, :nnocall]) |> DataFrame

# ╔═╡ 5ec60c75-033c-46d1-9cd1-28eb515c1da6
wt3mk = CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/modkit_mods/WT_SPK_mods.bed", delim='\t', header=[:ref_id, :pos, :end, :mod, :score, :strand, :start, :end2, :color, :nvalid, :ratio, :nmod, :ncanonical, :nother, :ndelete, :nfail, :ndiff, :nnocall]) |> DataFrame

# ╔═╡ dc4f92af-2a11-4d75-b0f0-bf49b88ea345
# wtmk = combine(groupby(filter(r -> r.nvalid >= 30, vcat(wt1mk, wt2mk, wt3mk)), [:ref_id, :pos, :end]),
# 			   :mod => unique => :mod,
# 			   :ratio => mean => :ratio)
begin
	wtmk = filter(r -> r.nvalid >= 30,
					  combine(groupby(vcat(wt1mk, wt2mk, wt3mk), [:ref_id, :pos]),
							  :mod => unique => :mod,
							  :nvalid => sum => :nvalid,
							  :nmod => sum => :nmod))
	wtmk[!, :ratio] = wtmk.nmod ./ wtmk.nvalid
	wtmk
end

# ╔═╡ e38e88da-a4bb-4bed-9402-d878af4b454e
ivt1mk = CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/modkit_mods/IVT_tailed_1_mods.bed", delim='\t', header=[:ref_id, :pos, :end, :mod, :score, :strand, :start, :end2, :color, :nvalid, :ratio, :nmod, :ncanonical, :nother, :ndelete, :nfail, :ndiff, :nnocall], skipto=2) |> DataFrame

# ╔═╡ c82eaa43-29f4-4a05-bfcb-8f8e30cd86cc
data(innerjoin(wt1mk[wt1mk.nvalid .> 20 .&& wt1mk.ratio .>= 30, :], wt2mk[wt2mk.nvalid .> 20 .&& wt2mk.ratio .>= 30, :], on=[:ref_id, :pos], makeunique=true)) *
	mapping(:ratio => "Ratio WT1", :ratio_1 => "Ratio WT2") *
	AlgebraOfGraphics.density() |> draw

# ╔═╡ 710c7efd-53fd-45a0-83bb-7571e78110ad
wt1mk[wt1mk.nvalid .>= 30, :mod] |> countmap

# ╔═╡ bff0c718-491d-4936-9b73-c95406488576
begin
	local df = leftjoin(an_sig_peaks_ivt, wt1mk,
						on=[:ref_id, :pos],
						makeunique=true)
	map(r -> (r.mod, r.mod_1), eachrow(df)) |> countmap |> sort
end

# ╔═╡ 3c99cb4e-1d7a-4e81-889c-231bb63fea07
begin
	local df = innerjoin(wt1mk[wt1mk.nvalid .>= 30, :], ivt[ivt.predicted .> 2, :],
						 on=[:ref_id, :pos],
						 makeunique=true)

	# data(df[(df.mod_ratio .* 100) ./ df.ratio .> 1.1, :]) *
	# 	mapping(:ratio => "Modkit ratio", :mod_ratio => "Nanocompore ratio") *
	# 	visual(Scatter, markersize=5) |> draw
	println(nrow(df))
	
	local c = round(cor(df.mod_ratio, df.ratio), digits = 2)
	
	(data(df) *
		mapping(:ratio => (x -> x/100) => "Modkit ratio", :mod_ratio => "Nanocompore ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 	 data(DataFrame(x = [0, 1], y = [0, 1])) *
		 mapping(:x, :y) *
		 visual(Lines; linestyle = :dash)
	) |> draw(axis = (; title = "Mod. ratio correlation between Modkit and Nanocompore: $c", xlabel = "Modkit mod rato", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 954a0a6c-9df4-40ed-afa0-3fe8fdc7cc1a
begin
	local df = innerjoin(wt1mk[wt1mk.nvalid .>= 30, :], ivt[ivt.predicted .> 2, :],
						 on=[:ref_id, :pos],
						 makeunique=true)

	df = df[(df.mod_ratio .* 100) ./ df.ratio .> 1.5, :]

	maybe_bad_modkit = df[(df.mod_ratio .* 100) ./ df.ratio .> 1.5, :]

	# data(df[(df.mod_ratio .* 100) ./ df.ratio .> 1.1, :]) *
	# 	mapping(:ratio => "Modkit ratio", :mod_ratio => "Nanocompore ratio") *
	# 	visual(Scatter, markersize=5) |> draw
	local c = round(cor(df.mod_ratio, df.ratio), digits = 2)
	
	(data(df) *
		mapping(:ratio => (x -> x/100) => "Modkit ratio", :mod_ratio => "Nanocompore ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 	 data(DataFrame(x = [0, 1], y = [0, 1])) *
		 mapping(:x, :y) *
		 visual(Lines; linestyle = :dash)
	) |> draw(axis = (; title = "Mod. ratio correlation between Modkit and Nanocompore: $c", xlabel = "Modkit mod rato", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 0527aada-3d52-4b73-9849-4e7c47c58d6f
begin
	local df = innerjoin(wt1mk, ivt[ivt.predicted .> 2, :],
						 on=[:ref_id, :pos],
						 makeunique=true)
	df.mod |> countmap
end

# ╔═╡ a4f502ac-a048-481c-b497-18725863e2b7
begin
	# local df = innerjoin(select(wt1mk[wt1mk.nvalid .>= 30, :], Not(:strand)), storm[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
	# 	  on=[:ref_id, :pos])
	# df = df[df.mod .== "a", :]
	# println(nrow(df))
	local df = innerjoin(ivt[ivt.predicted .> 2, :], glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.mod_ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.mod_ratio .- df.glori_mean_ratio)))
	
	(data(df) * mapping(:glori_mean_ratio, :mod_ratio => "IVT ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Nanocompore IVT and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Nanocompore mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 705c8dac-fcc0-4d65-af7a-20c0c78beb3f
begin
	local df = innerjoin(select(wt1mk[wt1mk.nvalid .>= 30, :], Not(:strand)), storm[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
		  on=[:ref_id, :pos])
	# df = df[df.mod .== "a", :]
	println(nrow(df))
	df = innerjoin(df, glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.ratio .- df.glori_mean_ratio)))
	
	(data(df) * mapping(:glori_mean_ratio, :ratio => (x->x/100) => "Modkit ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Modkit and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Modkit mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ bfe2b3ed-16fc-439f-87ff-25d74237c300
begin
	local f = Figure(size=(900, 400))
	
	local df = innerjoin(ivt[ivt.predicted .> 2, :], glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.mod_ratio), digits = 2)
	println(nrow(df))
	println("Nanocompore Pearson corr: ", c)
	println("Nanocomopre MAE: ", mean(abs.(df.mod_ratio .- df.glori_mean_ratio)))
	
	local plt1 = (data(df) * mapping(:glori_mean_ratio, :mod_ratio => "IVT ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash))
	draw!(f[1, 1], plt1, axis = (; title = "Nanocompore (WT/IVT) vs. GLORI", xlabel = "GLORI mod. ratio", ylabel = "Nanocompore mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

	
	local df = innerjoin(wtmk, storm[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
		  on=[:ref_id, :pos])
	# df = df[df.mod .== "a", :]
	println(nrow(df))
	df = innerjoin(df, glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.ratio), digits = 2)
	println(nrow(df))
	println("Modkit Pearson corr: ", c)
	println("Modkit MAE: ", mean(abs.(df.ratio .- df.glori_mean_ratio)))
	
	local plt2 = (data(df) * mapping(:glori_mean_ratio, :ratio => "Modkit ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash))
	draw!(f[1, 2], plt2, axis = (; title = "Modkit (WT) vs GLORI", xlabel = "GLORI mod. ratio", ylabel = "Modkit mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

	for (label, layout) in zip(["A", "B"],
							   [f[1, 1], f[1, 2]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 16,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
end

# ╔═╡ 07c66f25-2f0c-49db-be23-3e26b34be383
begin
	local df = innerjoin(select(wt1mk[wt1mk.nvalid .>= 30, :], Not(:strand)), storm[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
		  on=[:ref_id, :pos])

	df = innerjoin(df, ivt[ivt.predicted .> 2, [:ref_id, :pos]], on=[:ref_id, :pos])
	# df = df[df.mod .== "a", :]
	println(nrow(df))
	df = innerjoin(df, glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.ratio), digits = 2)
	println(nrow(df))
	println(c)
	
	(data(df) * mapping(:glori_mean_ratio, :ratio => (x->x/100) => "Modkit ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Modkit and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Modkit mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 66397d80-9c61-4ccb-8ec9-a161b3a9eb59
begin
	local df = innerjoin(select(wt1mk[wt1mk.nvalid .>= 30, :], Not(:strand)), storm[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
		  on=[:ref_id, :pos])

	df = antijoin(df, maybe_bad_modkit, on=[:ref_id, :pos])
	# df = df[df.mod .== "a", :]
	println(nrow(df))
	df = innerjoin(df, glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.ratio), digits = 2)
	println(nrow(df))
	println(c)
	
	(data(df) * mapping(:glori_mean_ratio, :ratio => (x->x/100) => "Modkit ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash)) |> draw(axis = (; title = "Mod. ratio correlation between Modkit and GLORI: $c", xlabel = "GLORI mean mod ratio", ylabel = "Modkit mod ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))
end

# ╔═╡ 98121a40-0e41-42a7-ac26-d7cae9a68273
md"""
## New comods
"""

# ╔═╡ b6d27714-b4cb-4ed9-99d1-1d361814bd2d


# ╔═╡ c3ff86aa-3c15-4564-9153-b4e75306ed54
an_sig_peaks_ivt[between.(an_sig_peaks_ivt.genomicPos, 1344933-5, 1344933+5), :]

# ╔═╡ e2a597a1-db29-46ba-83dc-514f5a8e8f9f
import CategoricalArrays

# ╔═╡ b1fbabea-5f72-4069-a90a-1461b8f32dbc
import MultivariateStats

# ╔═╡ ef9b0edc-6a46-4f04-a611-287c94545a56
import UMAP

# ╔═╡ 631f6203-2979-47d7-a16c-047694bab32d
function assoc_type(r00, r01, r10, r11)
	if abs(r00 + r11) > abs(r01 + r10)
		ifelse(r00 + r11 > 0, :comod, :excl)
	else
		ifelse(r01 + r10 > 0, :excl, :comod)
	end
end

# ╔═╡ 6d436275-9bce-4e53-a6c5-9e59acd12fda
begin
	# comods = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/comods_50perc_always_fisher.tsv"))
	comods = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/comod_50perc_contingency.tsv")) # this is also always fisher
	comods[!, :pos1] .+= 4
	comods[!, :pos2] .+= 4
	comods[!, :chi2_pvalue] = comods.pvalue
	# comods[!, :fisher_pvalue] = [p == "missing" ? missing : parse(Float64, p)
	# 							 for p in comods.fisher_pvalue]
	comods[!, :pvalue] = ifelse.(comods.fisher_pvalue .!== missing,
								 comods.fisher_pvalue,
								 comods.pvalue)
	comods[!, :pvalue] .+= 1e-300
	comods[!, :relation] = map(assoc_type, comods.residual00, comods.residual01, comods.residual10, comods.residual11)
	comods[!, :phi] .*= ifelse.(comods.relation .== :comod, 1, -1)
	comods[:, :qvalue] = MultipleTesting.adjust(comods.pvalue, MultipleTesting.BenjaminiHochberg())
	# comods = comods[.! isnan.(comods.pvalue), :]
	comods = comods[abs.(comods.pos1 .- comods.pos2) .> 8, :]
	comods

	# add genomic positions
	comods = leftjoin(comods, ivt[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
				      on = [:reference => :ref_id, :pos1 => :pos])
	rename!(comods, :genomicPos => :genomicPos1)
	comods = leftjoin(comods, ivt[:, [:ref_id, :pos, :genomicPos]],
				      on = [:reference => :ref_id, :pos2 => :pos])
	rename!(comods, :genomicPos => :genomicPos2)
	
	comods
end

# ╔═╡ 6e6dc490-df49-4a4e-ba3a-11667920280c
begin
	# local p1 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	# rename!(p1, :mod => :mod1)
	# local p2 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	# rename!(p2, :mod => :mod2)

	# an_sig_comod = leftjoin(leftjoin(sig_comod, p1,
	# 	 	 		        		 on = [:reference => :ref_id, :pos1 => :pos]),
	# 		     		    p2,
	# 		 				on = [:reference => :ref_id, :pos2 => :pos])

	# local mods = copy(unique(an_ivt[:, [:ref_id, :pos, :mod]]))
	local mods = copy(unique(an_sig_peaks_ivt[:, [:ref_id, :pos, :mod]]))

	mods = Dict(zip(mods.ref_id, mods.pos) .=> mods.mod)

	local df = copy(comods)
	df[!, :mod1] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos1 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))
	df[!, :mod2] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos2 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))

	an_comods = df
end

# ╔═╡ 5cdc22db-0de2-4e82-a31b-57569f9427cb
sig_comods = comods[comods.qvalue .< 0.05, :]

# ╔═╡ 9393db10-0fc1-42a3-8447-f963c2052290
begin
	# local p1 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	# rename!(p1, :mod => :mod1)
	# local p2 = copy(unique(an_peaks_ivt[:, [:ref_id, :pos, :mod]]))
	# rename!(p2, :mod => :mod2)

	# an_sig_comod = leftjoin(leftjoin(sig_comod, p1,
	# 	 	 		        		 on = [:reference => :ref_id, :pos1 => :pos]),
	# 		     		    p2,
	# 		 				on = [:reference => :ref_id, :pos2 => :pos])

	# local mods = copy(unique(an_ivt[:, [:ref_id, :pos, :mod]]))
	local mods = copy(unique(an_sig_peaks_ivt[:, [:ref_id, :pos, :mod]]))

	mods = Dict(zip(mods.ref_id, mods.pos) .=> mods.mod)

	local df = copy(sig_comods)
	df[!, :mod1] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos1 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))
	df[!, :mod2] = map(r -> coalesce(map(o -> get!(mods, (r.reference, r.pos2 + o), missing), [0, 1, -1, 2, -2, 3, -3, 4, -4])...),
					   eachrow(df))


	# Annotate if we're quite certain that the mods are not m6A
	# (not in GLORI and not a DRACH)
	local not_m6A = high_conf_not_m6A(sig_comods)
	df[!, :pos1_not_m6A] = first.(not_m6A)
	df[!, :pos2_not_m6A] = last.(not_m6A)

	an_sig_comods = df
end

# ╔═╡ 60a54f27-165d-4ff1-a526-4f9c957800ef
begin
	local df = an_sig_comods[an_sig_comods.mod1 .!== missing .&& an_sig_comods.mod2 .!== missing, :]
	df[!, :pair] = map(mods -> sort(map(string, collect(mods))), zip(df.mod1, df.mod2))
	df = df[in.(df.pair, Ref(Set([["m6A", "m7G"], ["Y", "m5C"], ["Y", "m6A"], ["m5C", "m6A"]]))), :]
	# filter(p -> p[2] > 100 && p[1][1] != p[1][2], countmap(df.pair))
	
	local m = Matrix(df[:, [:residual00, :residual01, :residual10, :residual11]]) |> transpose
	local proj = UMAP.umap(m; n_neighbors = 15, min_dist = 0.1) |> permutedims
	proj = DataFrame(proj, [:umap1, :umap2])
	proj[!, :pair] = map(a -> join(a, "-") , df.pair)
	data(proj) *
		mapping(:umap1, :umap2, color = :pair) *
		visual(Scatter, markersize = 4) |> draw
end

# ╔═╡ aceaef06-1e07-4059-a228-628410ffdf50
begin
	local df = copy(an_sig_comods)
	df[!, :type] = ifelse.(df.mod1 .=== df.mod2,
						   ifelse.(df.mod1 .!== missing, :missing, df.mod1),
						   :different)
	# df = df[df.type .!== :missing, :]

	local f = Figure(size = (800, 400))
	
	local p1 = (data(DataFrame(x=[-20, 40], y=[-20, 40])) * mapping(:x, :y) * visual(Lines; color = :red)
		+
	 data(df) * mapping(:residual00, :residual11) * visual(Scatter; markersize = 4))
	draw!(f[1, 1], p1; axis = (; xlabel = "Residual 0-0", ylabel = "Residual 1-1", aspect = 1))
	
	local p2 = (data(DataFrame(x=[-20, 40], y=[-20, 40])) * mapping(:x, :y) * visual(Lines; color = :red)
		+
	 data(df) * mapping(:residual01, :residual10) * visual(Scatter; markersize = 4))
	draw!(f[1, 2], p2; axis = (; xlabel = "Residual 0-1", ylabel = "Residual 1-0", aspect = 1))

	f
end

# ╔═╡ fd2de26f-49b8-4171-acfc-7758152e4af1
begin
	local r00 = CategoricalArrays.cut(sig_comods.residual00, [-Inf, -1, 1, Inf], labels=["strong-neg", "weak", "strong-pos"])
	local r01 = CategoricalArrays.cut(sig_comods.residual01, [-Inf, -1, 1, Inf], labels=["strong-neg", "weak", "strong-pos"])
	local r10 = CategoricalArrays.cut(sig_comods.residual10, [-Inf, -1, 1, Inf], labels=["strong-neg", "weak", "strong-pos"])
	local r11 = CategoricalArrays.cut(sig_comods.residual11, [-Inf, -1, 1, Inf], labels=["strong-neg", "weak", "strong-pos"])
	zip(r00, r01, r10, r11) |> countmap
end

# ╔═╡ a0bfd280-8ff7-4376-a00f-aee71ae3f6eb
begin
	local pcres = fit(MultivariateStats.PCA, Matrix(sig_comods[:, [:residual00, :residual01, :residual10, :residual11]]); maxoutdim = 2)
	data(DataFrame(pcres.proj, [:PC1, :PC2])) * mapping(:PC1, :PC2) * visual(Scatter; markersize = 3) |> draw
end

# ╔═╡ 8e63a621-9eb5-4c8b-a7cd-74ad038e9814
begin
	local df = copy(sig_comods)
	local f = Figure(size = (900, 600))

	local cols = [:residual00, :residual01, :residual10, :residual11]
	local plot = 0
	for (i, col) in enumerate(cols)
		for col2 in cols[i+1:end]
			xs = [-30, 30]
			if (col, col2) in [(:residual00, :residual11), (:residual01, :residual10)]
				ys = [-30, 30]
			else
				ys = [30, -30]
			end

			plt = (data(DataFrame(x=xs, y=ys)) * mapping(:x, :y) * visual(Lines; color = :red)
				+
			 data(df) * mapping(col, col2) * visual(Scatter; markersize = 4))
			j = div(plot, 3) + 1
			k = mod(plot, 3) + 1
			draw!(f[j, k], plt; axis = (; xlabel = string(col), ylabel = string(col2), aspect = 1))
			plot += 1
		end
	end
	f
end

# ╔═╡ ad08769c-6387-4f68-b87e-8a085e50a453
begin
	local df = copy(sig_comods)
	df[!, :high00] = df.residual00 .> 10 .&& df.residual11 .< 10
	df[df.high00, :]
	# data(df) * mapping(:high00, (:pos1, :pos2) => ((p1, p2) -> abs(p1 - p2))) * visual(BoxPlot; show_outliers = false) |> draw
end

# ╔═╡ e6c1302f-dbdc-4cce-b9ca-f5f09ca60eed
begin
	local df = copy(sig_comods)
	df[!, :high00] = df.residual00 .> 10 .&& df.residual11 .< 10
	data(df) * mapping(:high00, (:pos1, :pos2) => ((p1, p2) -> abs(p1 - p2))) * visual(BoxPlot; show_outliers = false) |> draw
end

# ╔═╡ 2b2f4216-41fa-4f7f-b707-35bbcec979df
function get_pair(mod1, mod2)
	join(sort([mod1, mod2]), "-")
end

# ╔═╡ 9ec908a7-d59d-477f-80fe-8765dde8019e
# ╠═╡ disabled = true
#=╠═╡
begin
	local df = copy(comods)
	# df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	data(df) *
		mapping(:phi, :pvalue => (p -> -log10(p))) *
		visual(Scatter; markersize = 3) |> draw(; axis = (; xlabel = "φ", ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", limits = ((-1, 1), (nothing, nothing))))
end
  ╠═╡ =#

# ╔═╡ 24b4949a-8947-4be4-af8b-b5979dda3c30
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure(size = (600, 450))
	local sig_color = "#800080" # "#051094"

	local df = copy(comods)
	df[!, :significant] = df.qvalue .< 0.01
	local plt = data(df) *
		mapping(:phi, :pvalue => (p -> -log10(p)), color = :significant) *
		visual(Scatter; markersize = 5)
	draw!(f[1, 1],
		  plt,
		  scales(Color = (; palette = [:grey, sig_color]));
		  axis = (; title = "All modification pairs",
				  	xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
					ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
					limits = ((-1, 1), (-2, 310))))
	f
end
  ╠═╡ =#

# ╔═╡ 2912aef8-60da-4122-9392-e7223512b108
begin
	local df = copy(an_sig_comods)
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	df = df[in.(df.pair, Ref(Set(["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G"]))), :]
	# df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	data(df) *
		mapping(:phi => "φ", :pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", col = :pair) *
		visual(Scatter; markersize = 3.5) |> draw()
end

# ╔═╡ 656e85b6-30e7-4943-b15d-f1f8654fee26
begin
	local df = an_sig_comods[an_sig_comods.mod1 .!== missing .&& an_sig_comods.mod2 .!== missing, :]
	df[!, :pair] = map(mods -> sort(map(string, collect(mods))), zip(df.mod1, df.mod2))
	filter(p -> p[2] > 100 && p[1][1] != p[1][2], countmap(df.pair))
end

# ╔═╡ 94240437-5cb9-474b-9204-7eadb5aaa0a3
# ╠═╡ disabled = true
#=╠═╡
(data(DataFrame(x = [0, 10], y = [0, 10])) * mapping(:x, :y) * visual(Lines, color = :red)
	+
 data(comods) * mapping(:pvalue => (p -> -log10(p -1e-300 + 1e-12)) => "pvalue", :emp_pvalue => (p -> -log10(p)) => "emp pvalue") * visual(Scatter; markersize = 3)) |> draw(; axis = (; xlabel = "-log10(pvalue)", ylabel = "-log10(emp pvalue)"))
  ╠═╡ =#

# ╔═╡ 72b3a744-e8aa-4cd8-a178-579aa6d564d6
begin
	local df = copy(an_sig_comods)
	df[!, :pair] = map((x, y) -> "$x-$y", df.mod1, df.mod2)
	# df[!, :ordered_pair] = map(get_pair, df.mod1, df.mod2)
	# sort!(df, :ordered_pair)
	local selected = map(r -> r.pair in ["m6A-m7G", "m5C-m6A", "Y-m6A", "Y-m5C", "m7G-m6A", "m6A-m5C", "m6A-Y", "m5C-Y"], eachrow(df))
	df[!, :distance] = abs.(df.pos1 .- df.pos2)
	
	# sort(combine(groupby(df[selected, :], [:mods, :relation]), nrow => :count))

	local f = Figure(size = (1000, 300))
	local plt = data(df[selected, :]) *
		mapping(:pair => "Modifications pair",
				:distance => "Distance",
				color=:relation => "Association",
				dodge=:relation) *
		visual(BoxPlot; show_outliers = false)
	draw!(f[1, 1], plt)
	f
end

# ╔═╡ dc29fb1e-82b7-4db1-8b18-7d63f63c3395
begin
	local df = copy(an_sig_comods)
	df[!, :pair] = map((x, y) -> "$x-$y", df.mod1, df.mod2)
	df[!, :ordered_pair] = map(get_pair, df.mod1, df.mod2)
	# sort!(df, :ordered_pair)
	local selected = map(r -> r.pair in ["m6A-m7G", "m5C-m6A", "Y-m6A", "Y-m5C", "m7G-m6A", "m6A-m5C", "m6A-Y", "m5C-Y"], eachrow(df))
	df[!, :distance] = abs.(df.pos1 .- df.pos2)
	
	# sort(combine(groupby(df[selected, :], [:mods, :relation]), nrow => :count))

	local f = Figure(size = (1000, 300))
	local plt = data(df[selected, :]) *
		mapping(:pair => "Modifications pair",
				:distance => "Distance",
				color=:relation => "Association",
				dodge=:relation,
			    layout=:ordered_pair) *
		visual(BoxPlot; show_outliers = false)
	draw!(f[1, 1], plt)
	f
end

# ╔═╡ 61a0ec53-8374-4399-9600-87061c8f2c61
begin
	local df = copy(an_sig_comods)
	df[!, :pair] = map((x, y) -> "$x-$y", df.mod1, df.mod2)
	df[!, :ordered_pair] = map(get_pair, df.mod1, df.mod2)
	df[!, :distance] = abs.(df.pos1 .- df.pos2)

	local fig = Figure(size = (640, 640))
	for (i, ordered_pair) in enumerate(["m5C-m6A", "Y-m6A", "m6A-m7G", "Nm-m6A"])
		subset = df[df.ordered_pair .== ordered_pair, :]
		counts = countmap(subset.pair)
		subset[!, :pair] = map(p -> "$p ($(counts[p]))", subset.pair)
		println(ordered_pair, " ", nrow(subset))
		# println(:comod, " ", subset[subset.relation .== :comod, :distance] |> mean)
		println(combine(groupby(subset, [:pair, :relation]), :distance => median) |> sort)
		plt = data(subset) *
			mapping(:pair, :distance, color = :relation, dodge = :relation) *
			visual(BoxPlot; show_outliers = false)
		draw!(fig[(i - 1) % 2 + 1, div(i - 1, 2) + 1], plt;
			  axis = (;
					  yticks = 0:100:1100,
					  aspect = 1,
					  limits = (nothing, (0, 1100))))
	end
	Legend(fig[3, 1:2],
		   [PolyElement(color=:blue),
			PolyElement(color=:orange)],
	       ["Co-occurring", "Mutually exclusive"],
		   orientation = :horizontal,
		   framevisible = false)

	fig
end

# ╔═╡ ffd6defb-70f5-4366-a42a-588f37d25c66
an_sig_comods[occursin.("RPL34", an_sig_comods.reference), :]

# ╔═╡ f1e9ce86-974d-41cc-8cc0-13a90f20fc8b
begin
	local df = copy(an_sig_comods)
	df[!, :gene] = map(x -> x[6], split.(df.reference, "|"))
	local df = combine(groupby(df, [:chr, :strand, :genomicPos1, :genomicPos2, :gene]),
					   :phi => minimum => :min_phi,
					   :phi => maximum => :max_phi,
					   :mod1 => first => :mod1,
					   :mod2 => first => :mod2,
					   :reference => Set => :references)
	df[!, :phi_diff] = df.max_phi .- df.min_phi
	# df[!, :transcript] = map(x -> x[5], split.(df.reference, "|"))
	df = sort(df, [:phi_diff], rev = true)
	# df = df[df.mod1 .!== missing .&& df.mod2 .!== missing .&& df.mod1 .!= df.mod2, :]
	# df |> CSV.write("/home/mzdravkov/comod_changing_association.tsv", delim = '\t')
	df
end

# ╔═╡ b84d3ad5-f11e-4738-af7c-09c92c8d9de0


# ╔═╡ 282cd81c-43b0-4122-8f8b-6b6d9ff14f35
an_sig_comods[an_sig_comods.genomicPos1 .=== 76066564 .&& an_sig_comods.genomicPos2 .=== 76066573, :]

# ╔═╡ de2c4249-900d-4fee-941d-5a865ba1b296
an_sig_peaks_ivt[an_sig_peaks_ivt.genomicPos .== 76066564 .|| an_sig_peaks_ivt.genomicPos .== 76066573, :]

# ╔═╡ 8a95d8d5-daa6-4efc-92d5-35b363660dd7
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure(size = (1200, 1450))
	local ga = f[1, 1] = GridLayout()
	local gb = f[2, 1] = GridLayout()
	local gc = f[3, 1] = GridLayout()
	


	local df = renamemods(an_sig_comods)

	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	df = df[in.(df.pair, Ref(Set(["m⁵C-m⁶A", "m⁵C-Ψ", "m⁶A-Ψ", "m⁶A-m⁷G"]))), :]
	# df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	local plt = data(df) *
		mapping(:phi => "φ", :pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", col = :pair) *
		visual(Scatter; markersize = 5)
	draw!(ga[1, 1:2], plt)


	local gpos1 = 76066564
	local gpos2 = 76066573
	local row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
	local row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
	local tx1 = split(row1.reference, "|")[1]
	local tx2 = split(row2.reference, "|")[1]
	local name1 = split(row1.reference, "|")[5]
	local name2 = split(row2.reference, "|")[5]
	
	local ax = Axis(gb[1, 1:2],
				    height = 110)
	plot_isoforms_model!(ax, "MDH"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => :blue, tx2 => :orange]), focus = (gpos1-150, gpos2+200), rename = Dict([tx1 => name1, tx2 => name2]))
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.(tx1, an_sig_comods.reference)
	local tx2_mask = occursin.(tx2, an_sig_comods.reference)
	local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	arc_plot_genomic2(gb[2, 1],
					  renamemods(an_sig_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
					  range=(row1.pos1-2000, row1.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name1,
					  highlight=(row1.pos1, row1.pos2),
					  spinecolor=:blue)
	arc_plot_genomic2(gb[3, 1],
					  renamemods(an_sig_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
					  range=(row2.pos1-2000, row2.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name2,
					  highlight=(row2.pos1, row2.pos2),
					  spinecolor=:orange)
	Colorbar(gb[2:3, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	# local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000304494.10|ENSG00000147889.18|OTTHUMG00000019686.7|OTTHUMT00000051915.1|CDKN2A-201|CDKN2A|978|protein_coding|")
	# # local f = Figure(size = (600, 600))
	

	# local df = DataFrame(Tables.table(probs[:, [1385, 1584]] .> 0.75))
	# dropmissing!(df)
	# local n2 = nrow(df)
	# println(contingency(df.Column1, df.Column2))
	# local ax1 = Axis(gc[1, 1],
	# 				 title="YWHAE-201",
	# 			     # aspect = 0.3,
	# 			     yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 				 xticks=(0:(n2/10):n2, ["$x%" for x in 0:10:100]),
	# 			     xlabel="")
	# println("n2: $n2")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax1, Array(sort(df)), colormap=["#dbdbdb", :black])

	# local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000571732.5|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000437356.1|YWHAE-208|YWHAE|1818|protein_coding|")
	# # local f = Figure(size = (600, 600))

	# local df = DataFrame(Tables.table(probs[:, [1416, 1615]] .> 0.75))
	# dropmissing!(df)
	# local n1 = nrow(df)
	# local ax2 = Axis(gc[2, 1],
	# 			 title="YWHAE-208",
	# 			 # aspect = 0.3,
	# 			 yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 			 xticks=(0:(n1/10):n1, ["$x%" for x in 0:10:100]),
	# 			 xlabel="Reads")
	# println("n1: $n1")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax2, Array(sort(df)), colormap=["#dbdbdb", :black])


	
	
	# xlims!(ax1, (0, max(n1, n2)))
	# xlims!(ax2, (0, max(n1, n2)))
	

	# colsize!(f.layout, 2, Relative(7/10))
	# colsize!(f.layout, 3, Relative(1/10))
	rowsize!(f.layout, 1, Relative(1.8/7))
	rowsize!(f.layout, 2, Relative(3.9/7))
	rowsize!(f.layout, 3, Relative(1.3/7))

	rowsize!(gb, 2, Relative(2/5))
	rowsize!(gb, 3, Relative(2/5))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	Legend(gb[4, 1],
		   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π])],
	       ["Mutually exclusive", "Co-modified"],
		   orientation = :horizontal,
		   framevisible = false)
	Legend(gc[3, 1],
		   [PolyElement(color=:black),
			PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	       ["Modified", "Not modified"],
		   orientation = :horizontal,
		   framevisible = false)


	for (label, layout) in zip(["A", "B", "C"],
							   [ga[1, 1:2], gb[1, 1], gc[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end
  ╠═╡ =#

# ╔═╡ 879eaae7-21a2-40e7-9c31-ab9603e4947b
function contingency(a, b)::Matrix{Integer}
    [sum(a .== 0 .&& b .== 0) sum(a .== 0 .&& b .== 1)
     sum(a .== 1 .&& b .== 0) sum(a .== 1 .&& b .== 1)]
end

# ╔═╡ 6df047c1-f3a7-44ac-a444-f470e17d49a4
contingency(Int.(comods.fisher_pvalue .< 0.05), Int.(comods.emp_pvalue .< 0.05))

# ╔═╡ 3403dfc2-27fd-439d-875f-b981b946d603
begin
	local cont = contingency(Int.(comods.fisher_pvalue .< 0.05), Int.(comods.emp_pvalue .< 0.05))

	
	# (sum(cont[:, 1]), sum(cont[:, 2]), sum(cont[1, :]), sum(cont[2, :]))
	local mcc = (cont[1, 1]*cont[2, 2] - cont[1, 2]*cont[2, 1])/prod(sqrt.(cont))
end

# ╔═╡ a88a543d-d11d-4941-b588-469b67a8f1cc
contingency(Int.(ifelse.(comods.fisher_pvalue .!== missing, comods.fisher_pvalue, comods.pvalue) .< 0.01), Int.(comods.emp_pvalue .< 0.01))

# ╔═╡ ea8b9460-5973-45c5-ac50-181ef50800da
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure(size = (1200, 1450))
	local ga = f[1, 1] = GridLayout()
	local gb = f[2, 1] = GridLayout()
	local gc = f[3, 1] = GridLayout()


	local df = renamemods(an_sig_comods)
	# df[!, :relation] = map(r -> let s = (r.state1, r.state2)
	# 								if (s in [(1, 1), (0, 0)] && r.observed > r.expected) || (s in [(0, 1), (1, 0)] && r.observed < r.expected)
	# 									:positive
	# 								elseif (s in [(1, 1), (0, 0)] && r.observed < r.expected) || (s in [(0, 1), (1, 0)] && r.observed > r.expected)
	# 									:negative
	# 								else
	# 									:weird
	# 								end
	# 							end,
	# 					   eachrow(df))
	# df[!, "mods"] = map(r -> join(sort([coalesce(r.mod1, ""), coalesce(r.mod2, "")]), "-"),
	# 					 eachrow(df))
	# println(sort(combine(groupby(df, :mods), nrow => :count), :count, rev=true))
	# df[!, :distance] = abs.(df.pos1 .- df.pos2)
	# # local selected = map(r -> r.mods in ["m⁶A-m⁷G", "m⁵C-m⁶A", "m⁶A-Ψ", "m⁵C-Ψ"], eachrow(df))
	# local common = filter(r -> r.count > 40, combine(groupby(df, :mods), nrow => :count)).mods |> Set
	# println(common)
	# local selected = [non_trivial_simple(r) && r.mod1 !== missing && r.mod2 !== missing && r.mods in common && r.mod1 != r.mod2 for r in eachrow(df)]
	
	# # sort(combine(groupby(df[selected, :], [:mods, :relation]), nrow => :count))

	# local plt = data(df[selected, :]) *
	# 	mapping(:mods => "Modification pair", color=:relation, dodge=:relation) *
	# 	frequency()
	# draw!(ga[1, 1], plt)
	
	# local plt = data(df[selected, :]) *
	# 	mapping(:mods => "Modifications pair",
	# 			:distance => "Distance",
	# 			color=:relation => "Association",
	# 			dodge=:relation) *
	# 	visual(BoxPlot)
	# draw!(ga[1, 2], plt)

	# local df = copy(an_sig_comods)
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	df = df[in.(df.pair, Ref(Set(["m⁵C-m⁶A", "m⁵C-Ψ", "m⁶A-Ψ", "m⁶A-m⁷G"]))), :]
	# df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	local plt = data(df) *
		mapping(:phi => "φ", :pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", col = :pair) *
		visual(Scatter; markersize = 5)
	draw!(ga[1, 1:2], plt)

	# local ax = Axis(gb[1, 1:2],
	# 			    height = 110)
	# plot_isoforms_model!(ax, "PRRC2B"; transcripts = Set(["ENST00000682501.1", "ENST00000684596.1", "ENST00000683519.1"]), colors = Dict(["ENST00000682501.1" => :blue, "ENST00000684596.1" => :orange]), focus = (131496129-150, 131496179+200), rename = Dict(["ENST00000682501.1" => "PRRC2B-201", "ENST00000684596.1" => "PRRC2B-208", "ENST00000683519.1" => "PRRC2B-206"]))
	
	# # local ax = Axis(gb[1, 1])
	# local tx1_mask = occursin.("ENST00000682501.1", an_sig_comods.reference)
	# local tx2_mask = occursin.("ENST00000684596.1", an_sig_comods.reference)
	# local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
	# 			  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	# arc_plot_genomic2(gb[2, 1],
	# 				  renamemods(an_sig_comods[tx1_mask, :]),
	# 				  an_sig_peaks_ivt[occursin.("ENST00000682501.1", an_sig_peaks_ivt.ref_id), :];
	# 				  range=(5000, 5300),
	# 				  grange=(131496129-150, 131496179+100),
	# 				  colorrange=colorrange,
	# 				  title="PRRC2B-201",
	# 				  highlight=(5089, 5139),
	# 				  spinecolor=:blue)
	# arc_plot_genomic2(gb[3, 1],
	# 				  renamemods(an_sig_comods[tx2_mask, :]),
	# 				  an_sig_peaks_ivt[occursin.("ENST00000684596.1", an_sig_peaks_ivt.ref_id), :];
	# 				  range=(0, 12000),
	# 				  grange=(131496129-150, 131496179+100),
	# 				  colorrange=colorrange,
	# 				  title="PRRC2B-208",
	# 				  highlight=(7211, 7261),
	# 				  spinecolor=:orange)
	# Colorbar(gb[2:3, 2], limits = colorrange, colormap = :viridis,
 #    vertical = true, label = "-log₁₀(𝑃-value)")
	local gpos1 = 76066564
	local gpos2 = 76066573
	local row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
	local row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
	local tx1 = split(row1.reference, "|")[1]
	local tx2 = split(row2.reference, "|")[1]
	local name1 = split(row1.reference, "|")[5]
	local name2 = split(row2.reference, "|")[5]
	
	local ax = Axis(gb[1, 1:2],
				    height = 110)
	plot_isoforms_model!(ax, "MDH"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => :blue, tx2 => :orange]), focus = (gpos1-150, gpos2+200), rename = Dict([tx1 => name1, tx2 => name2]))
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.(tx1, an_sig_comods.reference)
	local tx2_mask = occursin.(tx2, an_sig_comods.reference)
	local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	arc_plot_genomic2(gb[2, 1],
					  renamemods(an_sig_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
					  range=(row1.pos1-2000, row1.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name1,
					  highlight=(row1.pos1, row1.pos2),
					  spinecolor=:blue)
	arc_plot_genomic2(gb[3, 1],
					  renamemods(an_sig_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
					  range=(row2.pos1-2000, row2.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name2,
					  highlight=(row2.pos1, row2.pos2),
					  spinecolor=:orange)
	Colorbar(gb[2:3, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000264335.13|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000259354.4|YWHAE-201|YWHAE|2052|protein_coding|")
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [1385, 1584]] .> 0.75))
	dropmissing!(df)
	local n2 = nrow(df)
	println(contingency(df.Column1, df.Column2))
	# local ax1 = Axis(gc[1, 1],
	# 				 title="YWHAE-201",
	# 			     # aspect = 0.3,
	# 			     yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 				 xticks=(0:(n2/10):n2, ["$x%" for x in 0:10:100]),
	# 			     xlabel="")
	# println("n2: $n2")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax1, Array(sort(df)), colormap=["#dbdbdb", :black])

	# local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000571732.5|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000437356.1|YWHAE-208|YWHAE|1818|protein_coding|")
	# # local f = Figure(size = (600, 600))

	# local df = DataFrame(Tables.table(probs[:, [1416, 1615]] .> 0.75))
	# dropmissing!(df)
	# local n1 = nrow(df)
	# local ax2 = Axis(gc[2, 1],
	# 			 title="YWHAE-208",
	# 			 # aspect = 0.3,
	# 			 yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 			 xticks=(0:(n1/10):n1, ["$x%" for x in 0:10:100]),
	# 			 xlabel="Reads")
	# println("n1: $n1")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax2, Array(sort(df)), colormap=["#dbdbdb", :black])


	
	
	# xlims!(ax1, (0, max(n1, n2)))
	# xlims!(ax2, (0, max(n1, n2)))
	

	# colsize!(f.layout, 2, Relative(7/10))
	# colsize!(f.layout, 3, Relative(1/10))
	rowsize!(f.layout, 1, Relative(1.8/7))
	rowsize!(f.layout, 2, Relative(3.9/7))
	rowsize!(f.layout, 3, Relative(1.3/7))

	rowsize!(gb, 2, Relative(2/5))
	rowsize!(gb, 3, Relative(2/5))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	Legend(gb[4, 1],
		   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π])],
	       ["Mutually exclusive", "Co-modified"],
		   orientation = :horizontal,
		   framevisible = false)
	Legend(gc[3, 1],
		   [PolyElement(color=:black),
			PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	       ["Modified", "Not modified"],
		   orientation = :horizontal,
		   framevisible = false)


	for (label, layout) in zip(["A", "B", "C"],
							   [ga[1, 1:2], gb[1, 1], gc[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end
  ╠═╡ =#

# ╔═╡ 0919f200-3984-4cab-9ac0-bf6e982ad566
comods[occursin.("CDKN2A", comods.reference) .&& comods.genomicPos1 .== 21968192 .&& comods.genomicPos2 .== 21968114, :]

# ╔═╡ 19894766-0820-47a2-9a0b-9dab37ca5daa
maximum(ivt[occursin.("GNAS", ivt.ref_id), :genomicPos])

# ╔═╡ faf85657-eb0e-4b9e-86db-74a3cc6b1c04
an_sig_comods[occursin.("GNAS", an_sig_comods.reference), :]

# ╔═╡ c107996c-abeb-423f-90a0-6aa5125df163
an_sig_comods[an_sig_comods.genomicPos1 .== 58910358 .&& an_sig_comods.genomicPos2 .== 58910385, :]

# ╔═╡ e2557d68-4418-4266-9e6b-adc99f5e5d9c
sort(combine(groupby(an_sig_comods[occursin.("NOTCH3", an_sig_comods.reference), :], [:genomicPos1, :genomicPos2]), :phi => (phis -> maximum(phis) - minimum(phis)) => :phi_diff), :phi_diff, rev = true)

# ╔═╡ 50e26ae2-bc4b-4324-88f7-8588276cc44a
begin
	gtf = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/gencode.v41.annotation.gtf", delim = '\t', comment = "#", header = [:chr, :source, :feature, :start, :end, :score, :strand, :frame, :attribute]))
	gtf[!, :transcript_id] = map(a -> replace(split(split(a, "; ")[2], " ")[2], "\"" => ""), gtf.attribute)
	gtf
end

# ╔═╡ 213085a9-e9de-4404-80de-a26b2566caa8
gtf[occursin.("PRRC2B", gtf.attribute) .&& gtf.feature .== "exon", :]

# ╔═╡ 4e8b6a0d-40bc-418c-8520-fa9365f43a85
begin
	local df = copy(an_sig_comods)
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	local selected = map(r -> r.pair in ["m6A-m7G", "m5C-m6A", "Y-m6A", "Y-m5C"], eachrow(df))
	df = df[selected, :]
	sort(combine(groupby(df, [:mod1, :mod2]), nrow => :count))
end

# ╔═╡ 0fc62ea4-35bd-410d-a0b9-2ad09140b327
begin
	local df = copy(an_sig_comods)
	local nearfull_refs = filter(r -> r.pos_minimum <= 50, combine(groupby(ivt, :ref_id), :pos => minimum)).ref_id |> Set
	df = df[in.(df.reference, Ref(nearfull_refs)), :]
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	local selected = map(r -> r.pair in ["m6A-m7G", "m5C-m6A", "Y-m6A", "Y-m5C"], eachrow(df))
	df = df[selected, :]
	df[!, :transcript_id] = map(r -> split(r, "|")[1], df.reference)
	local ref = copy(gtf[gtf.feature .== "start_codon", [:transcript_id, :start, :end, :strand]])
	ref[!, :start_codon] = ifelse.(ref.strand .== "+", ref.start, ref.end)
	df = innerjoin(df, ref[:, [:transcript_id, :start_codon]], on = [:transcript_id], makeunique = true)
	local ref = copy(gtf[gtf.feature .== "stop_codon", [:transcript_id, :start, :end, :strand]])
	ref[!, :stop_codon] = ifelse.(ref.strand .== "+", ref.end, ref.start)
	df = innerjoin(df, ref[:, [:transcript_id, :stop_codon]], on = [:transcript_id], makeunique = true)

	# df = df[df.feature .!== missing, :]
	# get_start = row -> row.strand == "+" ? row.start : row.end
	# get_stop = row -> row.strand == "+" ? row.end : row.start

	# local starts = Dict(map(r -> (r.transcript_id, get_start(r)), eachrow(gtf[gtf.feature .== "start_codon", :])))
	# local stops = Dict(map(r -> (r.transcript_id, get_stop(r)), eachrow(gtf[gtf.feature .== "stop_codon", :])))

	get_region = (row, pos) -> if row.strand == "+"
		if pos < row.start_codon
			:utr5
		elseif pos >= row.start_codon && pos <= row.stop_codon
			:cds
		else
			:utr3
		end
	else
		if pos > row.start_codon
			:utr5
		elseif pos <= row.start_codon && pos >= row.stop_codon
			:cds
		else
			:utr3
		end
	end

	df[!, :region1] = map((row, pos) -> get_region(row, pos), eachrow(df), ifelse.(df.strand .== "+", df.genomicPos1, df.genomicPos2))
	df[!, :region2] = map((row, pos) -> get_region(row, pos), eachrow(df), ifelse.(df.strand .== "+", df.genomicPos2, df.genomicPos1))

	# sort(combine(groupby(df, [:region1, :region2]), nrow => :count), by = (r -> Dict(:utr5 => 1, :cds => 2, :utr3 => 3)[r]))
	df = sort(combine(groupby(df, [:pair, :region1, :region2]), nrow => :count))


	function get_matrix(df)
		m = zeros(3, 3)
		for (i, r1) in enumerate([:utr5, :cds, :utr3])
			for (j, r2) in enumerate([:utr5, :cds, :utr3])
				m[i, j] = sum(df[df.region1 .== r1 .&& df.region2 .== r2, :count])
			end
		end
		m
	end
	local f = Figure(size = (1000, 250))
	for (i, pair) in enumerate(unique(df.pair))
		matrix = get_matrix(df[df.pair .== pair, :])
		# matrix ./= sum(matrix)
		ax = Axis(f[1, i],
				  title = pair,
				  xlabel = "Modification 1",
				  ylabel = "Modification 2",
				  xticks = (1:3, ["5' UTR", "CDS", "3' UTR"]),
				  yticks = (1:3, ["5' UTR", "CDS", "3' UTR"]))
		heatmap!(ax, matrix, colormap = :reds)
	end
	f
end

# ╔═╡ e266cb67-ce59-49e5-8ce3-5ca8e5128ded
combine(groupby(filter(r -> r.mod1 !== missing && r.mod2 !== missing && r.mod1 != r.mod2 && r.mod1 in ["m6A", "m7G"] && r.mod2 in ["m6A", "m7G"], an_sig_comods), [:mod1, :mod2]), nrow => :count)

# ╔═╡ 9e4877fe-5c7c-493d-aa24-a5d22fbde4fd
gtf[gtf.feature .== "CDS", :]

# ╔═╡ b5febe6e-e9b5-46b8-8e09-700e28abe37f
gtf.feature |> unique

# ╔═╡ 8334730d-6e9e-46f8-a797-0d4730df3ab4
sum(combine(groupby(ivt, :ref_id), :pos => minimum).pos_minimum .< 50)

# ╔═╡ 33a269c7-9271-4f58-ab63-1a0264cd1965
data(sig_comods) *
	mapping((:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))) => "Distance (log-scale)", :phi => "← Mutually exclusive    φ    Co-occurring →         ") *
	visual(Scatter; markersize = 4) |>
	draw(; axis = (; title = "Effect size of modification pair association by distance",
				   xticks = (map(log10, [(10:10:100)..., (200:100:1000)..., (2000:1000:10000)..., 11000, 12000]), ["10nt", repeat([""], 8)..., "100nt", repeat([""], 8)..., "1,000nt", repeat([""], 8)..., "10,000nt", "", ""]),
				   limits = ((0.9, 4.2), (-1, 1))))

# ╔═╡ cd4f0a10-228a-473d-a8d3-6e7af9fe42da
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure(size = (1600, 1100))
	local trow = f[1, 1] = GridLayout()
	local brow = f[2, 1] = GridLayout()

	local sig_color = "#800080" # "#051094"


	local df = copy(comods)[1:100, :]
	df[!, :significant] = df.qvalue .< 0.01
	local plt = data(df) *
		mapping(:phi, :pvalue => (p -> -log10(p)), color = :significant) *
		visual(Scatter; markersize = 5)
	draw!(trow[1, 1],
		  plt,
		  scales(Color = (; palette = [:grey, sig_color]));
		  axis = (; #title = "All modification pairs",
				  	xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
					ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
					limits = ((-1, 1), (-2, 310))))


	# local df = renamemods(an_comods)
	# df[!, :significant] = df.qvalue .< 0.05

	# # local df = copy(an_sig_comods)
	# df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	# df = df[in.(df.pair, Ref(Set(["m⁵C-m⁶A"]))), :]
	# # df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# # df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	# local plt = data(df) *
	# 	mapping(:phi,
	# 			:pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
	# 		    color = :significant) *
	# 	visual(Scatter; markersize = 5)
	# draw!(trow[1, 2],
	# 	  plt,
	# 	  scales(Color = (; palette = [:grey, sig_color]));
	# 	  axis = (; title = "m⁵C-m⁶A",
	# 			    xlabel = "← Mutually exclusive    𝛗    Co-occurring →       ",
	# 				ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
	# 				limits = ((-1, 1), (-2, 310))))

	# === PLOT B VARIANT 1

	# local df = copy(an_sig_comods)
	# df[!, :type] = sign.(df.phi)
	# df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	# df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	# # df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	# df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	# df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	# df[!, :lfc] = log2.(df.count ./ df.count_1)
	# df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	# df[!, :conf] = 1.96 .* df.SE
	# df[!, :total] = df.count .+ df.count_1
	# df = sort(df[df.count .> 30 .&& df.count_1 .> 30, :], :total)
	# local ax = Axis(trow[1, 2],
	# 				xlabel = "logFC (co-occurring/mutually-exclusive)",
	# 				ylabel = "Modification pair",
	# 			    yticks = (1:nrow(df), map((p, t) -> "$p ($t)", df.pair, df.total)))
	# for (i, row) in enumerate(eachrow(df))
	# 	scatter!(ax, [row.lfc], [i], color = :black)
	# 	lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)
	# end
	# xlims!(ax, [0, 4])

	# === PLOT B VARIANT 2
	
	# local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	# local df = renamemods(an_sig_comods)
	
	# df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# # df[!, :significant] = df.qvalue .< 0.01
	# df[!, :mod1] = ifelse.(df.mod1 .=== missing, "unknown non-m⁶A", df.mod1)
	# df[!, :mod2] = ifelse.(df.mod2 .=== missing, "unknown non-m⁶A", df.mod2)
	# local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	# local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	# df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	# df[!, :type] = sign.(df.phi)
	# local counts = combine(groupby(df, :pair), nrow => :count)
	# counts = counts[counts.count .> 10, :]
	# df = innerjoin(df, counts, on = [:pair])

	# local df2 = copy(df)

	# local perm_results = []
	# for i in 1:10_000
	# 	# df2[!, :mod1] = shuffle(df2.mod1)
	# 	# df2[!, :mod2] = shuffle(df2.mod2)
	# 	df2[!, :pair] = shuffle(df2.pair)
	# 	x = combine(groupby(df2, [:pair, :type]), nrow => :count)
	# 	x = innerjoin(x[x.type .== 1, [:pair, :count]],
	# 				  x[x.type .== -1, [:pair, :count]],
	# 				  on = :pair,
	# 				  makeunique = true)
	# 	x[!, :lfc] = log2.(x.count ./ x.count_1)
	# 	push!(perm_results, x)
	# end
	# perm_results = vcat(perm_results...)
	# # perm_results = combine(groupby(perm_results, :pair),
	# # 					   :lfc => mean => :perm_mean_lfc,
	# # 					   :lfc => std => :perm_std_lfc,
	# # 					   :lfc => (lfc -> 1.96 * std(lfc)/sqrt(length(lfc))) => :perm_conf)
	# local baseline = mean(perm_results.lfc)
	
	# df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	# # df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	# df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	# df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	# df[!, :lfc] = log2.(df.count ./ df.count_1)
	# df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	# df[!, :conf] = 1.96 .* df.SE
	# df[!, :total] = df.count .+ df.count_1
	# df = sort(df[df.count .> 2 .&& df.count_1 .> 2, :], :total)

	# # df = sort(innerjoin(df, perm_results, on = [:pair]), :total)
	
	# local ax = Axis(trow[1, 2],
	# 				xlabel = "logFC (co-occurring/mutually-exclusive)",
	# 				ylabel = "Modification pair",
	# 			    yticks = (1:nrow(df), map((p, t) -> "$p\n($t)", df.pair, df.total)))
	# lines!(ax, [baseline, baseline], [0.5, nrow(df)+0.5], color = :orange, linestyle = :dash)
	# for (i, row) in enumerate(eachrow(df))
	# 	scatter!(ax, [row.lfc], [i], color = :black)
	# 	lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)

	# 	# scatter!(ax, [row.perm_mean_lfc], [i], color = :orange)
	# 	# lines!(ax, [row.perm_mean_lfc - row.perm_conf, row.perm_mean_lfc + row.perm_conf], [i, i], color = :orange)
	# end
	# xlims!(ax, [-1, 5])

	# PLOT B VARIANT 3

	local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	local df = renamemods(an_sig_comods)
	
	df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# df[!, :significant] = df.qvalue .< 0.01
	df[!, :mod1] = ifelse.(df.mod1 .=== missing, "Unclassified", df.mod1)
	df[!, :mod2] = ifelse.(df.mod2 .=== missing, "Unclassified", df.mod2)
	local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	df[!, :type] = sign.(df.phi)
	local counts = combine(groupby(df, :pair), nrow => :count)
	counts = counts[counts.count .> 10, :]
	df = innerjoin(df, counts, on = [:pair])
	df[!, :pair_label] = map((p, c) -> "$p\n($c)", df.pair, df.count)

	sort!(df, :count)

	local groups = combine(groupby(df, :pair),
						   :phi => (ps -> skipmissing(ps)) => :phis)
	groups = map(collect, groups.phis)
	
	println(HypothesisTests.KruskalWallisTest(groups...))
	
	local plt = data(df) *
		mapping(:pair_label => presorted => "Modification pair",
				:phi => "← Mutually exclusive    𝛗    Co-occurring →        ") *
		visual(Violin;
			   orientation = :horizontal,
			   side = :right,
			   color = "#FF33B5")
	draw!(trow[1, 2], plt; axis = (; limits = ((-1, 1), (0.7, length(groups) + 0.7))))


	local panel_c = trow[1, 3] = GridLayout()
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2

	local plt = data(df) *
		mapping((:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))) => "Distance",
				:phi => "← Mutually exclusive    𝛗    Co-occurring →         ",
			    color = :mean_entropy) *
		visual(Scatter; markersize = 5)
	draw!(panel_c[1, 1],
		  plt;
		  axis = (;
			# title = "Modification pairs' association strength by distance between them",
			# xticks = (1:0.1:4.2,
			# 		  ["10nt", repeat([""], 9)..., "100nt", repeat([""], 9)..., "1,000nt", repeat([""], 9)..., "10,000nt", "", ""]),
			# xticks = (map(log10, [(10:10:100)..., (200:100:1000)..., (2000:1000:10000)..., 11000, 12000]), ["10nt", repeat([""], 8)..., "100nt", repeat([""], 8)..., "1,000nt", repeat([""], 8)..., "10,000nt", "", ""]),
			# xticksize = [3.0, repeat([1], 8)..., 3, repeat([1], 8)..., 3, repeat([1], 8)..., 3, repeat([1], 2)...],
			# xminorticksvisible = true,
		 #    xminorgridvisible = true,
			# xminorticksize = 1.5,
			# xminorticks = IntervalsBetween(5),
			xticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
			xminorticksvisible = true,
		    xminorgridvisible = true,
			xminorticksize = 2,
			xminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]),
			limits = ((0.9, 4.2), (-1, 1))))
	Colorbar(panel_c[1, 2], limits = [0, 1], colormap = :viridis,
    vertical = true, label = "Information content")
	colsize!(panel_c, 1, Relative(5/6))

	

	local ax = Axis(brow[1, 1:3],
				    xtickformat = "{:,d}")# (vals -> map(format_with_commas, vals)))
				    # height = 110)
	# plot_isoforms_model!(ax, "PRRC2B"; transcripts = Set(["ENST00000682501.1", "ENST00000684596.1", "ENST00000683519.1"]), colors = Dict(["ENST00000682501.1" => :blue, "ENST00000684596.1" => :orange]), focus = (131496129-150, 131496179+200), rename = Dict(["ENST00000682501.1" => "PRRC2B-201", "ENST00000684596.1" => "PRRC2B-208", "ENST00000683519.1" => "PRRC2B-206"]))
	
	# # local ax = Axis(gb[1, 1])
	# local tx1_mask = occursin.("ENST00000682501.1", an_sig_comods.reference)
	# local tx2_mask = occursin.("ENST00000684596.1", an_sig_comods.reference)
	# local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
	# 			  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	# arc_plot_genomic2(brow[2, 1],
	# 				  renamemods(an_sig_comods[tx1_mask, :]),
	# 				  an_sig_peaks_ivt[occursin.("ENST00000682501.1", an_sig_peaks_ivt.ref_id), :];
	# 				  range=(5000, 5300),
	# 				  grange=(131496129-150, 131496179+100),
	# 				  colorrange=colorrange,
	# 				  title="PRRC2B-201",
	# 				  highlight=(5089, 5139),
	# 				  spinecolor=:blue)
	# arc_plot_genomic2(brow[3, 1],
	# 				  renamemods(an_sig_comods[tx2_mask, :]),
	# 				  an_sig_peaks_ivt[occursin.("ENST00000684596.1", an_sig_peaks_ivt.ref_id), :];
	# 				  range=(0, 12000),
	# 				  grange=(131496129-150, 131496179+100),
	# 				  colorrange=colorrange,
	# 				  title="PRRC2B-208",
	# 				  highlight=(7211, 7261),
	# 				  spinecolor=:orange)
	# Colorbar(brow[2:3, 2], limits = colorrange, colormap = :viridis,
 #    vertical = true, label = "-log₁₀(𝑃-value)")
	local gpos1 = 76066564
	local gpos2 = 76066573
	local row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
	local row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
	local tx1 = split(row1.reference, "|")[1]
	local tx2 = split(row2.reference, "|")[1]
	local name1 = split(row1.reference, "|")[5]
	local name2 = split(row2.reference, "|")[5]
	
	# local ax = Axis(gb[1, 1:2],
	# 			    height = 110)
	plot_isoforms_model!(ax, "MDH"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => :blue, tx2 => :orange]), focus = (gpos1-100, gpos2+100), rename = Dict([tx1 => name1, tx2 => name2]))
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.(tx1, an_sig_comods.reference)
	local tx2_mask = occursin.(tx2, an_sig_comods.reference)
	local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	local ax1, max_radius1 = arc_plot_genomic2(brow[2, 1],
					  renamemods(an_sig_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
					  range=(row1.pos1-2000, row1.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name1,
					  highlight=(row1.pos1, row1.pos2),
					  spinecolor=:blue)
	text!(ax1, 2, 8; text = "Co-occurring", align = (:left, :top))
	text!(ax1, 2, -8; text = "Mutually exclusive", align = (:left, :bottom))
	local ax2, max_radius2 = arc_plot_genomic2(brow[3, 1],
					  renamemods(an_sig_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
					  range=(row2.pos1-2000, row2.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name2,
					  highlight=(row2.pos1, row2.pos2),
					  spinecolor=:orange)
	text!(ax2, 2, 8; text = "Co-occurring", align = (:left, :top))
	text!(ax2, 2, -8; text = "Mutually exclusive", align = (:left, :bottom))
	local max_radius = max(max_radius1, max_radius2) * 2
	ylims!(ax1, -max_radius, max_radius)
	ylims!(ax2, -max_radius, max_radius)
	Colorbar(brow[2:3, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	local ref = row1.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [row1.pos1+1-4, row1.pos2+1-4]] .> 0.5))
	dropmissing!(df)
	local observed1 = contingency(df.Column1, df.Column2) .+ 1
	local res1 = HypothesisTests.ChisqTest(observed1)

	local ref = row2.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [row2.pos1+1-4, row2.pos2+1-4]] .> 0.5))
	dropmissing!(df)
	local observed2 = contingency(df.Column1, df.Column2) .+ 1
	local res2 = HypothesisTests.ChisqTest(observed2)


	local grid = brow[2:3, 3] = GridLayout()

	local matrices = [(Int.(round.(res1.expected; digits = 0)), observed1),
					  (Int.(round.(res2.expected; digits = 0)), observed2)]
	local titles = [name1, name2]
	for (i, (expected, observed)) in enumerate(matrices)
		expected = expected[2:-1:1, :]
		observed = observed[2:-1:1, :]
		
	    # Create an axis for this heatmap
	    ax = Axis(grid[i, 1],
				  title = titles[i],
				  aspect = 1,
				  xticks = (1:2, ["Unmodified", "Modified"]),
				  yticks = (1:2, ["Unmodified", "Modified"]),
				  yticklabelrotation = π/2)

		println(transpose(log2.(observed ./ expected)))
	    # Plot the heatmap
	    heatmap!(ax, transpose(log2.(observed ./ expected)),
	             colormap = :diverging_bwr_40_95_c42_n256,
				 colorrange = (-2.5, 2.5))
	             # colorrange = (0, max_count))
	
	    # Add row and column labels
		# if i < 1
		ax.xlabel = "$(format_with_commas(gpos1))\nm⁵C"
		ax.ylabel = "?\n$(format_with_commas(gpos2))"
		# else
		# 	ax.xlabel = "Pos. 7211 (m⁶A)"
		# 	ax.ylabel = "Pos. 7261 (m⁵C)"
		# end
	
	    # Annotate cells with values (optional)
	    for j in 1:2, k in 1:2
	        text!(ax, string(observed[j, k]) * "/" * string(expected[j, k]),
	              position = (k, j),
				  fontsize = 13.5,
	              align = (:center, :center),
	              color = :black)
	    end
	end

	# Legend(brow[4, 1],
	# 	   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π]),
	# 		LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π])],
	#        ["Mutually exclusive", "Co-occurring"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)


	colsize!(trow, 1, Relative(5/18))
	colsize!(trow, 2, Relative(5/18))
	colsize!(trow, 3, Relative(8/18))

	colsize!(brow, 1, Relative(11/13))
	colsize!(brow, 2, Relative(0.5/13))
	colsize!(brow, 3, Relative(1.5/13))
	# rowsize!(f.layout, 1, Relative(1.8/7))
	
	rowsize!(f.layout, 1, Relative(3/8))
	rowsize!(f.layout, 2, Relative(5/8))

	rowsize!(brow, 1, Relative(1/7))
	rowsize!(brow, 2, Relative(3/7))
	rowsize!(brow, 3, Relative(3/7))
	# rowsize!(brow, 4, Relative(1/10))
	# rowsize!(brow, 5, Relative(2/10))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	
	# Legend(gc[3, 1],
	# 	   [PolyElement(color=:black),
	# 		PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	#        ["Modified", "Not modified"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)


	for (label, layout) in zip(["A", "B", "C", "D"],
							   [trow[1, 1], trow[1, 2], trow[1, 3], brow[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end
  ╠═╡ =#

# ╔═╡ a5392b4a-bf22-4936-8426-476e85ccd92e
import LinearAlgebra

# ╔═╡ 755f8682-9890-4edb-bd0e-7c02f8acb3b2
(rand()-0.5)/10

# ╔═╡ 1613721e-9aa9-4898-aeff-68e90cd2a6bb
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

# ╔═╡ 6b71f0e3-288b-490a-9ab8-aeef32021708
begin

	# Example usage
	local points = rand(Point2f, 10)
	local labels = ["Label $i" for i in 1:10]
	
	local fig = Figure()
	local ax = Axis(fig[1, 1])
	
	scatter!(ax, points, color=:blue, markersize=15)
	
	local label_positions = place_labels_nonoverlapping(points, labels; offset=0.06, padding=0.05)
	
	for (p, lp, label) in zip(points, label_positions, labels)
		xalign = (lp[1] > p[1] ? :left : :right)
	    text!(ax, label, position=lp, align=(xalign, :center), color=:black)
	    lines!(ax, [p, lp], color=:gray, linewidth=0.5)
	end
	
	fig

end

# ╔═╡ ebeafdca-d72f-4f7e-88ae-93074150e041
combine(groupby(an_sig_comods, :reference), (df -> df[1:min(5, size(df, 1)), :]))

# ╔═╡ 5f9e87a6-e6b2-4a81-8f0a-0a5ba174b45c
md"""
# Final paper plot
"""

# ╔═╡ 854eef7a-e911-4f5b-8562-5bbd22e0c900
LABELSIZE = 18

# ╔═╡ d57a1fed-bd4a-485f-8aef-8cacfb7a294a
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

# ╔═╡ eddcd4e1-1ae0-40bb-8d0a-751fdd416caf
begin
	local f = Figure(size = (1200, 1450))
	local ga = f[1, 1] = GridLayout()
	local gb = f[2, 1] = GridLayout()
	local gc = f[3, 1] = GridLayout()


	local df = renamemods(an_sig_comods)

	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	df = df[in.(df.pair, Ref(Set(["m⁵C-m⁶A", "m⁵C-Ψ", "m⁶A-Ψ", "m⁶A-m⁷G"]))), :]
	# df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	local plt = data(df) *
		mapping(:phi => "φ", :pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", col = :pair) *
		visual(Scatter; markersize = 5)
	draw!(ga[1, 1:2], plt)

	local ax = Axis(gb[1, 1:2],
				    height = 110)
	plot_isoforms_model!(ax, "RPL34"; transcripts = Set(["ENST00000394665.5", "ENST00000502534.5"]), colors = Dict(["ENST00000394665.5" => :blue, "ENST00000502534.5" => :orange]), focus = (108625238-150, 108625248+200), rename = Dict(["ENST00000394665.5" => "RPL34-201", "ENST00000502534.5" => "RPL34-204"]))
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.("ENST00000394665.5", an_sig_comods.reference)
	local tx2_mask = occursin.("ENST00000502534.5", an_sig_comods.reference)
	local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	arc_plot_genomic2(gb[2, 1],
					  renamemods(an_sig_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.("ENST00000394665.5", an_sig_peaks_ivt.ref_id), :];
					  range=(512-50, 522+50),
					  grange=(108625238-100, 108625248+100),
					  colorrange=colorrange,
					  title="RPL34-201",
					  highlight=(0, 1000),
					  spinecolor=:blue)
	arc_plot_genomic2(gb[3, 1],
					  renamemods(an_sig_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.("ENST00000502534.5", an_sig_peaks_ivt.ref_id), :];
					  range=(421-100, 431+100),
					  grange=(108625238-100, 108625248+100),
					  colorrange=colorrange,
					  title="RPL34-204",
					  highlight=(0, 1000),
					  spinecolor=:orange)
	Colorbar(gb[2:3, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	# local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000304494.10|ENSG00000147889.18|OTTHUMG00000019686.7|OTTHUMT00000051915.1|CDKN2A-201|CDKN2A|978|protein_coding|")
	# # local f = Figure(size = (600, 600))
	

	# local df = DataFrame(Tables.table(probs[:, [1385, 1584]] .> 0.75))
	# dropmissing!(df)
	# local n2 = nrow(df)
	# println(contingency(df.Column1, df.Column2))
	# local ax1 = Axis(gc[1, 1],
	# 				 title="YWHAE-201",
	# 			     # aspect = 0.3,
	# 			     yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 				 xticks=(0:(n2/10):n2, ["$x%" for x in 0:10:100]),
	# 			     xlabel="")
	# println("n2: $n2")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax1, Array(sort(df)), colormap=["#dbdbdb", :black])

	# local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000571732.5|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000437356.1|YWHAE-208|YWHAE|1818|protein_coding|")
	# # local f = Figure(size = (600, 600))

	# local df = DataFrame(Tables.table(probs[:, [1416, 1615]] .> 0.75))
	# dropmissing!(df)
	# local n1 = nrow(df)
	# local ax2 = Axis(gc[2, 1],
	# 			 title="YWHAE-208",
	# 			 # aspect = 0.3,
	# 			 yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 			 xticks=(0:(n1/10):n1, ["$x%" for x in 0:10:100]),
	# 			 xlabel="Reads")
	# println("n1: $n1")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax2, Array(sort(df)), colormap=["#dbdbdb", :black])


	
	
	# xlims!(ax1, (0, max(n1, n2)))
	# xlims!(ax2, (0, max(n1, n2)))
	

	# colsize!(f.layout, 2, Relative(7/10))
	# colsize!(f.layout, 3, Relative(1/10))
	rowsize!(f.layout, 1, Relative(1.8/7))
	rowsize!(f.layout, 2, Relative(3.9/7))
	rowsize!(f.layout, 3, Relative(1.3/7))

	rowsize!(gb, 2, Relative(2/5))
	rowsize!(gb, 3, Relative(2/5))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	Legend(gb[4, 1],
		   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π])],
	       ["Mutually exclusive", "Co-modified"],
		   orientation = :horizontal,
		   framevisible = false)
	Legend(gc[3, 1],
		   [PolyElement(color=:black),
			PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	       ["Modified", "Not modified"],
		   orientation = :horizontal,
		   framevisible = false)


	for (label, layout) in zip(["A", "B", "C"],
							   [ga[1, 1:2], gb[1, 1], gc[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end

# ╔═╡ db26dd6d-23c0-495c-973d-b4829d05ff57
begin
	local f = Figure(size = (1000, 100))
	local ax = Axis(f[1, 1])
	plot_isoforms_model!(ax, "PRRC2B"; transcripts = Set(["ENST00000682501.1", "ENST00000684596.1", "ENST00000683519.1"]), colors = Dict(["ENST00000682501.1" => :blue, "ENST00000684596.1" => :orange]), focus = (131496129-150, 131496179+200), rename = Dict(["ENST00000682501.1" => "PRRC2B-201", "ENST00000684596.1" => "PRRC2B-208", "ENST00000683519.1" => "PRRC2B-206"]))
	f
end

# ╔═╡ dbf67725-f43f-49a7-8ef5-57b60acdbc67
begin
	local f = Figure(size = (1200, 1450))
	local ga = f[1, 1] = GridLayout()
	local gb = f[2, 1] = GridLayout()
	local gc = f[3, 1] = GridLayout()


	local df = renamemods(an_sig_comods)

	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	df = df[in.(df.pair, Ref(Set(["m⁵C-m⁶A", "m⁵C-Ψ", "m⁶A-Ψ", "m⁶A-m⁷G"]))), :]
	# df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	local plt = data(df) *
		mapping(:phi => "φ", :pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", col = :pair) *
		visual(Scatter; markersize = 5)
	draw!(ga[1, 1:2], plt)

	local gene = "GNAS"
	local tx1 = "ENST00000371095.7"
	local tx2 = "ENST00000371085.8"

	local ax = Axis(gb[1, 1:2],
				    height = 110)
	plot_isoforms_model!(ax, "GNAS"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => :blue, tx2 => :orange]), focus = (58891379, 58911187), rename = Dict([tx1 => "ENST00000371095.7", tx2 => "ENST00000371085.8"]))
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.(tx1, an_sig_comods.reference)
	local tx2_mask = occursin.(tx2, an_sig_comods.reference)
	local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	arc_plot_genomic2(gb[2, 1],
					  renamemods(an_sig_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
					  range=(300-30, 1800+30),
					  grange=(58910358-100, 58910385+100),
					  colorrange=colorrange,
					  title="CDKN2A-201",
					  # highlight=(536, 614),
					  spinecolor=:blue)
	arc_plot_genomic2(gb[3, 1],
					  renamemods(an_sig_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
					  range=(300-30, 1800+30),
					  grange=(58910358-100, 58910385+100),
					  colorrange=colorrange,
					  title="CDKN2A-214",
					  # highlight=(610, 688),
					  spinecolor=:orange)
	Colorbar(gb[2:3, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	# local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000304494.10|ENSG00000147889.18|OTTHUMG00000019686.7|OTTHUMT00000051915.1|CDKN2A-201|CDKN2A|978|protein_coding|")
	# # local f = Figure(size = (600, 600))
	

	# local df = DataFrame(Tables.table(probs[:, [1385, 1584]] .> 0.75))
	# dropmissing!(df)
	# local n2 = nrow(df)
	# println(contingency(df.Column1, df.Column2))
	# local ax1 = Axis(gc[1, 1],
	# 				 title="YWHAE-201",
	# 			     # aspect = 0.3,
	# 			     yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 				 xticks=(0:(n2/10):n2, ["$x%" for x in 0:10:100]),
	# 			     xlabel="")
	# println("n2: $n2")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax1, Array(sort(df)), colormap=["#dbdbdb", :black])

	# local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", "ENST00000571732.5|ENSG00000108953.17|OTTHUMG00000134316.7|OTTHUMT00000437356.1|YWHAE-208|YWHAE|1818|protein_coding|")
	# # local f = Figure(size = (600, 600))

	# local df = DataFrame(Tables.table(probs[:, [1416, 1615]] .> 0.75))
	# dropmissing!(df)
	# local n1 = nrow(df)
	# local ax2 = Axis(gc[2, 1],
	# 			 title="YWHAE-208",
	# 			 # aspect = 0.3,
	# 			 yticks=(1:2, ["m⁶A", "m⁵C"]),
	# 			 xticks=(0:(n1/10):n1, ["$x%" for x in 0:10:100]),
	# 			 xlabel="Reads")
	# println("n1: $n1")
	# println(reduce(*, map(sum, eachcol(df)) ./ nrow(df)) * nrow(df))
	# println(sum(df.Column1 .== 1 .&& df.Column2 .== 1))
	# heatmap!(ax2, Array(sort(df)), colormap=["#dbdbdb", :black])


	
	
	# xlims!(ax1, (0, max(n1, n2)))
	# xlims!(ax2, (0, max(n1, n2)))
	

	# colsize!(f.layout, 2, Relative(7/10))
	# colsize!(f.layout, 3, Relative(1/10))
	rowsize!(f.layout, 1, Relative(1.8/7))
	rowsize!(f.layout, 2, Relative(3.9/7))
	rowsize!(f.layout, 3, Relative(1.3/7))

	rowsize!(gb, 2, Relative(2/5))
	rowsize!(gb, 3, Relative(2/5))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	Legend(gb[4, 1],
		   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π]),
			LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π])],
	       ["Mutually exclusive", "Co-modified"],
		   orientation = :horizontal,
		   framevisible = false)
	Legend(gc[3, 1],
		   [PolyElement(color=:black),
			PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	       ["Modified", "Not modified"],
		   orientation = :horizontal,
		   framevisible = false)


	for (label, layout) in zip(["A", "B", "C"],
							   [ga[1, 1:2], gb[1, 1], gc[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end

# ╔═╡ 9bd1d228-e255-48ec-ba0f-07e5ac026581
LABELPAD = 0

# ╔═╡ 7dd9fe5b-c6d2-4e69-adf5-c7240f54e020
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1])

	local inset_ax = Axis(fig[1, 1],
						  width=Relative(0.5),
						  height=Relative(0.2),
						  halign=0.05,
						  valign=0.95,
						  backgroundcolor="#F0F0F0")
	hidedecorations!(inset_ax)
	local tx = "ENST00000234875.9"
	local name = "RPL22"

	plot_isoforms_model!(inset_ax, "ENST00000234875.9"; transcripts = Set([tx]), colors = Dict([tx => :black]), rename = Dict([tx => name]), focus=(0, 100))

	
	fig
end

# ╔═╡ bf1d3910-2680-40c4-95e9-20969fb5c1f8
sort(combine(groupby(sig_comods, [:reference]),
	    :phi => (ps -> maximum(ps) - minimum(ps)) => :phi_diff,
	    nrow => :count), :phi_diff, rev = true)

# ╔═╡ 1bcc41a6-aa15-4118-8a88-0a160ce1df5f


# ╔═╡ 0bf65cb8-d144-4384-b7c1-682c219a71ba
sig_comods

# ╔═╡ 50cc58bc-2b7a-4501-8405-5f178beddc44


# ╔═╡ 4db1d361-eb5c-45cd-9b4f-aaeb9b5528b6
filter(r -> r.x1, combine(groupby(sig_comods, [:reference]),
	    (df -> any(df.pvalue .< 1e-50 .&& df.phi .> 0) && any(df.pvalue .< 1e-20 .&& df.phi .< 0))))

# ╔═╡ dc7a2b77-4b3c-46d6-bc09-1a98dabf8c88
filter(r -> r.x1, combine(groupby(sig_comods, [:reference]),
	    (df -> any(df.pvalue .< 1e-30 .&& df.phi .> 0) && any(df.pvalue .< 1e-10 .&& df.phi .< 0))))

# ╔═╡ 795cbda2-5cb9-4e63-bc09-43b55288dd5c
sort(combine(groupby(sig_comods, [:reference]),
		[:pvalue, :phi] => ((pvalue, phi) -> sum(pvalue .< 1e-10 .&& phi .> 0.1)) => :npos,
	    [:pvalue, :phi] => ((pvalue, phi) -> sum(pvalue .< 1e-10 .&& phi .< -0.1)) => :nneg), [:nneg, :npos], rev = true)

# ╔═╡ 347a22fe-b7a4-422b-a229-d34369e64422
Associations.association(Associations.GaussianMI, [0, 0, 1, 1, 0, 0, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1])

# ╔═╡ d762ff2f-e90d-420b-9ae8-f5186316b5fa
function phi_coefficient(x::BitVector, y::BitVector)
    # Check vectors are of equal length
    length(x) == length(y) || throw(ArgumentError("Vectors must be of equal length"))

    # Calculate counts for 2x2 contingency table
    a = sum(x .& y)        # both true
    b = sum(x .& .!y)      # x true, y false
    c = sum(.!x .& y)      # x false, y true
    d = sum(.!x .& .!y)    # both false

    # Calculate phi
    phi_coefficient(a, b, c, d)
end

# ╔═╡ 8b238fa0-94b1-46a1-b6c0-fb3be10ecd6f
function phi_coefficient(a, b, c, d)
	numerator = a * d - b * c
    denominator = sqrt((a + b) * (c + d) * (a + c) * (b + d))
	if denominator == 0
        return NaN
    end

    return numerator / denominator
end

# ╔═╡ 6623763d-a9f8-492f-9f72-9054ec3d0ac4
begin
	local res = []
	for i in 1:1000
		local x = convert(Vector{Bool}, randn(100) .> 1.5)
		local y = convert(Vector{Bool}, randn(100) .> 1.5)		
		phi = phi_coefficient(x, y)
		println(sum(x), " ", sum(y), " ", phi)
		if !isnan(phi)
			push!(res, phi)
		end
	end
	mean(res)
end

# ╔═╡ adb29d24-b447-4301-9a8b-bd197b4401da
begin
	sig_comods_mi = copy(sig_comods)
	local disc = Associations.CodifyVariables(Associations.UniqueElements())
	local est = Associations.JointProbabilities(Associations.MIShannon(), disc)
	local mi = []
	for row in eachrow(sig_comods_mi)
		local x = vcat(zeros(row.a),
					   zeros(row.b),
					   ones(row.c),
					   ones(row.d))
		local y = vcat(zeros(row.a),
					   ones(row.b),
					   zeros(row.c),
					   ones(row.d))
		push!(mi, Associations.association(est, x, y))
	end
	sig_comods_mi[!, :mi] = mi
	sig_comods_mi
end

# ╔═╡ 03357cc7-cba4-459d-8233-e268c24f2d91
data(sig_comods_mi) *
	mapping(:mi,
			:pvalue => (p -> -log10(p)) => "-log10(p-value)",
		    color = :qvalue => (q -> q < 0.01) => "Significant") *
	visual(Scatter; markersize = 5) |> draw

# ╔═╡ de9108ea-7d68-42f5-a637-704866510e3a
data(sig_comods_mi) *
	mapping(:phi, :mi, color = :qvalue => (q -> -log10(q)) => "-log10(q-value)") *
	visual(Scatter; markersize = 5) |> draw

# ╔═╡ 3257ee0b-5fa8-4680-8c46-e81d515c4b05
sig_comods_mi[sig_comods_mi.phi .< 0.1 .&& sig_comods_mi.mi .> 0.9, :]

# ╔═╡ 1c110e43-7ef9-43fb-ad2d-48bcf710e461
data(sig_comods_mi) *
	mapping(:phi,
		    (:pos1, :pos2) => ((p1, p2) -> log10(abs(p2-p1))) => "log10(Distance)",
		    color = :mi => "Mutual information") *
	visual(Scatter; markersize = 5) |> draw

# ╔═╡ a8344213-28f9-4627-b4ba-1c336763bc9b
data(sig_comods_mi) *
	mapping(:phi,
		    (:pos1, :pos2) => ((p1, p2) -> log10(abs(p2-p1))) => "log10(Distance)",
		    color = (:a, :b, :c, :d) => ((a, b, c, d) -> 
				mean([(c+d)/(a+b+c+d), (b+d)/(a+b+c+d)])) => "Mean stoichiometry") *
	# mapping(:phi,
	# 	    (:pos1, :pos2) => ((p1, p2) -> log10(abs(p2-p1))) => "log10(Distance)",
	# 	    color = (:a, :b, :c, :d) => ((a, b, c, d) -> if a*d > b*c
	# 			1 - mean([(c+d)/(a+b+c+d), (b+d)/(a+b+c+d)])
	# 		else
	# 			mean([(c+d)/(a+b+c+d), (b+d)/(a+b+c+d)])
	# 		end) => "relevance") *
	visual(Scatter; markersize = 5) |> draw

# ╔═╡ c72aa1ad-ea8d-4f4f-afde-abdc65d00fbf
data(sig_comods_mi) *
	mapping((:a, :b, :c, :d) => ((a, b, c, d) -> (c+d)/(a+b+c+d)) => "stoich1",
			(:a, :b, :c, :d) => ((a, b, c, d) -> (b+d)/(a+b+c+d)) => "stoich2",
			color = :phi) *
	visual(Scatter; markersize = 5) |> draw(scales(Color = (; colormap = :diverging_bwr_20_95_c54_n256)))

# ╔═╡ 77ffbc65-469d-4ec6-ab6e-102d41b9afa4
data(sig_comods_mi) *
	mapping(:phi,
		    (:pos1, :pos2) => ((p1, p2) -> log10(abs(p2-p1))) => "log10(Distance)",
		    color = :pvalue => ((p) -> -log10(p)) => "-log10(p-value)") *
	visual(Scatter; markersize = 5) |> draw(scales(Color = (; colormap = :seaborn_crest_gradient)))

# ╔═╡ 406367df-98a1-4112-ac0d-7b579440d8be
data(shuffle(sig_comods_mi)) *
	mapping(:phi,
		    (:pos1, :pos2) => ((p1, p2) -> log10(abs(p2-p1))) => "log10(Distance)",
		    layout = (:pos1, :pos2, :genomicPos1, :genomicPos2) => ((pos1, pos2, genomicPos1, genomicPos2) -> abs(pos2 - pos1) == abs(genomicPos2 - genomicPos1)) => "Same exon") *
	visual(Scatter; markersize = 5) |> draw

# ╔═╡ 84e72546-516a-488f-9d6c-e7d006a5c427
[abs(r.pos2 - r.pos1) == abs(r.genomicPos2 - r.genomicPos1) for r in eachrow(sig_comods)] |> countmap

# ╔═╡ 0b09c0af-c716-48d9-8a63-7a8a9e534c60
begin
	# x, y = randn(1000), randn(1000)
	# disc = Associations.CodifyVariables(Associations.ValueBinning(3))
	# local x = [0, 0, 1, 1, 0, 0, 1, 1]
	# local y = [0, 0, 1, 1, 0, 0, 1, 1]
	# local x = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
	# local y = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
	local x = vcat(zeros(1290), zeros(968), ones(946), ones(1010))
	local y = vcat(zeros(1290), ones(968), zeros(946), ones(1010))
	# local def = Associations.ConditionalEntropyShannon(base = 2)
	# local disc = Associations.CodifyVariables(Associations.UniqueElements())
	# local est = Associations.JointProbabilities(Associations.MIShannon(), disc)
	local est = Associations.JointProbabilities(Associations.MIShannon(), Associations.UniqueElements())
	# local est = Associations.JointProbabilities(def, disc)
	local ent = (xs) -> let s = sum(xs) / length(xs)
		-(s * log2(s) + (1-s) * log2(1-s))
	end
	Associations.association(est, x, y), cor(x, y), mean([ent(x), ent(y)])

	
end

# ╔═╡ efee0b3a-16be-42b0-a9c5-1d1244c702fc
Associations.codify(Associations.UniqueElements(), [0, 0, 1, 1, 0, 0, 1, 1])

# ╔═╡ 9ee57741-1ed8-4b3d-a665-c2e782140da5
(100*100 - 0*0)/sqrt(100*100*100*100)

# ╔═╡ 7a6530d5-cd06-4e50-aba7-685e7517fa4e
phi_coefficient(Bool.([0, 1, 1, 1, 1, 1, 1, 1]), Bool.([0, 1, 1, 1, 1, 1, 1, 1]))

# ╔═╡ c81bea74-6d32-4a18-ba04-92cdc1afaa6a
let s = 0.05
	-(s * log2(s) + (1-s) * log2(1-s))
end

# ╔═╡ bd96030a-8920-41ac-b58e-f204deea0662
begin
	local m6a_m6a = filter(r -> r.mod1 === "m6A" && r.mod2 === "m6A", an_sig_comods)
	m6a_m6a[!, :category] .= "m⁶A-xm⁶A"
	m6a_m6a[!, :m6a_upstream] .= true
	local m6a_notm6a = filter(r -> (r.mod1 === "m6A" && r.pos2_not_m6A) || (r.mod2 === "m6A" && r.pos1_not_m6A), an_sig_comods)
	m6a_notm6a[!, :category] .= "m⁶A-not m⁶A"
	m6a_notm6a[!, :m6a_upstream] = m6a_notm6a.mod1 .=== "m6A"

	local df = vcat(m6a_m6a, m6a_notm6a)


	df[!, :tx_len] = map(x -> parse(Int, x[7]), split.(df.reference, "|"))
	df[!, :pos1_norm] = df.pos1 ./ df.tx_len
	df[!, :pos2_norm] = df.pos2 ./ df.tx_len

	df = stack(df, [:pos1_norm, :pos2_norm])

	data(df) *
		mapping(:value => "Transcript position", color = :variable, row = :category, col = :m6a_upstream) *
		AlgebraOfGraphics.density() |> draw

end

# ╔═╡ 1eab17ad-dbba-482c-8936-8968eae477ea
begin
	local m6a_m6a = filter(r -> r.mod1 === "m6A" && r.mod2 === "m6A", an_sig_comods)
	m6a_m6a[!, :category] .= "m⁶A-m⁶A"
	m6a_m6a[!, :m6a_upstream] .= true
	local m6a_notm6a = filter(r -> (r.mod1 === "m6A" && r.pos2_not_m6A) || (r.mod2 === "m6A" && r.pos1_not_m6A), an_sig_comods)
	m6a_notm6a[!, :category] .= "m⁶A-not m⁶A"
	m6a_notm6a[!, :m6a_upstream] = m6a_notm6a.mod1 .=== "m6A"

	local df = vcat(m6a_m6a, m6a_notm6a)

	df[!, :tx_len] = map(x -> parse(Int, x[7]), split.(df.reference, "|"))
	df[!, :distance] = abs.(df.pos1 .- df.pos2)

	data(df) *
		mapping(:tx_len, :distance, layout = :category) *
		visual(Scatter; markersize = 5) |> draw
end

# ╔═╡ af23935b-6ef5-4926-8ad8-b4e6d3b9fd4d
function to_metagene_pos(pos, transcript_id, gtf, start_metapos = 0.2, stop_metapos = 0.7)
	# annots = gtf[occursin.(transcript_id, gtf.attribute), :]
	# annots = gtf[gtf.transcript_id .== transcript_id, :]
	annots = gtf[(transcript_id = transcript_id,)]
	
	starts = annots[annots.feature .== "start_codon", :]
	if nrow(starts) == 0
		starts = annots[annots.feature .== "CDS", :]
	end
	if nrow(starts) == 0
		throw(error("No start codon in annotation for transcript $transcript_id"))
	end
	start = first(starts).start

	stops = annots[annots.feature .== "stop_codon", :]
	if nrow(stops) == 0
		stops = sort(annots[annots.feature .== "CDS", :], :start, rev = true)
	end
	if nrow(stops) == 0
		throw(error("No stop codon in annotation for transcript $transcript_id"))
	end
	stop = first(stops).start

	tx_start = minimum(annots.start)
	tx_end = maximum(annots.end)

	if pos < start
		start_metapos * (pos - tx_start) / (start - tx_start)
	elseif pos < stop
		start_metapos + (stop_metapos - start_metapos)*(pos - start)/(stop - start)
	else
		stop_metapos + (1 - stop_metapos)*(pos - stop)/(tx_end - stop)
	end
end

# ╔═╡ 3ba622e6-0d55-413f-92f6-4aa54144920c
begin
	local m6a_m6a = filter(r -> r.mod1 === "m6A" && r.mod2 === "m6A", an_sig_comods)
	m6a_m6a[!, :category] .= "m⁶A-xm⁶A"
	m6a_m6a[!, :m6a_upstream] .= true
	local m6a_notm6a = filter(r -> (r.mod1 === "m6A" && r.pos2_not_m6A) || (r.mod2 === "m6A" && r.pos1_not_m6A), an_sig_comods)
	m6a_notm6a[!, :category] .= "m⁶A-not m⁶A"
	m6a_notm6a[!, :m6a_upstream] = m6a_notm6a.mod1 .=== "m6A"

	local df = vcat(m6a_m6a, m6a_notm6a)

	df[!, :tx_len] = map(x -> parse(Int, x[7]), split.(df.reference, "|"))
	df[!, :pos1_norm] = df.pos1 ./ df.tx_len
	df[!, :pos2_norm] = df.pos2 ./ df.tx_len

	df = df[occursin.("protein_coding", df.reference), :]

	local annots = Dict(pairs(groupby(gtf, :transcript_id)))
	df[!, :transcript] = first.(split.(df.reference, "|"))

	local up_pos = map(min, df.genomicPos1, df.genomicPos2)
	local down_pos = map(max, df.genomicPos1, df.genomicPos2)
	df[!, :metapos1] = to_metagene_pos.(up_pos, df.transcript, Ref(annots))
	df[!, :metapos2] = to_metagene_pos.(down_pos, df.transcript, Ref(annots))

	df[!, :distance] = abs.(df.pos1 .- df.pos2)

	sort!(df, :metapos1)

	# println(occursin.("protein_coding", df.reference) |> sum, nrow(df))

	local fig = Figure(size = (1000, 600))
	
	local ax1 = Axis(fig[1, 1],
					 xticks = ([0.2, 0.7], ["", ""]),
					 ylabel = "Pairs",
					 title = "Co-occurring m⁶A-m⁶A")

	for (i, row) in enumerate(eachrow(df[df.category .== "m⁶A-xm⁶A" .&& df.relation .== :comod, :]))
		# transcript = split(row.reference, "|")[1]
		# metapos1 = to_metagene_pos(row.genomicPos1, transcript, annots)
		# metapos2 = to_metagene_pos(row.genomicPos2, transcript, annots)
		lines!(ax1, [row.metapos1, row.metapos2], [i, i], color = :black)
	end

	local ax2 = Axis(fig[1, 2],
					 xticks = ([0.2, 0.7], ["", ""]),
					 title = "Mutually exclusive m⁶A-m⁶A")

	for (i, row) in enumerate(eachrow(df[df.category .== "m⁶A-xm⁶A" .&& df.relation .== :excl, :]))
		# transcript = split(row.reference, "|")[1]
		# metapos1 = to_metagene_pos(row.genomicPos1, transcript, annots)
		# metapos2 = to_metagene_pos(row.genomicPos2, transcript, annots)
		lines!(ax2, [row.metapos1, row.metapos2], [i, i], color = :black)
	end

	local ax3 = Axis(fig[2, 1],
					 xticks = ([0.2, 0.7], ["Start codon", "Stop codon"]),
					 ylabel = "Pairs",
					 title = "Co-occurring m⁶A-not m⁶A")

	for (i, row) in enumerate(eachrow(df[df.category .== "m⁶A-not m⁶A" .&& df.relation .== :comod, :]))
		# transcript = split(row.reference, "|")[1]
		# metapos1 = to_metagene_pos(row.genomicPos1, transcript, annots)
		# metapos2 = to_metagene_pos(row.genomicPos2, transcript, annots)
		lines!(ax3, [row.metapos1, row.metapos2], [i, i], color = :black)
	end

	local ax4 = Axis(fig[2, 2],
					 xticks = ([0.2, 0.7], ["Start codon", "Stop codon"]),
					 title = "Mutually exclusive m⁶A-not m⁶A")

	for (i, row) in enumerate(eachrow(df[df.category .== "m⁶A-not m⁶A" .&& df.relation .== :excl, :]))
		# transcript = split(row.reference, "|")[1]
		# metapos1 = to_metagene_pos(row.genomicPos1, transcript, annots)
		# metapos2 = to_metagene_pos(row.genomicPos2, transcript, annots)
		lines!(ax4, [row.metapos1, row.metapos2], [i, i], color = :black)
	end

	ax5 = Axis(fig[1, 3], title = "m⁶A-m⁶A", xlabel = "Transcript length [log10]")
	density!(ax5, log10.(df[df.category .== "m⁶A-xm⁶A", :tx_len]), color = :black)
	ax6 = Axis(fig[2, 3], title = "m⁶A-not m⁶A", xlabel = "Transcript length [log10]")
	density!(ax6, log10.(df[df.category .== "m⁶A-not m⁶A", :tx_len]), color = :black)

	ax7 = Axis(fig[1, 4], title = "m⁶A-m⁶A", xlabel = "Distance [log10]", limits = ((0, 4.2), nothing))
	density!(ax7, log10.(df[df.category .== "m⁶A-xm⁶A", :distance]), color = :black)
	ax8 = Axis(fig[2, 4], title = "m⁶A-not m⁶A", xlabel = "Distance [log10]", limits = ((0, 4.2), nothing))
	density!(ax8, log10.(df[df.category .== "m⁶A-not m⁶A", :distance]), color = :black)
	
	fig
end

# ╔═╡ a94207bd-7666-41ab-a858-bd1170693a39
round(to_metagene_pos(130885683, "ENST00000318560.6", gtf); digits = 5)

# ╔═╡ 266de19f-4d79-45f7-8af2-071904b08b47
gtf[occursin.("ENST00000529421.5", gtf.transcript_id), :]

# ╔═╡ 9bdf330c-f5c1-4a42-9a30-caee4bd8bb3e
cancer_genes = DataFrame(CSV.File("/home/mzdravkov/cancerGeneList.tsv",
								  delim = '\t',
								  normalizenames = true))

# ╔═╡ 731c89c3-ca71-404d-8f80-cf3062eaec01
an_sig_comods[occursin.("GNAS", an_sig_comods.reference), :]

# ╔═╡ 03f335b8-97de-4b2b-b03d-000800ae6565
function plot_gene_comods(sig_comods, gtf, gene; axis = nothing, colormap = :seaborn_crest_gradient, colorrange = (0, 100), highclip = :black, min_w = 2, max_w = 8, xlimits = nothing)
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

	if !isnothing(xlimits)
		xlims!(ax, xlimits)
	end

	if isnothing(axis)
		fig
	else
		axis
	end
end

# ╔═╡ 7c0d7864-53ce-4fa1-bc8d-b87ce0eb2801
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

	pairs = renamemods2(pairs)

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
		
		# arc!(ax, Point2f(center, 0), radius, 0, orientation*π,
		# 	 color = pair.neglogp, #pair.phi,
		# 	 linewidth = pair.linewidth,
		# 	 alpha = alpha,
		# 	 colormap = colormap,
		# 	 colorrange = colorrange,#(minimum(pairs.phi), maximum(pairs.phi)),
		# 	 highclip = highclip,
		# 	 resolution = 10000)
		println(pair.phi)
	end

	local dots = data(mods) *
		mapping(:pos => (p -> p - shift), :y, color=:mod => "Modification", marker=:mod => "Modification") *
		visual(Scatter; markersize=dotsize)
	local grid = draw!(ax, dots)
	if legend && !isnothing(fig)
		legend!(fig[1, 1], grid; tellheight=false, tellwidth=false, halign=:center, valign=:bottom, orientation=:horizontal, framevisible=false, labelsize=fontsize, titlesize=fontsize)
	end
	println(mods)

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

# ╔═╡ d112cefd-e4f9-45ba-a96e-9f1fe7fbd353
plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, an_sig_comods), gtf, "ENST00000234875.9"; minheight = 1, maxheight = 10, dotsize=10)#, xlimits=(1170, 1850)) # RPL22

# ╔═╡ 1bc79106-8075-4b6f-aac8-17cc2f58ea88
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1])
	plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, an_sig_comods), gtf, "ENST00000234875.9"; minheight = 5, xlimits=(1170, 1850), legend=false, axis = ax) # RPL22
	# ylims!(ax, -100, 100)
	fig
end

# ╔═╡ e35dd878-06c3-4b3e-99f1-bd24705d8a6b
plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.1, an_sig_comods), gtf, "ENST00000370994.8"; minheight = 1000)#, xlimits=(1170, 1850)) # SERBP1

# ╔═╡ 6bfecb6c-b750-46d6-b5c1-1a0caf732809
plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.1, an_sig_comods), gtf, "ENST00000370994.8"; minheight = 1000, xlimits=(4000, 6000)) # SERBP1

# ╔═╡ 6c4e95d8-39fb-4e44-9778-778df71d9888
an_sig_comods

# ╔═╡ 4ac53947-763f-40bc-afed-6b3c05de042e
split(sig_comods.reference[1], "|")[7]

# ╔═╡ 9ee56964-d715-4faf-b740-6ce72e6b8e51
function plot_transcript_comods_elliptical_arcs(sig_comods, gene; axis = nothing, colormap = :seaborn_crest_gradient, colorrange = (0, 100), highclip = :black, min_w = 2, max_w = 8, xlimits = nothing, minheight = nothing)
	pairs = sig_comods[occursin.(gene, sig_comods.reference), :]
	# exons = gtf[gtf.feature .== "exon" .&& occursin.(gene, gtf.attribute), :]


	# shift = minimum(exons.start)
	# minx = minimum(exons.start)
	minx = 0
	maxx = parse(Int, split(sig_comods[1, :reference], "|")[7])

	if isnothing(axis)
		fig = Figure(size = (1200, 600))
		ax = Axis(fig[1, 1])
	else
		ax = axis
	end

	ticks = round.(range(minx, maxx; length = 4); digits = 0)
	ax.xticks = (ticks, [format_with_commas(Int(t)) for t in ticks])

	chr = pairs[1, :chr]
	# text!(ax, -(maxx-minx)/30, 0; text = chr, align = (:right, :center))
	ax.xlabel = chr
	ax.xlabelsize = LABELSIZE

	lines!(ax, [minx, maxx], [0, 0], color = :black, linewidth = 1)
	scatter!(ax, range(minx, maxx; length = 40), zeros(40);
			 marker = pairs[1, :strand] == "+" ? :rtriangle : :ltriangle,
			 color = :black)

	# for exon in eachrow(exons)
	# 	lines!(ax, [exon.start - shift, exon.end - shift], [0, 0], color = :black, linewidth = 15)
	# end
	lines!(ax, [minx, maxx], [0, 0], color=:black, linewidth=15)

	pairs[!, :neglogp] = -log10.(pairs.fisher_pvalue)
	# tx_data[‘neglogp’] = -np.log10(tx_data[‘p_value’])
    # min_w, max_w = 0.5, 6
    # w_norm = (pairs.neglogp .- minimum(pairs.neglogp)) ./
                 # (maximum(pairs.neglogp) - minimum(pairs.neglogp) + 1e-9)
	w_norm = pairs.phi
	pairs[!, :linewidth] = min_w .+ (max_w - min_w) .* w_norm

	for pair in eachrow(pairs)
		alpha = 0.25 + (1 - 0.25) * abs(pair.phi)
		left = min(pair.pos1, pair.pos2)
		right = max(pair.pos1, pair.pos2)
		distance = right - left
		radius = distance/2
		
		center = left + radius
		orientation = sign(pair.phi)

		h = isnothing(minheight) ? radius : max(minheight, radius)

		t = range(0, orientation*π, 100)
		x = radius * cos.(t)
		y = minimum(h) * sin.(t)

		lines!(ax, center.+x, y,
			   color = pair.neglogp, #pair.phi,
			   linewidth = pair.linewidth,
			   alpha = alpha,
			   colormap = colormap,
			   colorrange = colorrange,#(minimum(pairs.phi), maximum(pairs.phi)),
			   highclip = highclip)
		
		# arc!(ax, Point2f(center, 0), radius, 0, orientation*π,
		# 	 color = pair.neglogp, #pair.phi,
		# 	 linewidth = pair.linewidth,
		# 	 alpha = alpha,
		# 	 colormap = colormap,
		# 	 colorrange = colorrange,#(minimum(pairs.phi), maximum(pairs.phi)),
		# 	 highclip = highclip,
		# 	 resolution = 10000)
		println(pair.phi)
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

# ╔═╡ cfd23630-c1dd-4848-b6eb-9a93f7f59cde
plot_transcript_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, sig_comods), "ENST00000234875.9"; minheight=20) # RPL22

# ╔═╡ 585b9beb-7e21-45ab-9794-e95c55799e55
let
	# Define ellipse parameters
a, b = 2, 10 # semi-major and semi-minor axes
t = range(0, π, 100)  # parameter for a half-ellipse (arc)

# Parametric equations for an ellipse
x = a * cos.(t)
y = b * sin.(t)

# Plot the arc
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, x, y, color = :blue)
fig
end

# ╔═╡ da330153-49ee-41bd-82b3-43eb87c8ffa1
plot_gene_comods(sig_comods, gtf, "ENST00000371095.7") # GNAS

# ╔═╡ a1df4246-3313-44f5-9d81-bf66b27b2eae
an_sig_comods[occursin.("RPL22", an_sig_comods.reference), :]

# ╔═╡ c49f764e-4718-4678-869d-ac1dfdfb9492
sort(sig_comods[occursin.("ENST00000234875.9", sig_comods.reference), :], :phi)

# ╔═╡ 4d554305-bc31-4d9c-b459-5bde0665fbea


# ╔═╡ b4718644-4b13-4473-851d-2fe8bb0be509
sig_comods[sig_comods.pvalue .< 1e-30 .&& sig_comods.pos2 .- sig_comods.pos1 .> 300, :]

# ╔═╡ e948eff2-d2e1-4067-8fb0-36297a48770a
plot_gene_comods(sig_comods, gtf, "ENST00000253024.10") # TRIM28

# ╔═╡ c8e12743-31d4-47eb-a2ec-64ab133e6cbe
plot_gene_comods(sig_comods, gtf, "ENST00000254719.10") # RPA1

# ╔═╡ 5b73e0b5-321b-4433-90c6-63c788f0d578

plot_gene_comods(filter(r -> 6185600 < min(r.genomicPos1, r.genomicPos2) && max(r.genomicPos1, r.genomicPos2) < 6186640, sig_comods), gtf, "ENST00000234875.9", xlimits = (1000, 2000)) # RPL22

# ╔═╡ b7ee8e0a-7f2f-4f05-ac43-ceb5ba893715
plot_gene_comods(sig_comods, gtf, "ENST00000234875.9") # RPL22

# ╔═╡ 3abe560f-6acd-43b3-83d8-850a5ccf3345
plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, an_sig_comods), gtf, "ENST00000234875.9"; minheight = 1000, xlimits=(1170, 1850)) # RPL22

# ╔═╡ b51dbb01-ce9b-4076-aa78-2b5020503539
plot_gene_comods(sig_comods, gtf, "ENST00000479035.7") # RPL23

# ╔═╡ 7627b3c9-5b1a-45fe-a885-5e42266629a9
plot_gene_comods(sig_comods, gtf, "ENST00000479035.7"; xlimits = (2000, 2500)) # RPL23

# ╔═╡ d448709c-ec77-4ab4-8f8a-eeaddc5a390c
plot_gene_comods(sig_comods, gtf, "ENST00000359995.10") # SRSF2

# ╔═╡ 6a1eb59c-9249-4f3c-a66c-c44dbf803869
an_sig_comods[occursin.("MDH2", an_sig_comods.reference), :]

# ╔═╡ de14604b-2100-4c92-9297-733e832f4f4e
plot_gene_comods(sig_comods, gtf, "ENST00000315758.10") # MDH2

# ╔═╡ 099f01b5-f129-4159-af77-4e8f7aee07b7
an_sig_comods[occursin.("RPL34", an_sig_comods.reference), :]

# ╔═╡ 64c40f7a-82f1-4362-9456-4f93284faa88
plot_gene_comods(sig_comods, gtf, "ENST00000394668.2") # RPL34

# ╔═╡ 2a95d9b1-e129-4025-8a2e-26ab4bd86bfb
an_sig_comods[occursin.("HSP90AA1", an_sig_comods.reference), :]

# ╔═╡ c4575309-b6f0-4b8b-bb40-173c9b28fd4d
plot_gene_comods(sig_comods, gtf, "ENST00000334701.11") # HSP90AA1

# ╔═╡ eb668cda-a026-4677-84c5-ab8da1850669
sig_comods

# ╔═╡ 9ec3cf67-128b-43c5-b993-30bdab77ddee
begin
	local df = copy(sig_comods)
	df[!, :gene] = [x[6] for x in split.(df.reference, "|")]
	sort(combine(groupby(df, [:gene]), nrow => :count), :count, rev = true)
end

# ╔═╡ 1f56cac3-c7db-448e-9d5c-50e5e9f2b490
begin
	local df = copy(sig_comods)
	df[!, :gene] = [x[6] for x in split.(df.reference, "|")]
	df = sort(combine(groupby(df, [:gene]), nrow => :count), :count, rev = true)
	df[occursin.("RPL", df.gene), :]
end

# ╔═╡ 9d6c1ce1-17ce-4269-ae77-8f1e0dc67a77


# ╔═╡ 3563fda2-b3e1-4bba-831c-63a7bcfda780
begin
	local x = [1 2; 3 4]

	x[2:-1:1, :]
end

# ╔═╡ fc3e907d-3fc4-42ee-9bd6-79a79aa1ee4f
import HypothesisTests

# ╔═╡ 9836bc74-eff7-4d7a-853d-2398df762316
begin
	local m6a_m6a = filter(r -> r.mod1 === "m6A" && r.mod2 === "m6A", an_sig_comods)
	m6a_m6a[!, :category] .= "m6A-m6A"
	m6a_m6a[!, :m6a_upstream] .= true
	local m6a_notm6a = filter(r -> (r.mod1 === "m6A" && r.pos2_not_m6A) || (r.mod2 === "m6A" && r.pos1_not_m6A), an_sig_comods)
	m6a_notm6a[!, :category] .= "m6A-not m6A"
	m6a_notm6a[!, :m6a_upstream] = m6a_notm6a.mod1 .=== "m6A"

	local df = vcat(m6a_m6a, m6a_notm6a)

	df[!, :tx_len] = map(x -> parse(Int, x[7]), split.(df.reference, "|"))
	df[!, :pos1_norm] = df.pos1 ./ df.tx_len
	df[!, :pos2_norm] = df.pos2 ./ df.tx_len

	df = df[occursin.("protein_coding", df.reference), :]

	local annots = Dict(pairs(groupby(gtf, :transcript_id)))
	df[!, :transcript] = first.(split.(df.reference, "|"))

	local up_pos = map(min, df.genomicPos1, df.genomicPos2)
	local down_pos = map(max, df.genomicPos1, df.genomicPos2)
	df[!, :metapos1] = to_metagene_pos.(up_pos, df.transcript, Ref(annots))
	df[!, :metapos2] = to_metagene_pos.(down_pos, df.transcript, Ref(annots))

	df[!, :same] = map(r -> (r.metapos1 < 0.2 && r.metapos2 < 0.2) ||
						    (r.metapos1 >= 0.2 && r.metapos1 < 0.7 && r.metapos2 >= 0.2 && r.metapos2 < 0.7) ||
						    (r.metapos1 >= 0.7 && r.metapos2 >= 0.7), eachrow(df))

	local cont = [sum(df.category .== "m6A-m6A" .&& .! df.same) sum(df.category .== "m6A-m6A" .&& df.same);
	sum(df.category .== "m6A-not m6A" .&& .! df.same) sum(df.category .== "m6A-not m6A" .&& df.same)]
	HypothesisTests.FisherExactTest(cont[1, 1], cont[1, 2], cont[2, 1], cont[2, 2])
end

# ╔═╡ 69e63c65-1fe3-4106-b90b-d80978c9c609
begin
	local gpos1 = 76066564
	local gpos2 = 76066573
	local row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
	local row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
	local tx1 = split(row1.reference, "|")[1]
	local tx2 = split(row2.reference, "|")[1]
	local name1 = split(row1.reference, "|")[5]
	local name2 = split(row2.reference, "|")[5]
	
	local ref = row1.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	local df1 = DataFrame(Tables.table(probs[:, [row1.pos1+1-4, row1.pos2+1-4]] .> 0.5))
	dropmissing!(df1)
	local observed1 = contingency(df1.Column1, df1.Column2) .+ 1
	local res1 = HypothesisTests.ChisqTest(observed1)
	local data1 = map(r -> join(Int.(collect(r)), "-"), eachrow(df1))
	local counts1 = data1 |> countmap

	local ref = row2.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	local df2 = DataFrame(Tables.table(probs[:, [row2.pos1+1-4, row2.pos2+1-4]] .> 0.5))
	dropmissing!(df2)
	local observed2 = contingency(df2.Column1, df2.Column2) .+ 1
	local res2 = HypothesisTests.ChisqTest(observed2)
	local data2 = map(r -> join(Int.(collect(r)), "-"), eachrow(df2))
	local counts2 = data2 |> countmap

	
	local fig = Figure(size = (600, 600))
	local top = Axis(fig[1, 1],
					 xlabel = "",
					 xticks = (1:4, ["", "", "", ""]),
					 ylabel = "Pattern occurrences")
	xlims!(top, 0.5, 4.5)
	ylims!(top, 0, 2200)

	df1 = DataFrame(tx = name1, pattern = data1)
	df2 = DataFrame(tx = name2, pattern = data2)
	local df = vcat(df1, df2)
	df = combine(groupby(df, [:tx, :pattern]), nrow => :count)
	local barplt = data(df) *
		mapping(:pattern, :count, dodge = :tx, color = :tx) *
		visual(BarPlot; width = 0.75, dodge_gap = 0.05)
	draw!(top, barplt, scales(Color = (; palette = [:blue, :orange])))

	local expected = DataFrame(tx = [name1, name1, name1, name1,
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
	local barplt = data(expected) *
		mapping(:pattern, :count, dodge = :tx) *
		visual(BarPlot;
			   width = 0.75,
			   dodge_gap = 0.05,
			   color = (:white, 0),
			   strokewidth = 1.5,
			   strokecolor = :black)
	draw!(top, barplt)

	local order = ["0-0", "0-1", "1-0", "1-1"]
	local bottom = Axis(fig[2, 1],
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

	rowsize!(fig.layout, 1, Relative(5/6))
  	rowsize!(fig.layout, 2, Relative(1/6))

	Legend(
        fig[1, 1],
		[PolyElement(color=:blue),
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
	
	fig
end

# ╔═╡ c9d4449b-4d4b-401b-86c3-c52cbaba8060
begin
	local ref = "ENST00000682501.1|ENSG00000288701.1|OTTHUMG00000020827.6|OTTHUMT00000054751.1|PRRC2B-201|PRRC2B|9157|protein_coding|"
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [5089+1-4, 5139+1-4]] .>= 0.75))
	dropmissing!(df)
	local observed1 = contingency(df.Column1, df.Column2) .+ 1
	local res1 = HypothesisTests.ChisqTest(observed1)

	local ref = "ENST00000684596.1|ENSG00000288701.1|OTTHUMG00000020827.6|-|PRRC2B-208|PRRC2B|11253|protein_coding|"
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [7211+1-4, 7261+1-4]] .>= 0.75))
	dropmissing!(df)
	local observed2 = contingency(df.Column1, df.Column2) .+ 1
	local res2 = HypothesisTests.ChisqTest(observed2)


	local f = Figure(size = (1200, 300))
	local grid = f[1, 1:4] = GridLayout()

	local matrices = [round.(res1.expected; digits = 2), observed1, round.(res2.expected; digits = 2), observed2]
	local max_count = max(map(maximum, matrices)...)
	local titles = ["Expected PRRC2B-201", "Observed PRRC2B-201", "Expected PRRC2B-208", "Observed PRRC2B-208"]
	for (i, mat) in enumerate(matrices)
		mat = mat[2:-1:1, :]
	    # Create an axis for this heatmap
	    ax = Axis(grid[1, i],
				  title = titles[i],
				  xticks = (1:2, ["Unmodified", "Modified"]),
				  yticks = (1:2, ["Unmodified", "Modified"]),
				  yticklabelrotation = π/2)
	
	    # Plot the heatmap
	    heatmap!(ax, transpose(mat),
	             colormap = :YlOrRd_3		
,
	             colorrange = (0, max_count))
	
	    # Add row and column labels
		if i < 3
			ax.xlabel = "Pos. 5089 (m⁶A)"
			ax.ylabel = "Pos. 5139 (m⁵C)"
		else
			ax.xlabel = "Pos. 7211 (m⁶A)"
			ax.ylabel = "Pos. 7261 (m⁵C)"
		end
	
	    # Annotate cells with values (optional)
	    for j in 1:2, k in 1:2
	        text!(ax, string(mat[j, k]),
	              position = (k, j),
				  fontsize = 22,
	              align = (:center, :center),
	              color = :black)
	    end
	end
	
	f
end

# ╔═╡ 1a0a63a1-4cb6-4578-a32a-7c623670f0f0
begin
	local ref = "ENST00000682501.1|ENSG00000288701.1|OTTHUMG00000020827.6|OTTHUMT00000054751.1|PRRC2B-201|PRRC2B|9157|protein_coding|"
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [5089+1-4, 5139+1-4]] .>= 0.75))
	dropmissing!(df)
	local observed1 = contingency(df.Column1, df.Column2) .+ 1
	local res1 = HypothesisTests.ChisqTest(observed1)

	local ref = "ENST00000684596.1|ENSG00000288701.1|OTTHUMG00000020827.6|-|PRRC2B-208|PRRC2B|11253|protein_coding|"
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [7211+1-4, 7261+1-4]] .>= 0.75))
	dropmissing!(df)
	local observed2 = contingency(df.Column1, df.Column2) .+ 1
	local res2 = HypothesisTests.ChisqTest(observed2)


	local f = Figure(size = (1200, 300))
	local grid = f[1, 1:4] = GridLayout()

	local matrices = [(Int.(round.(res1.expected; digits = 0)), observed1),
					  (Int.(round.(res2.expected; digits = 0)), observed2)]
	local titles = ["PRRC2B-201", "PRRC2B-208"]
	for (i, (expected, observed)) in enumerate(matrices)
		expected = expected[2:-1:1, :]
		observed = observed[2:-1:1, :]
		
	    # Create an axis for this heatmap
	    ax = Axis(grid[1, i],
				  title = titles[i],
				  aspect = 1,
				  xticks = (1:2, ["Unmodified", "Modified"]),
				  yticks = (1:2, ["Unmodified", "Modified"]),
				  yticklabelrotation = π/2)

		println(transpose(log2.(observed ./ expected)))
	    # Plot the heatmap
	    heatmap!(ax, transpose(log2.(observed ./ expected)),
	             colormap = :diverging_bwr_40_95_c42_n256,
				 colorrange = (-2.5, 2.5))
	             # colorrange = (0, max_count))
	
	    # Add row and column labels
		if i < 3
			ax.xlabel = "Pos. 5089 (m⁶A)"
			ax.ylabel = "Pos. 5139 (m⁵C)"
		else
			ax.xlabel = "Pos. 7211 (m⁶A)"
			ax.ylabel = "Pos. 7261 (m⁵C)"
		end
	
	    # Annotate cells with values (optional)
	    for j in 1:2, k in 1:2
	        text!(ax, string(observed[j, k]) * "/" * string(expected[j, k]),
	              position = (k, j),
				  fontsize = 22,
	              align = (:center, :center),
	              color = :black)
	    end
	end
	
	f
end

# ╔═╡ 212eda1c-da46-4833-9f58-c98d437c450f
begin
	local df = copy(an_sig_comods)
	df[!, :type] = sign.(df.phi)
	df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	# df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	df[!, :lfc] = log2.(df.count ./ df.count_1)
	df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	df[!, :conf] = 1.96 .* df.SE
	df[!, :total] = df.count .+ df.count_1
	df = sort(df[df.count .> 10 .&& df.count_1 .> 10, :], :total)
	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "logFC (co-occurring/mutually-exclusive)",
					ylabel = "Modification pair",
				    yticks = (1:nrow(df), map((p, t) -> "$p ($t)", df.pair, df.total)))
	for (i, row) in enumerate(eachrow(df))
		scatter!(ax, [row.lfc], [i], color = :black)
		lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)
	end
	xlims!(ax, [0, 4])
	f
end

# ╔═╡ 435ffc34-5696-49ac-8205-6ab4a92ca3b2
begin
	local df = copy(an_sig_comods)
	
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	df[!, :type] = sign.(df.phi)

	local df2 = copy(df)

	local perm_results = []
	for i in 1:10_000
		# df2[!, :mod1] = shuffle(df2.mod1)
		# df2[!, :mod2] = shuffle(df2.mod2)
		df2[!, :pair] = shuffle(df2.pair)
		x = combine(groupby(df2, [:pair, :type]), nrow => :count)
		x = innerjoin(x[x.type .== 1, [:pair, :count]],
					  x[x.type .== -1, [:pair, :count]],
					  on = :pair,
					  makeunique = true)
		x[!, :lfc] = log2.(x.count ./ x.count_1)
		push!(perm_results, x)
	end
	perm_results = vcat(perm_results...)
	# perm_results = combine(groupby(perm_results, :pair),
	# 					   :lfc => mean => :perm_mean_lfc,
	# 					   :lfc => std => :perm_std_lfc,
	# 					   :lfc => (lfc -> 1.96 * std(lfc)/sqrt(length(lfc))) => :perm_conf)
	local baseline = mean(perm_results.lfc)
	
	df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	# df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	df[!, :lfc] = log2.(df.count ./ df.count_1)
	df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	df[!, :conf] = 1.96 .* df.SE
	df[!, :total] = df.count .+ df.count_1
	df = sort(df[df.count .> 5 .&& df.count_1 .> 5, :], :total)

	# df = sort(innerjoin(df, perm_results, on = [:pair]), :total)
	
	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "logFC (co-occurring/mutually-exclusive)",
					ylabel = "Modification pair",
				    yticks = (1:nrow(df), map((p, t) -> "$p ($t)", df.pair, df.total)))
	lines!(ax, [baseline, baseline], [0.5, nrow(df)+0.5], color = :orange, linestyle = :dash)
	for (i, row) in enumerate(eachrow(df))
		scatter!(ax, [row.lfc], [i], color = :black)
		lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)

		# scatter!(ax, [row.perm_mean_lfc], [i], color = :orange)
		# lines!(ax, [row.perm_mean_lfc - row.perm_conf, row.perm_mean_lfc + row.perm_conf], [i, i], color = :orange)
	end
	
	xlims!(ax, [0, 4])
	f
end

# ╔═╡ f5996482-be98-4751-bbd8-19679d03658c
innerjoin(sig_peaks_ivt,
		  unique(stack(comods[:, [:reference, :pos1, :pos2]], [:pos1, :pos2])[:, [:reference, :value]]),
		  on = [:ref_id => :reference, :pos => :value])

# ╔═╡ 6a5c262d-b493-43fe-9f78-7223c15734ed
begin
	local df = stack(sig_comods[:, [:reference, :pos1, :pos2]], [:pos1, :pos2])[:, [:reference, :value]]
	intersect(Set(zip(sig_peaks_ivt.ref_id, sig_peaks_ivt.pos)),
		  	  Set(zip(df.reference, df.value))) |> length
end

# ╔═╡ 91728250-037c-419b-a531-523b16ea9aea
begin
	
	local df = stack(sig_comods[:, [:reference, :pos1, :pos2]], [:pos1, :pos2])[:, [:reference, :value]]
	intersect(Set(zip(sig_peaks_ivt.ref_id, sig_peaks_ivt.pos)),
		  	  Set(zip(df.reference, df.value))) |> length
end

# ╔═╡ 4dabb682-d787-4b18-b78e-a3be8ee90431
sig_comods

# ╔═╡ 71d16997-a2e4-4cc5-94f9-aaddf5d90339
begin
	local ref = "ENST00000682501.1|ENSG00000288701.1|OTTHUMG00000020827.6|OTTHUMT00000054751.1|PRRC2B-201|PRRC2B|9157|protein_coding|"
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))

	local modpos = sig_peaks_ivt[sig_peaks_ivt.ref_id .== ref, :pos]
	# modpos = modpos[modpos .> 5000]
	

	local df = DataFrame(Tables.table(probs[:, modpos .+ 1 .- 4] .>= 0.75))
	df = rename!(df, map(string, modpos))
	df
	local m = Float64.(Matrix(dropmissing!(df))) |> transpose
	local proj = UMAP.umap(m; n_neighbors = 2, min_dist = 0.1) |> permutedims
	proj = DataFrame(proj, [:umap1, :umap2])
	# proj[!, :pair] = map(a -> join(a, "-") , df.pair)
	data(proj) *
		mapping(:umap1, :umap2) *
		visual(Scatter, markersize = 4) |> draw
end

# ╔═╡ dcce6567-b5ad-4e7a-ace7-81490ed559be
function entropy(row)
	-(row.mod_ratio * log2(row.mod_ratio + 1e-300) + (1-row.mod_ratio) * log2(1-row.mod_ratio  + 1e-300))
end

# ╔═╡ 0555b650-e9a2-48c6-939e-c146241e125c
begin
	local f = Figure(size = (1600, 1100))
	local trow = f[1, 1] = GridLayout()
	local brow = f[2, 1] = GridLayout()

	local sig_color = "#800080" # "#051094"


	local df = copy(comods)[1:100, :]
	df[!, :significant] = df.qvalue .< 0.01
	local plt = data(df) *
		mapping(:phi, :pvalue => (p -> -log10(p)), color = :significant) *
		visual(Scatter; markersize = 5)
	draw!(trow[1, 1],
		  plt,
		  scales(Color = (; palette = [:grey, sig_color]));
		  axis = (; #title = "All modification pairs",
				  	xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
					ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
					limits = ((-1, 1), (-2, 310))))


	# PLOT B
	local plt_b = trow[1, 2] = GridLayout()	

	# local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	# local df = renamemods(an_sig_comods)
	
	# df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# # df[!, :significant] = df.qvalue .< 0.01
	# df[!, :mod1] = ifelse.(df.mod1 .=== missing, "Unclassified", df.mod1)
	# df[!, :mod2] = ifelse.(df.mod2 .=== missing, "Unclassified", df.mod2)
	# local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	# local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	# df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	# df[!, :type] = sign.(df.phi)
	# local counts = combine(groupby(df, :pair), nrow => :count)
	# counts = counts[counts.count .> 40, :]
	# df = innerjoin(df, counts, on = [:pair])
	# df[!, :pair_label] = map((p, c) -> "$p ($c)", df.pair, df.count)

	# sort!(df, :count, rev = true)

	# local ordered_labels = collect(unique(df.pair_label))

	# local groups = combine(groupby(df, :pair),
	# 					   :phi => (ps -> skipmissing(ps)) => :phis)
	# groups = map(collect, groups.phis)
	
	# println(HypothesisTests.KruskalWallisTest(groups...))

	# df = df[:, [:pair_label, :phi]]

	local m6a_m6a = filter(r -> r.mod1 === "m6A" && r.mod2 === "m6A", an_sig_comods)
	local m6a_nonm6a = filter(r -> (r.mod1 === "m6A" && r.pos2_not_m6A) || (r.mod2 === "m6A" && r.pos1_not_m6A), an_sig_comods)

	local some_m6a = vcat(m6a_m6a, m6a_nonm6a)
	some_m6a[!, :pair] = map((a, b) -> join([a, b], "-"), some_m6a.mod1, some_m6a.mod2)
	local df2 = copy(some_m6a[:, [:pair, :phi]])

	local perm_results = []
	for i in 1:10_00#0
		df2[!, :pair] = shuffle(df2.pair)
		push!(perm_results, df2)
	end
	perm_results = vcat(perm_results...)
	# sort!(perm_results, :pair, by = p -> findfirst(ordered_labels .== p))

	# df[!, :category] .= :original
	# perm_results[!, :category] .= :permutations
	# local combined_data = vcat(df, perm_results)
	# combined_data = combined_data[.! occursin.("m⁶A-Unclassified", combined_data.pair_label), :]
	# println(combined_data.pair_label |> unique)

	# combined_data[!, :label] = map(r -> if r.category == :permutations
	# 	"Null"
	# elseif occursin.("m⁶A-m⁶A", r.pair_label)
	# 	"m⁶A-m⁶A"
	# else
	# 	"m⁶A-not m⁶A"
	# end, eachrow(combined_data))

	local ax = Axis(plt_b[1, 1],
				    xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ")
	# density!(ax,
	# 		 combined_data[combined_data.category .== :permutations, :phi],
	# 		 label = "Null",
	# 		 bandwidth = 0.02,
	# 		 boundary = (-1, 1),
	# 		 color = (:grey, 0.0),
	# 		 strokecolor = :grey,
	# 		 strokewidth = 2)
	density!(ax,
			 m6a_m6a.phi,
			 label = "m⁶A-m⁶A",
			 bandwidth = 0.01,
			 boundary = (-1, 1),
			 color = (:blue, 0.0),
			 strokecolor = :blue,
			 strokewidth = 2)
	density!(ax,
			 m6a_nonm6a.phi,
			 label = "m⁶A-not m⁶A",
			 bandwidth = 0.01,
			 boundary = (-1, 1),
			 color = (:orange, 0.0),
			 strokecolor = :orange,
			 strokewidth = 2)
	
	# local f = Figure()

	# local axes = []
	# for (i, pair) in enumerate(ordered_labels)
	# 	ax = Axis(plt_b[i, 1],
	# 			  xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
	# 			  yticks = 1:2:9,
	# 			  yticklabelsize = 9)
	# 	density!(ax,
	# 			 combined_data[combined_data.pair_label .== pair .&& combined_data.category .== :permutations, :phi],
	# 			 color = (:grey, 0.2),
	# 			 linestyle = :dash,
	# 			 strokewidth = 1.5,
	# 			 strokecolor = (:grey, 1))
	# 	density!(ax,
	# 			 combined_data[combined_data.pair_label .== pair .&& combined_data.category .== :original, :phi],
	# 			 color = ("#FF33B5", 0.2),
	# 			 linestyle = :solid,
	# 			 strokewidth = 1.5,
	# 			 strokecolor = ("#FF33B5", 1))
	# 	Legend(
 #        	plt_b[i, 1],
	# 		[MarkerElement(marker = 'x', markersize = 0)],
	# 		[pair],
	# 		position = (0.1, 0.1),
	#         tellheight = false,
	#         tellwidth = false,
	# 		halign = :left,
	# 		valign = :top,
	#         # padding = (2, 2, -5, 0),
	# 		patchsize = (0, 0),
	# 		framevisible = false,
	# 		backgroundcolor = (:white, 0),
	# 		labelfont = :bold,
	# 		labelsize = 10
 #    	)
	# 	xlims!(ax, (-1, 1))
	# 	ylims!(ax, (-1, 10))
	# 	if i < length(ordered_labels)
	# 		hidexdecorations!(ax, grid = false)
	# 	end
	# 	ax.ylabelrotation = 0
	# 	push!(axes, ax)	
	# end
	# linkyaxes!(axes...)
	# Label(plt_b[1:length(ordered_labels), 0], "pdf", rotation = pi/2)
	# rowgap!(plt_b, 0)


	# plt = data(combined_data) *
	# 	mapping(:phi => "← Mutually exclusive    𝛗    Co-occurring →        ",
	# 		    color = :label) *
	# 	visual(Density)
	# draw!(trow[1, 2], plt)

	

	# local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	# local df = renamemods(an_sig_comods)
	
	# df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# # df[!, :significant] = df.qvalue .< 0.01
	# df[!, :mod1] = ifelse.(df.mod1 .=== missing, "Unclassified", df.mod1)
	# df[!, :mod2] = ifelse.(df.mod2 .=== missing, "Unclassified", df.mod2)
	# local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	# local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	# df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	# df[!, :type] = sign.(df.phi)
	# local counts = combine(groupby(df, :pair), nrow => :count)
	# counts = counts[counts.count .> 10, :]
	# df = innerjoin(df, counts, on = [:pair])
	# df[!, :pair_label] = map((p, c) -> "$p\n($c)", df.pair, df.count)

	# sort!(df, :count)

	# local groups = combine(groupby(df, :pair),
	# 					   :phi => (ps -> skipmissing(ps)) => :phis)
	# groups = map(collect, groups.phis)
	
	# println(HypothesisTests.KruskalWallisTest(groups...))

	# local df2 = copy(df)
	# local perm_results = []
	# for i in 1:10_000
	# 	df2[!, :mod1] = shuffle(df2.mod1)
	# 	df2[!, :mod2] = shuffle(df2.mod2)
	# 	m6a_nonm6a = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# 	append!(perm_results, m6a_nonm6a.phi)
	# end

	
	# df = vcat(df[:, [:pair_label, :phi]], DataFrame(pair_label = "Null", phi = perm_results))
	
	# local plt = data(df) *
	# 	mapping(:pair_label => presorted => "Modification pair",
	# 			:phi => "← Mutually exclusive    𝛗    Co-occurring →        ") *
	# 	visual(Violin;
	# 		   orientation = :horizontal,
	# 		   side = :right,
	# 		   color = "#FF33B5")
	# draw!(trow[1, 2], plt; axis = (; limits = ((-1, 1), (0.7, length(groups) + 1 + 0.7))))


	# PLOT C

	local panel_c = trow[1, 3] = GridLayout()
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2

	df[!, :gene] = [x[6] for x in split.(df.reference, "|")]

	local df_cancer = innerjoin(df, cancer_genes,
			  		   		    on = [:gene => :Hugo_Symbol])
	df_cancer = df_cancer[abs.(df_cancer.phi) .>= 0.3 .&& df_cancer.mean_entropy .>= 0.75, :]
	println(size(df_cancer))

	local ax = Axis(panel_c[1, 1],
				    xlabel = "← Mutually exclusive    𝛗    Co-occurring →                    ",
				    ylabel = "Distance between associated modifications",
				    yticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
					yminorticksvisible = true,
				    yminorgridvisible = true,
					yminorticksize = 2,
					yminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]))
	local plt = data(df) *
		mapping(:phi,
				(:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))),
			    color = :mean_entropy) *
		visual(Scatter; markersize = 5)
	draw!(ax, #panel_c[1, 1],
		  plt)
		 #  axis = (;
			# yticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
			# yminorticksvisible = true,
		 #    yminorgridvisible = true,
			# yminorticksize = 2,
			# yminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]),
			# limits = ((-1, 1), (0.9, 4.2))))
	Colorbar(panel_c[1, 2], limits = [0, 1], colormap = :viridis,
    vertical = true, label = "Information content\n(mean stoichiometry entropy)")
	colsize!(panel_c, 1, Relative(5/6))

	local points = Point2f.(df_cancer.phi, log10.(abs.(df_cancer.pos1 .- df_cancer.pos2)))
	local labels = df_cancer.gene
	local label_positions = place_labels_nonoverlapping(points, labels)
	for (p, lp, label) in zip(points, label_positions, labels)
		xalign = (lp[1] > p[1] ? :left : :right)
	    text!(ax, label, position=lp, align=(xalign, :center), color=:black, fontsize=8)
	    lines!(ax, [p, lp], color=:gray, linewidth=0.5)
	end
	xlims!(ax, -1, 1.3)
	
	# local ax = Axis(panel_c[1, 1])
	# for row in eachrow(df_cancer)
	# 	text!(ax, row.phi, log10(abs(row.pos1 - row.pos2)); text = row.gene, fontsize = 8)
	# end
	
	# local plt2 = data(df_cancer) *
	# 	mapping(:phi => "← Mutually exclusive    𝛗    Co-occurring →        ",
	# 			(:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))) => "Distance between associated modifications",
	# 		    color = :mean_entropy) *
	# 	visual(Scatter; markersize = 10, marker = 'x')
	# draw!(panel_c[1, 1],
	# 	  plt2;
	# 	  axis = (;
	# 		yticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
	# 		yminorticksvisible = true,
	# 	    yminorgridvisible = true,
	# 		yminorticksize = 2,
	# 		yminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]),
	# 		limits = ((-1, 1), (0.9, 4.2))))

	

	local ax = Axis(brow[1, 1:3],
				    xtickformat = "{:,d}")# (vals -> map(format_with_commas, vals)))
			
	local gpos1 = 76066564
	local gpos2 = 76066573
	local row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
	local row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
	local tx1 = split(row1.reference, "|")[1]
	local tx2 = split(row2.reference, "|")[1]
	local name1 = split(row1.reference, "|")[5]
	local name2 = split(row2.reference, "|")[5]
	
	# local ax = Axis(gb[1, 1:2],
	# 			    height = 110)
	plot_isoforms_model!(ax, "MDH"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => :blue, tx2 => :orange]), focus = (gpos1-100, gpos2+100), rename = Dict([tx1 => name1, tx2 => name2]))
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.(tx1, an_sig_comods.reference)
	local tx2_mask = occursin.(tx2, an_sig_comods.reference)
	local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	local ax1, max_radius1 = arc_plot_genomic2(brow[2, 1],
					  renamemods(an_sig_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
					  range=(row1.pos1-2000, row1.pos2+2000),
					  grange=(gpos1-30, gpos2+30),
					  colorrange=colorrange,
					  title=name1,
					  highlight=(row1.pos1, row1.pos2),
					  spinecolor=:blue)
	text!(ax1, 2, 8; text = "Co-occurring", align = (:left, :top))
	text!(ax1, 2, -8; text = "Mutually exclusive", align = (:left, :bottom))
	local ax2, max_radius2 = arc_plot_genomic2(brow[3, 1],
					  renamemods(an_sig_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
					  range=(row2.pos1-2000, row2.pos2+2000),
					  grange=(gpos1-30, gpos2+30),
					  colorrange=colorrange,
					  title=name2,
					  highlight=(row2.pos1, row2.pos2),
					  spinecolor=:orange)
	text!(ax2, 2, 8; text = "Co-occurring", align = (:left, :top))
	text!(ax2, 2, -8; text = "Mutually exclusive", align = (:left, :bottom))
	local max_radius = max(max_radius1, max_radius2) * 2
	ylims!(ax1, -max_radius, max_radius)
	ylims!(ax2, -max_radius, max_radius)
	Colorbar(brow[2:3, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	local ref = row1.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	local df1 = DataFrame(Tables.table(probs[:, [row1.pos1+1-4, row1.pos2+1-4]] .> 0.5))
	dropmissing!(df1)
	local observed1 = contingency(df1.Column1, df1.Column2) .+ 1
	local res1 = HypothesisTests.ChisqTest(observed1)
	local data1 = map(r -> join(Int.(collect(r)), "-"), eachrow(df1))
	local counts1 = data1 |> countmap

	local ref = row2.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	local df2 = DataFrame(Tables.table(probs[:, [row2.pos1+1-4, row2.pos2+1-4]] .> 0.5))
	dropmissing!(df2)
	local observed2 = contingency(df2.Column1, df2.Column2) .+ 1
	local res2 = HypothesisTests.ChisqTest(observed2)
	local data2 = map(r -> join(Int.(collect(r)), "-"), eachrow(df2))
	local counts2 = data2 |> countmap

	
	# local fig = Figure(size = (600, 600))

	local grid = brow[2:3, 3] = GridLayout()	
	local top = Axis(grid[1, 1],
					 xlabel = "",
					 xticks = (1:4, ["", "", "", ""]),
					 ylabel = "Pattern occurrences")
	xlims!(top, 0.5, 4.5)
	ylims!(top, 0, 2350)

	df1 = DataFrame(tx = name1, pattern = data1)
	df2 = DataFrame(tx = name2, pattern = data2)
	local df = vcat(df1, df2)
	df = combine(groupby(df, [:tx, :pattern]), nrow => :count)
	local barplt = data(df) *
		mapping(:pattern, :count, dodge = :tx, color = :tx) *
		visual(BarPlot; width = 0.7, dodge_gap = 0.15)
	draw!(top, barplt, scales(Color = (; palette = [:blue, :orange])))

	local expected = DataFrame(tx = [name1, name1, name1, name1,
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
	local barplt = data(expected) *
		mapping(:pattern, :count, dodge = :tx) *
		visual(BarPlot;
			   width = 0.7,
			   dodge_gap = 0.15,
			   color = (:white, 0),
			   strokewidth = 1.5,
			   strokecolor = :black)
	draw!(top, barplt)

	local order = ["0-0", "0-1", "1-0", "1-1"]
	local bottom = Axis(grid[2, 1],
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
		[PolyElement(color=:blue),
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


	colsize!(trow, 1, Relative(6/18))
	colsize!(trow, 2, Relative(6/18))
	colsize!(trow, 3, Relative(6/18))

	colsize!(brow, 1, Relative(10/13))
	colsize!(brow, 2, Relative(0.25/13))
	colsize!(brow, 3, Relative(2.75/13))
	# rowsize!(f.layout, 1, Relative(1.8/7))
	
	rowsize!(f.layout, 1, Relative(3.5/8))
	rowsize!(f.layout, 2, Relative(4.5/8))

	rowsize!(brow, 1, Relative(1/7))
	rowsize!(brow, 2, Relative(3/7))
	rowsize!(brow, 3, Relative(3/7))
	# rowsize!(brow, 4, Relative(1/10))
	# rowsize!(brow, 5, Relative(2/10))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	
	# Legend(gc[3, 1],
	# 	   [PolyElement(color=:black),
	# 		PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	#        ["Modified", "Not modified"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)


	for (label, layout) in zip(["A", "B", "C", "D"],
							   [trow[1, 1], trow[1, 2], trow[1, 3], brow[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end

# ╔═╡ a7910931-bbdf-42c5-821f-6af9996933c7
begin
	local f = Figure(size = (1600, 1100))
	local trow = f[1, 1] = GridLayout()
	local brow = f[2, 1] = GridLayout()

	# local sig_color = "#800080" # "#051094"
	local sig_color = "#2E2585"
	local insig_color = "#DDDDDD"


	local df = copy(comods)[1:100, :]
	df[!, :significant] = df.qvalue .< 0.05
	local plt = data(df) *
		mapping(:phi, :pvalue => (p -> -log10(p)), color = :significant) *
		visual(Scatter; markersize = 5)
	draw!(trow[1, 1],
		  plt,
		  scales(Color = (; palette = [insig_color, sig_color]));
		  axis = (; #title = "All modification pairs",
				  	xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
					ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
				    xlabelsize = LABELSIZE,
				  	ylabelsize = LABELSIZE,
					limits = ((-1, 1), (-2, 310))))


	# # PLOT B
	# local plt_b = trow[1, 2] = GridLayout()	


	# local m6a_m6a = filter(r -> r.mod1 === "m6A" && r.mod2 === "m6A", an_sig_comods)
	# local m6a_nonm6a = filter(r -> (r.mod1 === "m6A" && r.pos2_not_m6A) || (r.mod2 === "m6A" && r.pos1_not_m6A), an_sig_comods)

	# local some_m6a = vcat(m6a_m6a, m6a_nonm6a)
	# some_m6a[!, :pair] = map((a, b) -> join([a, b], "-"), some_m6a.mod1, some_m6a.mod2)
	# local df2 = copy(some_m6a[:, [:pair, :phi]])

	# local perm_results = []
	# for i in 1:10_00#0
	# 	df2[!, :pair] = shuffle(df2.pair)
	# 	push!(perm_results, df2)
	# end
	# perm_results = vcat(perm_results...)
	

	# local ax = Axis(plt_b[1, 1],
	# 			    xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ")
	
	# density!(ax,
	# 		 m6a_m6a.phi,
	# 		 label = "m⁶A-m⁶A",
	# 		 bandwidth = 0.01,
	# 		 boundary = (-1, 1),
	# 		 color = (:blue, 0.0),
	# 		 strokecolor = :blue,
	# 		 strokewidth = 2)
	# density!(ax,
	# 		 m6a_nonm6a.phi,
	# 		 label = "m⁶A-not m⁶A",
	# 		 bandwidth = 0.01,
	# 		 boundary = (-1, 1),
	# 		 color = (:orange, 0.0),
	# 		 strokecolor = :orange,
	# 		 strokewidth = 2)
	

	# PANEL B

	local panel_b = brow[1:3, 1] = GridLayout()
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(filter(r -> abs(r.phi) > 0.1, sig_comods),
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2

	df[!, :gene] = [x[6] for x in split.(df.reference, "|")]

	local df_cancer = innerjoin(df, cancer_genes,
			  		   		    on = [:gene => :Hugo_Symbol])
	# df_cancer = df_cancer[abs.(df_cancer.phi) .>= 0.3 .&& df_cancer.mean_entropy .>= 0.75, :]
	println(size(df_cancer))

	df[!, :log10_distance] = log10.(abs.(df.pos1 .- df.pos2))
	df[!, :log10_distance_bin] = round.(df.log10_distance ./ 0.5; digits = 0)
	println("Mann U entropy ", HypothesisTests.MannWhitneyUTest(df[df.log10_distance .< log10(40), :mean_entropy], df[df.log10_distance .>= log10(40), :mean_entropy]))
	println("mean entropy ", mean(df[df.log10_distance .< log10(40), :mean_entropy]), " ", mean(df[df.log10_distance .>= log10(40), :mean_entropy]))
	println("Mann U phi ", HypothesisTests.MannWhitneyUTest(df[df.log10_distance .< log10(40), :phi], df[df.log10_distance .>= log10(40), :phi]))
	println("mean phi ", mean(df[df.log10_distance .< log10(40), :phi]), " ", mean(df[df.log10_distance .>= log10(40), :phi]))
	df[!, :combined_score] = 100 .* abs.(df.phi) .+ clamp.(-log10.(df.pvalue), 0, 100)
	local labels_per_bin = 5
	local labeled = combine(groupby(df, [:log10_distance_bin]),
							# df -> vcat(sort(df[df.phi .>= 0, :], :combined_score, rev = true)[1:(df[1, :log10_distance_bin] == 2 ? 2*labels_per_bin : labels_per_bin), :],
							# 		   sort(df[df.phi .< 0, :], :combined_score, rev = true)[1:(df[1, :log10_distance_bin] == 2 ? 2*labels_per_bin : labels_per_bin), :]))
						    df -> vcat(sort(df[df.phi .>= 0, :], :phi, by = abs, rev = true)[1:(df[1, :log10_distance_bin] <= 2 ? 0 : labels_per_bin), :],
									   sort(df[df.phi .< 0, :], :phi, by = abs, rev = true)[1:(df[1, :log10_distance_bin] == 2 ? labels_per_bin : labels_per_bin), :]))
	local high_significance = df[df.pvalue .< 1e-100 .&&
								 (.! (df.log10_distance_bin .<= 2 .&&
									 df.phi .> 0 .&&
		 							 df.phi .< 0.6)), :]
	labeled = unique(vcat(labeled, high_significance))
	# local srsf2 = df[df.gene .== "SRSF2", :]
	# labeled = vcat(labeled, srsf2[argmax(srsf2.phi):argmax(srsf2.phi), :])
	local cancer_geneset = Set(cancer_genes.Hugo_Symbol)
	labeled[!, :cancer] = [gene in cancer_geneset for gene in labeled.gene]
	# df_cancer[!, :cancer] .= true
	# labeled = leftjoin(labeled, df_cancer[:, [:gene, :cancer]],
	# 				   on = [:gene => :gene])

	

	local ax = Axis(panel_b[1, 1],
				    xlabel = "← Mutually exclusive    𝛗    Co-occurring →                     ",
				    ylabel = "Distance between associated modifications",
					xlabelsize = LABELSIZE,
				  	ylabelsize = LABELSIZE,
				    yticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
					yminorticksvisible = true,
				    yminorgridvisible = true,
					yminorticksize = 2,
					yminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]))
	local plt = data(df) *
		mapping(:phi,
				(:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))),
				color = :pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)") *
			    # color = :mean_entropy) *
		visual(Scatter; markersize = 5)
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
	Colorbar(panel_b[1, 2], limits = [0, 100], colormap = :seaborn_crest_gradient, #:viridis,
    vertical = true, label = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", labelsize = LABELSIZE, labelpadding = LABELPAD)
	colsize!(panel_b, 1, Relative(8/9))

	local points = Point2f.(labeled.phi, labeled.log10_distance)
	local labels = labeled.gene
	println(labeled)
	local label_positions = place_labels_nonoverlapping(points, labels; offset = 0.1, gene_offsets = Dict("SRSF2" => 0.45), height=0.075, width_factor=0.0475)
	for (p, lp, label, cancer) in zip(points, label_positions, labels, labeled.cancer)
		xalign = (lp[1] > p[1] ? :left : :right)
	    text!(ax, label, position=lp, align=(xalign, :center), color=:black, fontsize = 12, font = cancer ? :bold : :italic)
	    lines!(ax, [p, lp], color=:gray, linewidth=0.5)
	end
	xlims!(ax, -1.2, 1.55)


	# PANEL C

	local ax = Axis(trow[1, 2])

	# Plot SRSF2
	# plot_gene_comods(sig_comods, gtf, "ENST00000359995.10"; axis = ax, colormap = :seaborn_crest_gradient, highclip = :black)
	# text!(ax, 1600, -60; text = "SRSF2", align = (:center, :top))
	# Plot RPL22
	#plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, an_sig_comods), gtf, "ENST00000234875.9"; axis = ax, minheight=800, maxheight=2200, colormap = :seaborn_crest_gradient, highclip=:darkblue, fig=trow[1, 2])
	plot_gene_comods_elliptical_arcs(filter(r -> abs(r.phi) > 0.15, an_sig_comods), gtf, "ENST00000234875.9"; minheight = 5, maxheight = 150, dotsize=15, xlimits=(1170, 1820), legend=true, colormap = :seaborn_crest_gradient, highclip=:darkblue, fig=trow[1, 2], axis = ax, fontsize=20) # RPL22
	ylims!(ax, -150, 150)
	# local xs, xe = 1170-30, 1850+30
	# poly!(ax, Point2f[(xs, 2000), (xe, 2000), (xe, -2000), (xs, -2000)], color=(:white, 0), strokecolor=:red, strokewidth=1)
	# plot_gene_comods(sig_comods, gtf, "ENST00000234875.9"; axis = ax, colormap = :seaborn_crest_gradient, highclip = :black)
	# text!(ax, 0, -400; text = "RPL22", align = (:left, :top))
	text!(ax, 0, -30; text = "RPL22", align = (:left, :top))
	text!(ax, 1720, 130; text = "↑ Co-occurring", align = (:left, :top))
	text!(ax, 1720, -130; text = "↓ Mutually exclusive", align = (:left, :bottom))

	local inset_ax = Axis(trow[1, 2],
						  width=Relative(0.5),
						  height=Relative(0.2),
						  halign=0.05,
						  valign=0.95)
						  # backgroundcolor="#F0F0F0")
	hidedecorations!(inset_ax)
	local tx = "ENST00000234875.9"
	local name = "RPL22"
	local gene_start = minimum(gtf[gtf.transcript_id .== tx, :start])

	plot_isoforms_model!(inset_ax, "RPL22"; transcripts = Set([tx]), colors = Dict([tx => :black]), rename = Dict([tx => name]), focus=(gene_start+1170, gene_start+1820), fontsize=16, lpadding=2000)

	
	# lines!(ax, [2200, 3200], [1000, 1000], color = :black)
	Colorbar(trow[1, 3], limits = [0, 100], colormap = :seaborn_crest_gradient, #:viridis,
    vertical = true, label = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", labelsize = LABELSIZE, labelpadding = LABELPAD)
	# colsize!(trow, 3, Relative(1/15))


	
	# PLOT D


	local ax = Axis(brow[1, 2:4],
				    xtickformat = "{:,d}")# (vals -> map(format_with_commas, vals)))
			
	local gpos1 = 76066564
	local gpos2 = 76066573
	local row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
	local row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
	local tx1 = split(row1.reference, "|")[1]
	local tx2 = split(row2.reference, "|")[1]
	local name1 = split(row1.reference, "|")[5]
	local name2 = split(row2.reference, "|")[5]

	local blue = "#00aaff"
	
	# local ax = Axis(gb[1, 1:2],
	# 			    height = 110)
	plot_isoforms_model!(ax, "MDH"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => blue, tx2 => :orange]), focus = (gpos1-100, gpos2+100), rename = Dict([tx1 => name1, tx2 => name2]), lpadding=1650)
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.(tx1, an_sig_comods.reference)
	local tx2_mask = occursin.(tx2, an_sig_comods.reference)
	# local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
	# 			  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	local colorrange = (0, 100)
	local seq = tx_ref[row1.reference][(row1.pos1-10):(row1.pos2+10)]
	println(seq)
	local mod_symbols = Dict(["m⁵C" => :utriangle,
							  "m⁶A" => :cross,
							  "?" => :circle,
							  "Ψ" => :rect])
	local mod_colors = Dict(["m⁵C" => "#e69f00",
							 "m⁶A" => "#009e73",
							 "?" => "#0072b2",
							 "Ψ" => "#cc79a7"])
	local ax1, max_radius1 = arc_plot_genomic2(brow[2, 2],
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
	text!(ax1, 1, 13; text = "↑ Co-occurring", align = (:left, :top))
	text!(ax1, 1, -13; text = "↓ Mutually exclusive", align = (:left, :bottom))
	Legend(brow[2, 2],
    	   [MarkerElement(marker=mod_symbols["?"], markersize=16, color=mod_colors["?"]),
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
	local ax2, max_radius2 = arc_plot_genomic2(brow[3, 2],
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
	text!(ax2, 1, 13; text = "↑ Co-occurring", align = (:left, :top))
	text!(ax2, 1, -13; text = "↓ Mutually exclusive", align = (:left, :bottom))
	Legend(brow[3, 2],
    	   [MarkerElement(marker=mod_symbols["?"], markersize=16, color=mod_colors["?"]),
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
	local max_radius = max(max_radius1, max_radius2) * 2
	ylims!(ax1, -max_radius, max_radius)
	ylims!(ax2, -max_radius, max_radius)
	Colorbar(brow[2:3, 3], limits = colorrange, colormap = :seaborn_crest_gradient, #viridis
    vertical = true, label = "-log₁₀(𝑃-value)", labelsize = LABELSIZE, labelpadding = LABELPAD)
	

	local ref = row1.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	local df1 = DataFrame(Tables.table(probs[:, [row1.pos1+1-4, row1.pos2+1-4]] .> 0.5))
	dropmissing!(df1)
	local observed1 = contingency(df1.Column1, df1.Column2) .+ 1
	local res1 = HypothesisTests.ChisqTest(observed1)
	local data1 = map(r -> join(Int.(collect(r)), "-"), eachrow(df1))
	local counts1 = data1 |> countmap

	local ref = row2.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	local df2 = DataFrame(Tables.table(probs[:, [row2.pos1+1-4, row2.pos2+1-4]] .> 0.5))
	dropmissing!(df2)
	local observed2 = contingency(df2.Column1, df2.Column2) .+ 1
	local res2 = HypothesisTests.ChisqTest(observed2)
	local data2 = map(r -> join(Int.(collect(r)), "-"), eachrow(df2))
	local counts2 = data2 |> countmap

	
	# local fig = Figure(size = (600, 600))

	local grid = brow[2:3, 4] = GridLayout()	
	local top = Axis(grid[1, 1],
					 xlabel = "",
					 xticks = (1:4, ["", "", "", ""]),
					 ylabel = "Pattern occurrences",
					 ylabelsize = LABELSIZE)
	xlims!(top, 0.5, 4.5)
	ylims!(top, 0, 2350)

	df1 = DataFrame(tx = name1, pattern = data1)
	df2 = DataFrame(tx = name2, pattern = data2)
	local df = vcat(df1, df2)
	df = combine(groupby(df, [:tx, :pattern]), nrow => :count)
	local barplt = data(df) *
		mapping(:pattern, :count, dodge = :tx, color = :tx) *
		visual(BarPlot; width = 0.7, dodge_gap = 0.15)
	draw!(top, barplt, scales(Color = (; palette = [blue, :orange])))

	local expected = DataFrame(tx = [name1, name1, name1, name1,
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
	local barplt = data(expected) *
		mapping(:pattern, :count, dodge = :tx) *
		visual(BarPlot;
			   width = 0.7,
			   dodge_gap = 0.15,
			   color = (:white, 0),
			   strokewidth = 1.5,
			   strokecolor = :black)
	draw!(top, barplt)

	local order = ["0-0", "0-1", "1-0", "1-1"]
	local bottom = Axis(grid[2, 1],
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


	colsize!(trow, 1, Relative(5/18))
	# colsize!(trow, 2, Relative(6/18))
	# colsize!(trow, 3, Relative(6/18))

	colsize!(brow, 1, Relative(9/18))
	colsize!(brow, 2, Relative(5/18))
	colsize!(brow, 3, Relative(0.25/18))
	colsize!(brow, 4, Relative(3.75/18))
	# rowsize!(f.layout, 1, Relative(1.8/7))
	
	rowsize!(f.layout, 1, Relative(3/8))
	rowsize!(f.layout, 2, Relative(5/8))

	rowsize!(brow, 1, Relative(1/7))
	rowsize!(brow, 2, Relative(3/7))
	rowsize!(brow, 3, Relative(3/7))
	# rowsize!(brow, 4, Relative(1/10))
	# rowsize!(brow, 5, Relative(2/10))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	
	# Legend(gc[3, 1],
	# 	   [PolyElement(color=:black),
	# 		PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	#        ["Modified", "Not modified"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)


	for (label, layout) in zip(["A", "B", "C", "D", "E"],
							   [trow[1, 1], brow[1, 1], trow[1, 2], brow[1, 2], brow[2, 4]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end

# ╔═╡ 946b304a-0b43-4e3b-a9d4-550bbdf3e89e
begin
	local m6a_m6a = filter(r -> r.mod1 === "m6A" && r.mod2 === "m6A", an_sig_comods)
	m6a_m6a[!, :category] .= "m⁶A-m⁶A"
	local m6a_notm6a = filter(r -> (r.mod1 === "m6A" && r.pos2_not_m6A) || (r.mod2 === "m6A" && r.pos1_not_m6A), an_sig_comods)
	m6a_notm6a[!, :category] .= "m⁶A-not m⁶A"

	local df = vcat(m6a_m6a, m6a_notm6a)

	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(df,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2

	data(df) *
		mapping(:phi, :mean_entropy, row = :category) *
		histogram(bins = 100, normalization = :pdf) |> draw
		# visual(Scatter; markersize = 5) |> draw
end

# ╔═╡ d47b96bd-3b7e-49ca-a240-4accfa3027e7
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	# Makie.density(df.mean_entropy)
	df = df[df.mean_entropy .> 0.95, :]
	df[!, :gene] = [x[6] for x in split.(df.reference, "|")]

	local df_cancer = innerjoin(df, cancer_genes,
			  		   			on = [:gene => :Hugo_Symbol])

	sort(combine(groupby(df_cancer, [:gene, :_of_occurrence_within_resources_Column_J_P_]),
				 nrow => :count,
				 :phi => minimum => :min_phi,
				 :phi => maximum => :max_phi), :count, rev = true)
	
end

# ╔═╡ 575a5fd5-c5ed-4a05-9927-dc06e56817b9
begin
	local df = copy(sig_peaks_ivt)
	df[!, :entropy] = entropy.(eachrow(df))
	# df[df.ref_id .== "ENST00000222374.3|ENSG00000105767.3|OTTHUMG00000182735.3|OTTHUMT00000463352.2|CADM4-201|CADM4|2189|protein_coding|", :]
	df = sort(combine(groupby(df, :ref_id),
				 :entropy => (e -> mean(skipmissing(e))) => :mean_entropy,
				 nrow => :count), :mean_entropy, rev = true)
	data(df) * mapping(:count, :mean_entropy) * visual(Scatter; markersize = 4) |> draw
end

# ╔═╡ ded94e23-bdbb-4aa6-9dfd-1956ff3fc89a
begin
	local f = Figure(size = (800, 400))
	local df = innerjoin(combine(groupby(comods, :reference),
								 :phi => (p -> mean(skipmissing(p))) => :phi_mean,
					   			 :phi => (p -> std(skipmissing(p))) => :phi_std),
						 combine(groupby(sig_peaks_ivt, :ref_id), nrow => :count),
						 on = [:reference => :ref_id])


	local c = data(df) * mapping(:count) * histogram(bins = 100)
	draw!(Axis(f[1, 1], xlabel = "Mod. Count"), c)

	local m = data(df) * mapping(:phi_mean) * histogram(bins = 100)
	draw!(Axis(f[1, 2], xlabel = "mean(𝛗)"), m)

	local s = data(df) * mapping(:phi_std) * histogram(bins = 100)
	draw!(Axis(f[1, 3], xlabel = "σ(𝛗)"), s)

	local cm = data(df) * mapping(:count, :phi_mean) * visual(Scatter; markersize = 5)
	draw!(Axis(f[2, 1], xlabel = "Mod. Count", ylabel = "mean(𝛗)"), cm)
	
	local cs = data(df) * mapping(:count, :phi_std) * visual(Scatter; markersize = 5)
	draw!(Axis(f[2, 2], xlabel = "Mod. Count", ylabel = "σ(𝛗)"), cs)
	
	local ms = data(df) * mapping(:phi_mean, :phi_std) * visual(Scatter; markersize = 5)
	draw!(Axis(f[2, 3], xlabel = "mean(𝛗)", ylabel = "σ(𝛗)"), ms)

	f
end

# ╔═╡ e064a2e3-63dd-4b95-a905-e27e9d1fd09f
begin
	local df = innerjoin(comods,
			 		     combine(groupby(sig_peaks_ivt, :ref_id), nrow => :count),
					     on = [:reference => :ref_id])
end

# ╔═╡ 79a8f8c0-0578-4729-a13f-9245dfe67593
begin
	local c = contingency(comods.pvalue .< 0.01, comods.emp_pvalue .< 0.01)
	Tables.table(
	hcat(["chi 0"; "chi 1"], c); header = ["", "perm 0", "perm 1"])
end

# ╔═╡ 61fc008e-6f8c-4088-b70c-954e84626284
178144/(178144+26755)

# ╔═╡ 2bc6e9c6-1745-4701-ae0f-8c15eb2cee17
178144/(178144+1007)

# ╔═╡ 0a6a8645-8dbd-4605-b6e4-3fd54ef7d4a4
begin
	local t = 0.01 * median(comods.pvalue ./ comods.emp_pvalue)
	local c = contingency(comods.pvalue .< t, comods.emp_pvalue .< 0.01)
	Tables.table(
	hcat(["chi 0"; "chi 1"], c); header = ["", "perm 0", "perm 1"])
end

# ╔═╡ 52583e7c-20d1-42a9-a282-6c6940bb5f20
176621/(176621+21739)

# ╔═╡ 21f353bc-05b7-416e-84af-faec44839f8e
176621/(176621+2530)

# ╔═╡ cdc3d63a-446e-4958-84ab-fe7ffd33946c
begin
	local df = copy(comods)
	df[!, :significant] = df.qvalue .< 0.01
	df[!, :transcript_len] = map(t -> parse(Int, t[7]), split.(df.reference, "|"))
	df[!, :type] = map(t -> t[8], split.(df.reference, "|"))
	df = combine(groupby(df, :reference),
			     nrow => :num_pairs,
			     :significant => sum => :num_significant,
			     :transcript_len => first => :transcript_len,
				 :type => first => :type)
	df[!, :sig_perc] = df.num_significant ./ df.num_pairs
	df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))
	local df = innerjoin(df,
						 combine(groupby(sig_peaks_ivt, :ref_id), nrow => :mod_count),
						 on = [:reference => :ref_id])
	comod_summary = df
end

# ╔═╡ 7e10b591-32a0-453e-a8f5-91860a8ddc01
begin
	local df = combine(groupby(comod_summary[comod_summary.mod_count .> 3, :], :gene),
					   :sig_perc => mean => :mean_sig_perc)
	for row in eachrow(sort(df, :mean_sig_perc, rev = true))
		println(row.gene, "\t", row.mean_sig_perc)
	end
end

# ╔═╡ f2eff5a8-6f20-4577-89b5-ded89cc0bc08
data(filter(r -> r.type in ["protein_coding", "retained_intron", "nonsense_mediated_decay", "processed_transcript"], comod_summary)) * mapping(:type, :sig_perc => "% significant pairs") * visual(BoxPlot) |> draw(; axis = (; xticklabelrotation = π / 3))

# ╔═╡ 429172af-a05c-45ff-a1e6-70f70c88d158
comod_summary.type |> countmap

# ╔═╡ 0567fc6c-c45a-4cd3-96dd-acc870f35150
data(comod_summary) * mapping(:transcript_len => log10 => "log10(transcript length)", :sig_perc => "% significant pairs") * visual(Scatter) |> draw

# ╔═╡ 5c1fb6d1-7942-47b6-ab4c-c37bfedd158e
data(comod_summary) * mapping(:mod_count => "Num. mods", :sig_perc => "% significant pairs") * visual(Scatter; markersize = 5) |> draw

# ╔═╡ 0ad5996c-1b2a-4250-b857-b0ddfec3b388
map(println, comod_summary[comod_summary.mod_count .> 3 .&& comod_summary.sig_perc .> 0.5, :gene])

# ╔═╡ 144b4518-7fde-411d-86ea-56279dea4d03
begin
	local path = "/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/tail_lengths"
	local dfs = []
	for f in ["WT_1.tsv", "WT_2.tsv", "WT_SPK.tsv"]
		df = DataFrame(CSV.File("$path/$f", delim = '\t', header = ["read", "len"]))
		push!(dfs, df)
	end
	wt_tails = vcat(dfs...)
	wt_tails = Dict(wt_tails.read .=> wt_tails.len)
end

# ╔═╡ 5d8b869d-44ff-487e-a3ca-62787d706830
sig_comod.reference |> unique |> length

# ╔═╡ d74712d3-7fa1-4ff6-85ec-80aaa9eaa682
begin
	local refs = sig_comod.reference |> unique

	local co_tails = []
	local ex_tails = []
	local rlock = ReentrantLock();

	Threads.@threads for ref in refs #[1:1000]
		read_ids, probs = get_read_mods_with_ids("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)

		ref_co_tails = []
		ref_ex_tails = []

		pairs = sig_comods[sig_comods.reference .== ref .&& sig_comods.phi .> 0, :]
		for r in eachrow(pairs)
			preds = probs[:, [r.pos1 + 1 - 4, r.pos2 + 1 - 4]] .>= 0.75
			valid = preds[:, 1] .!== missing .&& preds[:, 2] .!== missing
			co = valid .&& sum.(eachrow(preds)) .!= 1
			ex = valid .&& .! co
			try
				mean_co_len = mean(skipmissing(map(r -> haskey(wt_tails, r) ? wt_tails[r] : missing, read_ids[co])))
				mean_ex_len = mean(skipmissing(map(r -> haskey(wt_tails, r) ? wt_tails[r] : missing, read_ids[ex])))
				# println(mean_co_len, " ", mean_ex_len)
				push!(ref_co_tails, mean_co_len)
				push!(ref_ex_tails, mean_ex_len)
			catch
			end
		end
		lock(rlock) do
			append!(co_tails, ref_co_tails)New
			append!(ex_tails, ref_ex_tails)
		end
	end
	local df = vcat(DataFrame(type = "Co-occurrence reads", len = co_tails), DataFrame(type = "Other reads", len = ex_tails))
	data(df) * mapping(:type, :len => "Tail length [nt]") * visual(BoxPlot) |> draw
end

# ╔═╡ 529bc9f9-3537-45bd-ae86-4b32079febf1
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	data(df) * mapping((:entropy1, :entropy2) => ((e1, e2) -> (e1 + e2)/2) => "Mean entropy",
					   :phi => "𝛗") * visual(Scatter; markersize = 5) |> draw	
end

# ╔═╡ 9e75dd8f-f3fa-4fc6-a654-c6699c6f09ee
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	df[!, :selected] = df.mean_entropy .+ abs.(df.phi) .> 1.3

	data(df) * mapping(:mean_entropy => "Mean entropy",
					   :phi => "𝛗",
					   color = :selected) * visual(Scatter; markersize = 5) |> draw	
end

# ╔═╡ 5682cadd-0a47-47be-afff-8902e92a6898


# ╔═╡ c2d014bf-553f-4a6b-bbfc-a96db1bfa11e
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	df[!, :selected] = df.mean_entropy .+ abs.(df.phi) .> 1.3
	df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))

	map(println, unique(df[df.selected, :gene]))
end

# ╔═╡ d10bc4aa-09b2-401d-860e-07adb32ff17e
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	df[!, :different_mods] = df.mod1 .!== missing .&& df.mod2 .!== missing .&& df.mod1 .!= df.mod2 #df.mean_entropy .+ abs.(df.phi) .> 1.3

	df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))
	map(println, unique(df[df.different_mods .|| df.phi .< 0, :gene]))

	data(df) * mapping(:mean_entropy => "Mean entropy",
					   :phi => "𝛗",
					   color = :different_mods) * visual(Scatter; markersize = 5) |> draw	
end

# ╔═╡ ab42f3bb-4473-4bf4-967e-85bb4cbb3532
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	# df[!, :selected] = df.mean_entropy .+ abs.(df.phi) .> 1.3
	# df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))

	df = innerjoin(df, comod_summary[:, [:reference, :mod_count]],
				   on = [:reference])

	df = sort(combine(groupby(df, :reference),
		   	  :mod_count => first => :mod_count,
		      [:mean_entropy, :phi] => ((es, ps) -> mean(es .* abs.(ps))) => :et_complexity), :et_complexity, rev = true)

	df[!, :selected] .= false
	for i in 3:5:180
		subset = (1:nrow(df))[df.mod_count .>= i .&& df.mod_count .< i + 10]
		top = first(subset, Int(round(length(subset)/10; digits = 0)))
		# top = sortperm(subset.et_complexity, rev = true)[1:min(10, nrow(subset))]
		# println(top)
		for k in top
			df[k, :selected] = true
		end
	end

	df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))
	map(r -> println(r.gene, "\t", r.max_interconnectedness),
		eachrow(combine(groupby(df[df.mod_count .> 2, :], :gene),
						:et_complexity => maximum => :max_interconnectedness)))
	# map(println, unique(df[df.selected, :gene]))

	# println(df[occursin.("CDK16", df.reference), :reference])
	
	# map(println, unique(df[df.selected, :gene]))
	data(df) * mapping(:mod_count, :et_complexity, color = :selected) * visual(Scatter; markersize = 5) |> draw
end

# ╔═╡ 387e4645-cc41-4184-9612-d242129b27e3
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	df[!, :relevance] = df.mean_entropy
	
	data(df) *
	mapping((:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))) => "Distance (log-scale)",
			:phi => "← Mutually exclusive    φ    Co-occurring →         ",
		    color = :relevance => "Mean entropy") *
	visual(Scatter; markersize = 4) |>
	draw(; axis = (; title = "Effect size of modification pair association by distance",
				   xticks = (1:0.1:4.2, ["10nt", repeat([""], 9)..., "100nt", repeat([""], 9)..., "1,000nt", repeat([""], 9)..., "10,000nt", "", ""]),
				   limits = ((0.9, 4.2), (-1, 1))))
end

# ╔═╡ adced47b-69d1-43f4-931c-15d80801991e
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))

	# df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	df[!, :size] = sqrt.(df.entropy1 .^ 2 .+ df.entropy2 .^ 2 .+ df.phi .^ 2)
	map(r -> println(r.gene, "\t", r.size), eachrow(sort(df, :size, rev = true)))
end

# ╔═╡ 29d6e7ac-3d56-4473-82b6-4d17fbda0e97
map(t -> println(t[6]), split.(unique(sig_comods.reference), "|"))

# ╔═╡ ec460cc7-a688-488f-ac70-43576329a4f7


# ╔═╡ 6d788695-bbba-443f-a8d2-ee7672fcb6cd
an_sig_comods

# ╔═╡ 922095ef-fcba-416b-92f8-aea0d1dbd961
begin
	local sig_color = "#800080" # "#051094"
	
	local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	local df = renamemods(an_comods)

	function select_m6A_nonm6A(r)
		if r.mod1 !== "m⁶A" && r.mod2 !== "m⁶A"
			return false
		elseif r.mod1 === "m⁶A"
			return !((r.chr, r.strand, r.genomicPos2) in glori_mods) && r.mod2 !== "m⁶A"
		else
			return !((r.chr, r.strand, r.genomicPos1) in glori_mods) && r.mod1 !== "m⁶A"
		end
	end
	
	df = filter(select_m6A_nonm6A, df)
	df[!, :significant] = df.qvalue .< 0.01
	df[!, :mod1] = ifelse.(df.mod1 .=== missing, "unknown", df.mod1)
	df[!, :mod2] = ifelse.(df.mod2 .=== missing, "unknown", df.mod2)
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	local counts = combine(groupby(df[df.significant, :], :pair), nrow => :count)
	counts = counts[counts.count .> 20, :]
	df = innerjoin(df, counts, on = [:pair])

	local f = Figure(size = (1000, 700))
	local plt = data(df) *
		mapping(:phi, :pvalue => (p -> -log10(p)), color = :significant, layout = :pair) *
		visual(Scatter; markersize = 5)
	draw!(f, plt, scales(Color = (; palette = [:lightgrey, sig_color]));
		  axis = (; xlabel = "← Mutually exclusive    𝛗    Co-occurring →       ", ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)", limits = ((-1, 1), (nothing, nothing))))
	f
end

# ╔═╡ 0e5025d6-3fc0-4380-9c1d-e1265ca4e58b
begin
	local f = Figure(size = (1600, 1100))
	local trow = f[1, 1] = GridLayout()
	local brow = f[2, 1] = GridLayout()

	local sig_color = "#800080" # "#051094"


	local df = copy(comods)[1:100, :]
	df[!, :significant] = df.qvalue .< 0.01
	local plt = data(df) *
		mapping(:phi, :pvalue => (p -> -log10(p)), color = :significant) *
		visual(Scatter; markersize = 5)
	draw!(trow[1, 1],
		  plt,
		  scales(Color = (; palette = [:grey, sig_color]));
		  axis = (; #title = "All modification pairs",
				  	xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
					ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
					limits = ((-1, 1), (-2, 310))))


	# local df = renamemods(an_comods)
	# df[!, :significant] = df.qvalue .< 0.05

	# # local df = copy(an_sig_comods)
	# df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	# df = df[in.(df.pair, Ref(Set(["m⁵C-m⁶A"]))), :]
	# # df[!, :type] = map(assoc_type, df.residual00, df.residual01, df.residual10, df.residual11)
	# # df[!, :phi] .*= ifelse.(df.type .== :comod, 1, -1)
	# local plt = data(df) *
	# 	mapping(:phi,
	# 			:pvalue => (p -> -log10(p)) => "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
	# 		    color = :significant) *
	# 	visual(Scatter; markersize = 5)
	# draw!(trow[1, 2],
	# 	  plt,
	# 	  scales(Color = (; palette = [:grey, sig_color]));
	# 	  axis = (; title = "m⁵C-m⁶A",
	# 			    xlabel = "← Mutually exclusive    𝛗    Co-occurring →       ",
	# 				ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
	# 				limits = ((-1, 1), (-2, 310))))

	# === PLOT B VARIANT 1

	# local df = copy(an_sig_comods)
	# df[!, :type] = sign.(df.phi)
	# df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	# df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	# # df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	# df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	# df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	# df[!, :lfc] = log2.(df.count ./ df.count_1)
	# df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	# df[!, :conf] = 1.96 .* df.SE
	# df[!, :total] = df.count .+ df.count_1
	# df = sort(df[df.count .> 30 .&& df.count_1 .> 30, :], :total)
	# local ax = Axis(trow[1, 2],
	# 				xlabel = "logFC (co-occurring/mutually-exclusive)",
	# 				ylabel = "Modification pair",
	# 			    yticks = (1:nrow(df), map((p, t) -> "$p ($t)", df.pair, df.total)))
	# for (i, row) in enumerate(eachrow(df))
	# 	scatter!(ax, [row.lfc], [i], color = :black)
	# 	lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)
	# end
	# xlims!(ax, [0, 4])

	# === PLOT B VARIANT 2
	
	# local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	# local df = renamemods(an_sig_comods)
	
	# df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# # df[!, :significant] = df.qvalue .< 0.01
	# df[!, :mod1] = ifelse.(df.mod1 .=== missing, "unknown non-m⁶A", df.mod1)
	# df[!, :mod2] = ifelse.(df.mod2 .=== missing, "unknown non-m⁶A", df.mod2)
	# local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	# local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	# df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	# df[!, :type] = sign.(df.phi)
	# local counts = combine(groupby(df, :pair), nrow => :count)
	# counts = counts[counts.count .> 10, :]
	# df = innerjoin(df, counts, on = [:pair])

	# local df2 = copy(df)

	# local perm_results = []
	# for i in 1:10_000
	# 	# df2[!, :mod1] = shuffle(df2.mod1)
	# 	# df2[!, :mod2] = shuffle(df2.mod2)
	# 	df2[!, :pair] = shuffle(df2.pair)
	# 	x = combine(groupby(df2, [:pair, :type]), nrow => :count)
	# 	x = innerjoin(x[x.type .== 1, [:pair, :count]],
	# 				  x[x.type .== -1, [:pair, :count]],
	# 				  on = :pair,
	# 				  makeunique = true)
	# 	x[!, :lfc] = log2.(x.count ./ x.count_1)
	# 	push!(perm_results, x)
	# end
	# perm_results = vcat(perm_results...)
	# # perm_results = combine(groupby(perm_results, :pair),
	# # 					   :lfc => mean => :perm_mean_lfc,
	# # 					   :lfc => std => :perm_std_lfc,
	# # 					   :lfc => (lfc -> 1.96 * std(lfc)/sqrt(length(lfc))) => :perm_conf)
	# local baseline = mean(perm_results.lfc)
	
	# df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	# # df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	# df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	# df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	# df[!, :lfc] = log2.(df.count ./ df.count_1)
	# df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	# df[!, :conf] = 1.96 .* df.SE
	# df[!, :total] = df.count .+ df.count_1
	# df = sort(df[df.count .> 2 .&& df.count_1 .> 2, :], :total)

	# # df = sort(innerjoin(df, perm_results, on = [:pair]), :total)
	
	# local ax = Axis(trow[1, 2],
	# 				xlabel = "logFC (co-occurring/mutually-exclusive)",
	# 				ylabel = "Modification pair",
	# 			    yticks = (1:nrow(df), map((p, t) -> "$p\n($t)", df.pair, df.total)))
	# lines!(ax, [baseline, baseline], [0.5, nrow(df)+0.5], color = :orange, linestyle = :dash)
	# for (i, row) in enumerate(eachrow(df))
	# 	scatter!(ax, [row.lfc], [i], color = :black)
	# 	lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)

	# 	# scatter!(ax, [row.perm_mean_lfc], [i], color = :orange)
	# 	# lines!(ax, [row.perm_mean_lfc - row.perm_conf, row.perm_mean_lfc + row.perm_conf], [i, i], color = :orange)
	# end
	# xlims!(ax, [-1, 5])

	# PLOT B VARIANT 3

	local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	local df = renamemods(an_sig_comods)
	
	df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# df[!, :significant] = df.qvalue .< 0.01
	df[!, :mod1] = ifelse.(df.mod1 .=== missing, "Unclassified", df.mod1)
	df[!, :mod2] = ifelse.(df.mod2 .=== missing, "Unclassified", df.mod2)
	local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	df[!, :type] = sign.(df.phi)
	local counts = combine(groupby(df, :pair), nrow => :count)
	counts = counts[counts.count .> 10, :]
	df = innerjoin(df, counts, on = [:pair])
	df[!, :pair_label] = map((p, c) -> "$p\n($c)", df.pair, df.count)

	sort!(df, :count)

	local groups = combine(groupby(df, :pair),
						   :phi => (ps -> skipmissing(ps)) => :phis)
	groups = map(collect, groups.phis)
	
	println(HypothesisTests.KruskalWallisTest(groups...))
	
	local plt = data(df) *
		mapping(:pair_label => presorted => "Modification pair",
				:phi => "← Mutually exclusive    𝛗    Co-occurring →        ") *
		visual(Violin;
			   orientation = :horizontal,
			   side = :right,
			   color = "#FF33B5")
	draw!(trow[1, 2], plt; axis = (; limits = ((-1, 1), (0.7, length(groups) + 0.7))))


	local panel_c = trow[1, 3] = GridLayout()
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2

	local plt = data(df) *
		mapping((:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))) => "Distance",
				:phi => "← Mutually exclusive    𝛗    Co-occurring →         ",
			    color = :mean_entropy) *
		visual(Scatter; markersize = 5)
	draw!(panel_c[1, 1],
		  plt;
		  axis = (;
			# title = "Modification pairs' association strength by distance between them",
			# xticks = (1:0.1:4.2,
			# 		  ["10nt", repeat([""], 9)..., "100nt", repeat([""], 9)..., "1,000nt", repeat([""], 9)..., "10,000nt", "", ""]),
			# xticks = (map(log10, [(10:10:100)..., (200:100:1000)..., (2000:1000:10000)..., 11000, 12000]), ["10nt", repeat([""], 8)..., "100nt", repeat([""], 8)..., "1,000nt", repeat([""], 8)..., "10,000nt", "", ""]),
			# xticksize = [3.0, repeat([1], 8)..., 3, repeat([1], 8)..., 3, repeat([1], 8)..., 3, repeat([1], 2)...],
			# xminorticksvisible = true,
		 #    xminorgridvisible = true,
			# xminorticksize = 1.5,
			# xminorticks = IntervalsBetween(5),
			xticks = (1:1:4, ["10nt", "100nt", "1,000nt", "10,000nt"]),
			xminorticksvisible = true,
		    xminorgridvisible = true,
			xminorticksize = 2,
			xminorticks = map(log10, [(20:10:90)..., (200:100:900)..., (2000:1000:9000)..., 11000, 12000]),
			limits = ((0.9, 4.2), (-1, 1))))
	Colorbar(panel_c[1, 2], limits = [0, 1], colormap = :viridis,
    vertical = true, label = "Information content")
	colsize!(panel_c, 1, Relative(5/6))

	

	local ax = Axis(brow[1, 1:3],
				    xtickformat = "{:,d}")# (vals -> map(format_with_commas, vals)))
				    # height = 110)
	# plot_isoforms_model!(ax, "PRRC2B"; transcripts = Set(["ENST00000682501.1", "ENST00000684596.1", "ENST00000683519.1"]), colors = Dict(["ENST00000682501.1" => :blue, "ENST00000684596.1" => :orange]), focus = (131496129-150, 131496179+200), rename = Dict(["ENST00000682501.1" => "PRRC2B-201", "ENST00000684596.1" => "PRRC2B-208", "ENST00000683519.1" => "PRRC2B-206"]))
	
	# # local ax = Axis(gb[1, 1])
	# local tx1_mask = occursin.("ENST00000682501.1", an_sig_comods.reference)
	# local tx2_mask = occursin.("ENST00000684596.1", an_sig_comods.reference)
	# local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
	# 			  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	# arc_plot_genomic2(brow[2, 1],
	# 				  renamemods(an_sig_comods[tx1_mask, :]),
	# 				  an_sig_peaks_ivt[occursin.("ENST00000682501.1", an_sig_peaks_ivt.ref_id), :];
	# 				  range=(5000, 5300),
	# 				  grange=(131496129-150, 131496179+100),
	# 				  colorrange=colorrange,
	# 				  title="PRRC2B-201",
	# 				  highlight=(5089, 5139),
	# 				  spinecolor=:blue)
	# arc_plot_genomic2(brow[3, 1],
	# 				  renamemods(an_sig_comods[tx2_mask, :]),
	# 				  an_sig_peaks_ivt[occursin.("ENST00000684596.1", an_sig_peaks_ivt.ref_id), :];
	# 				  range=(0, 12000),
	# 				  grange=(131496129-150, 131496179+100),
	# 				  colorrange=colorrange,
	# 				  title="PRRC2B-208",
	# 				  highlight=(7211, 7261),
	# 				  spinecolor=:orange)
	# Colorbar(brow[2:3, 2], limits = colorrange, colormap = :viridis,
 #    vertical = true, label = "-log₁₀(𝑃-value)")
	local gpos1 = 76066564
	local gpos2 = 76066573
	local row1 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][2, :]
	local row2 = an_sig_comods[an_sig_comods.genomicPos1 .=== gpos1 .&& an_sig_comods.genomicPos2 .=== gpos2, :][1, :]
	local tx1 = split(row1.reference, "|")[1]
	local tx2 = split(row2.reference, "|")[1]
	local name1 = split(row1.reference, "|")[5]
	local name2 = split(row2.reference, "|")[5]
	
	# local ax = Axis(gb[1, 1:2],
	# 			    height = 110)
	plot_isoforms_model!(ax, "MDH"; transcripts = Set([tx1, tx2]), colors = Dict([tx1 => :blue, tx2 => :orange]), focus = (gpos1-100, gpos2+100), rename = Dict([tx1 => name1, tx2 => name2]))
	
	# local ax = Axis(gb[1, 1])
	local tx1_mask = occursin.(tx1, an_sig_comods.reference)
	local tx2_mask = occursin.(tx2, an_sig_comods.reference)
	local colorrange = (-log10(maximum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])),
				  		-log10(minimum(an_sig_comods[tx1_mask .|| tx2_mask, :pvalue])))
	local ax1, max_radius1 = arc_plot_genomic2(brow[2, 1],
					  renamemods(an_sig_comods[tx1_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx1, an_sig_peaks_ivt.ref_id), :];
					  range=(row1.pos1-2000, row1.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name1,
					  highlight=(row1.pos1, row1.pos2),
					  spinecolor=:blue)
	text!(ax1, 2, 8; text = "Co-occurring", align = (:left, :top))
	text!(ax1, 2, -8; text = "Mutually exclusive", align = (:left, :bottom))
	local ax2, max_radius2 = arc_plot_genomic2(brow[3, 1],
					  renamemods(an_sig_comods[tx2_mask, :]),
					  an_sig_peaks_ivt[occursin.(tx2, an_sig_peaks_ivt.ref_id), :];
					  range=(row2.pos1-2000, row2.pos2+2000),
					  grange=(gpos1-100, gpos2+100),
					  colorrange=colorrange,
					  title=name2,
					  highlight=(row2.pos1, row2.pos2),
					  spinecolor=:orange)
	text!(ax2, 2, 8; text = "Co-occurring", align = (:left, :top))
	text!(ax2, 2, -8; text = "Mutually exclusive", align = (:left, :bottom))
	local max_radius = max(max_radius1, max_radius2) * 2
	ylims!(ax1, -max_radius, max_radius)
	ylims!(ax2, -max_radius, max_radius)
	Colorbar(brow[2:3, 2], limits = colorrange, colormap = :viridis,
    vertical = true, label = "-log₁₀(𝑃-value)")
	

	local ref = row1.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df1 = DataFrame(Tables.table(probs[:, [row1.pos1+1-4, row1.pos2+1-4]] .> 0.5))
	dropmissing!(df1)
	# local observed1 = contingency(df.Column1, df.Column2) .+ 1
	# local res1 = HypothesisTests.ChisqTest(observed1)

	local ref = row2.reference
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df2 = DataFrame(Tables.table(probs[:, [row2.pos1+1-4, row2.pos2+1-4]] .> 0.5))
	dropmissing!(df2)
	# local observed2 = contingency(df.Column1, df.Column2) .+ 1
	# local res2 = HypothesisTests.ChisqTest(observed2)


	local grid = brow[2:3, 3] = GridLayout()

	local titles = [name1, name2]
	for (i, df) in enumerate([df1, df2])
		local observed = contingency(df.Column1, df.Column2) .+ 1
		local res = HypothesisTests.ChisqTest(observed)
		expected = Int.(round.(res.expected; digits = 0))

		# expected = expected[2:-1:1, :]
		# observed = observed[2:-1:1, :]

		values = [expected[1, 1] observed[1, 1];
				  expected[1, 2] observed[1, 2];
				  expected[2, 1] observed[2, 1];
				  expected[2, 2] observed[2, 2]]
		values = values[end:-1:1, :]
		color = Float64.(values)
		color[:, :] ./= color[:, 1]
		color = log2.(color)

		ax = Axis(grid[i, 1],
				  title = titles[i])
		hidedecorations!(ax)
				  # aspect = 1,
				  # xticks = (1:2, ["Unmodified", "Modified"]),
				  # yticks = (1:2, ["Unmodified", "Modified"]),
				  # yticklabelrotation = π/2)

		table = fill(0.0, (5, 3))
		table[1:4, 2:3] = color
		
		# println(transpose(log2.(observed ./ expected)))
	    # Plot the heatmap
	    heatmap!(ax, transpose(table),
				 colormap = :RdBu_9,
				 colorrange = (-maximum(abs.(color))*1.5, maximum(abs.(color))*1.5))

		text!(ax, "Pattern", position = (1, 5), fontsize = 13.5, align = (:center, :center), color = :black, font = :bold)
		text!(ax, "Expected", position = (2, 5), fontsize = 13.5, align = (:center, :center), color = :black, font = :bold)
		text!(ax, "Observed", position = (3, 5), fontsize = 13.5, align = (:center, :center), color = :black, font = :bold)
		text!(ax, "0-0", position = (1, 4), fontsize = 13.5, align = (:center, :center), color = :black, font = :bold)
		text!(ax, "0-1", position = (1, 3), fontsize = 13.5, align = (:center, :center), color = :black, font = :bold)
		text!(ax, "1-0", position = (1, 2), fontsize = 13.5, align = (:center, :center), color = :black, font = :bold)
		text!(ax, "1-1", position = (1, 1), fontsize = 13.5, align = (:center, :center), color = :black, font = :bold)
		# # Annotate cells with values (optional)
	    for j in 1:4, k in 1:2
	        text!(ax, string(values[j, k]),
	              position = (k+1, j),
				  fontsize = 13.5,
	              align = (:center, :center),
	              color = :black)
	    end
	end
	

	# local matrices = [(Int.(round.(res1.expected; digits = 0)), observed1),
	# 				  (Int.(round.(res2.expected; digits = 0)), observed2)]
	# local titles = [name1, name2]
	# for (i, (expected, observed)) in enumerate(matrices)
	# 	expected = expected[2:-1:1, :]
	# 	observed = observed[2:-1:1, :]
		
	#     # Create an axis for this heatmap
	#     ax = Axis(grid[i, 1],
	# 			  title = titles[i],
	# 			  aspect = 1,
	# 			  xticks = (1:2, ["Unmodified", "Modified"]),
	# 			  yticks = (1:2, ["Unmodified", "Modified"]),
	# 			  yticklabelrotation = π/2)

	# 	println(transpose(log2.(observed ./ expected)))
	#     # Plot the heatmap
	#     heatmap!(ax, transpose(log2.(observed ./ expected)),
	#              colormap = :diverging_bwr_40_95_c42_n256,
	# 			 colorrange = (-2.5, 2.5))
	#              # colorrange = (0, max_count))
	
	#     # Add row and column labels
	# 	# if i < 1
	# 	ax.xlabel = "$(format_with_commas(gpos1))\nm⁵C"
	# 	ax.ylabel = "?\n$(format_with_commas(gpos2))"
	# 	# else
	# 	# 	ax.xlabel = "Pos. 7211 (m⁶A)"
	# 	# 	ax.ylabel = "Pos. 7261 (m⁵C)"
	# 	# end
	
	#     # Annotate cells with values (optional)
	#     for j in 1:2, k in 1:2
	#         text!(ax, string(observed[j, k]) * "/" * string(expected[j, k]),
	#               position = (k, j),
	# 			  fontsize = 13.5,
	#               align = (:center, :center),
	#               color = :black)
	#     end
	# end

	# Legend(brow[4, 1],
	# 	   [LineElement(points = [Point2f(cos(x)/2+0.3, sin(x) + 1) for x in 0:-0.1:-π]),
	# 		LineElement(points = [Point2f(cos(x)/2+0.3, sin(x)) for x in 0:0.1:π])],
	#        ["Mutually exclusive", "Co-occurring"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)


	colsize!(trow, 1, Relative(5/18))
	colsize!(trow, 2, Relative(5/18))
	colsize!(trow, 3, Relative(8/18))

	colsize!(brow, 1, Relative(10.5/13))
	colsize!(brow, 2, Relative(0.5/13))
	colsize!(brow, 3, Relative(2/13))
	# rowsize!(f.layout, 1, Relative(1.8/7))
	
	rowsize!(f.layout, 1, Relative(3/8))
	rowsize!(f.layout, 2, Relative(5/8))

	rowsize!(brow, 1, Relative(1/7))
	rowsize!(brow, 2, Relative(3/7))
	rowsize!(brow, 3, Relative(3/7))
	# rowsize!(brow, 4, Relative(1/10))
	# rowsize!(brow, 5, Relative(2/10))

	# Legend(ga[2, 1:2],
	# 	   [PolyElement(color=Makie.wong_colors()[1]),
	# 	    PolyElement(color=Makie.wong_colors()[2])],
	#        ["Negative association", "Positive association"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)
	
	# Legend(gc[3, 1],
	# 	   [PolyElement(color=:black),
	# 		PolyElement(color="#dbdbdb", strokecolor=:darkgrey, strokewidth=1.5)],
	#        ["Modified", "Not modified"],
	# 	   orientation = :horizontal,
	# 	   framevisible = false)


	for (label, layout) in zip(["A", "B", "C", "D"],
							   [trow[1, 1], trow[1, 2], trow[1, 3], brow[1, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
	
end

# ╔═╡ a5f61fab-b1b4-4834-8945-b63947f81776
begin
	local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	local df = renamemods(an_sig_comods)

	
	df = filter(select_m6A_nonm6A, df)
	df[!, :mod1] = ifelse.(df.mod1 .=== missing, "unknown", df.mod1)
	df[!, :mod2] = ifelse.(df.mod2 .=== missing, "unknown", df.mod2)
	df[!, :pair] = map(get_pair, df.mod1, df.mod2)
	local counts = combine(groupby(df, :pair), nrow => :count)
	counts = counts[counts.count .> 20, :]
	df = innerjoin(df, counts, on = [:pair])

	df[!, :distance] = abs.(df.pos1 - df.pos2)
	df[!, :relation] = ifelse.(df.relation .== :comod, "Co-occurrence", "Mutually exclusive")

	# local f = Figure(size = (1000, 300))
	local plt = data(df) *
		mapping(:pair => "Modifications pair",
				:distance => "Distance",
				color=:relation => "Association",
				dodge=:relation) *
		visual(BoxPlot; show_outliers = false) |> draw(; axis = (; width = 1000, height = 300))
	# draw!(f[1, 1], plt)
	# f
end

# ╔═╡ be6569c6-00d7-4201-b003-6b6852a3cd02
begin
	local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	local df = renamemods(an_sig_comods)
	
	df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# df[!, :significant] = df.qvalue .< 0.01
	df[!, :mod1] = ifelse.(df.mod1 .=== missing, "unknown non-m⁶A", df.mod1)
	df[!, :mod2] = ifelse.(df.mod2 .=== missing, "unknown non-m⁶A", df.mod2)
	local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	df[!, :type] = sign.(df.phi)
	local counts = combine(groupby(df, :pair), nrow => :count)
	counts = counts[counts.count .> 10, :]
	df = innerjoin(df, counts, on = [:pair])

	local df2 = copy(df)

	local perm_results = []
	for i in 1:10_000
		# df2[!, :mod1] = shuffle(df2.mod1)
		# df2[!, :mod2] = shuffle(df2.mod2)
		df2[!, :pair] = shuffle(df2.pair)
		x = combine(groupby(df2, [:pair, :type]), nrow => :count)
		x = innerjoin(x[x.type .== 1, [:pair, :count]],
					  x[x.type .== -1, [:pair, :count]],
					  on = :pair,
					  makeunique = true)
		x[!, :lfc] = log2.(x.count ./ x.count_1)
		push!(perm_results, x)
	end
	perm_results = vcat(perm_results...)
	# perm_results = combine(groupby(perm_results, :pair),
	# 					   :lfc => mean => :perm_mean_lfc,
	# 					   :lfc => std => :perm_std_lfc,
	# 					   :lfc => (lfc -> 1.96 * std(lfc)/sqrt(length(lfc))) => :perm_conf)
	local baseline = mean(perm_results.lfc)
	
	df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	# df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	df[!, :lfc] = log2.(df.count ./ df.count_1)
	df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	df[!, :conf] = 1.96 .* df.SE
	df[!, :total] = df.count .+ df.count_1
	df = sort(df[df.count .> 2 .&& df.count_1 .> 2, :], :total)

	# df = sort(innerjoin(df, perm_results, on = [:pair]), :total)
	
	local f = Figure()
	local ax = Axis(f[1, 1],
					xlabel = "logFC (co-occurring/mutually-exclusive)",
					ylabel = "Modification pair",
				    yticks = (1:nrow(df), map((p, t) -> "$p ($t)", df.pair, df.total)))
	lines!(ax, [baseline, baseline], [0.5, nrow(df)+0.5], color = :orange, linestyle = :dash)
	for (i, row) in enumerate(eachrow(df))
		scatter!(ax, [row.lfc], [i], color = :black)
		lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)

		# scatter!(ax, [row.perm_mean_lfc], [i], color = :orange)
		# lines!(ax, [row.perm_mean_lfc - row.perm_conf, row.perm_mean_lfc + row.perm_conf], [i, i], color = :orange)
	end
	
	xlims!(ax, [-1, 5])
	f
end

# ╔═╡ db0b4044-a42a-4c4e-9911-c24baa258818
# ╠═╡ disabled = true
#=╠═╡
begin
	local glori_mods = Set(zip(glori.chrom, glori.strand, glori.start))
	local df = renamemods(an_sig_comods)
	
	df = filter(r -> select_m6A_nonm6A(r) || (r.mod1 === "m⁶A" && r.mod2 === "m⁶A"), df)
	# df[!, :significant] = df.qvalue .< 0.01
	df[!, :mod1] = ifelse.(df.mod1 .=== missing, "Unclassified", df.mod1)
	df[!, :mod2] = ifelse.(df.mod2 .=== missing, "Unclassified", df.mod2)
	local m1 = ifelse.(df.mod1 .== "m⁶A", df.mod1, df.mod2)
	local m2 = ifelse.(df.mod1 .== "m⁶A", df.mod2, df.mod1) # map(get_pair, df.mod1, df.mod2)
	df[!, :pair] = map((a, b) -> join([a, b], "-"), m1, m2)
	df[!, :type] = sign.(df.phi)
	local counts = combine(groupby(df, :pair), nrow => :count)
	counts = counts[counts.count .> 40, :]
	df = innerjoin(df, counts, on = [:pair])
	df[!, :pair_label] = map((p, c) -> "$p ($c)", df.pair, df.count)

	sort!(df, :count, rev = true)

	local ordered_labels = collect(unique(df.pair_label))

	# local groups = combine(groupby(df, :pair),
	# local groups = combine(groupby(df[df.pair .!== "m⁶A-m⁶A" .&& df.pair .!== "m⁶A-Unclassified", :], :pair),
	local groups = combine(groupby(df[df.pair .!== "m⁶A-m⁶A" .&& df.pair .!== "m⁶A-Nm" .&& df.pair .!== "m⁶A-Unclassified", :], :pair),
						   :phi => (ps -> skipmissing(ps)) => :phis)
	groups = map(collect, groups.phis)
	
	println(HypothesisTests.KruskalWallisTest(groups...))

	df = df[:, [:pair_label, :phi]]

	local perm_results = []
	perm_ks_stats = Dict([pair => [] for pair in ordered_labels])
	local PERMS = 5_000
	for i in 1:PERMS
		df2 = copy(df)
		df2[!, :pair_label] = shuffle(df2.pair_label)
		for pair in ordered_labels
			if occursin("m⁶A-m⁶A", pair) || occursin("m⁶A-Unclassified", pair)
				continue
			end
			res = HypothesisTests.ApproximateTwoSampleKSTest(df2[df2.pair_label .== pair, :phi], df2[.! occursin.("m⁶A-m⁶A", df2.pair_label) .&& .! occursin.("m⁶A-Unclassified", df2.pair_label) .&& df2.pair_label .!= pair, :phi])
			ks_stat = sqrt(res.n_x*res.n_y/(res.n_x+res.n_y))*res.δ
			push!(perm_ks_stats[pair], ks_stat)
		end
		push!(perm_results, df2)
	end
	perm_results = vcat(perm_results...)
	sort!(perm_results, :pair_label, by = p -> findfirst(ordered_labels .== p))

	for pair in ordered_labels
		if occursin("m⁶A-m⁶A", pair) || occursin("m⁶A-Unclassified", pair)
			continue
		end
		res = HypothesisTests.ApproximateTwoSampleKSTest(df[df.pair_label .== pair, :phi], df[.! occursin.("m⁶A-m⁶A", df.pair_label) .&& .! occursin.("m⁶A-Unclassified", df.pair_label) .&& df.pair_label .!= pair, :phi])
		ks_stat = sqrt(res.n_x*res.n_y/(res.n_x+res.n_y))*res.δ
		emp_pval = sum(perm_ks_stats[pair] .> ks_stat) / PERMS
		println("$pair emp_pvalue = $emp_pval")
	end

	df[!, :category] .= :original
	perm_results[!, :category] .= :permutations
	local combined_data = vcat(df, perm_results)
	
	local f = Figure()

	local axes = []
	for (i, pair) in enumerate(ordered_labels)
		ax = Axis(f[i, 1],
				  xlabel = "← Mutually exclusive    𝛗    Co-occurring →        ",
				  yticks = 1:2:9,
				  yticklabelsize = 9)
		density!(ax,
				 combined_data[combined_data.pair_label .== pair .&& combined_data.category .== :permutations, :phi],
				 color = (:grey, 0.2),
				 linestyle = :dash,
				 strokewidth = 1.5,
				 strokecolor = (:grey, 1))
		density!(ax,
				 combined_data[combined_data.pair_label .== pair .&& combined_data.category .== :original, :phi],
				 color = ("#FF33B5", 0.2),
				 linestyle = :solid,
				 strokewidth = 1.5,
				 strokecolor = ("#FF33B5", 1))
		Legend(
        	f[i, 1],
			[MarkerElement(marker = 'x', markersize = 0)],
			[pair],
	        tellheight = false,
	        tellwidth = false,
			halign = :left,
			valign = :top,
	        margin = (10, 10, 10, 10),
			patchsize = (0, 0),
			framevisible = false,
			backgroundcolor = (:white, 0),
			labelfont = :bold
    	)
		xlims!(ax, (-1, 1))
		ylims!(ax, (-1, 10))
		if i < length(ordered_labels)
			hidexdecorations!(ax, grid = false)
		end
		ax.ylabelrotation = 0
		push!(axes, ax)	
	end
	linkyaxes!(axes...)
	Label(f[1:length(ordered_labels), 0], "pdf", rotation = pi/2)
	rowgap!(f.layout, 0)

	
	# local plt1 = data(perm_results) * mapping(:pair_label => presorted => "Modification pair", :phi => "← Mutually exclusive    𝛗    Co-occurring →        ") * visual(Violin; orientation = :horizontal, side = :right, color = (:grey, 0.5), show_median = false)



	
	# # df = vcat(df[:, [:pair_label, :phi]], DataFrame(pair_label = "Null", phi = perm_results))
	
	
	# local plt2 = data(df) * mapping(:pair_label => presorted => "Modification pair", :phi => "← Mutually exclusive    𝛗    Co-occurring →        ") * visual(Violin; orientation = :horizontal, side = :right, color = ("#FF33B5", 0.5), show_median = false)
	# draw!(f[1, 1], plt2; axis = (; limits = ((-1, 1), nothing)))
	# draw!(f[1, 1], plt1; axis = (; limits = ((-1, 1), nothing)))

	f

	
	# local df2 = copy(df)

	# local perm_results = []
	# for i in 1:10_000
	# 	# df2[!, :mod1] = shuffle(df2.mod1)
	# 	# df2[!, :mod2] = shuffle(df2.mod2)
	# 	df2[!, :pair] = shuffle(df2.pair)
	# 	x = combine(groupby(df2, [:pair, :type]), nrow => :count)
	# 	x = innerjoin(x[x.type .== 1, [:pair, :count]],
	# 				  x[x.type .== -1, [:pair, :count]],
	# 				  on = :pair,
	# 				  makeunique = true)
	# 	x[!, :lfc] = log2.(x.count ./ x.count_1)
	# 	push!(perm_results, x)
	# end
	# perm_results = vcat(perm_results...)
	# # perm_results = combine(groupby(perm_results, :pair),
	# # 					   :lfc => mean => :perm_mean_lfc,
	# # 					   :lfc => std => :perm_std_lfc,
	# # 					   :lfc => (lfc -> 1.96 * std(lfc)/sqrt(length(lfc))) => :perm_conf)
	# local baseline = mean(perm_results.lfc)
	
	# df = filter(r -> r.mod1 !== missing && r.mod2 !== missing, df)
	# # df = filter(r -> r.pair in ["m5C-m6A", "Y-m5C", "Y-m6A", "m6A-m7G", "m5C-m7G"], df)
	# df = sort(combine(groupby(df, [:pair, :type]), nrow => :count))

	# df = innerjoin(df[df.type .== 1, [:pair, :count]], df[df.type .== -1, [:pair, :count]], on = :pair, makeunique = true)
	# df[!, :lfc] = log2.(df.count ./ df.count_1)
	# df[!, :SE] = sqrt.(1 ./ df.count .+ 1 ./ df.count_1)
	# df[!, :conf] = 1.96 .* df.SE
	# df[!, :total] = df.count .+ df.count_1
	# df = sort(df[df.count .> 2 .&& df.count_1 .> 2, :], :total)

	# # df = sort(innerjoin(df, perm_results, on = [:pair]), :total)
	
	# local f = Figure()
	# local ax = Axis(f[1, 1],
	# 				xlabel = "logFC (co-occurring/mutually-exclusive)",
	# 				ylabel = "Modification pair",
	# 			    yticks = (1:nrow(df), map((p, t) -> "$p ($t)", df.pair, df.total)))
	# lines!(ax, [baseline, baseline], [0.5, nrow(df)+0.5], color = :orange, linestyle = :dash)
	# for (i, row) in enumerate(eachrow(df))
	# 	scatter!(ax, [row.lfc], [i], color = :black)
	# 	lines!(ax, [row.lfc - row.conf, row.lfc + row.conf], [i, i], color = :black)

	# 	# scatter!(ax, [row.perm_mean_lfc], [i], color = :orange)
	# 	# lines!(ax, [row.perm_mean_lfc - row.perm_conf, row.perm_mean_lfc + row.perm_conf], [i, i], color = :orange)
	# end
	
	# xlims!(ax, [-1, 5])
	# f
end
  ╠═╡ =#

# ╔═╡ 33c0f5a5-1eab-4857-b6e5-dc1625a05c61
begin
	local df = copy(an_sig_comods)
	df = df[df.mod1 .!== missing .&&
			df.mod2 .!== missing .&&
			(df.mod1 .== "m6A" .||
		 	 df.mod2 .== "m6A") .&&
			.! (df.mod1 .== "m6A" .&& df.mod2 .== "m6A"), :]
	df[!, :pair] = map((a, b) -> join(sort([a, b]), "-"), df.mod1, df.mod2)
	df = df[:, [:pair, :phi]]
	for pair in unique(df.pair)
		res = HypothesisTests.ApproximatePermutationTest(df[df.pair .== pair, :phi],
															   df[df.pair .!= pair, :phi],
															   median,
															   10000)
		println("$pair emp_pvalue = $(HypothesisTests.pvalue(res))")
	end
	println("")

	df = df[in.(df.pair, Ref(["m5C-m6A", "Y-m6A", "m6A-m7G"])), :]

	local groups = combine(groupby(df, :pair),
						   :phi => (ps -> skipmissing(ps)) => :phis)
	CSV.write("/home/mzdravkov/group_phis.csv", groups)
	groups = map(collect, groups.phis)
	local adtest = HypothesisTests.KSampleADTest(groups...)
	local ad_stat = (adtest.A²k - adtest.k + 1) / adtest.σ
	println(adtest)

	
	local sim_adtest = HypothesisTests.KSampleADTest(groups..., nsim = 20_000)
	println(sim_adtest)

	# local perm_stats = []
	# local PERM = 20_000
	# for i in 1:PERM
	# 	x = HypothesisTests.KSampleADTest(groups...)
	# 	stat = (x.A²k - x.k + 1) / x.σ
	# 	push!(perm_stats, stat)
	# end
	# println("Emp pvalue = ", sum(perm_stats .> ad_stat) / PERM)
	df[df.pair .== "Nm-m6A", :]
end

# ╔═╡ 54f04bdf-240e-4c55-83a1-7539e8679b95


# ╔═╡ a5c05aef-84f5-4a58-a663-0478c4fdcc35
begin
	local df = innerjoin(combine(groupby(ivt, :ref_id),
							     :pos => minimum => :min_pos,
							     :pos => maximum => :max_pos),
						 comod_summary,
						 on = [:ref_id => :reference])
	well_covered_refs = df[(df.max_pos .- df.min_pos) ./ df.transcript_len .> 0.9, :ref_id]


	data(df[in.(df.ref_id, Ref(well_covered_refs)), :]) * mapping(:transcript_len => log10 => "log10(transcript length)", :sig_perc => "% significant pairs") * visual(Scatter) |> draw
end

# ╔═╡ 4470fa12-f8cb-4f98-9c33-95310b7dacda
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	# df[!, :selected] = df.mean_entropy .+ abs.(df.phi) .> 1.3
	# df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))

	df = innerjoin(df, comod_summary[:, [:reference, :mod_count, :num_pairs, :num_significant]],
				   on = [:reference])

	df = sort(combine(groupby(df, :reference),
		   	  :mod_count => first => :mod_count,
			  :num_pairs => first => :num_pairs,
			  :num_significant => first => :num_significant,
		      [:mean_entropy, :phi] => ((es, ps) -> mean(es .* abs.(ps))) => :et_complexity), :et_complexity, rev = true)

	df[!, :selected] .= false
	for i in 3:5:180
		subset = (1:nrow(df))[df.mod_count .>= i .&& df.mod_count .< i + 10]
		top = first(subset, Int(round(length(subset)/10; digits = 0)))
		# top = sortperm(subset.et_complexity, rev = true)[1:min(10, nrow(subset))]
		# println(top)
		for k in top
			df[k, :selected] = true
		end
	end

	df[!, :gene] = map(t -> t[6], split.(df.reference, "|"))
	map(r -> println(r.gene, "\t", r.max_interconnectedness),
		eachrow(combine(groupby(df[df.mod_count .> 2, :], :gene),
						:et_complexity => maximum => :max_interconnectedness)))
	# map(println, unique(df[df.selected, :gene]))

	# println(df[occursin.("CDK16", df.reference), :reference])
	
	# map(println, unique(df[df.selected, :gene]))
	data(df[in.(df.reference, Ref(well_covered_refs)), :]) *
		mapping(:mod_count,
			(:num_pairs, :num_significant) => ((n, s) -> s/n) => "% of pairs is significant", color = :et_complexity) *
		visual(Scatter; markersize = 5) |> draw
end

# ╔═╡ 13765151-447d-41c8-afd1-8150d600deae
comod_summary[occursin.("CDK16", comod_summary.reference), :]

# ╔═╡ 2ec7f031-3e5d-4f1e-92ce-1a8f9f5387a8
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :entropy1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :entropy1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :entropy2] = mods.entropy1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :entropy2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[!, :mean_entropy] = (df.entropy1 .+ df.entropy2) ./ 2
	df[!, :relevance] = df.mean_entropy
	
	data(df[in.(df.reference, Ref(well_covered_refs)), :]) *
	mapping((:pos1, :pos2) => ((p1, p2) -> log10(abs(p1 - p2))) => "Distance (log-scale)",
			:phi => "← Mutually exclusive    φ    Co-occurring →         ",
		    color = :relevance => "Mean entropy") *
	visual(Scatter; markersize = 4) |>
	draw(; axis = (; title = "Effect size of modification pair association by distance",
				   xticks = (1:0.1:4.2, ["10nt", repeat([""], 9)..., "100nt", repeat([""], 9)..., "1,000nt", repeat([""], 9)..., "10,000nt", "", ""]),
				   limits = ((0.9, 4.2), (-1, 1))))
end

# ╔═╡ 3e03a99c-f2fa-4dd8-a10a-46a3a886a393
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :ratio1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :ratio1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :ratio2] = mods.ratio1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :ratio2]],
				   on = [:reference => :ref_id, :pos2 => :pos])
	data(df) * mapping(:ratio1, :ratio2,
					   markersize = :pvalue => (p -> -log10(p)),
					   color = :pvalue => (p -> -log10(p))) * visual(Scatter) |> draw
end

# ╔═╡ a84865e1-6f00-4fa2-89ce-254492b7abd6
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :ratio1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :ratio1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :ratio2] = mods.ratio1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :ratio2]],
				   on = [:reference => :ref_id, :pos2 => :pos])
	data(df) * mapping(:ratio1, :ratio2,
					   color = :phi) * visual(Scatter; alpha = 0.7, markersize = 5) |> draw(scales(Color = (; colormap = :diverging_bwr_20_95_c54_n256)))
end

# ╔═╡ 7abde529-0206-4cff-b1b6-cfa23b2a966b
begin
	local mods = copy(sig_peaks_ivt)
	mods[!, :ratio1] = entropy.(eachrow(mods))
	local df = innerjoin(an_sig_comods,
						 mods[:, [:ref_id, :pos, :ratio1]],
						 on = [:reference => :ref_id, :pos1 => :pos])
	mods[!, :ratio2] = mods.ratio1
	df = innerjoin(df,
				   mods[:, [:ref_id, :pos, :ratio2]],
				   on = [:reference => :ref_id, :pos2 => :pos])

	df[df.ratio1 .> 0.8 .&& df.ratio2 .> 0.8 .&& df.phi .< 0, :]
end

# ╔═╡ c568ef3d-a8f3-49dc-8bd3-76a5233f5442
begin
	local ref = "ENST00000228306.8|ENSG00000089157.16|OTTHUMG00000169317.10|-|RPLP0-201|RPLP0|1257|protein_coding|"
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [1241+1-4, 1251+1-4]] .>= 0.75))
	dropmissing!(df)
	local observed1 = contingency(df.Column1, df.Column2) .+ 1
	local res1 = HypothesisTests.ChisqTest(observed1)

	local ref = "ENST00000228306.8|ENSG00000089157.16|OTTHUMG00000169317.10|-|RPLP0-201|RPLP0|1257|protein_coding|"
	local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)
	# local f = Figure(size = (600, 600))
	

	local df = DataFrame(Tables.table(probs[:, [1241+1-4, 1251+1-4]] .>= 0.75))
	dropmissing!(df)
	local observed2 = contingency(df.Column1, df.Column2) .+ 1
	local res2 = HypothesisTests.ChisqTest(observed2)


	local f = Figure(size = (1200, 300))
	local grid = f[1, 1:4] = GridLayout()

	local matrices = [(Int.(round.(res1.expected; digits = 0)), observed1),
					  (Int.(round.(res2.expected; digits = 0)), observed2)]
	local titles = ["RPLP0-201", "RPLP0-201"]
	for (i, (expected, observed)) in enumerate(matrices)
		expected = expected[2:-1:1, :]
		observed = observed[2:-1:1, :]
		
	    # Create an axis for this heatmap
	    ax = Axis(grid[1, i],
				  title = titles[i],
				  aspect = 1,
				  xticks = (1:2, ["Unmodified", "Modified"]),
				  yticks = (1:2, ["Unmodified", "Modified"]),
				  yticklabelrotation = π/2)

		println(transpose(log2.(observed ./ expected)))
	    # Plot the heatmap
	    heatmap!(ax, transpose(log2.(observed ./ expected)),
	             colormap = :diverging_bwr_40_95_c42_n256,
				 colorrange = (-2.5, 2.5))
	             # colorrange = (0, max_count))
	
	    # Add row and column labels
		if i < 3
			ax.xlabel = "Pos. 5089 (m⁶A)"
			ax.ylabel = "Pos. 5139 (m⁵C)"
		else
			ax.xlabel = "Pos. 7211 (m⁶A)"
			ax.ylabel = "Pos. 7261 (m⁵C)"
		end
	
	    # Annotate cells with values (optional)
	    for j in 1:2, k in 1:2
	        text!(ax, string(observed[j, k]) * "/" * string(expected[j, k]),
	              position = (k, j),
				  fontsize = 22,
	              align = (:center, :center),
	              color = :black)
	    end
	end
	
	f
end

# ╔═╡ 21db476e-e901-478f-b5d6-d1c75a0017e3
begin
	ivt
end

# ╔═╡ 0cdbf519-0c60-48f9-b32a-37d6f3ac667c
begin
	local f = Figure(size=(900, 400))
	
	local df = innerjoin(ivt[ivt.predicted .> 2, :], glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.mod_ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.mod_ratio .- df.glori_mean_ratio)))
	
	local plt1 = (data(df) * mapping(:glori_mean_ratio, :mod_ratio => "IVT ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash))
	draw!(f[1, 1], plt1, axis = (; title = "Nanocompore (WT/IVT) vs. GLORI", xlabel = "GLORI mod. ratio", ylabel = "Nanocompore mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

	
	local df = innerjoin(wtmk, ivt[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
		  on=[:ref_id, :pos])
	# df = df[df.mod .== "a", :]
	println(nrow(df))
	df = innerjoin(df, glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	# df = innerjoin(df, ivt[ivt.predicted .> 2, :], on = [:ref_id, :pos, :chr, :strand, :genomicPos])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.ratio .- df.glori_mean_ratio)))
	
	local plt2 = (data(df) * mapping(:glori_mean_ratio, :ratio => "Modkit ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash))
	draw!(f[1, 2], plt2, axis = (; title = "Modkit (WT) vs GLORI", xlabel = "GLORI mod. ratio", ylabel = "Modkit mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

	for (label, layout) in zip(["A", "B"],
							   [f[1, 1], f[1, 2]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 16,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
end

# ╔═╡ 04db972c-b08b-44c1-be9d-388484b133e5
begin
	local f = Figure(size=(500, 500))
	
	# local df = innerjoin(ivt[ivt.predicted .> 2, :], glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	# df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2

	local glori_ratios = Dict(map(r -> (r.Chr, r.Strand, r.Sites) => (r.Ratio + r.Ratio_1)/2,
								  eachrow(glori_raw)))

	local df = copy(sig_peaks_ivt)
	df[!, :glori_mod_ratio] = map(
		r -> let ratios = [get(glori_ratios, (r.chr, r.strand, r.genomicPos + offset), missing)
			  			   for offset in [0, 1, -1, 2, -2, 3, -3, 4, -4]]
			existing = collect(skipmissing(ratios))
			length(existing) > 0 ? first(existing) : 0.0
		end,
		eachrow(df))

	df = df[df.glori_mod_ratio .> 0, :]
	
	local c = round(cor(df.glori_mod_ratio, df.mod_ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.mod_ratio .- df.glori_mod_ratio)))
	
	local plt1 = (data(df) * mapping(:glori_mod_ratio, :mod_ratio => "IVT ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash, color = :red))
	draw!(f[1, 1], plt1, axis = (; title = "Nanocompore (WT/IVT) vs. GLORI", aspect = 1, xlabel = "GLORI mod. ratio", ylabel = "Nanocompore mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))


	f
end

# ╔═╡ 372dcaef-18a1-4a1d-b292-6d2a285e4d83
begin
	local f = Figure(size=(500, 500))
	
	# local df = innerjoin(ivt[ivt.predicted .> 2, :], glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	# df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2

	local glori_ratios = Dict(map(r -> (r.Chr, r.Strand, r.Sites) => (r.Ratio + r.Ratio_1)/2,
								  eachrow(glori_raw)))

	local df = copy(sig_peaks_ivt)
	# df[!, :glori_mod_ratio] = 
	local closest_m6A = map(
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
	
	local c = round(cor(df.glori_mod_ratio, df.mod_ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.mod_ratio .- df.glori_mod_ratio)))
	
	local plt1 = data(df) * mapping(:glori_mod_ratio, :mod_ratio => "IVT ratio", color = :distance => "Distance to the GLORI site") * visual(Scatter, markersize = 5, alpha = 0.5) |> 
 # data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash, color = :red))
	# draw!(f[1, 1], plt1, 
	draw(; axis=(; title = "Nanocompore (WT/IVT) vs. GLORI", aspect = 1, xlabel = "GLORI mod. ratio", ylabel = "Nanocomporeco mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1), figure=(size=(550, 500),))#, colorbar=(size=5,))
	# resize_to_layout!(plt1)
	


	# colsize!(plt1.layout, 2, 0.1)
	# save("/home/mzdravkov/test.png", plt1)
	# plt1
end

# ╔═╡ 5bf138b3-7907-4fc6-96e7-d3a8bccc78af
begin
	local f = Figure(size=(500, 500))
	
	# local df = innerjoin(ivt[ivt.predicted .> 2, :], glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	# df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2

	local glori_ratios = Dict(map(r -> (r.Chr, r.Strand, r.Sites) => (r.Ratio + r.Ratio_1)/2,
								  eachrow(glori_raw)))

	# local df = copy(an_sig_peaks_ivt[in.(an_sig_peaks_ivt.mod, Ref(Set(["m6A", "m5C", "m7G", "Y", "Nm", missing]))), :])
	local df = copy(an_sig_peaks_ivt[in.(an_sig_peaks_ivt.mod, Ref(Set(["m6A"]))), :])
	# df[!, :glori_mod_ratio] = 
	local closest_m6A = map(
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
	
	local c = round(cor(df.glori_mod_ratio, df.mod_ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.mod_ratio .- df.glori_mod_ratio)))
	
	local plt1 = data(df) * mapping(:glori_mod_ratio, :mod_ratio => "IVT ratio", color = :category => "Mode of inferring") * visual(Scatter, markersize = 6, alpha = 0.9) |> 
 # data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash, color = :red))
	# draw!(f[1, 1], plt1, 
		draw(; axis = (; title = "Nanocompore (WT/IVT) vs. GLORI", aspect = 1, xlabel = "GLORI mod. ratio", ylabel = "Nanocompore mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))


	# f
end

# ╔═╡ bc7932b4-80e1-4f05-981e-0c869c0dce85
glori_raw

# ╔═╡ d45c9d47-68f2-44d7-82ad-ddcb4877263f
ivt

# ╔═╡ a7eec8b0-9a15-42c0-aa0c-95675a21dcf7
wtmk

# ╔═╡ faa56f3f-6891-423c-89e3-95099dcc149c
# ╠═╡ disabled = true
#=╠═╡
begin
	local refs = unique(sig_peaks_ivt.ref_id)

	local results = []
	local rlock = ReentrantLock();

	Threads.@threads for ref in refs
		local probs = get_read_mods("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_fix_mod_clust_inferring_read_level/out_sampComp_sql.db", ref)

		mods = sig_peaks_ivt[sig_peaks_ivt.ref_id .== ref, :pos] .+ 1 .- 4
		ref_res = []
		for pos in mods
			nvalid = sum(probs[:, pos] .!== missing)
			nmod = sum(skipmissing(probs[:, pos]) .>= 0.75)
			push!(ref_res, (ref, pos - 1 + 4, nvalid, nmod))
		end
		lock(rlock) do
			append!(results, ref_res)
		end
	end
	ivt_mods_75perc = DataFrame(ref = map(r -> r[1], results),
								pos = map(r -> r[2], results),
							    nvalid = map(r -> r[3], results),
							    nmod = map(r -> r[4], results))
	ivt_mods_75perc[!, :ratio] = ivt_mods_75perc.nmod ./ ivt_mods_75perc.nvalid
	ivt_mods_75perc
end
  ╠═╡ =#

# ╔═╡ 5efd1815-277b-4444-803d-b48336b0abd7
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure(size=(900, 400))


	local df = innerjoin(ivt_mods_75perc, ivt[:, [:ref_id, :pos, :chr, :strand, :genomicPos]], on = [:ref => :ref_id, :pos])
	df = innerjoin(df, glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.ratio .- df.glori_mean_ratio)))
	
	local plt1 = (data(df) * mapping(:glori_mean_ratio, :ratio => "IVT ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash))
	draw!(f[1, 1], plt1, axis = (; title = "Nanocompore (WT/IVT) vs. GLORI", xlabel = "GLORI mod. ratio", ylabel = "Nanocompore mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

	
	local df = innerjoin(wtmk, ivt[:, [:ref_id, :pos, :chr, :strand, :genomicPos]],
		  on=[:ref_id, :pos])
	# df = df[df.mod .== "a", :]
	println(nrow(df))
	df = innerjoin(df, glori_raw, on=[:chr => :Chr, :strand => :Strand, :genomicPos => :Sites])
	# df = innerjoin(df, ivt[ivt.predicted .> 2, :], on = [:ref_id, :pos, :chr, :strand, :genomicPos])
	df[:, :glori_mean_ratio] = (df.Ratio .+ df.Ratio_1)./2
	
	local c = round(cor(df.glori_mean_ratio, df.ratio), digits = 2)
	println(nrow(df))
	println(c)
	println(mean(abs.(df.ratio .- df.glori_mean_ratio)))
	
	local plt2 = (data(df) * mapping(:glori_mean_ratio, :ratio => "Modkit ratio") * visual(Scatter, markersize = 5, alpha = 0.5) +
 data(DataFrame(x = [0, 1], y = [0, 1])) * mapping(:x, :y) * visual(Lines; linestyle = :dash))
	draw!(f[1, 2], plt2, axis = (; title = "Modkit (WT) vs GLORI", xlabel = "GLORI mod. ratio", ylabel = "Modkit mod. ratio", xticks = 0:0.1:1, yticks = 0:0.1:1))

	for (label, layout) in zip(["A", "B"],
							   [f[1, 1], f[1, 2]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 16,
	        font = :bold,
	        padding = (0, 5, 5, 0),
	        halign = :right)
	end

	f
end
  ╠═╡ =#

# ╔═╡ da76c242-e945-4336-9aa9-f2aa267fb975
begin
	local df = copy(an_sig_comods)
	df[!, :gene] = [r[6] for r in split.(df.reference, "|")]
	df[!, :transcript] = [r[5] for r in split.(df.reference, "|")]
	df[!, :distance] = abs.(df.pos1 .- df.pos2)
	df = sort(df[:, [:gene, :transcript, :pos1, :pos2, :distance, :mod1, :mod2, :pvalue, :qvalue, :phi, :chr, :strand, :genomicPos1, :genomicPos2]], [:gene, :pos1, :transcript, :pos2])
	df |> CSV.write("/home/mzdravkov/comod_table.tsv", delim = '\t')
	df
end

# ╔═╡ 7ef78147-f3a9-4a73-a179-6d9b302629dd
md"""
## Percentage of detected mods with an association
"""

# ╔═╡ 7e5f1c7b-1957-453c-8aaa-90719e9ff215
vcat(map(tuple, sig_comods.reference, sig_comods.pos1),
	 map(tuple, sig_comods.reference, sig_comods.pos2)) |> Set |> length

# ╔═╡ a7269b72-6398-4016-ab2d-28cb6a024431
14129/nrow(sig_peaks_ivt)

# ╔═╡ fff284bb-8469-4441-ba77-f20bb296fc79
10630 / nrow(sig_peaks_ivt)

# ╔═╡ a03581d9-edf1-4932-9029-5336baa4b1d6
begin
	local refs = sig_peaks_ivt.ref_id |> unique
	local min_dist = 1000
	for ref in refs
		positions = sig_peaks_ivt[sig_peaks_ivt.ref_id .== ref, :pos]
		for (i, pos1) in enumerate(positions)
			for j in i+1:length(positions)
				pos2 = positions[j]
				min_dist = min(min_dist, abs(pos1 - pos2))
			end
		end
	end
	min_dist
end

# ╔═╡ ede4a4d6-206b-4f24-9485-417c2fa982cc
function dunn(x, y)
	nx = length(x)
	ny = length(y)
	N = nx + ny
	ranks = tiedrank(vcat(x, y))
	xranks = ranks[1:nx]
	yranks = ranks[nx+1:end]
	(sum(xranks) - sum(yranks)) / ((N*(N+1)/12)*(1/nx + 1/ny))
end

# ╔═╡ 6050ae2f-92c9-427e-8058-3a61d8f4e8a3
ivt[ivt.genomicPos .== 76066564 .|| ivt.genomicPos .== 76066573, :]

# ╔═╡ 479fb2d0-bae4-4a3f-b0a7-03f94ae442cb
comods

# ╔═╡ 1007c9b7-2aaa-4cdc-9cfe-4efc14642ef7
sum(sig_comods.pos2 .- sig_comods.pos1 .> 1000) / nrow(sig_comods)

# ╔═╡ d8795ad6-415d-46e0-8abb-e9c8f7bf926d
sum(sig_comods.pos2 .- sig_comods.pos1 .< 40) / nrow(sig_comods)

# ╔═╡ 2ee422ad-994f-42ab-9127-f77f157bebbd
nrow(sig_comods)*0.1

# ╔═╡ 8eb74e1d-a6e7-4b13-85c2-8b29ca971496
sort(sig_comods.pos2 .- sig_comods.pos1)

# ╔═╡ 3f622f6d-9828-4742-b040-e3bfd4749479
sort(an_sig_comods[occursin.("ENST00000234875.9", an_sig_comods.reference), :], :phi, by=abs, rev=true)

# ╔═╡ Cell order:
# ╠═41818e42-3cb8-4d48-93d9-b53e0eea7873
# ╠═164f8adc-3bba-11f0-3c64-19ee3bf9097e
# ╠═93db1a5d-4de2-482b-a4e2-540db2059689
# ╠═6e8e2f5d-d998-4a4c-a14f-dc63e247b8a7
# ╠═e6362147-b298-4b52-a2bd-b38dee879b49
# ╠═bc77ac69-225a-4f81-b219-15683bb5eb4a
# ╠═82dba028-3fbe-44e2-9e8f-377f738363a7
# ╠═c88509a9-5a1d-4895-9ae7-5174cc199772
# ╠═dedfd9b6-9b9e-4149-834a-547cbac99873
# ╠═1795bc42-bdc8-4685-bb52-39271581bc06
# ╠═74f8fb7a-e501-490a-a530-224f51d6312c
# ╠═7a783418-ca3a-4ac7-a8a1-8968a5bc9096
# ╠═88dac2fb-7062-4e75-aaaf-c2596d90b299
# ╠═14c823fa-6412-4f38-88c2-66dfa4bc4842
# ╠═5a2ad3e8-9f14-49ad-a78b-6fdd4be05220
# ╠═8445b562-cb03-4a7e-911b-5dc44778ee5c
# ╠═8b930a4c-cef4-49e6-8a76-d59d4df8a349
# ╠═a6d02ee4-2e77-4520-b140-36fb2150b1b0
# ╠═b89d56f0-5104-4b0d-a1ae-6090d262db25
# ╠═566437e3-1abe-49bf-a262-e9a4c702432e
# ╠═b5d8eafa-9ea4-4c59-bef8-3070126bb095
# ╠═c9fabd63-2131-41fc-a50a-bda27319bc6d
# ╠═43a5d6b6-652d-403c-b4fd-945ba2b80ead
# ╠═5abecd8c-e06c-477c-8938-545f64591385
# ╠═b0a1b932-7f0b-40ac-a2d6-e4867d2b54da
# ╠═7ca021f1-532f-4a3c-98b6-0ff08c27d01e
# ╠═33a21a05-dac1-40e0-97a4-8dae37fa8bfc
# ╠═9232e724-d8df-4970-931c-f1c11f0ce859
# ╠═82a9954d-2c23-46bc-a0bb-c2fc078eabb1
# ╠═16109b94-06e5-4285-9333-5b8a98761ef7
# ╠═57a4e1e9-59cb-4180-b735-9da13feb64a6
# ╠═fbd57196-5c5b-49a9-a5f5-a39d25a33c87
# ╠═d41911f6-ab79-4a7f-b245-76ef3cfcd1e6
# ╠═5aaee16e-1121-4bf7-b4c6-60941b62414c
# ╠═4016b836-3669-4098-9725-d3708755db40
# ╠═9539bce6-7042-4396-b4f0-1badf268b068
# ╠═f3dba5a1-0352-4334-ab79-89b73172a289
# ╠═5e1553f2-370d-4529-b9af-16764a2018d5
# ╠═390466c1-9b31-460b-82cd-676d85ca41da
# ╠═9f605ad2-6d05-49d5-a743-2dbf960b1a51
# ╠═4cb861eb-f76b-4871-a783-554c85642e42
# ╠═fe4155bb-d187-4113-9e85-dfc3fedee4ad
# ╠═081adcae-8bfe-4c6f-9a17-1ddd8e9dd87e
# ╠═3a91b67e-8efa-4e84-a519-6426594ea30e
# ╠═a4d697a0-71b6-4e01-879e-6a91d9060b8d
# ╠═23ce8398-5f79-47d4-afe9-c46818011b9a
# ╠═25f512fd-e5f8-40d7-a7eb-2d04e18b884b
# ╠═6669794a-cebf-46f7-8d1b-5bb91ccb7280
# ╠═c71344c3-4faa-4849-a35a-6f511eaf4f34
# ╠═16fff5ab-601e-42d7-9f74-1670690733bb
# ╠═86179fb3-d3d0-4f5b-a8d3-14fb43dac834
# ╠═5fde98ae-bdf2-4dcb-90a0-f5017b97262a
# ╠═d90664a5-8e7d-469b-8b05-b103d829f672
# ╠═ab03ca5a-3186-493e-9e3e-13367a85e624
# ╠═cefdac37-7b63-4c8c-ad92-41ca18b78452
# ╠═140e4f5f-3c19-4c94-8829-4a46a95eb6cd
# ╠═9513c89b-dd9e-40dd-91be-c480438e3f07
# ╠═a3f63803-5803-4942-8dce-b1b90b86f6a1
# ╠═56ed77ee-fcfb-498a-ae5c-f9e77ceb971a
# ╠═e429bac0-0add-4e15-8a7a-6ab22a5c9943
# ╠═2b618c14-b478-4d5b-abfc-19bbeb073fde
# ╠═ea86f47c-ed99-4108-bb5c-c9ee225dc7ef
# ╠═7fcfeb9f-5e35-4408-8cfb-b24c5f52c6a2
# ╠═eb2f9c7c-b5bb-493d-8cc1-abdc452e646d
# ╠═ed261b2a-459a-440e-99ad-91afdeda2952
# ╠═dc85be9a-f0af-4b8b-924e-758b3ce45a7a
# ╠═bfe35e43-60fd-4ba9-9f05-f5a7f868e011
# ╠═a81dc14a-fd54-417b-bcd7-84615a62aa37
# ╠═a2dbec3c-961f-4904-961d-be2d9d538720
# ╠═7c1bce27-3c08-4b0a-9364-9e992a793a97
# ╠═bf514f88-b2f4-4a35-8c4a-08fec382e74c
# ╠═c1002495-1c24-450c-9585-1aeb7d588bb5
# ╠═823b1914-ab79-499a-98ec-348f24b0107a
# ╠═80836ca5-66c5-4f9c-a815-9c57230f6948
# ╠═66983a7b-8a6a-4b9f-91ba-bb40b5b5ebce
# ╠═237e9e39-95c8-40a1-b8d3-306a7b98ce0f
# ╠═cf152b80-2eea-4ed4-8d7f-d0d30a1aec2f
# ╠═6c7f7c71-ed77-46b1-97a7-0d8b9c8d758e
# ╠═29c659f1-4e5a-4c65-aa53-33dea7d730f3
# ╠═d3fc7047-56a6-4dce-96d3-61189556299e
# ╠═1d820491-8699-4324-b113-b4d7ac0a819a
# ╠═da1c3021-eccf-40bd-8988-4a0d068eb0ba
# ╠═bc42854a-ecbe-46c2-b289-114b12cf6fe9
# ╠═3de4cb79-b2c6-4c12-9e8c-8a78a3af4457
# ╠═be74ab35-49fe-401e-9829-6dbc8d90f37f
# ╠═97711651-3a07-4c6b-b749-572dddb4a145
# ╠═d8feca5f-d4e8-4e78-8d87-4f0f7a4755c4
# ╠═ca2f63ea-677c-42e7-a736-b12a280cad22
# ╠═f28e2414-efe5-4a73-91e3-79ed437d34e5
# ╠═cf8c009c-a109-4121-8e39-b3e96cbd4092
# ╠═a7fe57e3-cfcc-4973-810e-98edacaf6fdb
# ╠═39ad0f2b-299e-4b35-9ecf-8cc83b55cd3f
# ╠═e113b1dc-a755-438a-8c08-6747cde5afde
# ╠═414f6a39-e6c8-4dbe-9f35-41f83a065040
# ╠═98ecad14-b670-4ba6-bb55-8adbef018a95
# ╠═435b9d3a-1eb7-4755-8778-05ce6f7def8c
# ╠═1aa2a4e9-fb40-4573-982b-cb6a0a8c2b59
# ╠═f203a977-1521-460d-969d-796e5e4fecd2
# ╠═2825eb3b-12a2-44ec-9146-06d95381dcbe
# ╠═87144656-e38e-4e91-aa08-33a9c7b790e8
# ╠═f064b8c4-d0ad-4e00-8296-4bbed1316256
# ╠═134194d3-c29f-4eb3-878d-a3c24892dd6c
# ╠═d2c86575-c50f-4069-a48b-cced8084a84b
# ╠═70c7949f-915c-4ded-96fd-bf158d9a8109
# ╠═02599f81-c0c5-453a-a278-4dfa7104a001
# ╠═5542defd-747d-462e-abf2-6894ddb54081
# ╠═ac85532f-1b1d-4571-a187-e4304882ba27
# ╠═47ff4422-b91d-4748-9937-b0ddc8e9ecd5
# ╠═9fc1f07c-48e1-420b-9335-4bfa031e86cf
# ╠═fc0494df-e2fa-4a36-a0e0-dc6de09350cb
# ╠═fd8b9aeb-d76b-431d-b75f-ef86bfa4c7de
# ╠═89aa206f-14bd-485d-ae1c-d7b6cabdfee6
# ╠═5c70bdac-40a0-4378-9d00-1bc8966da523
# ╠═694914ad-e7d0-4cbf-a9db-55fb49e393ac
# ╠═968ac33d-b86d-46cc-820f-c951f61a8e24
# ╠═3b5d190e-5b56-4b7c-b389-68bbbb81f836
# ╠═a617f69a-669d-499d-a40c-ef8f8b747e49
# ╠═dd97de04-b552-4d11-a9a4-f46e0f5d670c
# ╠═8ad73f0f-7e4f-4848-ae37-4e5780a1b0da
# ╠═31ef15f2-8cf6-4711-82b8-59f34d2e043e
# ╠═82a98275-f6d6-4ed1-9aa7-dd8ab59d9370
# ╠═ec42c786-4e6d-41e0-a40b-8694f4f99590
# ╠═61e39473-08e6-45fb-bfd2-545001a2ecef
# ╠═dd6d235f-b102-46b0-9ea7-71857a4db0ba
# ╠═1513ad3e-b7c9-4a14-8374-27d79eb6630e
# ╠═40eaa862-f1ee-4fc1-ad89-5897f2938f0d
# ╠═3e577c30-ae50-4deb-a361-83912f0904fc
# ╠═619bda6e-9209-4634-88bb-a3fbff14f74b
# ╠═301d30fe-4485-47a7-97e5-d093d2b8feff
# ╠═73ddf0c4-5961-4c81-b5d6-833cd0c06f0d
# ╠═e3665e85-ba56-4929-bd10-20c949e951e7
# ╠═79736277-c65c-4ec5-94b0-2268b38bc095
# ╠═f58d900d-a5f3-40eb-8533-95542d92de30
# ╠═a5a081d7-2ea0-4f7f-afd8-4e2f70ff3345
# ╠═c90a46c3-e3ca-4088-83af-1a1368559f20
# ╠═38cc2f30-0ece-4dd7-ae64-c5b3a28b99fd
# ╠═8c658c04-0948-4d33-83fe-5d9756f43552
# ╠═dd730961-b1f9-4b2b-af69-84b2ac751117
# ╠═5ed770d3-da06-4f86-b8fc-3425ddc7bdfc
# ╠═5be0dc53-9039-4324-b1e2-e24c344f069c
# ╠═46068dad-f73f-4660-9c09-2fbc85be0a09
# ╠═6ba634cb-e770-44ba-a044-578029cfd8ec
# ╠═d4443059-dade-4a9a-aeef-75408913ce16
# ╠═10c9acd5-47d0-4f2a-a619-44fd903c4b36
# ╠═26d0dc17-cd92-47d4-bec0-33ad54c002cd
# ╠═c2cf87ff-d11c-4bb3-b44d-6a16f4040c88
# ╠═fd549de3-8c1f-4358-8ec0-48fdd85d0bd1
# ╠═cbe49401-5a36-4f40-9ece-4cccfc89f877
# ╠═b37f9ca0-5a53-4716-9589-9b57268f72b2
# ╠═67638f93-2a44-4772-9c35-2e728cebc277
# ╠═c99dbb5d-77da-4c2e-9b7a-874f5ead6574
# ╠═409ec9ef-22c6-4192-a53c-382d059ff299
# ╠═664f5ea5-378a-45cc-a6cf-cab779b7f933
# ╠═d1c2e993-4ba4-46a5-840c-38bfbf397858
# ╠═cb81e3ba-5eea-4169-a4ff-5b7433e8b171
# ╠═07f53334-5428-4890-b02a-9e3c5de7b22c
# ╠═40bbff56-4264-4d9a-af0e-59485c923f92
# ╠═0ee4e9af-4b1e-41fb-8746-12ff51d20e27
# ╠═9503a9e1-8c7e-4dfe-aeed-015d7ec1095c
# ╠═9fbe17fb-5b2a-4cda-a70d-1598de0a6ab0
# ╠═58e326c7-e5b5-4f97-884e-7321a8084139
# ╠═e30b68c9-9aaa-44f7-86d9-c8e1ab9254dd
# ╠═3dc430dc-a218-42bf-90da-8ba8ba6dd56a
# ╠═46b1d7d3-2e9e-4d19-91d8-07150e981ef0
# ╠═809b0c5c-1538-4bbe-bb4d-293f6058b2a6
# ╠═ab71927a-b85c-4834-a20a-cebe1f734b7b
# ╠═3a9274bf-3af7-46d0-880e-1ea0e8b421de
# ╠═b5423ce9-cad4-48ec-88ec-33d12319d478
# ╠═5229ac7e-5a3c-4339-8c06-ce4e07f379c3
# ╠═25e2f5f3-0612-48ea-83fc-3c4e85eb793d
# ╠═beb418fd-d39b-4c5e-bf9d-98e936b9c819
# ╠═de7c2cc7-7931-4290-97be-39b24634d755
# ╠═088e365f-05f2-429c-9fa2-763ddf30b82b
# ╠═ea219179-f224-4aa3-bb1d-fdf494742363
# ╠═bcc0fec3-9925-4c52-8411-a59583f9171e
# ╠═b117d0ed-b44b-46a6-b9ea-204e1d5195e9
# ╠═5408e59e-53ae-4ef9-b9ae-e86be3119575
# ╠═72faa038-d71c-46f9-8110-d68b7633f2e6
# ╠═9b4cf8f2-84ae-4aba-981a-7ef0e669a453
# ╠═d7fb857c-0871-4be2-8322-e36ab95e6f19
# ╠═05ff34ad-5cfc-48a4-bfc9-4d4c7cfda4d1
# ╠═171a33c2-e532-4167-82a6-27d830402104
# ╠═657e7156-57ef-4f64-a4ba-13af4bb2784e
# ╠═89114c37-3d83-4e4c-a46f-32cd04e85bf0
# ╠═d6b27bbe-824e-492b-9aea-dd5041f4917f
# ╠═e0bfd78c-25fa-4752-a4d7-660aa4b6cd66
# ╠═b6a72db1-ec4d-4cc1-8174-c6179316dec7
# ╠═f2829ecb-cbdf-4e0d-953c-0072f03546ac
# ╠═74279cfc-ab7a-4e9a-b048-ad17a378c20f
# ╠═5b5358c8-f35c-473f-9e92-56bdb4dbc6a7
# ╠═9031f497-26d2-4b3b-ba34-3438df409fd9
# ╠═ad11f743-1eee-4298-87bf-620a0f1cd1ec
# ╠═f793c68f-731e-426f-ad6e-ef5419f16370
# ╠═cc022b63-4a04-4142-976c-02756047076d
# ╠═a6b2a25c-6fb2-4dd0-a7f1-ea16ceb63a97
# ╠═69798f9b-f9a8-41f7-a2ff-8f333a6ed728
# ╠═741a21f9-ea94-4620-83a6-a43dd5fc300d
# ╠═e3f624cc-1c45-461a-ac19-bed551b1ba4e
# ╠═5e8c4ba1-9f8d-40df-b02e-cfeec4bf4d57
# ╠═e8cee36d-6886-4e27-806e-b4df1dadff7d
# ╠═59a166ff-402f-44cd-91d4-1f58ca676e33
# ╠═25c85ec0-af9e-4f74-bb9f-df7325604100
# ╠═4cea1683-7f74-42dd-a0d5-5d66eca41253
# ╠═58a03eb1-4847-4f48-a86e-9da1233023a3
# ╠═56b54a7c-fb57-4aa0-9e84-45ade0d06b5c
# ╠═52a1db73-e957-41fd-a6dc-f483b874e9f2
# ╠═fbd6502e-efa4-4e52-81db-3bb998314b14
# ╠═bf177ece-0e74-45f2-8bb0-2b5a2ac2ddbe
# ╠═8d2aafde-d745-4a0d-bcdb-d82a7877292f
# ╠═3d10cafd-f87c-40a7-ba15-f64e924cc686
# ╠═28779e98-1a0e-481f-bb10-c6973bd98976
# ╠═40173831-c379-4e62-90fc-f7f2e55ff8c6
# ╠═95a3d0c3-09c8-490d-ae5c-993d8ec009e8
# ╠═7bede161-f145-44f4-b6de-9fa00a72be52
# ╠═a927ecd1-b616-4479-88b3-a070e1cc95b8
# ╠═25c51988-445b-4253-9763-22e1dd8ee210
# ╠═4444b683-1571-4df8-9841-0e28df365dec
# ╠═757fa776-c741-4eda-9e85-c1110231fcd4
# ╠═b8437ca7-b885-4dad-a72d-dc0b9c9a6be5
# ╠═90c789e2-999c-4e6c-8d49-0926ca080e71
# ╠═f1795d88-1e1f-4aa7-b245-b2f2cef7db17
# ╠═967bcd01-03a6-4cc5-ac07-b997b65d8a8f
# ╠═322240e8-a405-4ccf-8d94-3c52aa60eb26
# ╠═e698c85d-ab27-4dda-9349-cf679b22a0f5
# ╠═4415b732-d854-4eb9-9ae6-ed872b087c3e
# ╠═8dd89e2e-78e3-4d2a-b2f8-c9d29e949752
# ╠═c817df72-596c-41be-9d55-80dcbe73c485
# ╠═76fbfe22-1223-454d-8f6c-0505cdd6c8c4
# ╠═e33af351-2a7f-45fc-be7e-b36467b3e324
# ╠═b2fabd16-f406-4b70-b5ab-343ff87b6a05
# ╠═7b353ba1-cc33-405b-962e-5b7f82d1323b
# ╠═cb6f8399-addd-4095-98bf-915f1b38ee4d
# ╠═3cbff0f7-eb89-49e5-88e4-43ec57ab87b3
# ╠═46ce6ac5-c128-46df-8d38-90776f42861b
# ╠═ca061644-6cf0-4547-bb0d-c21642f76f3f
# ╠═fa354baf-9fcc-4022-828a-3458e4cdae0c
# ╠═68870781-fd02-4df4-b531-1ae27f53c0fd
# ╠═fbea4e5f-4c3a-4492-97cb-d2eeeb4d1c63
# ╠═be8d5828-880e-41dd-9195-95908bd3e829
# ╠═d888fa77-b67a-4654-8556-2c3dc2cd9f44
# ╠═8def7cc6-83ec-406c-b3d8-24c138d78199
# ╠═167ea40a-7ab7-42a3-af60-8f2b9367bf5f
# ╠═b747c8a5-6cb8-43ae-8c77-1c41260dbe0f
# ╠═451907e7-eb22-4d8f-aa4e-44d23708861c
# ╠═dd689831-a388-44df-9eed-97add91d37bf
# ╠═d636baaf-b65f-4c5e-a9ed-775a7f1c9909
# ╠═134e34a4-e467-49e6-9c39-12123401d1b7
# ╠═19e82456-11a8-441c-b753-208aba90af10
# ╠═1bf15795-87e8-4cf7-92b1-9db7b102450a
# ╠═d9be422c-3630-4c86-859a-3df1cd725405
# ╠═33cc218c-f02f-4542-ba7a-b6075d7e4357
# ╠═10a0495d-b861-4959-a5c2-1d74b6a0f5c6
# ╠═0e7cfe6c-a700-4ec1-b8dc-a67f246782b3
# ╠═3b6c3b40-d54e-41c2-afac-6bf8234ae648
# ╠═f97950ff-8c2b-4fc4-99e2-51c29bf7b489
# ╠═b9ee5aa2-ced4-4f71-b43a-bbfd6d5afef9
# ╠═bd94a0c6-0e84-4699-b449-36e45b2c2c23
# ╠═b79d6323-c6be-4b18-ae93-1a1ce48a1db3
# ╠═ca38417c-1d5c-45c6-816c-d2d6898af7db
# ╠═692af519-1ac2-483c-86c5-f1801d8426bd
# ╠═76d3d1c4-0e01-4360-bc3d-bc274ff1c39c
# ╠═3bfb00e8-f963-4c03-825f-71c4343074fa
# ╠═1c622f87-b763-4f8b-a253-fd38a1fd2830
# ╠═45a5a8f7-a1bb-434c-b882-3b24d3efa412
# ╠═96482819-8af4-4ce7-bfda-2e20ac05a66c
# ╠═6a70233a-4754-49e8-bd81-292d364b4031
# ╠═445fe109-2361-4bea-bd03-74f3f12dca65
# ╠═f7a4a1b1-5970-4e1c-a37a-5db901c1115c
# ╠═4605495f-6df1-44fb-8a33-6e84bebe83df
# ╠═90818cac-2f00-4fc8-9225-321eeac107e7
# ╠═38f8fbed-20a9-46ab-92e9-edf3cb693bfa
# ╠═f84bfc2a-c037-4627-b0c6-4f47d4f45a1b
# ╠═71261e27-8c1b-47e0-8fdd-68fc1cecee39
# ╠═d836b201-4413-412a-8ec2-be15d23504c1
# ╠═10a15475-c736-48f7-a4a4-cf972134a0f0
# ╠═a0798ff2-c72f-4de5-b31a-3c8a3bacf117
# ╠═007e6893-0718-41c0-8ca0-03fd367f6cc7
# ╠═49252d74-9a99-4795-9527-3e7eccf20eed
# ╠═e9788e02-f17f-4bd8-b344-caa18863e92e
# ╠═9ba8c503-f59e-4d6a-b97c-9826279d744a
# ╠═3266d3d6-64a2-4936-91ac-059a1ccf76c4
# ╠═da65e2fe-54a8-4ca6-b599-411ee2a3e9d4
# ╠═7212801b-ee05-4699-b4fb-c02377710e10
# ╠═aa699be6-0fa9-468d-847c-b7b3b7f91c72
# ╠═6a423afd-67f4-4152-9aa8-7d69c4262f63
# ╠═c02b568d-e607-4729-9cef-1d189ae22bc1
# ╠═7aa12ba0-e6e7-4c36-9620-7ed443013d40
# ╠═a11ab145-e6d2-4452-8df8-6976f5f67453
# ╠═046608c2-5c79-4651-bcfe-1b8b965e63c2
# ╠═e2bf25a2-b5c2-48be-95b3-9ee2905b7b35
# ╠═5a1dfded-7f69-4fae-ad5c-9242686026d8
# ╠═77dbbadb-2588-4db7-8e84-5b375dbec467
# ╠═ce3be394-4d52-40c6-9df2-4cad1ca55fd7
# ╠═90cffe20-514b-4bac-8e6f-04e30bc80f67
# ╠═9f184944-6093-487f-a095-98a8572af620
# ╠═ff570f16-c8fb-4e0d-8851-a0fb93b23645
# ╠═b23610f5-240b-4385-955b-216b3221d1d8
# ╠═119b9407-c1cd-4eff-a0a6-993b9b05f50b
# ╠═f208fdf1-e949-468a-a987-d798baa939ca
# ╠═434ef740-5386-4a64-9ebe-8211e40b0c2b
# ╠═f9dbfebc-9c8e-451a-a46c-03b548b630de
# ╠═46aee8ff-83fc-4e88-95d7-20a1c076d3b2
# ╠═6a44e37f-9925-44a5-afcd-6f3379796e86
# ╠═3fd1b51b-3a5c-4d36-9cf2-1281c0e4483c
# ╠═2357c8ef-6b46-4f17-97da-cdb99668ca85
# ╠═29283d8b-c564-4f0a-a6a2-6d198d320162
# ╠═452d1d7b-afae-4bf7-b0a4-88474eaf947c
# ╠═95b209bc-4fbb-4713-b1b0-8f6b4b1a14b4
# ╠═f652d459-817e-47c2-8456-3de6b5f89f9b
# ╠═17b3496f-d5bd-4782-9ec9-65723d9ced58
# ╠═a8ba8f4e-96ae-49d6-91c4-52d7f839f58a
# ╠═dbef3c64-4989-4bdd-916c-045bbeb03263
# ╠═c2ca3807-116d-4251-b864-7bbf7e63d41f
# ╠═f3b640b8-ab69-4ce0-b597-d667170f9cc2
# ╠═43603d8b-82dd-4670-a999-51802c6173df
# ╠═d5a93508-e559-4767-a078-54f9bb8fb3bf
# ╠═4ad4579a-4b64-460b-9ec4-a80f86ee391a
# ╠═446f9868-3f80-4200-97a3-a351de532be2
# ╠═93e58616-9b19-42f4-84ab-da73c10561e4
# ╠═5948af6a-d5df-46a5-a832-ad33d3feb8c5
# ╠═820fb1ae-e3fe-4aa8-8163-e72c4748f9be
# ╠═5f66681b-fda7-48fe-8771-d0d44564df45
# ╠═75c4e402-c255-49d5-a1c9-df2ffde1a578
# ╠═eed7c85f-5ec5-41ef-9bd9-250dd07137db
# ╠═6628a74b-4f8c-45eb-85f5-16681b2c849e
# ╠═3ef8cb41-d3de-4f38-a93e-2d0211b15263
# ╠═e86d7b59-449e-4e8f-8ad2-1c708ea15076
# ╠═6412e0da-abe4-44fb-9e06-7cce1835d686
# ╠═2d9caf6d-154d-408e-8e4d-aef4f1c3ff78
# ╠═1a8aac68-4308-4700-ad2d-8e6989757416
# ╠═4b5089dd-74aa-4dbf-8e26-bdaecc0b5b36
# ╠═7361c452-d2c9-4793-855e-055911702d8a
# ╠═ae549f94-6e59-48f2-aedd-5703b7301202
# ╠═0a46f8ed-2a03-4b04-98be-5d85473343ea
# ╠═3c36cc71-db84-4f59-81fe-e4dd7761f72b
# ╠═e9d02b89-5a52-424f-a371-d0ca248bf316
# ╠═913d87d1-6d0a-4db3-9e46-f8a10d247dd5
# ╠═faddf0dd-f931-4379-a13a-021bd5b23cae
# ╠═c4a7d77d-68cd-40ef-a168-6dbf9dc410b8
# ╠═d3567161-530b-476a-b339-a0b2d415da7e
# ╠═ff2cea1e-4bb4-4cd2-a1f4-a68cfeb6efde
# ╠═26299f99-3851-4130-84b6-eaac2c59fbd8
# ╠═da893d5b-1fb4-41f1-a8c9-16c09196739a
# ╠═f474351b-6b02-4402-a722-9bd03ac6b1b9
# ╠═78ac5d01-2ded-41c5-9b6b-64dfc94a2141
# ╠═69b7aa3a-7d0f-4af9-8c71-f1479d60225b
# ╠═577ba306-497a-4cfc-8ea4-bdbc8afb2076
# ╠═f621f4d2-eb36-4733-84a1-9a77f3d88cb9
# ╠═31aa6d10-3ba3-4af0-9458-2c6a60aa01d1
# ╠═552238b4-26f1-4a3c-9706-65410a8ea4d1
# ╠═966e302e-5ab9-4806-ac72-244eb4be2aae
# ╠═eb32157e-74aa-4e61-9b4c-3b57bfa1cb0a
# ╠═734aebf9-6da9-4e6e-bc79-2905bcd75e7a
# ╠═a6b67788-7967-47fc-9d47-bc4e7c202f5a
# ╠═9a545506-f35d-48b8-aeff-ca24d8ab6442
# ╠═3b9fd69d-ca22-4f05-a563-3492ddd66285
# ╠═4167bc84-b881-4eb0-b140-aa35c0acd365
# ╠═c92f80c5-6f47-47c3-92d8-8d1cc9fd6c42
# ╠═50443a47-3600-4733-8b77-1d56717b00b4
# ╠═3060c2aa-078b-43e3-8451-b17642f470eb
# ╠═a65930bd-58ca-4326-a04d-0416614b7abb
# ╠═5b380673-9c93-4e6f-a158-efe469ab8c52
# ╠═e265dcee-67b7-4238-987a-19242ffe9f03
# ╠═d3dd2610-f200-4b22-bfb3-6d28df96de13
# ╠═e6e36c0d-f522-4959-a979-20573ffda8c2
# ╠═55243de1-4f6e-4184-a1c3-40b3e24b1a1b
# ╠═ee927604-3622-4dd0-8f02-b6838e036fc1
# ╠═2254d304-dc0a-4dae-87d7-f3cd18ebe152
# ╠═4e1e1f43-de2e-4584-ac04-31bcfd7036ad
# ╠═5fbffc42-e89e-47d8-b493-e297ad66ca0b
# ╠═8201e87f-19c7-4b61-bfb8-46ed75bfedcc
# ╠═3dd015d7-3c68-4b7b-9331-d6bd51c3702d
# ╠═e9c15162-c3d4-449f-88db-0280b5874a69
# ╠═70267a1f-34da-49c6-a09f-1c214d0dfecb
# ╠═c15eca41-5fda-4869-a298-4991fa4ad72d
# ╠═a0e6818d-4d4d-4644-84ab-a4c45c9dcb26
# ╠═36477658-b2f5-4ff4-b57a-72367ed61819
# ╠═4593f1eb-c752-4043-98d9-1a00954d738e
# ╠═0f4380ab-a77f-4413-823f-c16044efd721
# ╠═91d20be7-a7a4-4120-8756-e8e5656ec804
# ╠═0e7d33ba-f19d-45b3-8606-0fe9eb1845c3
# ╠═bcd1df63-e0b6-467b-9858-de5a15a61765
# ╠═05558267-7c22-4480-99f7-847b02a273cd
# ╠═ce1658a1-67fd-4743-90f3-9dc914d857be
# ╠═a1052731-3d45-42e5-b851-9dc0a8cf8177
# ╠═749d91db-cebd-43f2-a22b-248fd9b8587d
# ╠═730a04cc-b8db-406a-94ad-f951baf08adf
# ╠═8805231b-7e9f-4dec-a999-ac8e750a49b9
# ╠═a6b2b943-c350-4d99-9124-cfaaed1ddd85
# ╠═9cb352a6-ae66-49a5-8104-9dfb5efeaf33
# ╠═e75d091d-a8e2-47c6-b6a6-1564bb0e196c
# ╠═215df8c8-90a0-4edb-8a19-158c4d2ecafe
# ╠═e8d917c2-269f-4fa8-8aa5-98761db2168f
# ╠═99f3cb7a-6da6-46d0-a21a-81142b26d81f
# ╠═e66bccc8-3650-45cc-93a3-b131d28aacb4
# ╠═5ec60c75-033c-46d1-9cd1-28eb515c1da6
# ╠═dc4f92af-2a11-4d75-b0f0-bf49b88ea345
# ╠═e38e88da-a4bb-4bed-9402-d878af4b454e
# ╠═c82eaa43-29f4-4a05-bfcb-8f8e30cd86cc
# ╠═710c7efd-53fd-45a0-83bb-7571e78110ad
# ╠═bff0c718-491d-4936-9b73-c95406488576
# ╠═3c99cb4e-1d7a-4e81-889c-231bb63fea07
# ╠═954a0a6c-9df4-40ed-afa0-3fe8fdc7cc1a
# ╠═0527aada-3d52-4b73-9849-4e7c47c58d6f
# ╠═a4f502ac-a048-481c-b497-18725863e2b7
# ╠═705c8dac-fcc0-4d65-af7a-20c0c78beb3f
# ╠═bfe2b3ed-16fc-439f-87ff-25d74237c300
# ╠═07c66f25-2f0c-49db-be23-3e26b34be383
# ╠═66397d80-9c61-4ccb-8ec9-a161b3a9eb59
# ╠═98121a40-0e41-42a7-ac26-d7cae9a68273
# ╠═6d436275-9bce-4e53-a6c5-9e59acd12fda
# ╠═6df047c1-f3a7-44ac-a444-f470e17d49a4
# ╠═3403dfc2-27fd-439d-875f-b981b946d603
# ╠═a88a543d-d11d-4941-b588-469b67a8f1cc
# ╠═6e6dc490-df49-4a4e-ba3a-11667920280c
# ╠═5cdc22db-0de2-4e82-a31b-57569f9427cb
# ╠═9393db10-0fc1-42a3-8447-f963c2052290
# ╠═b6d27714-b4cb-4ed9-99d1-1d361814bd2d
# ╠═c3ff86aa-3c15-4564-9153-b4e75306ed54
# ╠═f4100e9f-9456-4663-9328-485adc82baf2
# ╠═e2a597a1-db29-46ba-83dc-514f5a8e8f9f
# ╠═705da010-749c-4841-9259-c323ce711464
# ╠═b1fbabea-5f72-4069-a90a-1461b8f32dbc
# ╠═fd2de26f-49b8-4171-acfc-7758152e4af1
# ╠═a0bfd280-8ff7-4376-a00f-aee71ae3f6eb
# ╠═a731b6b7-bc7e-4e92-a54d-d0ccc8ff678a
# ╠═ef9b0edc-6a46-4f04-a611-287c94545a56
# ╠═60a54f27-165d-4ff1-a526-4f9c957800ef
# ╠═aceaef06-1e07-4059-a228-628410ffdf50
# ╠═8e63a621-9eb5-4c8b-a7cd-74ad038e9814
# ╠═ad08769c-6387-4f68-b87e-8a085e50a453
# ╠═e6c1302f-dbdc-4cce-b9ca-f5f09ca60eed
# ╠═631f6203-2979-47d7-a16c-047694bab32d
# ╠═2b2f4216-41fa-4f7f-b707-35bbcec979df
# ╠═9ec908a7-d59d-477f-80fe-8765dde8019e
# ╠═24b4949a-8947-4be4-af8b-b5979dda3c30
# ╠═2912aef8-60da-4122-9392-e7223512b108
# ╠═656e85b6-30e7-4943-b15d-f1f8654fee26
# ╠═94240437-5cb9-474b-9204-7eadb5aaa0a3
# ╠═72b3a744-e8aa-4cd8-a178-579aa6d564d6
# ╠═dc29fb1e-82b7-4db1-8b18-7d63f63c3395
# ╠═61a0ec53-8374-4399-9600-87061c8f2c61
# ╠═ffd6defb-70f5-4366-a42a-588f37d25c66
# ╠═f1e9ce86-974d-41cc-8cc0-13a90f20fc8b
# ╠═b84d3ad5-f11e-4738-af7c-09c92c8d9de0
# ╠═282cd81c-43b0-4122-8f8b-6b6d9ff14f35
# ╠═de2c4249-900d-4fee-941d-5a865ba1b296
# ╠═8a95d8d5-daa6-4efc-92d5-35b363660dd7
# ╠═eddcd4e1-1ae0-40bb-8d0a-751fdd416caf
# ╠═879eaae7-21a2-40e7-9c31-ab9603e4947b
# ╠═d57a1fed-bd4a-485f-8aef-8cacfb7a294a
# ╠═213085a9-e9de-4404-80de-a26b2566caa8
# ╠═db26dd6d-23c0-495c-973d-b4829d05ff57
# ╠═ea8b9460-5973-45c5-ac50-181ef50800da
# ╠═0919f200-3984-4cab-9ac0-bf6e982ad566
# ╠═19894766-0820-47a2-9a0b-9dab37ca5daa
# ╠═dbf67725-f43f-49a7-8ef5-57b60acdbc67
# ╠═faf85657-eb0e-4b9e-86db-74a3cc6b1c04
# ╠═c107996c-abeb-423f-90a0-6aa5125df163
# ╠═e2557d68-4418-4266-9e6b-adc99f5e5d9c
# ╠═50e26ae2-bc4b-4324-88f7-8588276cc44a
# ╠═4e8b6a0d-40bc-418c-8520-fa9365f43a85
# ╠═0fc62ea4-35bd-410d-a0b9-2ad09140b327
# ╠═e266cb67-ce59-49e5-8ce3-5ca8e5128ded
# ╠═9e4877fe-5c7c-493d-aa24-a5d22fbde4fd
# ╠═b5febe6e-e9b5-46b8-8e09-700e28abe37f
# ╠═8334730d-6e9e-46f8-a797-0d4730df3ab4
# ╠═33a269c7-9271-4f58-ab63-1a0264cd1965
# ╠═cd4f0a10-228a-473d-a8d3-6e7af9fe42da
# ╠═0e5025d6-3fc0-4380-9c1d-e1265ca4e58b
# ╠═a5392b4a-bf22-4936-8426-476e85ccd92e
# ╠═755f8682-9890-4edb-bd0e-7c02f8acb3b2
# ╠═1613721e-9aa9-4898-aeff-68e90cd2a6bb
# ╠═6b71f0e3-288b-490a-9ab8-aeef32021708
# ╠═0555b650-e9a2-48c6-939e-c146241e125c
# ╠═ebeafdca-d72f-4f7e-88ae-93074150e041
# ╠═5f9e87a6-e6b2-4a81-8f0a-0a5ba174b45c
# ╠═854eef7a-e911-4f5b-8562-5bbd22e0c900
# ╠═9bd1d228-e255-48ec-ba0f-07e5ac026581
# ╠═a7910931-bbdf-42c5-821f-6af9996933c7
# ╠═7dd9fe5b-c6d2-4e69-adf5-c7240f54e020
# ╠═d112cefd-e4f9-45ba-a96e-9f1fe7fbd353
# ╠═1bc79106-8075-4b6f-aac8-17cc2f58ea88
# ╠═bf1d3910-2680-40c4-95e9-20969fb5c1f8
# ╠═1bcc41a6-aa15-4118-8a88-0a160ce1df5f
# ╠═0bf65cb8-d144-4384-b7c1-682c219a71ba
# ╠═50cc58bc-2b7a-4501-8405-5f178beddc44
# ╠═4db1d361-eb5c-45cd-9b4f-aaeb9b5528b6
# ╠═dc7a2b77-4b3c-46d6-bc09-1a98dabf8c88
# ╠═795cbda2-5cb9-4e63-bc09-43b55288dd5c
# ╠═e35dd878-06c3-4b3e-99f1-bd24705d8a6b
# ╠═6bfecb6c-b750-46d6-b5c1-1a0caf732809
# ╠═347a22fe-b7a4-422b-a229-d34369e64422
# ╠═d762ff2f-e90d-420b-9ae8-f5186316b5fa
# ╠═8b238fa0-94b1-46a1-b6c0-fb3be10ecd6f
# ╠═6623763d-a9f8-492f-9f72-9054ec3d0ac4
# ╠═adb29d24-b447-4301-9a8b-bd197b4401da
# ╠═03357cc7-cba4-459d-8233-e268c24f2d91
# ╠═de9108ea-7d68-42f5-a637-704866510e3a
# ╠═3257ee0b-5fa8-4680-8c46-e81d515c4b05
# ╠═1c110e43-7ef9-43fb-ad2d-48bcf710e461
# ╠═a8344213-28f9-4627-b4ba-1c336763bc9b
# ╠═c72aa1ad-ea8d-4f4f-afde-abdc65d00fbf
# ╠═77ffbc65-469d-4ec6-ab6e-102d41b9afa4
# ╠═406367df-98a1-4112-ac0d-7b579440d8be
# ╠═84e72546-516a-488f-9d6c-e7d006a5c427
# ╠═0b09c0af-c716-48d9-8a63-7a8a9e534c60
# ╠═efee0b3a-16be-42b0-a9c5-1d1244c702fc
# ╠═9ee57741-1ed8-4b3d-a665-c2e782140da5
# ╠═7a6530d5-cd06-4e50-aba7-685e7517fa4e
# ╠═c81bea74-6d32-4a18-ba04-92cdc1afaa6a
# ╠═946b304a-0b43-4e3b-a9d4-550bbdf3e89e
# ╠═bd96030a-8920-41ac-b58e-f204deea0662
# ╠═1eab17ad-dbba-482c-8936-8968eae477ea
# ╠═3ba622e6-0d55-413f-92f6-4aa54144920c
# ╠═9836bc74-eff7-4d7a-853d-2398df762316
# ╠═af23935b-6ef5-4926-8ad8-b4e6d3b9fd4d
# ╠═a94207bd-7666-41ab-a858-bd1170693a39
# ╠═266de19f-4d79-45f7-8af2-071904b08b47
# ╠═9bdf330c-f5c1-4a42-9a30-caee4bd8bb3e
# ╠═d47b96bd-3b7e-49ca-a240-4accfa3027e7
# ╠═731c89c3-ca71-404d-8f80-cf3062eaec01
# ╠═03f335b8-97de-4b2b-b03d-000800ae6565
# ╠═7c0d7864-53ce-4fa1-bc8d-b87ce0eb2801
# ╠═6c4e95d8-39fb-4e44-9778-778df71d9888
# ╠═4ac53947-763f-40bc-afed-6b3c05de042e
# ╠═9ee56964-d715-4faf-b740-6ce72e6b8e51
# ╠═cfd23630-c1dd-4848-b6eb-9a93f7f59cde
# ╠═585b9beb-7e21-45ab-9794-e95c55799e55
# ╠═da330153-49ee-41bd-82b3-43eb87c8ffa1
# ╠═a1df4246-3313-44f5-9d81-bf66b27b2eae
# ╠═c49f764e-4718-4678-869d-ac1dfdfb9492
# ╠═4d554305-bc31-4d9c-b459-5bde0665fbea
# ╠═b4718644-4b13-4473-851d-2fe8bb0be509
# ╠═e948eff2-d2e1-4067-8fb0-36297a48770a
# ╠═c8e12743-31d4-47eb-a2ec-64ab133e6cbe
# ╠═5b73e0b5-321b-4433-90c6-63c788f0d578
# ╠═b7ee8e0a-7f2f-4f05-ac43-ceb5ba893715
# ╠═3abe560f-6acd-43b3-83d8-850a5ccf3345
# ╠═b51dbb01-ce9b-4076-aa78-2b5020503539
# ╠═7627b3c9-5b1a-45fe-a885-5e42266629a9
# ╠═d448709c-ec77-4ab4-8f8a-eeaddc5a390c
# ╠═6a1eb59c-9249-4f3c-a66c-c44dbf803869
# ╠═de14604b-2100-4c92-9297-733e832f4f4e
# ╠═099f01b5-f129-4159-af77-4e8f7aee07b7
# ╠═64c40f7a-82f1-4362-9456-4f93284faa88
# ╠═2a95d9b1-e129-4025-8a2e-26ab4bd86bfb
# ╠═c4575309-b6f0-4b8b-bb40-173c9b28fd4d
# ╠═eb668cda-a026-4677-84c5-ab8da1850669
# ╠═9ec3cf67-128b-43c5-b993-30bdab77ddee
# ╠═1f56cac3-c7db-448e-9d5c-50e5e9f2b490
# ╠═69e63c65-1fe3-4106-b90b-d80978c9c609
# ╠═9d6c1ce1-17ce-4269-ae77-8f1e0dc67a77
# ╠═3563fda2-b3e1-4bba-831c-63a7bcfda780
# ╠═fc3e907d-3fc4-42ee-9bd6-79a79aa1ee4f
# ╠═c9d4449b-4d4b-401b-86c3-c52cbaba8060
# ╠═1a0a63a1-4cb6-4578-a32a-7c623670f0f0
# ╠═212eda1c-da46-4833-9f58-c98d437c450f
# ╠═435ffc34-5696-49ac-8205-6ab4a92ca3b2
# ╠═f5996482-be98-4751-bbd8-19679d03658c
# ╠═6a5c262d-b493-43fe-9f78-7223c15734ed
# ╠═91728250-037c-419b-a531-523b16ea9aea
# ╠═4dabb682-d787-4b18-b78e-a3be8ee90431
# ╠═71d16997-a2e4-4cc5-94f9-aaddf5d90339
# ╠═dcce6567-b5ad-4e7a-ace7-81490ed559be
# ╠═575a5fd5-c5ed-4a05-9927-dc06e56817b9
# ╠═ded94e23-bdbb-4aa6-9dfd-1956ff3fc89a
# ╠═e064a2e3-63dd-4b95-a905-e27e9d1fd09f
# ╠═79a8f8c0-0578-4729-a13f-9245dfe67593
# ╠═61fc008e-6f8c-4088-b70c-954e84626284
# ╠═2bc6e9c6-1745-4701-ae0f-8c15eb2cee17
# ╠═0a6a8645-8dbd-4605-b6e4-3fd54ef7d4a4
# ╠═52583e7c-20d1-42a9-a282-6c6940bb5f20
# ╠═21f353bc-05b7-416e-84af-faec44839f8e
# ╠═cdc3d63a-446e-4958-84ab-fe7ffd33946c
# ╠═7e10b591-32a0-453e-a8f5-91860a8ddc01
# ╠═f2eff5a8-6f20-4577-89b5-ded89cc0bc08
# ╠═429172af-a05c-45ff-a1e6-70f70c88d158
# ╠═0567fc6c-c45a-4cd3-96dd-acc870f35150
# ╠═5c1fb6d1-7942-47b6-ab4c-c37bfedd158e
# ╠═0ad5996c-1b2a-4250-b857-b0ddfec3b388
# ╠═144b4518-7fde-411d-86ea-56279dea4d03
# ╠═5d8b869d-44ff-487e-a3ca-62787d706830
# ╠═d74712d3-7fa1-4ff6-85ec-80aaa9eaa682
# ╠═529bc9f9-3537-45bd-ae86-4b32079febf1
# ╠═9e75dd8f-f3fa-4fc6-a654-c6699c6f09ee
# ╠═5682cadd-0a47-47be-afff-8902e92a6898
# ╠═c2d014bf-553f-4a6b-bbfc-a96db1bfa11e
# ╠═d10bc4aa-09b2-401d-860e-07adb32ff17e
# ╠═ab42f3bb-4473-4bf4-967e-85bb4cbb3532
# ╠═387e4645-cc41-4184-9612-d242129b27e3
# ╠═adced47b-69d1-43f4-931c-15d80801991e
# ╠═29d6e7ac-3d56-4473-82b6-4d17fbda0e97
# ╠═ec460cc7-a688-488f-ac70-43576329a4f7
# ╠═6d788695-bbba-443f-a8d2-ee7672fcb6cd
# ╠═922095ef-fcba-416b-92f8-aea0d1dbd961
# ╠═a5f61fab-b1b4-4834-8945-b63947f81776
# ╠═be6569c6-00d7-4201-b003-6b6852a3cd02
# ╠═db0b4044-a42a-4c4e-9911-c24baa258818
# ╠═33c0f5a5-1eab-4857-b6e5-dc1625a05c61
# ╠═54f04bdf-240e-4c55-83a1-7539e8679b95
# ╠═a5c05aef-84f5-4a58-a663-0478c4fdcc35
# ╠═4470fa12-f8cb-4f98-9c33-95310b7dacda
# ╠═13765151-447d-41c8-afd1-8150d600deae
# ╠═2ec7f031-3e5d-4f1e-92ce-1a8f9f5387a8
# ╠═3e03a99c-f2fa-4dd8-a10a-46a3a886a393
# ╠═a84865e1-6f00-4fa2-89ce-254492b7abd6
# ╠═7abde529-0206-4cff-b1b6-cfa23b2a966b
# ╠═c568ef3d-a8f3-49dc-8bd3-76a5233f5442
# ╠═21db476e-e901-478f-b5d6-d1c75a0017e3
# ╠═0cdbf519-0c60-48f9-b32a-37d6f3ac667c
# ╠═04db972c-b08b-44c1-be9d-388484b133e5
# ╠═372dcaef-18a1-4a1d-b292-6d2a285e4d83
# ╠═5bf138b3-7907-4fc6-96e7-d3a8bccc78af
# ╠═bc7932b4-80e1-4f05-981e-0c869c0dce85
# ╠═d45c9d47-68f2-44d7-82ad-ddcb4877263f
# ╠═a7eec8b0-9a15-42c0-aa0c-95675a21dcf7
# ╠═faa56f3f-6891-423c-89e3-95099dcc149c
# ╠═5efd1815-277b-4444-803d-b48336b0abd7
# ╠═da76c242-e945-4336-9aa9-f2aa267fb975
# ╠═7ef78147-f3a9-4a73-a179-6d9b302629dd
# ╠═7e5f1c7b-1957-453c-8aaa-90719e9ff215
# ╠═a7269b72-6398-4016-ab2d-28cb6a024431
# ╠═fff284bb-8469-4441-ba77-f20bb296fc79
# ╠═a03581d9-edf1-4932-9029-5336baa4b1d6
# ╠═ede4a4d6-206b-4f24-9485-417c2fa982cc
# ╠═6050ae2f-92c9-427e-8058-3a61d8f4e8a3
# ╠═479fb2d0-bae4-4a3f-b0a7-03f94ae442cb
# ╠═1007c9b7-2aaa-4cdc-9cfe-4efc14642ef7
# ╠═d8795ad6-415d-46e0-8abb-e9c8f7bf926d
# ╠═2ee422ad-994f-42ab-9127-f77f157bebbd
# ╠═8eb74e1d-a6e7-4b13-85c2-8b29ca971496
# ╠═3f622f6d-9828-4742-b040-e3bfd4749479
