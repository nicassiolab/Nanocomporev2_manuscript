### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ bee0055d-8805-47e7-b995-f74d059dd9d4
using Pkg

# ╔═╡ 1c8593c2-a355-4016-a7b5-9dc58570226a
Pkg.add("Peaks")

# ╔═╡ f7b34149-9459-47fc-837b-43a581fc86ed
Pkg.add("Typst_jll")

# ╔═╡ f86b30f1-3a4d-4e7e-94f6-67ed90b5b239
Pkg.add("ConfSets")

# ╔═╡ 8c95186a-1d0c-4c50-943f-80fe43a83fa5
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

# ╔═╡ f519ebe0-e962-46f7-a7e1-1e9a0a741ddc
using Peaks

# ╔═╡ 5da72439-e4f1-4900-bcf4-9a4e977a64b0
using AlgebraOfGraphics

# ╔═╡ 9b2b3288-467b-4c5a-8972-4900b7182a17
using Typst_jll

# ╔═╡ a3ec7589-3e97-4238-91d3-a4bc286fb486
using ConfSets

# ╔═╡ 034a73ca-3f8e-11f0-0aea-77f9aa7561f4
include("lib.jl")

# ╔═╡ 6e0b2772-54ac-4bba-8a26-647c7b32f520
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

# ╔═╡ c0ccdc56-02cf-46b3-b721-dca81ab3f8e6
# ╠═╡ disabled = true
#=╠═╡
rmbase_psu = DataFrame(CSV.File("/home/mzdravkov/rmbase_mods/human.hg38.Pseudo.result.col29.bed",
				   header = ["chr", "start", "end", "mod_id", "score", "strand", "mod_type", "support_num", "support_list", "support_list_sub", "pub_list", "cell_list", "seq_type_list", "gene_id", "transcript_id", "gene_name", "gene_type", "region", "seq", "motif_score", "conserved_sites_list", "snoRNA_detail_info", "snoRNA_guide_site", "snoRNA_name_list", "snoRNA_database", "writer_loc", "writer_id", "writer_name_list", "source"]))
  ╠═╡ =#

# ╔═╡ 1e15d868-6cf1-4607-bc5c-567980ed2a9e
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ ac1c85bc-a140-4ccd-812f-4d3885eccf14
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 98a50473-a972-4739-8701-1ad562712dcc
begin
	ref002 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv"))
	ref002 = combine(groupby(ref002, ["chr", "strand", "genomicPos"]), :modified => maximum => :modified)
end

# ╔═╡ 5a8124d1-9d98-4639-bcbf-ffe378de6b7e
begin
	ref004 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx_GLORI_176k.tsv"))
	ref004 = combine(groupby(ref004, ["chr", "strand", "genomicPos"]), :modified => maximum => :modified)
end

# ╔═╡ e23c227e-29f2-4757-8b16-4e5c5c07ebd9
rna004 = annotate_results(
	read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_fix_mod_clust_inferring/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=true),
	ref004)

# ╔═╡ cd61aa21-9ad5-4dc2-ba26-9685f86186d4
rna002 = annotate_results(
	read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
	    genomic_collapse=true),
	ref002)

# ╔═╡ 5d08f4fc-52af-4e5a-ba7f-de5a3b8ef00b
rna004_p12 = annotate_results(
	read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_p12_motor_fix_gap_predictions_logdwell/out_nanocompore_results_gx.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=true),
	ref004)

# ╔═╡ ca0c4038-4104-4fb1-b9e5-53b282ad1f6a
rna002_p10 = annotate_results(
	read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_p10_motor_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
	    genomic_collapse=true),
	ref002)

# ╔═╡ 72a07892-065a-4781-ac1c-7fd4c60538d2
md"""
## RNA002
"""

# ╔═╡ a117edb7-04fc-4b6c-8b44-04a695d23117
rna002_tx = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
	    genomic_collapse=false)

# ╔═╡ fa1a7e46-84ab-4c56-a8bb-5f8dab4d05f9
rna002_p10_tx = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_p10_motor_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
	    genomic_collapse=false)

# ╔═╡ 29d93abb-89d1-407f-9d1e-cfc2f80853b8
sum(rna002_tx.GMM_chi2_pvalue .=== missing)

# ╔═╡ 952c3e5f-775d-4bf1-b4d6-99c4fd437aa6
sum(rna002_tx.GMM_chi2_qvalue .< 0.01)

# ╔═╡ 97574ea9-2c75-4161-9e7e-d3a90f5eac96
sum(rna002_p10_tx.GMM_chi2_pvalue .=== missing)

# ╔═╡ c4d62836-34ea-4b5b-8eec-40b74ee5c882
sum(rna002_p10_tx.GMM_chi2_qvalue .< 0.01)

# ╔═╡ 32d4efb9-0ec6-4499-a47d-0ebc12ec68a1
begin
	all_positions_002 = combine(groupby(rna002_tx, :ref_id),
							:pos => (p -> minimum(p):maximum(p)) => :pos)
	all_positions_002 = leftjoin(all_positions_002, rna002_tx,
			 				 on = [:ref_id, :pos])
	all_positions_002 = sort(all_positions_002, [:ref_id, :pos])
end

# ╔═╡ 9c116561-0239-497a-abda-7236159f34f0
begin
	local fn = function(i, d)
		zip(-25:25,
			crosscor(-log10.(coalesce.(i, 1)),
					 -log10.(coalesce.(d, 1)),
					 -25:25))
	end
	local valid_refs = filter(r -> r.count > 50,
							  combine(groupby(all_positions_002, :ref_id), nrow => :count)).ref_id |> Set
	local valid_mask = map(r -> in(r, valid_refs),
						   all_positions_002.ref_id)
	ref_crosscorrs_002 = combine(groupby(all_positions_002[valid_mask, :], :ref_id),
							 [:KS_intensity_pvalue, :KS_dwell_pvalue] => fn => [:offset, :xcorr])
end

# ╔═╡ cace108a-ea0a-4815-9590-4618ef8d114c
mean_xcorrs_002 = combine(groupby(ref_crosscorrs_002, :offset),
						  :xcorr => (vs -> mean(filter(v -> !isnan(v), vs))) => :mean_xcorr)

# ╔═╡ 88599c98-1f77-496c-9273-215ac8f25c01
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "RNA002: STORM/WT intensity/dwell cross-correlation",
				    xlabel = "Offset",
				    ylabel = "Mean cross-correlation",
				    yticks = 0:0.01:0.13)
	ylims!(ax, (0, 0.13))
	scatter!(ax, mean_xcorrs_002.offset, mean_xcorrs_002.mean_xcorr)
	f
end

# ╔═╡ 336c7723-a9c5-42f0-9f7a-f4381385ddce
# ╠═╡ disabled = true
#=╠═╡
begin
	cors_002 = []
	local t = copy(rna002_tx)
	for i in -25:25
		
		t[!, "motor_pos"] = disallowmissing(t.pos) .+ i
		tj = innerjoin(t, t,
					   on=[:ref_id => :ref_id,
					   	   :motor_pos => :pos],
					   makeunique=true)
		push!(cors_002, cor(tj[tj.KS_intensity_qvalue .!== missing .&& tj.KS_intensity_qvalue .<= 0.01, "KS_intensity_pvalue"],
						tj[tj.KS_intensity_qvalue .!== missing .&& tj.KS_intensity_qvalue .<= 0.01, "KS_dwell_pvalue_1"]))
	end
	# scatter(-15:15, cors)
end
  ╠═╡ =#

# ╔═╡ 63ea2754-bd01-4a94-a703-cf7906e72131
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "RNA002: Correlation of intensity/dwell for significant (m6A) sites",
				    xlabel = "Offset",
				    ylabel = "Correlation",
				    yticks = -0.1:0.01:0.16)
	ylims!(ax, (-0.01, 0.16))
	scatter!(ax, -25:25, cors_002)
	f
end
  ╠═╡ =#

# ╔═╡ ca66bb62-13e5-4572-a545-d99c47ec9bb5
md"""
## RNA004
"""

# ╔═╡ d9832607-e012-4d6a-a0b8-cdbd615c4b77
rna004_tx = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_fix_mod_clust_inferring/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false)

# ╔═╡ 03e999dd-a88d-40e6-b780-af26f51af72b
rna004_p12_tx = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_p12_motor_fix_gap_predictions_logdwell/out_nanocompore_results_gx.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false)

# ╔═╡ 3ea52b6d-3896-4ac9-a180-180ed0da3d3f
begin
	all_positions = combine(groupby(rna004_tx, :ref_id),
							:pos => (p -> minimum(p):maximum(p)) => :pos)
	all_positions = leftjoin(all_positions, rna004_tx,
			 				 on = [:ref_id, :pos])
	all_positions = sort(all_positions, [:ref_id, :pos])
end

# ╔═╡ 252600bf-4ed1-4978-a673-301d14090eea
begin
	local fn = function(i, d)
		zip(-25:25,
			crosscor(-log10.(coalesce.(i, 1)),
					 -log10.(coalesce.(d, 1)),
					 -25:25))
	end
	local valid_refs = filter(r -> r.count > 50,
							  combine(groupby(all_positions, :ref_id), nrow => :count)).ref_id |> Set
	local valid_mask = map(r -> in(r, valid_refs),
						   all_positions.ref_id)
	ref_crosscorrs = combine(groupby(all_positions[valid_mask, :], :ref_id),
	    [:KS_intensity_pvalue, :KS_dwell_pvalue] => fn => [:offset, :xcorr])
end

# ╔═╡ a5a923da-22fc-4798-9265-0614ce977466
mean_xcorrs = combine(groupby(ref_crosscorrs, :offset),
					  :xcorr => (vs -> mean(filter(v -> !isnan(v), vs))) => :mean_xcorr)


# ╔═╡ 5ee87cfe-d49d-4173-873c-4303bec3f1e1
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "RNA004: STORM/WT intensity/dwell cross-correlation",
				    xlabel = "Offset",
				    ylabel = "Mean cross-correlation",
				    yticks = 0:0.02:0.25)
	ylims!(ax, (0, 0.25))
	scatter!(ax, mean_xcorrs.offset, mean_xcorrs.mean_xcorr)
	f
end

# ╔═╡ e88ffa49-1220-4795-98e2-f422e4ab616f
function expand_labels(df, pos_col, left, right)
	N = left + right + 1
	expanded_df = repeat(df, N)
	pos_offsets = sort(repeat(-left:right, size(df, 1)))
	expanded_df[:, pos_col] += pos_offsets
	expanded_df
end

# ╔═╡ 05df07a5-ed51-4ddf-9c24-abbb0564a342
begin
	local significant = copy(rna004_tx[rna004_tx.KS_intensity_qvalue .!== missing .&& rna004_tx.KS_intensity_qvalue .< 0.01, [:ref_id, :pos]])
	significant[:, :site] = significant.pos
	local expanded = sort(expand_labels(significant, :pos, 25, 25),
						  [:ref_id, :pos])
	local max_len = map(r -> parse(Int, r[7]), split.(expanded.ref_id, "|"))
	expanded = expanded[expanded.pos .< max_len .&& expanded.pos .>= 0, :]
	expanded = unique(expanded)
	local subset = leftjoin(expanded, rna004_tx,
			 				on = [:ref_id, :pos])
	subset = sort(subset, [:ref_id, :pos])
	local fn = function(i, d)
		zip(-25:25,
			crosscor(-log10.(coalesce.(i, 1)),
					 -log10.(coalesce.(d, 1)),
					 -25:25))
	end
	local valid_sites = filter(r -> r.count >= 50,
							  combine(groupby(subset, [:ref_id, :site]), nrow => :count))
	valid_sites = Set([(r.ref_id, r.site) for r in eachrow(valid_sites)])
	local valid_mask = map(r -> in((r.ref_id, r.site), valid_sites),
						     eachrow(subset))
	ref_crosscorrs_sig = combine(groupby(subset[valid_mask, :], [:ref_id, :site]),
								 [:KS_intensity_pvalue, :KS_dwell_pvalue] => fn => [:offset, :xcorr])
end

# ╔═╡ 3a1cecf8-f3e2-48dc-8034-b93ec7ac6033
begin
	local significant = copy(rna004_tx[rna004_tx.KS_intensity_qvalue .!== missing .&& rna004_tx.KS_intensity_qvalue .< 0.01, [:ref_id, :pos]])
	significant[:, :site] = significant.pos
	local expanded = sort(expand_labels(significant, :pos, 25, 25),
						  [:ref_id, :pos])
	local max_len = map(r -> parse(Int, r[7]), split.(expanded.ref_id, "|"))
	expanded = expanded[expanded.pos .< max_len .&& expanded.pos .>= 0, :]
	expanded = unique(expanded)
	local subset = leftjoin(expanded, rna004_tx,
			 				on = [:ref_id, :pos])
	subset_df = sort(subset, [:ref_id, :pos])
end

# ╔═╡ d0f2466a-b371-45fd-8552-82c0ea2b1b9f
begin
	local tx = unique(subset_df.ref_id)[65]
	local df = subset_df[subset_df.ref_id .== tx, :]
	local site = df.site[1]
	df = df[df.site .== site, :]
	local f = Figure()
	local ax1 = Axis(f[1, 1])
	local ax2 = Axis(f[2, 1])
	scatter!(ax1, -25:25, -log10.(df.KS_intensity_pvalue))
	scatter!(ax2, -25:25, -log10.(df.KS_dwell_pvalue))
	f
end

# ╔═╡ 3b9902ae-8c77-4680-ad2b-3cfde857ec4f
mean_xcorrs_sig = combine(groupby(ref_crosscorrs_sig, :offset),
						  :xcorr => (vs -> mean(filter(v -> !isnan(v), vs))) => :mean_xcorr)

# ╔═╡ 4dc3c11f-99ac-4d45-988e-9c6e1122cf63
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Cross-correlation for significant (m6A) sites",
				    xlabel = "Offset",
				    ylabel = "Mean cross-correlation")
	scatter!(ax, mean_xcorrs_sig.offset, mean_xcorrs_sig.mean_xcorr)
	f
end

# ╔═╡ 3d04269e-ba19-4c0e-bd41-4c65e06dd9f0
# ╠═╡ disabled = true
#=╠═╡
begin
	cors = []
	local t = copy(rna004_tx)
	for i in -25:25
		
		t[!, "motor_pos"] = disallowmissing(t.pos) .+ i
		tj = innerjoin(t, t,
					   on=[:ref_id => :ref_id,
					   	   :motor_pos => :pos],
					   makeunique=true)
		push!(cors, cor(tj[tj.KS_intensity_qvalue .!== missing .&& tj.KS_intensity_qvalue .<= 0.01, "KS_intensity_pvalue"],
						tj[tj.KS_intensity_qvalue .!== missing .&& tj.KS_intensity_qvalue .<= 0.01, "KS_dwell_pvalue_1"]))
	end
	# scatter(-15:15, cors)
end
  ╠═╡ =#

# ╔═╡ eaa5e268-4704-47dc-b5f9-70bc28ccf737
#=╠═╡
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "RNA004: Correlation of intensity/dwell for significant (m6A) sites",
				    xlabel = "Offset",
				    ylabel = "Correlation",
				    yticks = -0.1:0.01:0.16)
	ylims!(ax, (-0.01, 0.16))
	scatter!(ax, -25:25, cors)
	f
end
  ╠═╡ =#

# ╔═╡ fa551c76-a994-4b2d-a299-d87ffc4818f2
md"""
## IVT
"""

# ╔═╡ 5d414a5b-633b-44c4-bdbd-300a7ac66629
ivt = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_gmm_gpu_v0.2.8/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false)

# ╔═╡ 18932f3b-af22-424b-97e1-a9a6411e6d78
ivt_p12 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_p12_motor_gmm_gpu_v0.2.8_fixed_motor_assignent//out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=false)

# ╔═╡ 173b84e8-c0ff-4aa7-8416-3c20fd5273a2
begin
	all_positions_ivt = combine(groupby(ivt, :ref_id),
							:pos => (p -> minimum(p):maximum(p)) => :pos)
	all_positions_ivt = leftjoin(all_positions_ivt, ivt,
			 				 on = [:ref_id, :pos])
	all_positions_ivt = sort(all_positions_ivt, [:ref_id, :pos])
end

# ╔═╡ 67fe3eff-4f4a-4948-8bb9-f74b03ab24e3
begin
	local fn = function(i, d)
		zip(-25:25,
			crosscor(-log10.(coalesce.(i, 1)),
					 -log10.(coalesce.(d, 1)),
					 -25:25))
	end
	local valid_refs = filter(r -> r.count > 50,
							  combine(groupby(all_positions_ivt, :ref_id), nrow => :count)).ref_id |> Set
	local valid_mask = map(r -> in(r, valid_refs),
						   all_positions_ivt.ref_id)
	ref_crosscorrs_ivt = combine(groupby(all_positions_ivt[valid_mask, :], :ref_id),
								 [:KS_intensity_pvalue, :KS_dwell_pvalue] => fn => [:offset, :xcorr])
end

# ╔═╡ d5336fe2-b5fb-4c2d-891b-67459a1af7ed
mean_xcorrs_ivt = combine(groupby(ref_crosscorrs_ivt, :offset),
					      :xcorr => (vs -> mean(filter(v -> !isnan(v), vs))) => :mean_xcorr)

# ╔═╡ 9b248593-f585-4031-8f13-b21c50fb2874
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "IVT intensity/dwell cross-correlation",
				    xlabel = "Offset",
				    ylabel = "Mean cross-correlation")
	scatter!(ax, mean_xcorrs_ivt.offset, mean_xcorrs_ivt.mean_xcorr)
	f
end

# ╔═╡ a7d34031-3843-4b69-8956-de98e6a02bd4
begin
	local significant = copy(ivt[ivt.KS_intensity_qvalue .!== missing .&& ivt.KS_intensity_qvalue .< 0.01, [:ref_id, :pos]])
	significant[:, :site] = significant.pos
	local expanded = sort(expand_labels(significant, :pos, 25, 25),
						  [:ref_id, :pos])
	local max_len = map(r -> parse(Int, r[7]), split.(expanded.ref_id, "|"))
	expanded = expanded[expanded.pos .< max_len .&& expanded.pos .>= 0, :]
	expanded = unique(expanded)
	local subset = leftjoin(expanded, ivt,
			 				on = [:ref_id, :pos])
	subset = sort(subset, [:ref_id, :pos])
	local fn = function(i, d)
		zip(-25:25,
			crosscor(-log10.(coalesce.(i, 1)),
					 -log10.(coalesce.(d, 1)),
					 -25:25))
	end
	local valid_sites = filter(r -> r.count >= 50,
							  combine(groupby(subset, [:ref_id, :site]), nrow => :count))
	valid_sites = Set([(r.ref_id, r.site) for r in eachrow(valid_sites)])
	local valid_mask = map(r -> in((r.ref_id, r.site), valid_sites),
						     eachrow(subset))
	ref_crosscorrs_sig_ivt = combine(groupby(subset[valid_mask, :], [:ref_id, :site]),
								 [:KS_intensity_pvalue, :KS_dwell_pvalue] => fn => [:offset, :xcorr])
end

# ╔═╡ cbc3cd90-9c58-4e7f-be58-b5520221f0a9
mean_xcorrs_sig_ivt = combine(groupby(ref_crosscorrs_sig_ivt, :offset),
						      :xcorr => (vs -> mean(filter(v -> !isnan(v), vs))) => :mean_xcorr)

# ╔═╡ 627d783e-7b95-4cb7-a25c-45182d0a54dc
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Cross-correlation for significant sites",
				    xlabel = "Offset",
				    ylabel = "Mean cross-correlation")
	scatter!(ax, mean_xcorrs_sig_ivt.offset, mean_xcorrs_sig_ivt.mean_xcorr)
	f
end

# ╔═╡ c0ce3531-414c-4080-90aa-6f7dd28848dd
begin
	local f = Figure()
	local ax1 = Axis(f[1, 1])
	local ax2 = Axis(f[2, 1])
	local ax3 = Axis(f[3, 1])
	local ax4 = Axis(f[1, 2])
	local ax5 = Axis(f[2, 2])
	local ax6 = Axis(f[3, 2])
	
	local i = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 5, 7, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	scatter!(ax1, i)
	local d = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 5, 7, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	scatter!(ax2, d)
	scatter!(ax3, -25:25, crosscor(i, d, -25:25))

	scatter!(ax4, i)
	local d = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 5, 7, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	scatter!(ax5, d)
	scatter!(ax6, -25:25, crosscor(i, d, -25:25))
	f
end

# ╔═╡ 97e6502b-6a72-4b7e-926a-f58524e94c89
md"""
## Plot
"""

# ╔═╡ 589da79b-73c1-446e-9253-c1a020c69a0f
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "STORM/WT intensity/dwell crosscorr.",
				    xlabel = "Offset",
				    ylabel = "Mean cross-correlation")
	lines!(ax,
		   mean_xcorrs_002.offset,
		   mean_xcorrs_002.mean_xcorr,
		   label = "RNA002",
		   color = COL_RNA002)
	scatter!(ax,
			 mean_xcorrs_002.offset,
			 mean_xcorrs_002.mean_xcorr,
			 markersize = 6,
			 color = COL_RNA002)
	lines!(ax,
		   mean_xcorrs.offset,
		   mean_xcorrs.mean_xcorr,
		   label = "RNA004",
		   color = COL_RNA004)
	scatter!(ax,
			 mean_xcorrs.offset,
			 mean_xcorrs.mean_xcorr,
			 markersize = 6,
			 color = COL_RNA004)
	axislegend(ax)
	local ax2 = Axis(f[1, 2],
					 title = "IVT/WT intensity/dwell crosscorr.",
					 xlabel = "Offset",
					 ylabel = "Mean cross-correlation")
	lines!(ax2,
		   mean_xcorrs_ivt.offset,
		   mean_xcorrs_ivt.mean_xcorr,
		   label = "RNA004",
		   color = COL_RNA004)
	scatter!(ax2,
			 mean_xcorrs_ivt.offset,
			 mean_xcorrs_ivt.mean_xcorr,
			 markersize = 6,
			 color = COL_RNA004)
	axislegend(a2)
	linkaxes!(ax, ax2)
	f
end

# ╔═╡ 76fb3894-6d04-4c97-b7cc-1918cc0be5a0
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
			    title = "RNA002: motor dwell effect",
			    aspect = 1,
		        xlabel = "Recall",
		        ylabel = "Precision")
	
	local a, p, r = auprc(rna002.predicted,
	                	  Int.(rna002.modified),
	                	  Set([1]))
	lines!(ax, r, p,
	       label = "Base\nAUC=$(round(a; digits = 2))",
	       color = :black)
	local evaluation = roc(rna002.modified, rna002.predicted .>= -log10(0.01))
	local p001 = precision(evaluation)
	local r001 = recall(evaluation)
	scatter!(ax, [r001], [p001], color = :black)
	
	local a, p, r = auprc(rna002_p12.predicted,
	                	  Int.(rna002_p12.modified),
	                	  Set([1]))
	lines!(ax, r, p,
	       label = "Motor\nAUC=$(round(a; digits = 2))",
	       color = :red,
		   linestyle = :dash)
	local evaluation = roc(rna002_p12.modified, rna002_p12.predicted .>= -log10(0.01))
	local p001 = precision(evaluation)
	local r001 = recall(evaluation)
	scatter!(ax, [r001], [p001], color = :red)
	axislegend(ax)
	f
end

# ╔═╡ afaa1214-30d1-4cdd-ae35-4e56b5783a01
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
			    title = "RNA004: motor dwell effect",
			    aspect = 1,
		        xlabel = "Recall",
		        ylabel = "Precision")
	
	local a, p, r = auprc(rna004.predicted,
	                	  Int.(rna004.modified),
	                	  Set([1]))
	lines!(ax, r, p,
	       label = "Base\nAUC=$(round(a; digits = 2))",
	       color = :black)
	local evaluation = roc(rna004.modified, rna004.predicted .>= -log10(0.01))
	local p001 = precision(evaluation)
	local r001 = recall(evaluation)
	scatter!(ax, [r001], [p001], color = :black)
	
	local a, p, r = auprc(rna004_p12.predicted,
	                	  Int.(rna004_p12.modified),
	                	  Set([1]))
	lines!(ax, r, p,
	       label = "Motor\nAUC=$(round(a; digits = 2))",
	       color = :red,
		   linestyle = :dash)
	local evaluation = roc(rna004_p12.modified, rna004_p12.predicted .>= -log10(0.01))
	local p001 = precision(evaluation)
	local r001 = recall(evaluation)
	scatter!(ax, [r001], [p001], color = :red)
	axislegend(ax)
	f
end

# ╔═╡ 5fc3bfd3-b429-4f6a-aa39-8eca8c3a62a4
rna004_tx

# ╔═╡ 67db0de7-c7eb-4dce-b06e-9a2832dab5bb
md"""
## Test specific site
"""

# ╔═╡ cc119196-64be-457d-9515-5efccc36251c
mod_ref = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/combined_mod_ref.tsv", delim = '\t'))

# ╔═╡ a65a2f01-6f6f-4627-8f40-42095a410251
begin
	local df = innerjoin(ivt, mod_ref[:, [:chr, :strand, :pos, :mod, :source]],
		  				 on = [:chr, :strand, :genomicPos => :pos])
	df[df.mod .== "Y", :]
end

# ╔═╡ 32d711b6-9c46-4f10-8b53-2dca3669ac7b
begin
	local df = innerjoin(innerjoin(ivt, ivt_p12,
								   on = [:ref_id, :pos],
								   makeunique = true),
						 mod_ref,
						 on = [:chr, :strand, :genomicPos => :pos])

	local pseUs = df[df.mod .== "Y", :]

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Reference pseU sites",
				    xlabel = "No motor dwell: -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)",
				    ylabel = "Motor dwell: -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)")
	scatter!(ax, pseUs.predicted, pseUs.predicted_1,
			 markersize = 6,
			 color = :black)
	# pseUs[:, :change] = abs.(log2.(pseUs.predicted_1 ./ (pseUs.predicted .+ 1e-300)))
	# sort(pseUs, :change, rev = true)

	f
end

# ╔═╡ 1e30107f-49c3-4479-be52-347d4ebb6ac0
begin
	local df = innerjoin(innerjoin(ivt, ivt_p12,
								   on = [:ref_id, :pos],
								   makeunique = true),
						 mod_ref,
						 on = [:chr, :strand, :genomicPos => :pos])

	local pseUs = df[df.mod .== "Y", :]

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Reference pseU sites",
					xlabel = "|LORmotor| - |LOR|",
				    ylabel = "-log10(𝑄-𝑣𝑎𝑙𝑢𝑒 motor) - -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)",
				    xticks = -3:1:3)

	local color = map(r -> if r.predicted >= 2 && r.predicted_1 >= 2
						       :black
						   elseif r.predicted <= 2 && r.predicted_1 >= 2
							   :orange
						   elseif r.predicted >= 2 && r.predicted_1 <= 2
							   :blue
						   else
							   :grey
						   end, eachrow(pseUs))
	
	scatter!(ax, abs.(pseUs.GMM_LOR_1) .- abs.(pseUs.GMM_LOR),
			 -log10.(pseUs.GMM_chi2_pvalue_1) .- -log10.(pseUs.GMM_chi2_pvalue),
			 markersize = 6,
			 color = color,
			 alpha = 0.75)

	local colorcounts = countmap(color)

	Legend(f[1, 2],
		   [PolyElement(color = :orange),
			PolyElement(color = :blue),
			PolyElement(color = :black),
			PolyElement(color = :grey)],
	       ["Becomes significant ($(colorcounts[:orange]))", "Becomes insignificant ($(colorcounts[:blue]))", "Remains significant ($(colorcounts[:black]))", "Remains insignificant ($(colorcounts[:grey]))"],
		   orientation = :vertical,
		   framevisible = false)
	# pseUs[:, :change] = abs.(log2.(pseUs.predicted_1 ./ (pseUs.predicted .+ 1e-300)))
	# sort(pseUs, :change, rev = true)

	f
end

# ╔═╡ fc6da8ec-a2b3-4786-ae9a-978c28d66872
begin
	local df = innerjoin(innerjoin(ivt, ivt_p12,
								   on = [:ref_id, :pos],
								   makeunique = true),
						 mod_ref,
						 on = [:chr, :strand, :genomicPos => :pos])

	local pseUs = df[df.mod .== "Y", :]
	pseUs[:, :change] = abs.(log2.((pseUs.predicted_1 .+ 1e-300) ./ (pseUs.predicted .+ 1e-300)))
	sort(pseUs, :change, rev = true)
end

# ╔═╡ 6dd6aac1-058f-4269-81ad-c5b727229fe3


# ╔═╡ 852fc9dd-0ba1-4c95-baea-19846e5403dc
begin
	local df = innerjoin(innerjoin(ivt, ivt_p12,
								   on = [:ref_id, :pos],
								   makeunique = true),
						 mod_ref,
						 on = [:chr, :strand, :genomicPos => :pos])

	local pseUs = df[df.mod .== "Y", :]

	# local f = Figure()
	# local ax = Axis(f[1, 1],
	# 				title = "Reference pseU sites",
	# 			    xlabel = "No motor dwell: -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)",
	# 			    ylabel = "Motor dwell: -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)")
	# scatter!(ax, pseUs.predicted, pseUs.predicted_1,
	# 		 markersize = 6,
	# 		 color = :black)
	pseUs[:, :change] = log2.(pseUs.predicted_1 ./ (pseUs.predicted .+ 1e-300))
	sort(pseUs, :change, rev = true)
end

# ╔═╡ 5c60b99d-ee4c-4310-b7ac-0c6856141cf8
begin
	local df = innerjoin(innerjoin(ivt, ivt_p12,
								   on = [:ref_id, :pos],
								   makeunique = true),
						 mod_ref,
						 on = [:chr, :strand, :genomicPos => :pos])

	local pseUs = df[df.mod .== "m7G", :]

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Reference m7G sites",
					xlabel = "|LORmotor| - |LOR|",
				    ylabel = "-log10(𝑄-𝑣𝑎𝑙𝑢𝑒 motor) - -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)",
				    xticks = -3:1:3)

	local color = map(r -> if r.predicted >= 2 && r.predicted_1 >= 2
						       :black
						   elseif r.predicted <= 2 && r.predicted_1 >= 2
							   :orange
						   elseif r.predicted >= 2 && r.predicted_1 <= 2
							   :blue
						   else
							   :grey
						   end, eachrow(pseUs))
	
	scatter!(ax, abs.(pseUs.GMM_LOR_1) .- abs.(pseUs.GMM_LOR),
			 -log10.(pseUs.GMM_chi2_pvalue_1) .- -log10.(pseUs.GMM_chi2_pvalue),
			 markersize = 6,
			 color = color,
			 alpha = 0.75)

	local colorcounts = countmap(color)

	Legend(f[1, 2],
		   [PolyElement(color = :orange),
			PolyElement(color = :blue),
			PolyElement(color = :black),
			PolyElement(color = :grey)],
	       ["Becomes significant ($(colorcounts[:orange]))", "Becomes insignificant ($(colorcounts[:blue]))", "Remains significant ($(colorcounts[:black]))", "Remains insignificant ($(colorcounts[:grey]))"],
		   orientation = :vertical,
		   framevisible = false)
	# pseUs[:, :change] = abs.(log2.(pseUs.predicted_1 ./ (pseUs.predicted .+ 1e-300)))
	# sort(pseUs, :change, rev = true)

	f
end

# ╔═╡ 368b97ad-6667-4fd1-89c5-c10890f39bce
begin
	local df = innerjoin(innerjoin(ivt, ivt_p12,
								   on = [:ref_id, :pos],
								   makeunique = true),
						 mod_ref,
						 on = [:chr, :strand, :genomicPos => :pos])

	local pseUs = df[df.mod .== "m1A", :]

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Reference m1A sites",
					xlabel = "|LORmotor| - |LOR|",
				    ylabel = "-log10(𝑄-𝑣𝑎𝑙𝑢𝑒 motor) - -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)",
				    xticks = -3:1:3)

	local color = map(r -> if r.predicted >= 2 && r.predicted_1 >= 2
						       :black
						   elseif r.predicted <= 2 && r.predicted_1 >= 2
							   :orange
						   elseif r.predicted >= 2 && r.predicted_1 <= 2
							   :blue
						   else
							   :grey
						   end, eachrow(pseUs))
	
	scatter!(ax, abs.(pseUs.GMM_LOR_1) .- abs.(pseUs.GMM_LOR),
			 -log10.(pseUs.GMM_chi2_pvalue_1) .- -log10.(pseUs.GMM_chi2_pvalue),
			 markersize = 6,
			 color = color,
			 alpha = 0.75)

	local colorcounts = countmap(color)

	Legend(f[1, 2],
		   [PolyElement(color = :orange),
			PolyElement(color = :blue),
			PolyElement(color = :black),
			PolyElement(color = :grey)],
	       ["Becomes significant ($(get(colorcounts, :orange, 0)))", "Becomes insignificant ($(colorcounts[:blue]))", "Remains significant ($(colorcounts[:black]))", "Remains insignificant ($(colorcounts[:grey]))"],
		   orientation = :vertical,
		   framevisible = false)
	# pseUs[:, :change] = abs.(log2.(pseUs.predicted_1 ./ (pseUs.predicted .+ 1e-300)))
	# sort(pseUs, :change, rev = true)

	f
end

# ╔═╡ f23a1b62-7a57-4a58-b830-e3e0a124681a
begin
	local df = innerjoin(innerjoin(ivt, ivt_p12,
								   on = [:ref_id, :pos],
								   makeunique = true),
						 mod_ref,
						 on = [:chr, :strand, :genomicPos => :pos])

	local pseUs = df[df.mod .== "Nm", :]

	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "Reference Nm sites",
					xlabel = "|LORmotor| - |LOR|",
				    ylabel = "-log10(𝑄-𝑣𝑎𝑙𝑢𝑒 motor) - -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)",
				    xticks = -3:1:3)

	local color = map(r -> if r.predicted >= 2 && r.predicted_1 >= 2
						       :black
						   elseif r.predicted <= 2 && r.predicted_1 >= 2
							   :orange
						   elseif r.predicted >= 2 && r.predicted_1 <= 2
							   :blue
						   else
							   :grey
						   end, eachrow(pseUs))
	
	scatter!(ax, abs.(pseUs.GMM_LOR_1) .- abs.(pseUs.GMM_LOR),
			 -log10.(pseUs.GMM_chi2_pvalue_1) .- -log10.(pseUs.GMM_chi2_pvalue),
			 markersize = 6,
			 color = color,
			 alpha = 0.75)

	local colorcounts = countmap(color)

	Legend(f[1, 2],
		   [PolyElement(color = :orange),
			PolyElement(color = :blue),
			PolyElement(color = :black),
			PolyElement(color = :grey)],
	       ["Becomes significant ($(colorcounts[:orange]))", "Becomes insignificant ($(colorcounts[:blue]))", "Remains significant ($(colorcounts[:black]))", "Remains insignificant ($(colorcounts[:grey]))"],
		   orientation = :vertical,
		   framevisible = false)
	# pseUs[:, :change] = abs.(log2.(pseUs.predicted_1 ./ (pseUs.predicted .+ 1e-300)))
	# sort(pseUs, :change, rev = true)

	f
end

# ╔═╡ 77e469bb-673e-4aa1-8034-1607c738f903
ivt_p12.GMM_3D |> countmap

# ╔═╡ 784cd890-11af-4656-a014-093223f54919
filter(r-> r.mod == "Y", eachrow(innerjoin(ivt_p12, mod_ref,
		  on = [:chr, :strand, :genomicPos => :pos]))).GMM_3D |> countmap

# ╔═╡ 9f03d039-c4de-4d6f-a05e-d1faaaa2dc64
805/(805+178)

# ╔═╡ 1407a831-0d6f-42a4-bb7d-93ee5435098e
antijoin(rna002[rna002.modified .== 0 .&& rna002.predicted .> 2, :],
		 rna002_p10[rna002_p10.modified .== 0 .&& rna002_p10.predicted .> 2, :],
		 on = [:chr, :strand, :genomicPos])

# ╔═╡ 03414a73-8f17-4c96-bb84-daf43119aafc
rna002[rna002.modified .== 0 .&& rna002.predicted .> 2, :]

# ╔═╡ d3ea3b64-373a-4410-a17e-825deaf7501d
rna002_p10[rna002_p10.modified .== 0 .&& rna002_p10.predicted .> 2, :]

# ╔═╡ 7096fcfb-739e-4775-a4f5-1238d97f4a20
rna002[rna002.genomicPos .== 47663400, :]

# ╔═╡ 686af1c3-3e1d-48ca-8efb-e8ba0a818a69
rna002_p10[rna002_p10.genomicPos .== 47663400, :]

# ╔═╡ bd6ad199-ab74-45ae-b4d4-5d65d03753e1
rna002_tx[rna002_tx.genomicPos .== 47663400, :]

# ╔═╡ 34644635-e88d-47aa-9441-8e7583240bea
ref002[ref002.chr .== "chr6" .&& ref002.genomicPos .> 143504128 .&& ref002.genomicPos .< 143504138, :]

# ╔═╡ 5cdf3882-4651-415e-9699-481a7dbe0324
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

# ╔═╡ 13abd3a3-2f85-4c23-9383-15d3df668ad4
rna002_peaks = peaks(rna002_tx, 4)

# ╔═╡ fc50841f-43d9-414b-8fa7-b95b511046d4
function peak_annotate(df, ref, peak_radius)
	df_peaks = peaks(df, peak_radius)
	collapsed_peaks = combine(groupby(df_peaks, [:chr, :strand, :genomicPos]),
							  :predicted => maximum => :predicted)
	result = leftjoin(ref, collapsed_peaks,
			 		  on = [:chr, :strand, :genomicPos])
	result[!, :predicted] = coalesce.(result.predicted, 0)
	sort(result, [:chr, :strand, :genomicPos])
end

# ╔═╡ 2f86ca19-0c1b-44e4-92cb-4b269f3fa12b
peakannot_002 = peak_annotate(rna002_tx, ref002, 4)

# ╔═╡ 7016d0c0-3720-4b50-a01e-6ae04cb013f3
sort(filter(v -> v > 0, peakannot_002.predicted))

# ╔═╡ fe258560-13ce-45f8-96a9-431081463201
peakannot_002_p10 = peak_annotate(rna002_p10_tx, ref002, 4)

# ╔═╡ e1d21630-577f-4888-a0b4-3c75f75a6e67
peakannot_004 = peak_annotate(rna004_tx, ref004, 4)

# ╔═╡ 98dbea9c-6611-4a53-a53b-c57784208fc6
peakannot_004_p12 = peak_annotate(rna004_p12_tx, ref004, 4)

# ╔═╡ 9d3a986e-f40a-4214-9267-dcaff4b013b7
rna002_p10_peaks = peaks(rna002_p10_tx, 4)

# ╔═╡ 081a4b47-79ab-41ec-b688-1a37a5e093a1
function annotate_peaks(peaks, glori)
	local df = copy(peaks)
	df[:, :peak_id] = 1:nrow(df)
	df[:, :start] = max.(df.genomicPos .- 4, 0)
	df[:, :end] = df.genomicPos .+ 4
	local ranges = df[:, [:peak_id, :chr, :strand]]
	ranges[:, :pos] = range.(df.start, df.end)
	local mod_overlaps = innerjoin(flatten(ranges, :pos), glori,
			  					   on = [:chr, :strand, :pos])
	rename!(mod_overlaps, :pos => :mod_pos)
	annot_peaks = leftjoin(df, mod_overlaps,
			 			   on = [:peak_id, :chr, :strand])
	annot_peaks[!, :modified] = coalesce.(annot_peaks.modified, 0)
	annot_peaks[!, :predicted] = coalesce.(annot_peaks.predicted, 0)
	annot_peaks
end

# ╔═╡ c6877e60-d3ff-4122-bbc7-d83da7173fce
begin
	glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k.bed"))
	# glori[!, "pos"] = glori.start .+ 2
	glori[:, :pos] = glori.start
	glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
	glori[!, "modified"] .= 1
	rename!(glori, :chrom => :chr)
	glori
end

# ╔═╡ 5ae2462c-a387-45cc-af66-3a474952d251
rna002_annot_peaks = annotate_peaks(rna002_peaks, glori[:, [:chr, :strand, :pos, :modified]])

# ╔═╡ 5733debc-825e-49a1-8a1b-3c286f71c935
rna002_p10_annot_peaks = annotate_peaks(rna002_p10_peaks, glori[:, [:chr, :strand, :pos, :modified]])

# ╔═╡ 0ab9322b-175a-49c6-973a-369f6bb975e7
rna004_peaks = peaks(rna004_tx, 4)

# ╔═╡ 26ae5e1b-29ed-4a13-8a51-a7bf4337a9f1
rna004_annot_peaks = annotate_peaks(rna004_peaks, glori[:, [:chr, :strand, :pos, :modified]])

# ╔═╡ fb4a5c1b-4a4c-41f1-9ed6-69bf589684b3
rna004_p12_peaks = peaks(rna004_p12_tx, 4)

# ╔═╡ 823797fe-09a0-4de2-b3fa-1bf310a77ab1
rna004_p12_annot_peaks = annotate_peaks(rna004_p12_peaks, glori[:, [:chr, :strand, :pos, :modified]])

# ╔═╡ bb89baaf-d5cf-4e35-982a-83d7cab2e87e
function plot_peak_dist_hist!(fig, gridpos, res, reference, peaks, motor_peaks; title = "", radius = 15, col_palette = [:black, :red], style_palette = [:solid, :dash], titlesize = 28, linewidth = 2, maxy=nothing)
	local ref = innerjoin(reference, res[:, [:chr, :strand, :genomicPos]],
						  on = [:chr, :strand, :genomicPos])
	local mods = unique(ref[ref.modified .== 1, :])

	local peaks2d = peaks[peaks.predicted .> 2, :]

	local distances2d = map(r -> begin
			distances = peaks2d[peaks2d.chr .== r.chr .&& peaks2d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	local peaks3d = motor_peaks[motor_peaks.predicted .> 2, :]

	local distances3d = map(r -> begin
			distances = peaks3d[peaks3d.chr .== r.chr .&& peaks3d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	distances2d = distances2d[abs.(distances2d) .<= radius]
	distances3d = distances3d[abs.(distances3d) .<= radius]

	local distances = vcat(DataFrame(distance = distances2d, type = "MD−"),
						   DataFrame(distance = distances3d, type = "MD+"))

	
	plt = data(distances) *
		  mapping(:distance => "Distance",
				  color = :type  => presorted => "GMM test",
				  linestyle = :type => presorted => "GMM test") * 
		  AlgebraOfGraphics.histogram(Stairs;
									  normalization = :pdf,
									  bins = -radius:(radius+1)) *
		  visual(; linewidth = linewidth)
	# |> draw(; axis = (; title = "RNA002", width = 800, height = 500, xticks = (-14.5:15.5, [string(t) for t in -15:15])))

	# local fig = Figure(size = (800, 500))
	# local gridpos = fig[1, 1]

	local f = draw!(gridpos,
					plt,
					scales(Color = (; palette = col_palette),
						   LineStyle = (; palette = style_palette));
					axis = (; title = title,
							  aspect=1,
							  titlesize = titlesize,
							  xticks = ((-radius + 0.5):(radius+0.5),
										[string(t) for t in -radius:radius]),
				   			  ylabel = "PDF",
						      limits=(nothing, (0, maxy))))
	
	legend!(gridpos, f;
			tellwidth = false,
			halign = :right,
			valign = :top,
			margin = (10, 10, 10, 10),
			patchsize = (40, 40))
	
	# fig
end

# ╔═╡ 191652bf-c223-44f5-991e-9bc49e4d2fa2
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure()
	plot_peak_dist_hist!(f, f[1, 1], rna002_tx, ref002, peakannot_002, peakannot_002_p10; title = "RNA002", radius = 10)
	f
end
  ╠═╡ =#

# ╔═╡ ef3b9161-a00d-480b-823f-6752906f7b85


# ╔═╡ 9daa78c5-d8be-4a43-8d43-364e16dec13e
md"""
## Plot
"""

# ╔═╡ 498343e1-7d59-4688-9075-55501a7115ca
begin
	rna002_cors_sig = []
	local t = copy(rna002_tx)
	for i in -25:25
		t[!, "motor_pos"] = disallowmissing(t.pos) .+ i
		tj = innerjoin(t[t.KS_intensity_qvalue .<= 0.01, :], t,
					   on=[:ref_id => :ref_id,
					   	   :motor_pos => :pos],
					   makeunique=true)
		push!(rna002_cors_sig,
			  cor(tj.KS_intensity_pvalue,
				  tj.KS_dwell_pvalue_1))
	end
end

# ╔═╡ bfcad6a4-aad1-4615-ad1f-ed90935bb8ea
begin
	rna004_cors_sig = []
	local t = copy(rna004_tx)
	for i in -25:25
		t[!, "motor_pos"] = disallowmissing(t.pos) .+ i
		tj = innerjoin(t[t.KS_intensity_qvalue .!== missing .&& t.KS_intensity_qvalue .<= 0.01, :], t,
					   on=[:ref_id => :ref_id,
					   	   :motor_pos => :pos],
					   makeunique=true)
		tj = tj[tj.KS_dwell_pvalue_1 .!== missing, :]
		push!(rna004_cors_sig,
			  cor(tj.KS_intensity_pvalue,
				  tj.KS_dwell_pvalue_1))
	end
end

# ╔═╡ 84a489f9-080e-422f-8eb7-588b89aafb89
function trapezoidal_rule(x, y)
    auc = 0.0
    for i in 2:length(x)
        auc += (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2
    end
    return auc
end

# ╔═╡ e6d9f53d-7ad3-4068-af7d-02ac384f1161
function auc(recall, precision)
    trapezoidal_rule(recall, precision)
end

# ╔═╡ 4625bce2-f17f-46b8-a674-4283ae0a3ae5
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

# ╔═╡ 4a673cb0-f5a5-46a5-8680-6724c51128a6
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
				    title = "RNA002: single base resolution",
				    aspect = 1)
	plot_prc!(ax,
			 peakannot_002.modified,
		 	 peakannot_002.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_002_p10.modified,
		 	 peakannot_002_p10.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_p10.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax, L"\textbf{GMM test}")
	
	f
end

# ╔═╡ 696b8bee-81d7-47d2-8e55-08f490341195
begin
	local peakannot_002_binned = copy(peakannot_002)
	peakannot_002_binned[!, :bin] = div.(peakannot_002.genomicPos, BIN_SIZE)
	peakannot_002_binned = combine(groupby(peakannot_002_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_002_p10_binned = copy(peakannot_002_p10)
	peakannot_002_p10_binned[!, :bin] = div.(peakannot_002_p10.genomicPos, BIN_SIZE)
	peakannot_002_p10_binned = combine(groupby(peakannot_002_p10_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	local f = Figure()
	local ax = Axis(f[1, 1],
				    title = "RNA002: binned",
				    aspect = 1)
	plot_prc!(ax,
			 peakannot_002_binned.modified,
		 	 peakannot_002_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_binned.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_002_p10_binned.modified,
		 	 peakannot_002_p10_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_p10_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax, L"\textbf{GMM test}")
	
	f
end

# ╔═╡ 9c838f8e-fe3a-4ee1-923b-146a165ff16f
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
				    title = "RNA004: single base resolution",
				    aspect = 1)
	plot_prc!(ax,
			 peakannot_004.modified,
		 	 peakannot_004.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_004_p12.modified,
		 	 peakannot_004_p12.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_p12.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax, L"\textbf{GMM test}")
	
	f
end

# ╔═╡ d0bd49e7-0e0e-44c0-a7cb-919f9c47ad6e
begin
	local peakannot_004_binned = copy(peakannot_004)
	peakannot_004_binned[!, :bin] = div.(peakannot_004.genomicPos, BIN_SIZE)
	peakannot_004_binned = combine(groupby(peakannot_004_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_004_p12_binned = copy(peakannot_004_p12)
	peakannot_004_p12_binned[!, :bin] = div.(peakannot_004_p12.genomicPos, BIN_SIZE)
	peakannot_004_p12_binned = combine(groupby(peakannot_004_p12_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	local f = Figure()
	local ax = Axis(f[1, 1],
				    title = "RNA004: binned",
				    aspect = 1)
	plot_prc!(ax,
			 peakannot_004_binned.modified,
		 	 peakannot_004_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_binned.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_004_p12_binned.modified,
		 	 peakannot_004_p12_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_p12_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax, L"\textbf{GMM test}")
	
	f
end

# ╔═╡ a82dddab-3df0-4c56-b3a2-944058a7228a
begin
	local peakannot_002_binned = copy(peakannot_002)
	peakannot_002_binned[!, :bin] = div.(peakannot_002.genomicPos, BIN_SIZE)
	peakannot_002_binned = combine(groupby(peakannot_002_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_002_p10_binned = copy(peakannot_002_p10)
	peakannot_002_p10_binned[!, :bin] = div.(peakannot_002_p10.genomicPos, BIN_SIZE)
	peakannot_002_p10_binned = combine(groupby(peakannot_002_p10_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	local f = Figure(size=(900, 450))
	local ax = Axis(f[1, 1],
				    title = "RNA002",
				    aspect = 1)
	plot_prc!(ax,
			 peakannot_002_binned.modified,
		 	 peakannot_002_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_binned.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_002_p10_binned.modified,
		 	 peakannot_002_p10_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_p10_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax, "GMM test")
	
	local peakannot_004_binned = copy(peakannot_004)
	peakannot_004_binned[!, :bin] = div.(peakannot_004.genomicPos, BIN_SIZE)
	peakannot_004_binned = combine(groupby(peakannot_004_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_004_p12_binned = copy(peakannot_004_p12)
	peakannot_004_p12_binned[!, :bin] = div.(peakannot_004_p12.genomicPos, BIN_SIZE)
	peakannot_004_p12_binned = combine(groupby(peakannot_004_p12_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	local ax2 = Axis(f[1, 2],
				     title = "RNA004",
				     aspect = 1)
	plot_prc!(ax2,
			 peakannot_004_binned.modified,
		 	 peakannot_004_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_binned.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax2,
			 peakannot_004_p12_binned.modified,
		 	 peakannot_004_p12_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_p12_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax2, "GMM test")

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

# ╔═╡ ea40bf7e-7626-44f3-bf22-e2a03897dc7b
# ╠═╡ disabled = true
#=╠═╡
begin
	local f = Figure(size = (2300, 1100), fontsize = 18)
	local ga = f[1, 1] = GridLayout()
	local gb = f[1, 2] = GridLayout()

	local cor_maxy = 0.185

	local ax = Axis(ga[1, 1],
					title = "RNA002 cross-correlation",
				    xlabel = "Offset",
				    ylabel = "Mean correlation",
				    yticks = -0.02:0.02:cor_maxy)
	ylims!(ax, (-0.02, cor_maxy))
	scatter!(ax, mean_xcorrs_002.offset, mean_xcorrs_002.mean_xcorr, markersize = 7, color = :black)
	lines!(ax, mean_xcorrs_002.offset, mean_xcorrs_002.mean_xcorr, color = :black)

	local ax = Axis(ga[1, 2],
					title = "RNA004 cross-correlation",
				    xlabel = "Offset",
				    ylabel = "Mean correlation",
				    yticks = -0.02:0.02:cor_maxy)
	ylims!(ax, (-0.02, cor_maxy))
	scatter!(ax, mean_xcorrs.offset, mean_xcorrs.mean_xcorr, markersize = 7, color = :black)
	lines!(ax, mean_xcorrs.offset, mean_xcorrs.mean_xcorr, color = :black)


	local ax = Axis(ga[2, 1],
					title = "RNA002 modified sites",
				    xlabel = "Offset",
				    ylabel = "Correlation",
				    yticks = -0.02:0.02:cor_maxy)
	ylims!(ax, (-0.02, cor_maxy))
	scatter!(ax, -25:25, rna002_cors_sig, markersize = 7, color = :black)
	lines!(ax, -25:25, rna002_cors_sig, color = :black)

	
	local ax = Axis(ga[2, 2],
					title = "RNA004 modified sites",
				    xlabel = "Offset",
				    ylabel = "Correlation",
				    yticks = -0.02:0.02:cor_maxy)
	ylims!(ax, (-0.02, cor_maxy))
	scatter!(ax, -25:25, rna004_cors_sig, markersize = 7, color = :black)
	lines!(ax, -25:25, rna004_cors_sig, color = :black)

	local ax = Axis(gb[1, 1],
			    title = "RNA002: binned",
			    aspect = 1,
		        xlabel = "Recall",
		        ylabel = "Precision")

	local peakannot_002_binned = copy(peakannot_002)
	peakannot_002_binned[!, :bin] = div.(peakannot_002.genomicPos, BIN_SIZE)
	peakannot_002_binned = combine(groupby(peakannot_002_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_002_p10_binned = copy(peakannot_002_p10)
	peakannot_002_p10_binned[!, :bin] = div.(peakannot_002_p10.genomicPos, BIN_SIZE)
	peakannot_002_p10_binned = combine(groupby(peakannot_002_p10_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	plot_prc!(ax,
			 peakannot_002_binned.modified,
		 	 peakannot_002_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_binned.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_002_p10_binned.modified,
		 	 peakannot_002_p10_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_p10_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax)

	local ax = Axis(gb[1, 2],
			    title = "RNA004: binned",
			    aspect = 1,
		        xlabel = "Recall",
		        ylabel = "Precision")
	
	local peakannot_004_binned = copy(peakannot_004)
	peakannot_004_binned[!, :bin] = div.(peakannot_004.genomicPos, BIN_SIZE)
	peakannot_004_binned = combine(groupby(peakannot_004_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_004_p12_binned = copy(peakannot_004_p12)
	peakannot_004_p12_binned[!, :bin] = div.(peakannot_004_p12.genomicPos, BIN_SIZE)
	peakannot_004_p12_binned = combine(groupby(peakannot_004_p12_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	plot_prc!(ax,
			 peakannot_004_binned.modified,
		 	 peakannot_004_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_binned.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_004_p12_binned.modified,
		 	 peakannot_004_p12_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_p12_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax)
	
	# local a, p, r = auprc(rna002_annot_peaks.predicted,
	#                 	  Int.(rna002_annot_peaks.modified),
	#                 	  Set([1]))

	# local rocs = roc(rna002_annot_peaks.modified,
	# 		     	 rna002_annot_peaks.predicted,
	# 				 nrow(rna002_annot_peaks))
	# local precisions = vcat([0], map(precision, rocs), [1])
	# local recalls = vcat([1], map(recall, rocs), [0])
	# local area = auc(recalls, precisions)
	
	# lines!(ax, recalls, precisions,
	#        label = "Base\nAUC=$(round(area; digits = 2))",
	#        color = :black)
	# local evaluation = roc(rna002_annot_peaks.modified, rna002_annot_peaks.predicted .>= -log10(0.01))
	# local p001 = precision(evaluation)
	# local r001 = recall(evaluation)
	# scatter!(ax, [r001], [p001], color = :black)
	
	local ax = Axis(gb[2, 1],
				    title = "RNA002: single-base resolution",
				    aspect = 1,
		        	xlabel = "Recall",
		        	ylabel = "Precision")
	plot_prc!(ax,
			 peakannot_002.modified,
		 	 peakannot_002.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002.predicted));
			 linestyle = :solid,
				 label = "Default")
	plot_prc!(ax,
			 peakannot_002_p10.modified,
		 	 peakannot_002_p10.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_p10.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax)

	local ax = Axis(gb[2, 2],
				    title = "RNA004: single-base resolution",
				    aspect = 1,
		        	xlabel = "Recall",
		        	ylabel = "Precision")
	plot_prc!(ax,
			 peakannot_004.modified,
		 	 peakannot_004.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004.predicted));
			 linestyle = :solid,
			 label = "Default")
	plot_prc!(ax,
			 peakannot_004_p12.modified,
		 	 peakannot_004_p12.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_p12.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = :red)

	axislegend(ax)

	plot_peak_dist_hist!(f, ga[3, 1], rna002_tx, ref002, peakannot_002, peakannot_002_p10; title = "RNA002: closest peak distance", radius = 10)
	plot_peak_dist_hist!(f, ga[3, 2], rna004_tx, ref004, peakannot_004, peakannot_004_p12; title = "RNA004: closest peak distance", radius = 10)

	# colsize!(f.layout, 1, Fixed(640))
	# colsize!(f.layout, 2, Fixed(640))
	# colsize!(f.layout, 3, Fixed(450))
	# colsize!(f.layout, 4, Fixed(450))

	
	# local a, p, r = auprc(rna002_p10_annot_peaks.predicted,
	#                 	  Int.(rna002_p10_annot_peaks.modified),
	#                 	  Set([1]))

	# local rocs = roc(rna002_p10_annot_peaks.modified,
	# 		     	 rna002_p10_annot_peaks.predicted,
	# 				 nrow(rna002_p10_annot_peaks))
	# local precisions = vcat([0], map(precision, rocs), [1])
	# local recalls = vcat([1], map(recall, rocs), [0])
	# local area = auc(recalls, precisions)
	
	# lines!(ax, recalls, precisions,
	#        label = "Motor\nAUC=$(round(area; digits = 2))",
	#        color = :red,
	# 	   linestyle = :dash)
	# local evaluation = roc(rna002_p10_annot_peaks.modified, rna002_p10_annot_peaks.predicted .>= -log10(0.01))
	# local p001 = precision(evaluation)
	# local r001 = recall(evaluation)
	# scatter!(ax, [r001], [p001], color = :red)
	# axislegend(ax)


	# local ax = Axis(f[2, 4],
	# 		    title = "RNA004: motor dwell effect",
	# 		    aspect = 1,
	# 	        xlabel = "Recall",
	# 	        ylabel = "Precision")
	
	# # local a, p, r = auprc(rna004.predicted,
	# #                 	  Int.(rna004.modified),
	# #                 	  Set([1]))
	# local rocs = roc(rna004_annot_peaks.modified,
	# 		     	 rna004_annot_peaks.predicted,
	# 				 nrow(rna004_annot_peaks))
	# local precisions = vcat([0], map(precision, rocs), [1])
	# local recalls = vcat([1], map(recall, rocs), [0])
	# local area = auc(recalls, precisions)
	
	# lines!(ax, recalls, precisions,
	#        label = "Base\nAUC=$(round(area; digits = 2))",
	#        color = :black)
	# local evaluation = roc(rna004_annot_peaks.modified, rna004_annot_peaks.predicted .>= -log10(0.01))
	# local p001 = precision(evaluation)
	# local r001 = recall(evaluation)
	# scatter!(ax, [r001], [p001], color = :black)
	
	# local a, p, r = auprc(rna004_p12.predicted,
	#                 	  Int.(rna004_p12.modified),
	#                 	  Set([1]))
	# local rocs = roc(rna004_p12_annot_peaks.modified,
	# 		     	 rna004_p12_annot_peaks.predicted,
	# 				 nrow(rna004_p12_annot_peaks))
	# local precisions = vcat([0], map(precision, rocs), [1])
	# local recalls = vcat([1], map(recall, rocs), [0])
	# local area = auc(recalls, precisions)
	
	# lines!(ax, recalls, precisions,
	#        label = "Motor\nAUC=$(round(area; digits = 2))",
	#        color = :red,
	# 	   linestyle = :dash)
	# local evaluation = roc(rna004_p12_annot_peaks.modified, rna004_p12_annot_peaks.predicted .>= -log10(0.01))
	# local p001 = precision(evaluation)
	# local r001 = recall(evaluation)
	# scatter!(ax, [r001], [p001], color = :red)
	# axislegend(ax)


	# local df = innerjoin(innerjoin(ivt, ivt_p12,
	# 							   on = [:ref_id, :pos],
	# 							   makeunique = true),
	# 					 mod_ref,
	# 					 on = [:chr, :strand, :genomicPos => :pos])

	# local pseUs = df[df.mod .== "Y", :]

	# local ax = Axis(f[2, 4],
	# 				title = "RNA004: Reference pseU sites",
	# 				xlabel = "|LORmotor| - |LOR|",
	# 			    ylabel = "-log10(𝑄-𝑣𝑎𝑙𝑢𝑒 motor) - -log10(𝑄-𝑣𝑎𝑙𝑢𝑒)",
	# 			    xticks = -3:1:3)

	# local color = map(r -> if r.predicted >= 2 && r.predicted_1 >= 2
	# 					       :black
	# 					   elseif r.predicted <= 2 && r.predicted_1 >= 2
	# 						   :orange
	# 					   elseif r.predicted >= 2 && r.predicted_1 <= 2
	# 						   :blue
	# 					   else
	# 						   :grey
	# 					   end, eachrow(pseUs))
	
	# scatter!(ax, abs.(pseUs.GMM_LOR_1) .- abs.(pseUs.GMM_LOR),
	# 		 -log10.(pseUs.GMM_chi2_pvalue_1) .- -log10.(pseUs.GMM_chi2_pvalue),
	# 		 markersize = 6,
	# 		 color = color,
	# 		 alpha = 0.75)

	# local colorcounts = countmap(color)

	# Legend(f[1, 4],
	# 	   [PolyElement(color = :orange),
	# 		PolyElement(color = :blue),
	# 		PolyElement(color = :black),
	# 		PolyElement(color = :grey)],
	#        ["Becomes significant ($(colorcounts[:orange]))", "Becomes insignificant ($(colorcounts[:blue]))", "Remains significant ($(colorcounts[:black]))", "Remains insignificant ($(colorcounts[:grey]))"],
	# 	   orientation = :vertical,
	# 	   framevisible = false)

	f
	# data(mean_xcorrs_002) * mapping(:offset => "Offset", :mean_xcorr => "Mean correlation") * visual(Scatter) |> draw
end
  ╠═╡ =#

# ╔═╡ 552ac57a-fa1b-4de5-b686-8dbb56952a6e
inch = 96

# ╔═╡ 88663e95-2291-408c-9114-22cd225309a0
pt = 4/3

# ╔═╡ 4ab1f43e-bc3f-4865-8638-f53426ebdb7f
cm = inch / 2.54

# ╔═╡ 99c5d277-0490-44bf-83a4-85e34d8daec6
Makie.to_font("Blackchancery")


# ╔═╡ 3eb319a8-f26a-4508-a749-2d539c021600
L"\textbf{KS(intensity) and KS(dwell) correlation}"

# ╔═╡ 1069343c-a14b-42ae-bb1d-a2c5ae2d06c5
L"\textbf{RNA002: KS}_{\large{intensity}}-\textbf{KS}_{dwell}\textbf{ correlation}" 

# ╔═╡ 64ed1f82-51f0-448f-8e54-a8565898774c
begin
	local COL_ALL_SITES = :black # "#7f7f7f22"
	local COL_MOD_SITES = "#00A0A9" #"#FFD700"
	local COL_2D = "#004488"  # "#3E8241" # "#17becf"
	local COL_3D = "#D75F00" # "#6948A3" # "#bcbd22"
	
	local f = Figure(size = (1150, 1850),
					 # fonts = (; regular = "fonts/cmunorm.ttf", bold = "fonts/cmunso.ttf"),
					 fontsize = 26)
	# local f = Figure(size = (9cm, 7cm), fontsize = 12pt)
	# local ga = f[1, 1] = GridLayout()
	# local gb = f[1, 2] = GridLayout()

	local cor_maxy = 0.21
	local markersize = 9

	local ax = Axis(f[1, 1],
				    title="RNA002",
				    titlesize=34)
	hidedecorations!(ax)
	hidespines!(ax)
	local ax = Axis(f[1, 2],
				    title="RNA004",
				    titlesize=34)
	hidedecorations!(ax)
	hidespines!(ax)


	local ax = Axis(f[2, 1],
					# title = L"\textbf{KS}_{intensity}-\textbf{KS}_{dwell}\textbf{ correlation}",
					titlesize=20,
					aspect=1,
				    xlabel = "Offset",
				    ylabel = "Correlation",
				    yticks = 0:0.02:cor_maxy)
	ylims!(ax, (-0.01, cor_maxy))

	local mean_xcorrs_002_ = mean_xcorrs_002[abs.(mean_xcorrs_002.offset) .< 21, :]
	local rna002_cors_sig_ = rna002_cors_sig[6:46]
	scatter!(ax, mean_xcorrs_002_.offset, mean_xcorrs_002_.mean_xcorr, markersize = markersize, color = COL_ALL_SITES, label = "All sites")
	lines!(ax, mean_xcorrs_002_.offset, mean_xcorrs_002_.mean_xcorr, color = COL_ALL_SITES, linewidth = 2.5)

	scatter!(ax, -20:20, rna002_cors_sig_, markersize = markersize, color = COL_MOD_SITES, label = "Modified sites")
	lines!(ax, -20:20, rna002_cors_sig_, color = COL_MOD_SITES, linestyle = :dash, linewidth = 2.5)

	# axislegend(ax)
	local legend_elem_all = [LineElement(color = COL_ALL_SITES, linestyle = :solid),
          					 MarkerElement(color = COL_ALL_SITES, marker = :circle, 								 markersize = 15,
         				     strokecolor = COL_ALL_SITES)]
	local legend_elem_mod = [LineElement(color = COL_MOD_SITES, linestyle = :dash),
          					 MarkerElement(color = COL_MOD_SITES, marker = :circle, 								 markersize = 15,
         				     strokecolor = COL_MOD_SITES)]
	Legend(f[2, 1],
		   [legend_elem_all, legend_elem_mod],
		   ["All sites", "Modified sites"],
		   patchsize = (40, 40),
		   rowgap = 5,
		   tellwidth = false,
		   halign = :right,
		   valign = :top,
		   margin = (10, 10, 10, 10))
	
	local ax = Axis(f[2, 2],
					# title = L"\textbf{KS}_{intensity}-\textbf{KS}_{dwell}\textbf{ correlation}",
					titlesize=20,
					aspect=1,
				    xlabel = "Offset",
				    ylabel = "Correlation",
				    yticks = 0:0.02:cor_maxy)
	ylims!(ax, (-0.01, cor_maxy))
	
	local mean_xcorrs_ = mean_xcorrs[abs.(mean_xcorrs.offset) .< 21, :]
	local rna004_cors_sig_ = rna004_cors_sig[6:46]
	scatter!(ax, mean_xcorrs_.offset, mean_xcorrs_.mean_xcorr, markersize = markersize, color = COL_ALL_SITES, label = "All sites")
	lines!(ax, mean_xcorrs_.offset, mean_xcorrs_.mean_xcorr, color = COL_ALL_SITES, linewidth = 2.5)
	
	scatter!(ax, -20:20, rna004_cors_sig_, markersize = markersize, color = COL_MOD_SITES, label = "Modified sites")
	lines!(ax, -20:20, rna004_cors_sig_, color = COL_MOD_SITES, linestyle = :dash, linewidth = 2.5)

	# axislegend(ax)
	Legend(f[2, 2],
		   [legend_elem_all, legend_elem_mod],
		   ["All sites", "Modified sites"],
		   patchsize = (40, 40),
		   rowgap = 5,
		   tellwidth = false,
		   halign = :right,
		   valign = :top,
		   margin = (10, 10, 10, 10))
	
	
	local ax = Axis(f[3, 1],
				    # title = L"\textbf{Single base resolution}",
					titlesize=20,
				    aspect = 1,
		        	xlabel = "Recall",
		        	ylabel = "Precision")
	plot_prc!(ax,
			  peakannot_002.modified,
		 	  peakannot_002.predicted,
		 	  sort(filter(v -> v > 0, peakannot_002.predicted));
			  linestyle = :solid,
			  linewidth = 2.5,
			  color = COL_2D,
			  label = "MD−")
	plot_prc!(ax,
			  peakannot_002_p10.modified,
		 	  peakannot_002_p10.predicted,
		 	  sort(filter(v -> v > 0, peakannot_002_p10.predicted));
			  label = "MD+",
			  linestyle = :dash,
			  linewidth = 2.5,
			  color = COL_3D)

	axislegend(ax, "GMM test", patchsize = (40, 40))

	local ax = Axis(f[3, 2],
				    # title = L"\textbf{Single base resolution}",
					titlesize=20,
				    aspect = 1,
		        	xlabel = "Recall",
		        	ylabel = "Precision")
	plot_prc!(ax,
			  peakannot_004.modified,
		 	  peakannot_004.predicted,
		 	  sort(filter(v -> v > 0, peakannot_004.predicted));
			  linestyle = :solid,
			  linewidth = 2.5,
			  color = COL_2D,
			  label = "MD−")
	plot_prc!(ax,
			  peakannot_004_p12.modified,
		 	  peakannot_004_p12.predicted,
		 	  sort(filter(v -> v > 0, peakannot_004_p12.predicted));
			  label = "MD+",
			  linestyle = :dash,
			  linewidth = 2.5,
			  color = COL_3D)

	axislegend(ax, "GMM test", patchsize = (40, 40))

	# plot_peak_dist_hist!(f, f[4, 1], rna002_tx, ref002, peakannot_002, peakannot_002_p10; title = L"\textbf{Closest peak distance}", radius = 5, col_palette = [COL_2D, COL_3D], titlesize = 28, linewidth = 2.5, maxy=0.85)
	# plot_peak_dist_hist!(f, f[4, 2], rna004_tx, ref004, peakannot_004, peakannot_004_p12; title = L"\textbf{Closest peak distance}", radius = 5, col_palette = [COL_2D, COL_3D], titlesize = 28, linewidth = 2.5, maxy=0.85)
	plot_peak_dist_hist!(f, f[4, 1], rna002_tx, ref002, peakannot_002, peakannot_002_p10; radius = 5, col_palette = [COL_2D, COL_3D], titlesize = 28, linewidth = 2.5, maxy=0.85)
	plot_peak_dist_hist!(f, f[4, 2], rna004_tx, ref004, peakannot_004, peakannot_004_p12; radius = 5, col_palette = [COL_2D, COL_3D], titlesize = 28, linewidth = 2.5, maxy=0.85)

	# colsize!(f.layout, 1, Relative(0.4))
	# colsize!(f.layout, 2, Relative(0.25))
	# colsize!(f.layout, 3, Relative(0.35))
	# colsize!(f.layout, 4, Fixed(450))
	rowsize!(f.layout, 1, Fixed(20))

	for (label, layout) in zip(["A", "B", "C", "D", "E", "F"],
							   [f[2, 1], f[2, 2], f[3, 1], f[3, 2], f[4, 1], f[4, 2]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 26,
	        font = :bold,
	        padding = (50, 30, 5, 0),
	        halign = :right)
	end

	f
end

# ╔═╡ a5d05c3c-6dfd-414e-be41-5d032ece5561
begin
	local COL_2D = "#004488"  # "#3E8241" # "#17becf"
	local COL_3D = "#D75F00" # "#6948A3" # "#bcbd22"

	local f = Figure(size = (1200, 600))
	
	local peakannot_002_binned = copy(peakannot_002)
	peakannot_002_binned[!, :bin] = div.(peakannot_002.genomicPos, BIN_SIZE)
	peakannot_002_binned = combine(groupby(peakannot_002_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_002_p10_binned = copy(peakannot_002_p10)
	peakannot_002_p10_binned[!, :bin] = div.(peakannot_002_p10.genomicPos, BIN_SIZE)
	peakannot_002_p10_binned = combine(groupby(peakannot_002_p10_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	local ax = Axis(f[1, 1],
				    title = L"\textbf{RNA002: binned}",
				    aspect = 1)
	plot_prc!(ax,
			 peakannot_002_binned.modified,
		 	 peakannot_002_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_binned.predicted));
			 linestyle = :solid,
			 label = "Default",
			 color = COL_2D)
	plot_prc!(ax,
			 peakannot_002_p10_binned.modified,
		 	 peakannot_002_p10_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_002_p10_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = COL_3D)

	axislegend(ax, L"\textbf{GMM test}")

	local peakannot_004_binned = copy(peakannot_004)
	peakannot_004_binned[!, :bin] = div.(peakannot_004.genomicPos, BIN_SIZE)
	peakannot_004_binned = combine(groupby(peakannot_004_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	local peakannot_004_p12_binned = copy(peakannot_004_p12)
	peakannot_004_p12_binned[!, :bin] = div.(peakannot_004_p12.genomicPos, BIN_SIZE)
	peakannot_004_p12_binned = combine(groupby(peakannot_004_p12_binned, [:chr, :strand, :bin]),
								   :modified => maximum => :modified,
								   :predicted => maximum =>  :predicted)
	
	local ax = Axis(f[1, 2],
				    title = L"\textbf{RNA004: binned}",
				    aspect = 1)
	plot_prc!(ax,
			 peakannot_004_binned.modified,
		 	 peakannot_004_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_binned.predicted));
			 linestyle = :solid,
			 label = "Default",
			 color = COL_2D)
	plot_prc!(ax,
			 peakannot_004_p12_binned.modified,
		 	 peakannot_004_p12_binned.predicted,
		 	 sort(filter(v -> v > 0, peakannot_004_p12_binned.predicted));
			 label = "Motor",
			 linestyle = :dash,
			 color = COL_3D)

	axislegend(ax, L"\textbf{GMM test}")
	
	f
end

# ╔═╡ ce755b22-4dc5-4734-be8d-bf50264cf9aa
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
					title = "RNA002 modified sites",
				 	xlabel = "Offset",
				    ylabel = "Correlation")
	local cors = []
	local t = copy(rna002_tx)
	for i in -25:25
		t[!, "motor_pos"] = disallowmissing(t.pos) .+ i
		tj = innerjoin(t[t.KS_intensity_qvalue .<= 0.01, :], t,
					   on=[:ref_id => :ref_id,
					   	   :motor_pos => :pos],
					   makeunique=true)
		push!(cors, cor(tj.KS_intensity_pvalue,
						tj.KS_dwell_pvalue_1))
	end
	scatter!(ax, -25:25, cors)
	f
end

# ╔═╡ fc756d1b-fbe7-4ecb-bde7-66ee3eda7dc7
begin
	local df2d = leftjoin(rna002_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df2d[!, :modified] = coalesce.(df2d.modified, 0)

	local df3d = leftjoin(rna002_p10_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df3d[!, :modified] = coalesce.(df3d.modified, 0)

	innerjoin(df3d[df3d.modified .== 1 .&& df3d.predicted .> 2, :],
			  df2d[df2d.modified .== 1 .&& df2d.predicted .< 2, :],
			  on = [:ref_id, :pos],
			  makeunique = true)
end

# ╔═╡ 8f5527b9-4bea-462d-bce9-e66a6432fcc4
begin
	local ref = "ENST00000544848.3|ENSG00000246705.5|OTTHUMG00000168736.2|OTTHUMT00000400845.1|H2AJ-203|H2AJ|620|protein_coding|"

	local pos = 162
	local radius = 25

	local df2d = leftjoin(rna002_tx[rna002_tx.ref_id .== ref, :], ref002,
					      on = [:chr, :strand, :genomicPos])
	df2d[!, :modified] = coalesce.(df2d.modified, 0)

	local df3d = leftjoin(rna002_p10_tx[rna002_p10_tx.ref_id .== ref, :], ref002,
					      on = [:chr, :strand, :genomicPos])
	df3d[!, :modified] = coalesce.(df3d.modified, 0)

	df2d = df2d[df2d.pos .> pos - radius .&& df2d.pos .< pos + radius , [:pos, :modified, :GMM_chi2_pvalue]]
	df2d[!, :type] .= "2D"
	df3d = df3d[df3d.pos .> pos - radius .&& df3d.pos .< pos + radius , [:pos, :modified, :GMM_chi2_pvalue]]
	df3d[!, :type] .= "3D"

	local df = vcat(df2d, df3d)
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :log10_pval] = -log10.(df.GMM_chi2_pvalue)

	(
		data(df) * mapping(:pos => "Position", :log10_pval => "-log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", color = :type) * (visual(Lines) + visual(Scatter, colormap = [:black, :red]))
		+
	 	data(DataFrame(x = [162, 162], y = [-0, 20])) * mapping(:x, :y) * visual(Lines, color = :red, linestyle = :dash)
	)|> draw
end

# ╔═╡ d04d0bc3-c8ee-4894-8f3b-e7d383fdd834
begin
	local radius = 20

	local df2d = leftjoin(rna002_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df2d[!, :modified] = coalesce.(df2d.modified, 0)
	local mods = unique(df2d[df2d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df2d = leftjoin(ranges, df2d,
			 		on = [:ref_id, :range => :pos])
	df2d[!, :offset] = df2d.range .- df2d.pos
	df2d[!, :type] .= "2D"


	local df3d = leftjoin(rna002_p10_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df3d[!, :modified] = coalesce.(df3d.modified, 0)
	local mods = unique(df3d[df3d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df3d = leftjoin(ranges, df3d,
			 		on = [:ref_id, :range => :pos])
	df3d[!, :offset] = df3d.range .- df3d.pos
	df3d[!, :type] .= "3D"
	

	local df = vcat(df2d[:, [:offset, :GMM_chi2_pvalue, :type]],
					df3d[:, [:offset, :GMM_chi2_pvalue, :type]])
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :log10_pval] = -log10.(df.GMM_chi2_pvalue .+ eps())

	df = combine(groupby(df, [:offset, :type]),
				 :log10_pval => mean => :log10_pval)
	colors = ["#FC7808", "#8C00EC"]
	data(df) * mapping(:offset => "Offset", :log10_pval => "Mean -log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", color = :type)  * (visual(Lines) + visual(Scatter)) |> draw(scales(Color = (; palette = [:tomato, :teal])))
end

# ╔═╡ f440817e-5e8e-4a59-b4fb-02b32b331463
begin
	local radius = 20

	local df2d = leftjoin(rna002_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df2d[!, :modified] = coalesce.(df2d.modified, 0)
	local mods = unique(df2d[df2d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df2d = leftjoin(ranges, df2d,
			 		on = [:ref_id, :range => :pos])
	df2d[!, :offset] = df2d.range .- df2d.pos
	df2d[!, :type] .= "2D"


	local df3d = leftjoin(rna002_p10_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df3d[!, :modified] = coalesce.(df3d.modified, 0)
	local mods = unique(df3d[df3d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df3d = leftjoin(ranges, df3d,
			 		on = [:ref_id, :range => :pos])
	df3d[!, :offset] = df3d.range .- df3d.pos
	df3d[!, :type] .= "3D"
	

	local df = vcat(df2d[:, [:offset, :predicted, :type]],
					df3d[:, [:offset, :predicted, :type]])
	# df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :predicted] = coalesce.(df.predicted, 0) .> 2
	# df[!, :log10_pval] = -log10.(df.GMM_chi2_pvalue .+ eps())

	df = combine(groupby(df, [:offset, :type]),
				 :predicted => mean => :predicted)
	
	data(df) * mapping(:offset => "Offset", :predicted => "Mean significant", color = :type)  * (visual(Lines) + visual(Scatter)) |> draw
end

# ╔═╡ 81a1bc5f-450e-4b66-9cfb-7ebc9459e65f
confusmat(2, peakannot_002.modified .+ 1, Int.(peakannot_002.predicted .> 2) .+ 1)

# ╔═╡ 473afb06-f4fa-463b-a332-39ed43f76c02
10174/(10174+902)

# ╔═╡ 4e6873d1-dda9-4ec4-a0e4-8ce452214c02
confusmat(2, peakannot_002_p10.modified .+ 1, Int.(peakannot_002_p10.predicted .> 2) .+ 1)

# ╔═╡ e6a6fcd7-0625-4dda-b2e5-5f3889c39ec0
9206/(9206+868)

# ╔═╡ f99cdd17-0c1d-41e8-9838-05093cdd071f
sum(rna002_tx.GMM_chi2_pvalue .=== missing)/nrow(rna002_tx), sum(rna002_p10_tx.GMM_chi2_pvalue .=== missing)/nrow(rna002_p10_tx)

# ╔═╡ b16ede32-1f2e-4d85-9ead-453365dab4f2
begin
	local radius = 20

	local df2d = leftjoin(rna002_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df2d[!, :modified] = coalesce.(df2d.modified, 0)
	local mods = unique(df2d[df2d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df2d = leftjoin(ranges, df2d,
			 		on = [:ref_id, :range => :pos])
	df2d[!, :offset] = df2d.range .- df2d.pos
	df2d[!, :type] .= "2D"


	local df3d = leftjoin(rna002_p10_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df3d[!, :modified] = coalesce.(df3d.modified, 0)
	local mods = unique(df3d[df3d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df3d = leftjoin(ranges, df3d,
			 		on = [:ref_id, :range => :pos])
	df3d[!, :offset] = df3d.range .- df3d.pos
	df3d[!, :type] .= "3D"

	local df = innerjoin(df2d, df3d,
			  on = [:ref_id, :pos, :offset],
			  makeunique = true)

	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1-eps())
	df[!, :log10_pval] = -log10.(clamp.(df.GMM_chi2_pvalue, eps(), 1 - eps()))
	df[!, :GMM_chi2_pvalue_1] = coalesce.(df.GMM_chi2_pvalue_1, 1-eps())
	df[!, :log10_pval_1] = -log10.(clamp.(df.GMM_chi2_pvalue_1, eps(), 1-eps()))	

	df[!, :ratio] = log2.(df.log10_pval_1 ./ df.log10_pval)

	df = combine(groupby(df, [:offset]),
				 :ratio => mean => :ratio)
	

	# local df = vcat(df2d[:, [:offset, :GMM_chi2_pvalue, :type]],
	# 				df3d[:, [:offset, :GMM_chi2_pvalue, :type]])
	# df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	# df[!, :log10_pval] = -log10.(df.GMM_chi2_pvalue .+ eps())

	# df = combine(groupby(df, [:offset, :type]),
	# 			 :log10_pval => mean => :log10_pval)
	
	data(df) * mapping(:offset => "Offset", :ratio => "Mean ratio")  * (visual(Lines) + visual(Scatter)) |> draw
end

# ╔═╡ ba062430-7a52-4a49-94b2-b0935a918eff
ConfSets.confint([2,8,9,7,6,5,7,9,0,7,4,4,4,6], 0.05, "clt")

# ╔═╡ d0ebdaf8-06e3-4922-9a40-b11c345e71b5
begin
	local radius = 20

	local df2d = leftjoin(rna002_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df2d[!, :modified] = coalesce.(df2d.modified, 0)
	local mods = unique(df2d[df2d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df2d = leftjoin(ranges, df2d,
			 		on = [:ref_id, :range => :pos])
	df2d[!, :offset] = df2d.range .- df2d.pos
	df2d[!, :Type] .= "2D"


	local df3d = leftjoin(rna002_p10_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df3d[!, :modified] = coalesce.(df3d.modified, 0)
	local mods = unique(df3d[df3d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df3d = leftjoin(ranges, df3d,
			 		on = [:ref_id, :range => :pos])
	df3d[!, :offset] = df3d.range .- df3d.pos
	df3d[!, :Type] .= "3D"
	

	local df = vcat(df2d[:, [:offset, :GMM_chi2_pvalue, :Type]],
					df3d[:, [:offset, :GMM_chi2_pvalue, :Type]])
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :log10_pval] = -log10.(df.GMM_chi2_pvalue .+ eps())

	# df = combine(groupby(df, [:offset, :type]),
	# 			 :log10_pval => mean => :log10_pval,
	# 			 :log10_pval => (v -> ConfSets.confint(collect(v), 0.05, "clt")) => [:lower, :upper])
	
	# data(df) * mapping(:offset => "Offset", :log10_pval => "Mean -log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", color = :type)  * (visual(Scatter)) |> draw

	local plt = data(df) * mapping(:offset => "Offset from m⁶A", :log10_pval => "-log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", color = :Type, dodge = :Type) * visual(BoxPlot, show_outliers = false, dodge_gap = 0.15, gap = 0.7, width = 3)

	local fig = Figure()
	local gridpos = fig[1, 1]

	local f = draw!(gridpos, plt)
	legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))
	
	fig
end

# ╔═╡ 8b71fdc6-a527-4248-948c-c95c408604cb
begin
	local radius = 20

	local df2d = leftjoin(rna002_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df2d[!, :modified] = coalesce.(df2d.modified, 0)
	local mods = unique(df2d[df2d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df2d = leftjoin(ranges, df2d,
			 		on = [:ref_id, :range => :pos])
	df2d[!, :offset] = df2d.range .- df2d.pos
	df2d[!, :Type] .= "2D"

	df2d[!, :GMM_chi2_pvalue] = coalesce.(df2d.GMM_chi2_pvalue, 1)
	df2d[!, :log10_pval] = -log10.(clamp.(df2d.GMM_chi2_pvalue, eps(), 1-eps()))


	local df3d = leftjoin(rna002_p10_tx, ref002,
					      on = [:chr, :strand, :genomicPos])
	df3d[!, :modified] = coalesce.(df3d.modified, 0)
	local mods = unique(df3d[df3d.modified .== 1, [:ref_id, :pos]])
	mods[!, :range] = map(p -> (p-radius):(p+radius), mods.pos)
	local ranges = flatten(mods[:, [:ref_id, :range, :pos]], :range)
	df3d = leftjoin(ranges, df3d,
			 		on = [:ref_id, :range => :pos])
	df3d[!, :offset] = df3d.range .- df3d.pos
	df3d[!, :Type] .= "3D"
	
	df3d[!, :GMM_chi2_pvalue] = coalesce.(df3d.GMM_chi2_pvalue, 1)
	df3d[!, :log10_pval] = -log10.(clamp.(df3d.GMM_chi2_pvalue, eps(), 1-eps()))


	df2d = innerjoin(df2d, df2d[df2d.pos .== df2d.range, :],
					 on = [:ref_id, :pos],
					 makeunique = true)
	df2d[!, :ratio] = log2.(df2d.log10_pval ./ df2d.log10_pval_1)

	df3d = innerjoin(df3d, df3d[df3d.pos .== df3d.range, :],
					 on = [:ref_id, :pos],
					 makeunique = true)
	df3d[!, :ratio] = log2.(df3d.log10_pval ./ df3d.log10_pval_1)


	df2d

	
	local df = vcat(df2d[:, [:offset, :ratio, :Type]],
				    df3d[:, [:offset, :ratio, :Type]])
	

	# # df = combine(groupby(df, [:offset, :type]),
	# # 			 :log10_pval => mean => :log10_pval,
	# # 			 :log10_pval => (v -> ConfSets.confint(collect(v), 0.05, "clt")) => [:lower, :upper])
	
	# # data(df) * mapping(:offset => "Offset", :log10_pval => "Mean -log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", color = :type)  * (visual(Scatter)) |> draw

	local plt = data(df) * mapping(:offset => "Offset from m⁶A", :ratio => "Ratio to m6A prediction", color = :Type, dodge = :Type) * visual(BoxPlot, show_outliers = false, dodge_gap = 0.15, gap = 1.5, width = 2)

	local fig = Figure()
	local gridpos = fig[1, 1]

	local f = draw!(gridpos, plt)
	legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))
	
	fig
end

# ╔═╡ 6d9fb141-e38a-4fe2-822f-58053304ded5
# ╠═╡ disabled = true
#=╠═╡
begin
	local df = innerjoin(rna002_tx, rna002_p10_tx,
			  			 on = [:ref_id, :pos],
						 makeunique = true)

	df = innerjoin(df, ref002[:, [:chr, :strand, :genomicPos, :modified]],
				   on = [:chr, :strand, :genomicPos])
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :GMM_chi2_pvalue_1] = coalesce.(df.GMM_chi2_pvalue_1, 1)
	df[!, :log10_pval_2d] = -log10.(clamp.(df.GMM_chi2_pvalue, 1e-306, 1))
	df[!, :log10_pval_3d] = -log10.(clamp.(df.GMM_chi2_pvalue_1, 1e-306, 1))
	df = df[df.GMM_chi2_pvalue .!= 1 .|| df.GMM_chi2_pvalue .!= 1, :]
	data(df) * mapping(:log10_pval_2d, :log10_pval_3d, color = :modified) * visual(Scatter, markersize = 5) |> draw
end
  ╠═╡ =#

# ╔═╡ a247b841-ecbc-4bee-a1c8-ec86099ac9dd
# ╠═╡ disabled = true
#=╠═╡
begin
	local df = innerjoin(rna002_tx, rna002_p10_tx,
			  			 on = [:ref_id, :pos],
						 makeunique = true)

	df = innerjoin(df, ref002[:, [:chr, :strand, :genomicPos, :modified]],
				   on = [:chr, :strand, :genomicPos])
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :GMM_chi2_pvalue_1] = coalesce.(df.GMM_chi2_pvalue_1, 1)
	df[!, :log10_pval_2d] = -log10.(clamp.(df.GMM_chi2_pvalue, 1e-306, 1))
	df[!, :log10_pval_3d] = -log10.(clamp.(df.GMM_chi2_pvalue_1, 1e-306, 1))
	df = df[df.GMM_chi2_pvalue .!= 1 .|| df.GMM_chi2_pvalue .!= 1, :]
	df[!, :modified] = df.modified .== 1
	data(df) * mapping(:log10_pval_2d, :log10_pval_3d, color = :modified) * linear() |> draw
end
  ╠═╡ =#

# ╔═╡ 8efc4f06-d436-45b1-955c-914bf599eb3c
# ╠═╡ disabled = true
#=╠═╡
begin
	local df = innerjoin(rna002_tx, rna002_p10_tx,
			  			 on = [:ref_id, :pos],
						 makeunique = true)

	df = innerjoin(df, ref002[:, [:chr, :strand, :genomicPos, :modified]],
				   on = [:chr, :strand, :genomicPos])
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :GMM_chi2_pvalue_1] = coalesce.(df.GMM_chi2_pvalue_1, 1)
	df[!, :log10_pval_2d] = -log10.(clamp.(df.GMM_chi2_pvalue, 1e-306, 1))
	df[!, :log10_pval_3d] = -log10.(clamp.(df.GMM_chi2_pvalue_1, 1e-306, 1))
	df = df[df.GMM_chi2_pvalue .!= 1 .|| df.GMM_chi2_pvalue .!= 1, :]
	df[!, :modified] = df.modified .== 1
	data(df) * mapping(:log10_pval_2d, :log10_pval_3d, color = :modified) * smooth(npoints = 2000) |> draw
end
  ╠═╡ =#

# ╔═╡ f2dbafe4-b755-4c8c-a2dc-d19272b93c86
# ╠═╡ disabled = true
#=╠═╡
begin
	local df = innerjoin(rna002_tx, rna002_p10_tx,
			  			 on = [:ref_id, :pos],
						 makeunique = true)

	df = innerjoin(df, ref002[:, [:chr, :strand, :genomicPos, :modified]],
				   on = [:chr, :strand, :genomicPos])
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :GMM_chi2_pvalue_1] = coalesce.(df.GMM_chi2_pvalue_1, 1)
	df[!, :log10_pval_2d] = -log10.(clamp.(df.GMM_chi2_pvalue, 1e-306, 1))
	df[!, :log10_pval_3d] = -log10.(clamp.(df.GMM_chi2_pvalue_1, 1e-306, 1))
	df = df[df.GMM_chi2_pvalue .!= 1 .|| df.GMM_chi2_pvalue .!= 1, :]
	df[!, :modified] = df.modified .== 1
	local map_legend_labels = modified -> modified ? "GLORI+" : "GLORI-"
	local plt = data(df) * mapping(:log10_pval_2d => "2D: -log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", :log10_pval_3d => "3D: -log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", color = :modified => map_legend_labels => "Site") * (visual(Scatter, markersize = 5, alpha = 0.5) + smooth(npoints = 2000))
	# |> draw(; axis = (; title = "RNA002"))

	local fig = Figure(size = (500, 500))
	local gridpos = fig[1, 1]

	local f = draw!(gridpos, plt; axis = (; title = "RNA002", aspect = 1))
	legend!(gridpos, f; tellwidth=false, halign=:left, valign=:top, margin=(10, 10, 10, 10))
	
	fig
end
  ╠═╡ =#

# ╔═╡ 3de9eae6-3865-48e6-8d8e-027c54024b02
# ╠═╡ disabled = true
#=╠═╡
begin
	local df = innerjoin(rna004_tx, rna004_p12_tx,
			  			 on = [:ref_id, :pos],
						 makeunique = true)

	df = innerjoin(df, ref004[:, [:chr, :strand, :genomicPos, :modified]],
				   on = [:chr, :strand, :genomicPos])
	df[!, :GMM_chi2_pvalue] = coalesce.(df.GMM_chi2_pvalue, 1)
	df[!, :GMM_chi2_pvalue_1] = coalesce.(df.GMM_chi2_pvalue_1, 1)
	df[!, :log10_pval_2d] = -log10.(clamp.(df.GMM_chi2_pvalue, 1e-306, 1))
	df[!, :log10_pval_3d] = -log10.(clamp.(df.GMM_chi2_pvalue_1, 1e-306, 1))
	df = df[df.GMM_chi2_pvalue .!= 1 .|| df.GMM_chi2_pvalue .!= 1, :]
	df[!, :modified] = df.modified .== 1
	local map_legend_labels = modified -> modified ? "GLORI+" : "GLORI-"
	local plt = data(df) * mapping(:log10_pval_2d => "2D: -log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", :log10_pval_3d => "3D: -log₁₀(GMM 𝑃-𝑣𝑎𝑙𝑢𝑒)", color = :modified => map_legend_labels => "Site") * (visual(Scatter, markersize = 5, alpha = 0.5) + smooth(npoints = 2000))
	# |> draw(; axis = (; title = "RNA002"))

	local fig = Figure(size = (500, 500))
	local gridpos = fig[1, 1]

	local f = draw!(gridpos, plt; axis = (; title = "RNA004", aspect = 1))
	legend!(gridpos, f; tellwidth=false, halign=:left, valign=:top, margin=(10, 10, 10, 10))
	
	fig
end
  ╠═╡ =#

# ╔═╡ 10a23508-008e-40a5-843b-bd4a1ae849c8
# ╠═╡ disabled = true
#=╠═╡
begin
	local radius = 15
	local ref = innerjoin(ref002, rna002_tx[:, [:chr, :strand, :genomicPos]],
						  on = [:chr, :strand, :genomicPos])
	local mods = unique(ref[ref.modified .== 1, :])

	local peaks2d = peakannot_002[peakannot_002.predicted .> 2, :]

	local distances2d = map(r -> begin
			distances = peaks2d[peaks2d.chr .== r.chr .&& peaks2d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	local peaks3d = peakannot_002_p10[peakannot_002_p10.predicted .> 2, :]

	local distances3d = map(r -> begin
			distances = peaks3d[peaks3d.chr .== r.chr .&& peaks3d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	distances2d = distances2d[abs.(distances2d) .<= radius]
	distances3d = distances3d[abs.(distances3d) .<= radius]

	local distances = vcat(DataFrame(distance = distances2d, type = "2D"),
						   DataFrame(distance = distances3d, type = "3D"))

	data(distances) * mapping(:distance, color = :type, dodge = :type) * histogram() |> draw

	# mods[!, :offsetpos] = map(p -> (p-radius):(p+radius), mods.genomicPos)
	# local ranges = flatten(mods[:, [:chr, :strand, :genomicPos, :offsetpos]], :offsetpos)
	# local df2d = leftjoin(ranges, peakannot_002,
	# 		 			  on = [:chr, :strand, :offsetpos => :genomicPos])
	# df2d = df2d[df2d.modified .!== missing, :]
	# df2d[!, :predicted] .= 0
	# df2d
end
  ╠═╡ =#

# ╔═╡ f66b0014-35b4-4d3c-bed0-5eded1c27d47
# ╠═╡ disabled = true
#=╠═╡
begin
	local radius = 15
	local ref = innerjoin(ref002, rna002_tx[:, [:chr, :strand, :genomicPos]],
						  on = [:chr, :strand, :genomicPos])
	local mods = unique(ref[ref.modified .== 1, :])

	local peaks2d = peakannot_002[peakannot_002.predicted .> 2, :]

	local distances2d = map(r -> begin
			distances = peaks2d[peaks2d.chr .== r.chr .&& peaks2d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	local peaks3d = peakannot_002_p10[peakannot_002_p10.predicted .> 2, :]

	local distances3d = map(r -> begin
			distances = peaks3d[peaks3d.chr .== r.chr .&& peaks3d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	distances2d = distances2d[abs.(distances2d) .<= radius]
	distances3d = distances3d[abs.(distances3d) .<= radius]

	local distances = vcat(DataFrame(distance = distances2d, type = "2D"),
						   DataFrame(distance = distances3d, type = "3D"))

	
	local plt = data(distances) * mapping(:distance => "Distance to the closest peak", color = :type => "GMM test") * AlgebraOfGraphics.histogram(Stairs; normalization = :pdf, bins = -radius:(radius+1))
	# |> draw(; axis = (; title = "RNA002", width = 800, height = 500, xticks = (-14.5:15.5, [string(t) for t in -15:15])))

	local fig = Figure(size = (800, 500),
					   fontsize = 18,
					   fonts = (; regular = "fonts/cmunorm.ttf", bold = "fonts/cmunso.ttf"),)
	local gridpos = fig[1, 1]

	local f = draw!(gridpos, plt; axis = (; title = "RNA002", xticks = (-14.5:15.5, [string(t) for t in -15:15])))
	legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))
	
	fig
end
  ╠═╡ =#

# ╔═╡ ea7a956e-fc58-4c06-9e32-a3921220bd6a
Makie.to_font("fonts/cmunorm.ttf")

# ╔═╡ 17972770-7786-464e-8bd8-709852c548d9
# ╠═╡ disabled = true
#=╠═╡
begin
	local radius = 15
	local ref = innerjoin(ref004, rna004_tx[:, [:chr, :strand, :genomicPos]],
						  on = [:chr, :strand, :genomicPos])
	local mods = unique(ref[ref.modified .== 1, :])

	local peaks2d = peakannot_004[peakannot_004.predicted .> 2, :]

	local distances2d = map(r -> begin
			distances = peaks2d[peaks2d.chr .== r.chr .&& peaks2d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	local peaks3d = peakannot_004_p12[peakannot_004_p12.predicted .> 2, :]

	local distances3d = map(r -> begin
			distances = peaks3d[peaks3d.chr .== r.chr .&& peaks3d.strand .== r.strand, :genomicPos] .- r.genomicPos
			push!(distances, 1000)
			distances[argmin(abs.(distances))]
		end,
		eachrow(mods))

	distances2d = distances2d[abs.(distances2d) .<= radius]
	distances3d = distances3d[abs.(distances3d) .<= radius]

	local distances = vcat(DataFrame(distance = distances2d, type = "2D"),
						   DataFrame(distance = distances3d, type = "3D"))

	local plt = data(distances) * mapping(:distance => "Distance to the closest peak", color = :type => "GMM test") * AlgebraOfGraphics.histogram(Stairs; normalization = :pdf, bins = -radius:(radius+1))
	# |> draw(; axis = (; width = 800, height = 500, xticks = (-14.5:15.5, [string(t) for t in -15:15])))

	local fig = Figure(size = (800, 500))
	local gridpos = fig[1, 1]

	local f = draw!(gridpos, plt; axis = (; title = "RNA004", 	xticks = (-14.5:15.5, [string(t) for t in -15:15])))
	legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))
	
	fig
end
  ╠═╡ =#

# ╔═╡ 173461ca-2322-4677-9d3a-8d49dbd8aa4e
findmin(abs, [-1, 5, -10])

# ╔═╡ bf32966d-3022-40ab-b1c8-1b874ac9095e
sort(rna002_tx.GMM_chi2_pvalue)

# ╔═╡ Cell order:
# ╠═034a73ca-3f8e-11f0-0aea-77f9aa7561f4
# ╠═8c95186a-1d0c-4c50-943f-80fe43a83fa5
# ╠═bee0055d-8805-47e7-b995-f74d059dd9d4
# ╠═6e0b2772-54ac-4bba-8a26-647c7b32f520
# ╠═c0ccdc56-02cf-46b3-b721-dca81ab3f8e6
# ╠═1e15d868-6cf1-4607-bc5c-567980ed2a9e
# ╠═ac1c85bc-a140-4ccd-812f-4d3885eccf14
# ╠═98a50473-a972-4739-8701-1ad562712dcc
# ╠═5a8124d1-9d98-4639-bcbf-ffe378de6b7e
# ╠═e23c227e-29f2-4757-8b16-4e5c5c07ebd9
# ╠═cd61aa21-9ad5-4dc2-ba26-9685f86186d4
# ╠═5d08f4fc-52af-4e5a-ba7f-de5a3b8ef00b
# ╠═ca0c4038-4104-4fb1-b9e5-53b282ad1f6a
# ╠═72a07892-065a-4781-ac1c-7fd4c60538d2
# ╠═a117edb7-04fc-4b6c-8b44-04a695d23117
# ╠═fa1a7e46-84ab-4c56-a8bb-5f8dab4d05f9
# ╠═29d93abb-89d1-407f-9d1e-cfc2f80853b8
# ╠═952c3e5f-775d-4bf1-b4d6-99c4fd437aa6
# ╠═97574ea9-2c75-4161-9e7e-d3a90f5eac96
# ╠═c4d62836-34ea-4b5b-8eec-40b74ee5c882
# ╠═32d4efb9-0ec6-4499-a47d-0ebc12ec68a1
# ╠═9c116561-0239-497a-abda-7236159f34f0
# ╠═cace108a-ea0a-4815-9590-4618ef8d114c
# ╠═88599c98-1f77-496c-9273-215ac8f25c01
# ╠═336c7723-a9c5-42f0-9f7a-f4381385ddce
# ╠═63ea2754-bd01-4a94-a703-cf7906e72131
# ╠═ca66bb62-13e5-4572-a545-d99c47ec9bb5
# ╠═d9832607-e012-4d6a-a0b8-cdbd615c4b77
# ╠═03e999dd-a88d-40e6-b780-af26f51af72b
# ╠═3ea52b6d-3896-4ac9-a180-180ed0da3d3f
# ╠═252600bf-4ed1-4978-a673-301d14090eea
# ╠═a5a923da-22fc-4798-9265-0614ce977466
# ╠═5ee87cfe-d49d-4173-873c-4303bec3f1e1
# ╠═e88ffa49-1220-4795-98e2-f422e4ab616f
# ╠═05df07a5-ed51-4ddf-9c24-abbb0564a342
# ╠═3a1cecf8-f3e2-48dc-8034-b93ec7ac6033
# ╠═d0f2466a-b371-45fd-8552-82c0ea2b1b9f
# ╠═3b9902ae-8c77-4680-ad2b-3cfde857ec4f
# ╠═4dc3c11f-99ac-4d45-988e-9c6e1122cf63
# ╠═3d04269e-ba19-4c0e-bd41-4c65e06dd9f0
# ╠═eaa5e268-4704-47dc-b5f9-70bc28ccf737
# ╠═fa551c76-a994-4b2d-a299-d87ffc4818f2
# ╠═5d414a5b-633b-44c4-bdbd-300a7ac66629
# ╠═18932f3b-af22-424b-97e1-a9a6411e6d78
# ╠═173b84e8-c0ff-4aa7-8416-3c20fd5273a2
# ╠═67fe3eff-4f4a-4948-8bb9-f74b03ab24e3
# ╠═d5336fe2-b5fb-4c2d-891b-67459a1af7ed
# ╠═9b248593-f585-4031-8f13-b21c50fb2874
# ╠═a7d34031-3843-4b69-8956-de98e6a02bd4
# ╠═cbc3cd90-9c58-4e7f-be58-b5520221f0a9
# ╠═627d783e-7b95-4cb7-a25c-45182d0a54dc
# ╠═c0ce3531-414c-4080-90aa-6f7dd28848dd
# ╠═97e6502b-6a72-4b7e-926a-f58524e94c89
# ╠═589da79b-73c1-446e-9253-c1a020c69a0f
# ╠═76fb3894-6d04-4c97-b7cc-1918cc0be5a0
# ╠═afaa1214-30d1-4cdd-ae35-4e56b5783a01
# ╠═5fc3bfd3-b429-4f6a-aa39-8eca8c3a62a4
# ╠═67db0de7-c7eb-4dce-b06e-9a2832dab5bb
# ╠═cc119196-64be-457d-9515-5efccc36251c
# ╠═a65a2f01-6f6f-4627-8f40-42095a410251
# ╠═32d711b6-9c46-4f10-8b53-2dca3669ac7b
# ╠═1e30107f-49c3-4479-be52-347d4ebb6ac0
# ╟─fc6da8ec-a2b3-4786-ae9a-978c28d66872
# ╠═6dd6aac1-058f-4269-81ad-c5b727229fe3
# ╠═852fc9dd-0ba1-4c95-baea-19846e5403dc
# ╠═5c60b99d-ee4c-4310-b7ac-0c6856141cf8
# ╠═368b97ad-6667-4fd1-89c5-c10890f39bce
# ╠═f23a1b62-7a57-4a58-b830-e3e0a124681a
# ╠═77e469bb-673e-4aa1-8034-1607c738f903
# ╠═784cd890-11af-4656-a014-093223f54919
# ╠═9f03d039-c4de-4d6f-a05e-d1faaaa2dc64
# ╠═1407a831-0d6f-42a4-bb7d-93ee5435098e
# ╠═03414a73-8f17-4c96-bb84-daf43119aafc
# ╠═d3ea3b64-373a-4410-a17e-825deaf7501d
# ╠═7096fcfb-739e-4775-a4f5-1238d97f4a20
# ╠═686af1c3-3e1d-48ca-8efb-e8ba0a818a69
# ╠═bd6ad199-ab74-45ae-b4d4-5d65d03753e1
# ╠═34644635-e88d-47aa-9441-8e7583240bea
# ╠═1c8593c2-a355-4016-a7b5-9dc58570226a
# ╠═f519ebe0-e962-46f7-a7e1-1e9a0a741ddc
# ╠═5cdf3882-4651-415e-9699-481a7dbe0324
# ╠═13abd3a3-2f85-4c23-9383-15d3df668ad4
# ╠═fc50841f-43d9-414b-8fa7-b95b511046d4
# ╠═4a673cb0-f5a5-46a5-8680-6724c51128a6
# ╠═4625bce2-f17f-46b8-a674-4283ae0a3ae5
# ╠═696b8bee-81d7-47d2-8e55-08f490341195
# ╠═9c838f8e-fe3a-4ee1-923b-146a165ff16f
# ╠═d0bd49e7-0e0e-44c0-a7cb-919f9c47ad6e
# ╠═a82dddab-3df0-4c56-b3a2-944058a7228a
# ╠═7016d0c0-3720-4b50-a01e-6ae04cb013f3
# ╠═2f86ca19-0c1b-44e4-92cb-4b269f3fa12b
# ╠═fe258560-13ce-45f8-96a9-431081463201
# ╠═e1d21630-577f-4888-a0b4-3c75f75a6e67
# ╠═98dbea9c-6611-4a53-a53b-c57784208fc6
# ╠═9d3a986e-f40a-4214-9267-dcaff4b013b7
# ╠═081a4b47-79ab-41ec-b688-1a37a5e093a1
# ╠═c6877e60-d3ff-4122-bbc7-d83da7173fce
# ╠═5ae2462c-a387-45cc-af66-3a474952d251
# ╠═5733debc-825e-49a1-8a1b-3c286f71c935
# ╠═0ab9322b-175a-49c6-973a-369f6bb975e7
# ╠═26ae5e1b-29ed-4a13-8a51-a7bf4337a9f1
# ╠═fb4a5c1b-4a4c-41f1-9ed6-69bf589684b3
# ╠═823797fe-09a0-4de2-b3fa-1bf310a77ab1
# ╠═bb89baaf-d5cf-4e35-982a-83d7cab2e87e
# ╠═191652bf-c223-44f5-991e-9bc49e4d2fa2
# ╠═ef3b9161-a00d-480b-823f-6752906f7b85
# ╟─9daa78c5-d8be-4a43-8d43-364e16dec13e
# ╠═498343e1-7d59-4688-9075-55501a7115ca
# ╠═bfcad6a4-aad1-4615-ad1f-ed90935bb8ea
# ╠═5da72439-e4f1-4900-bcf4-9a4e977a64b0
# ╠═84a489f9-080e-422f-8eb7-588b89aafb89
# ╠═e6d9f53d-7ad3-4068-af7d-02ac384f1161
# ╠═ea40bf7e-7626-44f3-bf22-e2a03897dc7b
# ╠═f7b34149-9459-47fc-837b-43a581fc86ed
# ╠═9b2b3288-467b-4c5a-8972-4900b7182a17
# ╠═552ac57a-fa1b-4de5-b686-8dbb56952a6e
# ╠═88663e95-2291-408c-9114-22cd225309a0
# ╠═4ab1f43e-bc3f-4865-8638-f53426ebdb7f
# ╠═99c5d277-0490-44bf-83a4-85e34d8daec6
# ╠═3eb319a8-f26a-4508-a749-2d539c021600
# ╠═1069343c-a14b-42ae-bb1d-a2c5ae2d06c5
# ╠═64ed1f82-51f0-448f-8e54-a8565898774c
# ╠═a5d05c3c-6dfd-414e-be41-5d032ece5561
# ╠═ce755b22-4dc5-4734-be8d-bf50264cf9aa
# ╠═fc756d1b-fbe7-4ecb-bde7-66ee3eda7dc7
# ╠═8f5527b9-4bea-462d-bce9-e66a6432fcc4
# ╠═d04d0bc3-c8ee-4894-8f3b-e7d383fdd834
# ╠═f440817e-5e8e-4a59-b4fb-02b32b331463
# ╠═81a1bc5f-450e-4b66-9cfb-7ebc9459e65f
# ╠═473afb06-f4fa-463b-a332-39ed43f76c02
# ╠═4e6873d1-dda9-4ec4-a0e4-8ce452214c02
# ╠═e6a6fcd7-0625-4dda-b2e5-5f3889c39ec0
# ╠═f99cdd17-0c1d-41e8-9838-05093cdd071f
# ╠═b16ede32-1f2e-4d85-9ead-453365dab4f2
# ╠═f86b30f1-3a4d-4e7e-94f6-67ed90b5b239
# ╠═a3ec7589-3e97-4238-91d3-a4bc286fb486
# ╠═ba062430-7a52-4a49-94b2-b0935a918eff
# ╠═d0ebdaf8-06e3-4922-9a40-b11c345e71b5
# ╠═8b71fdc6-a527-4248-948c-c95c408604cb
# ╠═6d9fb141-e38a-4fe2-822f-58053304ded5
# ╠═a247b841-ecbc-4bee-a1c8-ec86099ac9dd
# ╠═8efc4f06-d436-45b1-955c-914bf599eb3c
# ╠═f2dbafe4-b755-4c8c-a2dc-d19272b93c86
# ╠═3de9eae6-3865-48e6-8d8e-027c54024b02
# ╠═10a23508-008e-40a5-843b-bd4a1ae849c8
# ╠═f66b0014-35b4-4d3c-bed0-5eded1c27d47
# ╠═ea7a956e-fc58-4c06-9e32-a3921220bd6a
# ╠═17972770-7786-464e-8bd8-709852c548d9
# ╠═173461ca-2322-4677-9d3a-8d49dbd8aa4e
# ╠═bf32966d-3022-40ab-b1c8-1b874ac9095e
