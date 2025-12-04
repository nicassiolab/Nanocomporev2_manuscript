### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 78c86cd4-1dce-4afb-bd82-633875018901
begin
	using DataFrames
	using CSV
	using MLBase
	using CairoMakie
	using Random
	using StatsBase
	using Printf
end

# ╔═╡ 71667c1d-cf1a-43a6-b3a6-ffafd87fc069
using AlgebraOfGraphics

# ╔═╡ 11b96b0e-33e1-11f0-23a0-f998c7e38e3f
include("lib.jl")

# ╔═╡ b0e757dd-c8ce-45f6-b448-e8543168c7bf
import MultipleTesting

# ╔═╡ b7f89781-6332-4942-8d5b-7d95e0550b63
LOR_THRESHOLD = 0.8

# ╔═╡ 6608e7a5-fde0-4f0c-9a33-a89c0fe0b858
begin
	# glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/high_confidence_GLORI_sites.bed"))
	# glori[!, "pos"] = glori.start .+ 2
	glori = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/GLORI_intersection_176k_with_ratio.bed"))
	glori[:, :pos] = glori.start
	glori = rename!(glori, :chr => :chrom)
	
	glori[!, "bin"] = map(p -> div(p, BIN_SIZE), glori.pos)
	glori[!, "modified"] .= 1
	glori
end

# ╔═╡ 05ce3923-41da-4ab4-9ae5-860e0531f8f8
begin
	binned_glori = unique(glori[!, ["chrom", "strand", "bin", "modified"]])
	binned_glori = combine(groupby(glori, [:chrom, :strand, :bin, :modified]),
						   :ratio => maximum => :ratio)
	binned_glori = rename!(binned_glori, :chrom => :chr)
end

# ╔═╡ 672662f5-68a9-4c48-a3b8-45a9fb5683d0
begin
	ref002 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/WT_STORM_reference_cov_30_gx_GLORI.tsv"))
	ref002 = combine(groupby(ref002, ["chr", "strand", "genomicPos"]), :modified => maximum => :modified)
end

# ╔═╡ 7d536c8b-d7c7-42cb-a086-2d76b3748380
begin
	ref004 = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/WT_STORM_reference_cov_30_gx_GLORI.tsv"))
	ref004 = combine(groupby(ref004, ["chr", "strand", "genomicPos"]), :modified => maximum => :modified)
end


# ╔═╡ 52c91271-3317-427e-9472-28ff7e9a46c9
rna002 = annotate_results(
	read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
	    genomic_collapse=true; LOR_threshold = LOR_THRESHOLD),
	ref002)

# ╔═╡ ac6bc167-4e00-44df-b32a-397eb8f2300f
binned_rna002 = sort(annotate_binned_results(bin(rna002, BIN_SIZE), 
											 	binned_glori),
					 [:chr, :strand, :bin])

# ╔═╡ d5035c67-6f2e-4135-a9d1-bcd84a0f2c27
rna004 = annotate_results(
	read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_test_v2_release/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
	    genomic_collapse=true; LOR_threshold = LOR_THRESHOLD),
	ref004)

# ╔═╡ f56b4ce6-47d5-48a4-8da8-4e068b2e4da9
binned_rna004 = annotate_binned_results(bin(rna004, BIN_SIZE),
										binned_glori)

# ╔═╡ b6d5e7db-6e65-4120-901d-bcfce89b1537
begin
	rna002_100 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_depth_100_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
		genomic_collapse=true;
		LOR_threshold = LOR_THRESHOLD)
	binned_rna002_100 = annotate_binned_results(bin(rna002_100, BIN_SIZE),
										        binned_glori)
	:ok
end

# ╔═╡ 708e963c-f42c-415e-8d6b-f6836c878d87
begin
	rna002_500 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_depth_500_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
		genomic_collapse=true;
		LOR_threshold = LOR_THRESHOLD)
	binned_rna002_500 = annotate_binned_results(bin(rna002_500, BIN_SIZE),
										    	binned_glori)
	:ok
end

# ╔═╡ c87569b4-4ffd-44c2-9216-5193ddbcc736
begin
	rna002_1000 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_depth_1000_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=2,
		genomic_collapse=true;
		LOR_threshold = LOR_THRESHOLD)
	binned_rna002_1000 = annotate_binned_results(bin(rna002_1000, BIN_SIZE),
				  						     	 binned_glori)
	:ok
end

# ╔═╡ e2cb039e-8dd7-4983-965a-f52b7b718029
begin
	rna004_100 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_100_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=true;
		LOR_threshold = LOR_THRESHOLD)
	binned_rna004_100 = annotate_binned_results(bin(rna004_100, BIN_SIZE),
										        binned_glori)
	:ok
end

# ╔═╡ 749b06ac-a43b-4d1f-ba66-6a9b11ce5038
begin
	rna004_500 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_500_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=true;
		LOR_threshold = LOR_THRESHOLD)
	binned_rna004_500 = annotate_binned_results(bin(rna004_500, BIN_SIZE),
										        binned_glori)
	:ok
end

# ╔═╡ 8367cd69-1b11-4ac4-8b98-a3de8c5ad664
begin
	rna004_1000 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_STORM_eventalign_depth_1000_v2.0.0/out_nanocompore_results.tsv",
		"GMM_chi2_qvalue",
		shift=4,
		genomic_collapse=true;
		LOR_threshold = LOR_THRESHOLD)
	binned_rna004_1000 = annotate_binned_results(bin(rna004_1000, BIN_SIZE),
										         binned_glori)
	:ok
end

# ╔═╡ a4751773-afb8-4518-a639-e5c9aa1e1765
common_binned_100 = innerjoin(binned_rna002_100, binned_rna004_100,
							  on = [:chr, :strand, :bin],
							  makeunique = true)

# ╔═╡ 9ea0a845-e839-411c-99ad-051979516cb5
common_binned_500 = innerjoin(binned_rna002_500, binned_rna004_500,
							  on = [:chr, :strand, :bin],
							  makeunique = true)

# ╔═╡ 725002c7-b937-4240-bbb0-22bd28972e90
common_binned_1000 = innerjoin(binned_rna002_1000, binned_rna004_1000,
							   on = [:chr, :strand, :bin],
							   makeunique = true)

# ╔═╡ 35650746-09b3-4529-ab73-28769d178437
# ╠═╡ disabled = true
#=╠═╡
begin
	local v1_results_part1 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1/partial_first_run_results_with_KS_gx.tsv",
		"GMM_pvalue",
		shift = 2,
		genomic_collapse = false;
		LOR_threshold = 0)
	local v1_results_part2 = read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1_part2/outnanocompore_results_gx.tsv",
		"GMM_pvalue",
		shift = 2,
		genomic_collapse = false;
		LOR_threshold = 0)

	local cols = [x for x in intersect(Set(names(v1_results_part1)), Set(names(v1_results_part2)))]
	
	# v1_results = annotate_results(vcat(v1_results_part1, v1_results_part2),
	# 							  ref002)
	v1_results = vcat(v1_results_part1[:, cols], v1_results_part2[:, cols])

	# Correct for multiple testing
	local present = v1_results.GMM_pvalue .!== missing
	local qvals = MultipleTesting.adjust(
		disallowmissing(v1_results[present, :].GMM_pvalue),
		MultipleTesting.BenjaminiHochberg())
	v1_results[!, :GMM_qvalue] = copy(v1_results.GMM_pvalue)
	v1_results[present, :GMM_qvalue] .= qvals

	lor_corrected_significance = ifelse.(abs.(v1_results[!, :GMM_LOR]) .>= LOR_THRESHOLD, v1_results[!, :GMM_qvalue], 1.0)
	v1_results[!, :predicted] = -log10.(lor_corrected_significance)
	
	v1_results
end
  ╠═╡ =#

# ╔═╡ 3ec16ac1-dbec-4ce9-af4d-a913bdc68222
begin
	local v1_results_part1 = DataFrame(CSV.File(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1/partial_first_run_results_with_KS_gx.tsv",
		delim = '\t'))

	local v1_results_part2 = DataFrame(CSV.File(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1_part2/out_nanocompore_results_manual_db_export2.tsv",
		delim = '\t'))

	local cols = [x for x in intersect(Set(names(v1_results_part1)), Set(names(v1_results_part2)))]
	
	v1_results = vcat(v1_results_part1[:, cols], v1_results_part2[:, cols])
	
end

# ╔═╡ cc423cb8-7e32-4b76-b3e6-9fab041f2d1e
# rna002_tx = read_results(
# 		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv",
# 		"GMM_chi2_qvalue",
# 		shift=2,
# 		genomic_collapse=false;
# 	    LOR_threshold = LOR_THRESHOLD)

begin
	rna002_tx = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_eventalign_v2.0.0/out_nanocompore_results.tsv", delim = '\t'))
	# rna002_tx = rna002_tx[rna002_tx.GMM_chi2_pvalue .!== missing, :]
	rna002_tx[:, :GMM_chi2_pvalue] = ifelse.(rna002_tx.GMM_chi2_pvalue .=== missing,
											 1,
											 rna002_tx.GMM_chi2_pvalue)
	rna002_tx[:, :GMM_LOR] = ifelse.(rna002_tx.GMM_LOR .=== missing,
									 0,
									 rna002_tx.GMM_LOR)
	rna002_tx
end

# ╔═╡ a50423b1-d847-4429-9bb9-9d8a7523a459
v1_v2_common = innerjoin(v1_results, rna002_tx,
		  				 on = [:ref_id, :pos],
		  				 renamecols = "_v1" => "_v2")

# ╔═╡ 884239d8-5613-44bc-a14c-348b3e31d3c3
begin
	local common = copy(v1_v2_common[v1_v2_common.GMM_pvalue_v1 .!== missing .&& v1_v2_common.GMM_chi2_pvalue_v2 .!== missing, :])
	local pvals1 = -log10.(clamp.(common.GMM_pvalue_v1, 1e-300, 1))
	local pvals2 = -log10.(clamp.(common.GMM_chi2_pvalue_v2, 1e-300, 1))
	pvals1 = ifelse.(abs.(common.GMM_LOR_v1) .< 0.5, 0, pvals1)
	pvals2 = ifelse.(abs.(common.GMM_LOR_v2) .< 0.5, 0, pvals2)

	println(cor(pvals1, pvals2))
	
	# common = common[.! (pvals1 .== 0 .&& pvals2 .== 0), :]
	# local selected = pvals1 .!= 0 .|| pvals2 .!= 0
	
	scatter(pvals1,
			pvals2,
			markersize = 5,
			color = :black,
		    axis = (;
					aspect = 1,
				    xlabel = "v1: -log₁₀(𝑃-value)",
				    ylabel = "v2: -log₁₀(𝑃-value)"))
	# data(DataFrame(v1 = pvals1, v2 = pvals2)) * mapping(:v1 => "v1: -log₁₀(𝑃-value)", :v2 => "v2: -log₁₀(𝑃-value)") * visual(Scatter, markersize = 5) |> draw
end

# ╔═╡ 866b5d33-c2a7-4400-845e-50575053be71
begin
	local common = copy(v1_v2_common)
	common[:, :bin] = map(pos -> div(pos, BIN_SIZE), common.genomicPos_v2)
	common[:, :predicted_v1] = -log10.(ifelse.(abs.(common.GMM_LOR_v1) .< 0.5,
									   		   1,
									   		   common.GMM_pvalue_v1))
	common[:, :predicted_v2] = -log10.(ifelse.(abs.(common.GMM_LOR_v2) .< 0.5,
									   		   1,
									   		   common.GMM_chi2_pvalue_v2))
	col_selector = (pred, col) -> col[argmax(pred)]
	binned_v1_v2_common = combine(groupby(common, [:chr_v2, :strand_v2, :bin]),
		    [:predicted_v1, :GMM_LOR_v1] => col_selector => :GMM_LOR_v1,
		    [:predicted_v1, :GMM_pvalue_v1] => col_selector => :GMM_pvalue_v1,
		    [:predicted_v2, :GMM_LOR_v2] => col_selector => :GMM_LOR_v2,
		    [:predicted_v2, :GMM_chi2_pvalue_v2] => col_selector => :GMM_pvalue_v2)
	binned_v1_v2_common = rename!(binned_v1_v2_common, :chr_v2 => :chr, :strand_v2 
	=> :strand)

	binned_v1_v2_common[:, :GMM_pvalue_v1] = clamp.(binned_v1_v2_common.GMM_pvalue_v1, 1e-300, 1)
	binned_v1_v2_common[:, :GMM_pvalue_v2] = clamp.(binned_v1_v2_common.GMM_pvalue_v2, 1e-300, 1)

	binned_v1_v2_common = leftjoin(binned_v1_v2_common, binned_glori,
								   on = [:chr, :strand, :bin])
	binned_v1_v2_common[:, :modified] = coalesce.(binned_v1_v2_common.modified, 0)
	binned_v1_v2_common[:, :ratio] = coalesce.(binned_v1_v2_common.ratio, 0)
	binned_v1_v2_common = dropmissing(binned_v1_v2_common, disallowmissing = true)
	binned_v1_v2_common[:, :GMM_qvalue_v1] = MultipleTesting.adjust(
		binned_v1_v2_common.GMM_pvalue_v1,
		MultipleTesting.BenjaminiHochberg())
	binned_v1_v2_common[:, :GMM_qvalue_v2] = MultipleTesting.adjust(
		binned_v1_v2_common.GMM_pvalue_v2,
		MultipleTesting.BenjaminiHochberg())
	binned_v1_v2_common[:, :predicted_v1] = -log10.(binned_v1_v2_common.GMM_qvalue_v1)
	binned_v1_v2_common[:, :predicted_v2] = -log10.(binned_v1_v2_common.GMM_qvalue_v2)
	binned_v1_v2_common
end

# ╔═╡ 9193fd14-d45c-4534-a029-75ab8171ea40
begin
	local rocs1 = roc(binned_v1_v2_common.modified,
			     binned_v1_v2_common.predicted_v1, nrow(binned_v1_v2_common))
	local precisions1 = vcat([0], map(precision, rocs1), [1])
	local recalls1 = vcat([1], map(recall, rocs1), [0])
	local rocs2 = roc(binned_v1_v2_common.modified,
			     binned_v1_v2_common.predicted_v2, nrow(binned_v1_v2_common))
	local precisions2 = vcat([0], map(precision, rocs2), [1])
	local recalls2 = vcat([1], map(recall, rocs2), [0])

	local df1 = DataFrame(precision = precisions1, recall = recalls1)
	df1[:, :version] .= "v1"
	local df2 = DataFrame(precision = precisions2, recall = recalls2)
	df2[:, :version] .= "v2"

	local df = vcat(df1, df2)

	local p001_v1 = precision(roc(binned_v1_v2_common.modified,
			     binned_v1_v2_common.predicted_v1 .>= -log10(0.01)))
	local r001_v1 = recall(roc(binned_v1_v2_common.modified,
	    		  binned_v1_v2_common.predicted_v1 .>= -log10(0.01)))

	local p001_v2 = precision(roc(binned_v1_v2_common.modified,
			     binned_v1_v2_common.predicted_v2 .>= -log10(0.01)))
	local r001_v2 = recall(roc(binned_v1_v2_common.modified,
	    		  binned_v1_v2_common.predicted_v2 .>= -log10(0.01)))

	
	
	(
		data(df) * mapping(:recall => "Recall", :precision => "Precision", color = :version) * visual(Lines)
		+
		data(DataFrame(precision = [p001_v1, p001_v2], recall = [r001_v1, r001_v2], version = ["v1", "v2"])) * mapping(:recall => "Recall", :precision => "Precision", color = :version) * visual(Scatter)
	) |> draw(scales(Color = (; palette = [:black, COL_RNA002])); axis = (; aspect = 1))
end

# ╔═╡ 7a74c048-770c-4782-9c19-ab800f66b7e1
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
			  	    aspect = 1)


	local a, p, r = auprc(binned_v1_v2_common.predicted_v1,
	      		          Int.(binned_v1_v2_common.modified),
	            	      Set([1]))
	lines!(ax, r, p,
	       label = "v1 AUC=$(round(a; digits = 2))",
	       color = :black)
	local p001 = precision(roc(binned_v1_v2_common.modified,
			     binned_v1_v2_common.predicted_v1 .>= -log10(0.01)))
	local r001 = recall(roc(binned_v1_v2_common.modified,
	    		  binned_v1_v2_common.predicted_v1 .>= -log10(0.01)))
	# scatter!(ax, [r001], [p001], color = :black)

	local a, p, r = auprc(binned_v1_v2_common.predicted_v2,
	      		          Int.(binned_v1_v2_common.modified),
	            	      Set([1]))
	lines!(ax, r, p,
	       label = "v2 AUC=$(round(a; digits = 2))",
	       color = COL_RNA002)
	local p001 = precision(roc(binned_v1_v2_common.modified,
			     binned_v1_v2_common.predicted_v2 .>= -log10(0.01)))
	local r001 = recall(roc(binned_v1_v2_common.modified,
	    		  binned_v1_v2_common.predicted_v2 .>= -log10(0.01)))
	# scatter!(ax, [r001], [p001], color = COL_RNA002)
	
	axislegend(ax, "Version", position = :rt, labelsize = 12)
	
	f
end

# ╔═╡ 96a61659-2925-44ce-ae5a-a2e0c6f7b0f4
binned_v1_results = annotate_binned_results(
	bin(innerjoin(v1_results, v1_v2_common[:, [:ref_id, :pos]],
				  on = [:ref_id, :pos]), BIN_SIZE),
	binned_glori)

# ╔═╡ 9054c2f0-7050-4bb9-bd10-da57ec7fa69c


# ╔═╡ 9b8efe40-93b3-4dc9-8a01-ac06eae5024e
v1_refs = vcat(read_results(
					"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1/partial_first_run_results_with_KS_gx.tsv",
					"GMM_pvalue",
					shift = 2,
					genomic_collapse = false;
					LOR_threshold = LOR_THRESHOLD).ref_id,
			   read_results(
					"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1/partial_first_run_results_with_KS_gx.tsv",
					"GMM_pvalue",
					shift = 2,
					genomic_collapse = false;
				    LOR_threshold = LOR_THRESHOLD).ref_id) |> Set

# ╔═╡ 99b9d173-618d-4fd2-9bda-89d6e991711b
v1_positions = vcat(
   map(r -> (r.ref_id, r.pos),
	   eachrow(read_results(
		"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1/partial_first_run_results_with_KS_gx.tsv",
		"GMM_pvalue",
		shift = 2,
		genomic_collapse = false;
		LOR_threshold = LOR_THRESHOLD)[:, [:ref_id, :pos]])),
   map(r -> (r.ref_id, r.pos),
	   eachrow(read_results(
			"/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1/partial_first_run_results_with_KS_gx.tsv",
			"GMM_pvalue",
			shift = 2,
			genomic_collapse = false;
		    LOR_threshold = LOR_THRESHOLD)[:, [:ref_id, :pos]]))) |> Set

# ╔═╡ 095e0c0f-2af2-4dc1-b852-bb554d1d8af1


# ╔═╡ dc4ee9d0-7c6a-49fd-aa99-1fe3b5cb3f87
begin
	v2_matched = filter(r -> (r.ref_id, r.pos) in v1_positions, rna002_tx)
	local col_selector = (predicted, col) -> col[findmax(predicted)[2]]
	v2_matched = combine(groupby(v2_matched, ["chr", "strand", "genomicPos"]),
        	       :predicted => maximum => :predicted,
        	       [:predicted, :GMM_LOR] => col_selector => :GMM_LOR,
        	       :predicted_raw => maximum => :predicted_raw,
        	       [:predicted_raw, :GMM_LOR] => col_selector => :GMM_LOR_raw,
        	       [:predicted_raw, :mean_cov] => col_selector => :mean_cov)
end

# ╔═╡ d891c628-9785-400b-bdc4-bbc9831850f9
binned_v2_matched = annotate_binned_results(
	bin(v2_matched, BIN_SIZE),
	binned_glori)

# ╔═╡ 0e079cb8-3014-4149-ab6d-38f4ed3c256d
v1_v2_matched = innerjoin(binned_v1_results, binned_v2_matched,
						  on = [:chr, :strand, :bin],
						  makeunique = true,
						  renamecols = :_v1 => :_v2)

# ╔═╡ 0697fd84-432e-401b-97c3-dc56b0cccbff
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
			  	    aspect = 1)

	local a, p, r = auprc(v1_v2_matched.predicted_v1,
	      		          Int.(v1_v2_matched.modified_v1),
	            	      Set([1]))
	lines!(ax, r, p,
	       label = "v1 AUC=$(round(a; digits = 2))",
	       color = :black)
	local p001 = precision(roc(v1_v2_matched.modified_v1,
			     v1_v2_matched.predicted_v1 .>= -log10(0.01)))
	local r001 = recall(roc(v1_v2_matched.modified_v1,
	    		  v1_v2_matched.predicted_v1 .>= -log10(0.01)))
	scatter!(ax, [r001], [p001], color = :black)

	local a, p, r = auprc(v1_v2_matched.predicted_v2,
	      		          Int.(v1_v2_matched.modified_v2),
	            	      Set([1]))
	lines!(ax, r, p,
	       label = "v2 AUC=$(round(a; digits = 2))",
	       color = COL_RNA002)
	local p001 = precision(roc(v1_v2_matched.modified_v2,
			     v1_v2_matched.predicted_v2 .>= -log10(0.01)))
	local r001 = recall(roc(v1_v2_matched.modified_v2,
	    		  v1_v2_matched.predicted_v2 .>= -log10(0.01)))
	scatter!(ax, [r001], [p001], color = COL_RNA002)
	
	axislegend(ax, "Version", position = :rt, labelsize = 12)
	
	f
end

# ╔═╡ e48d049b-4297-4c1e-9162-85f947d1dc20
begin
	local f = Figure()
	local ax = Axis(f[1, 1],
			  aspect = 1)

	local df = leftjoin(v1_v2_common, glori,
						on = [:chr_v2 => :chrom,
							  :strand_v2 => :strand,
							  :genomicPos_v2 => :pos])
	df[:, :modified] = ifelse.(df.modified .=== missing, 0, df.modified)
	df = df[df.GMM_pvalue_v1 .!== missing .&& df.GMM_chi2_pvalue_v2 .!== missing, :]

	local a, p, r = auprc(-log10.(ifelse.(abs.(df.GMM_LOR_v1) .< 0.5,
								  		  1,
								  		  df.GMM_pvalue_v1)),
	      		          Int.(df.modified),
	            	      Set([1]))
	lines!(ax, r, p,
	       label = "v1 AUC=$(round(a; digits = 2))",
	       color = :black)

	local a, p, r = auprc(-log10.(ifelse.(abs.(df.GMM_LOR_v2) .< 0.5,
										  1,
										  df.GMM_chi2_pvalue_v2)),
	      		          Int.(df.modified),
	            	      Set([1]))
	lines!(ax, r, p,
	       label = "v2 AUC=$(round(a; digits = 2))",
	       color = COL_RNA002)
	f
	# df
end

# ╔═╡ 91c41c2e-5ae1-4cb2-86b9-d4fb953570d1
leftjoin(v1_v2_common, glori,
		 on = [:chr_v1 => :chrom, :strand_v1 => :strand, :genomicPos_v1 => :pos])

# ╔═╡ 2f1d642b-d646-416d-9550-bbe30f2044a1
scatter(v1_v2_matched.predicted_v1, v1_v2_matched.predicted_v2, markersize = 5)

# ╔═╡ 53df2546-878c-4027-bb5b-6f8a4af5bb5f
cor(v1_v2_matched.predicted_v1, v1_v2_matched.predicted_v2)

# ╔═╡ b94447fc-762a-43b8-9e0c-987064a62618
begin
	local rna002_fps = binned_rna002[binned_rna002.modified .== 0 .&& binned_rna002.predicted .>= -log10(0.01), :]
	local rna004_fps = binned_rna004[binned_rna004.modified .== 0 .&& binned_rna004.predicted .>= -log10(0.01), :]
	
	local common = size(innerjoin(rna002_fps, rna004_fps,
	                              on = [:chr, :strand, :bin],
	                              makeunique = true), 1)
	local rna002_only = size(antijoin(rna002_fps, rna004_fps,
	                                  on = [:chr, :strand, :bin],
	                                  makeunique = true), 1)
	local rna004_only = size(antijoin(rna004_fps, rna002_fps,
	                                  on = [:chr, :strand, :bin],
	                                  makeunique = true), 1)

	(rna002_only, common, rna004_only)
end

# ╔═╡ 31b0f531-c921-46b2-aada-87098aa16d97
begin
	f = Figure(size = (1400, 800))

	ga = f[1, 1:2] = GridLayout()
	gb = f[1, 3:8] = GridLayout()
	# gx = f[2, 1] = GridLayout()
	gc = f[2, 1:2] = GridLayout()
	gd = f[2, 3:6] = GridLayout()
	ge = f[2, 7:8] = GridLayout()

	
	# ax_a1 = Axis(ga[1, 1],
	# 	       aspect = 1,
	# 	       title = "Precision/recall for v1 and v2 (RNA002)",
	# 	       xlabel = "Recall",
	# 	       ylabel = "Precision")
	# xlims!(ax_a1, [0, 1])
	# ylims!(ax_a1, [0, 1])

	# local a, p, r = auprc(v1_v2_matched.predicted_v1,
	#       		          Int.(v1_v2_matched.modified_v1),
	#             	      Set([1]))
	# lines!(ax_a1, r, p,
	#        label = "v1 AUC=$(round(a; digits = 2))",
	#        color = :black)
	# local p001 = precision(roc(v1_v2_matched.modified_v1,
	# 		     v1_v2_matched.predicted_v1 .>= -log10(0.01)))
	# local r001 = recall(roc(v1_v2_matched.modified_v1,
	#     		  v1_v2_matched.predicted_v1 .>= -log10(0.01)))
	# scatter!(ax_a1, [r001], [p001], color = :black)

	# local a, p, r = auprc(v1_v2_matched.predicted_v2,
	#       		          Int.(v1_v2_matched.modified_v2),
	#             	      Set([1]))
	# lines!(ax_a1, r, p,
	#        label = "v2 AUC=$(round(a; digits = 2))",
	#        color = :red)
	# local p001 = precision(roc(v1_v2_matched.modified_v2,
	# 		     v1_v2_matched.predicted_v2 .>= -log10(0.01)))
	# local r001 = recall(roc(v1_v2_matched.modified_v2,
	#     		  v1_v2_matched.predicted_v2 .>= -log10(0.01)))
	# scatter!(ax_a1, [r001], [p001], color = :red)

	# axislegend(ax_a1, "Version", position = :lb, labelsize = 12)

	# ==== TOP ROW =====

	# Top left (correlation between v1 and v2)

	ax_a = Axis(ga[1, 1:2],
		       aspect = 1,
		       title = "Correlation between v1 and v2",
		       xlabel = "v1: -log₁₀(𝑃-value)",
		       ylabel = "v2: -log₁₀(𝑃-value)")

	scatter!(ax_a,
			 v1_v2_matched.predicted_v1,
			 v1_v2_matched.predicted_v2,
			 markersize = 5,
			 color = :black)

	# Top right (performance evaluation barplots)

	local ax_collapse = Axis(gb[1, 1],
							 aspect = 1/2.6,
							 title = "Preprocessing\ntime",
							 #titlesize = 13,
							 ylabel = "Execution time (s) / transcript",
							 xticks = ([1, 2], ["v1", "v2"]),
							  	 yticks = 0:0.5:2.5)
	ylims!(ax_collapse, (-0.1, 2.7))
	local ax_comp = Axis(gb[1, 2],
						 aspect = 1/2.6,
						 title = "Analysis\ntime",
						 ylabel = "Execution time (s) / transcript",
						 xticks = ([1, 2], ["v1", "v2"]),
						 yticks = 0:5:30)
	ylims!(ax_comp, (-1, 31))
	local ax_space = Axis(gb[1, 3],
						  aspect = 1/2.6,
						  title = "Preprocessing\ndisk space",
						  ylabel = "Disk space (GB) / 1 million reads",
						  xticks = ([1, 2], ["v1", "v2"]),
						  yticks = 10:10:60)
	ylims!(ax_space, (-1, 67))
	local ax_space2 = Axis(gb[1, 4],
					 	   aspect = 1/2.6,
						   title = "Analysis\ndisk space",
						   ylabel = "Disk space (MB) / 1,000 tested sites",
						   xticks = ([1, 2], ["v1", "v2"]),
						   yticks = 0:5:20)
	ylims!(ax_space2, (-0.05, 22.2))

	local t1 = (35*60+3)/1000
	local t2 = (7*60)/1000
	barplot!(ax_collapse,
			 #[(35*60+3)/60, (7*60)/60],
			 [t1, t2],
			 color = :black,
			 bar_labels = ["$(@sprintf("%.2f", t1))s",
						   "$(@sprintf("%.2f", t2))s"])
	
	local t1 = (59*60+57)/136
	local t2 = (3*60+21)/136
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
	

	# # == Total counts, second row left
	# local ax_totals = Axis(gx[1, 1],
	# 					 aspect = 1/2.6,
	# 					 title = "All tested positions",
	# 					 ylabel = "Count",
	# 					 xticks = ([1, 2], ["RNA002", "RNA004"]))
	# 					 # yticks = 0:5:30)

	local roc002 = roc(binned_rna002.modified .== 1,
			     	   binned_rna002.predicted .>= -log10(0.01))
	local roc004 = roc(binned_rna004.modified .== 1,
			     	   binned_rna004.predicted .>= -log10(0.01))
	
	# # ylims!(ax_totals, (-1, 31))
	# local df = DataFrame(
	#      kit = [1, 1, 1, 1, 2, 2, 2, 2],
	#      class = [0, 1, 2, 3, 0, 1, 2, 3],
	#      number = [true_negative(roc002),
	# 		       false_negative(roc002),
	# 		       true_positive(roc002),
	# 		       false_positive(roc002),
	# 		       true_negative(roc004),
	# 			   false_negative(roc004),
	# 		       true_positive(roc004),
	# 			   false_positive(roc004)])
	
	# local colormap = Dict(0 => :black,
	# 					  1 => :black,
	# 					  2 => :black,
	# 					  3 => :black)
	
	# barplot!(ax_totals, df.kit, df.number,
	#          stack = df.class,
	#          color = map(c -> colormap[c], df.class),
	#          direction = :y,
	#          bar_labels = [(@sprintf "%d" true_negative(roc002)),
	#                        (@sprintf "%d" false_negative(roc002)),
	#                        (@sprintf "%d" true_positive(roc002)),
	# 					   (@sprintf "%d" false_positive(roc002)),
	# 					   (@sprintf "%d" true_negative(roc004)),
	#                        (@sprintf "%d" false_negative(roc004)),
	#                        (@sprintf "%d" true_positive(roc004)),
	# 					   (@sprintf "%d" false_positive(roc004))],
	#          #label_position = :end, #[:start, :start, :end, :start, :center, :end],
	# 		 label_position = :center,
	# 		 label_size = 0,
	# 		 gap = 0.4,
	#          # label_align = :lb,
	#          #label_offset = [-160, -70, 5, -135, -80, 5],
	#          label_color = :black)
	

	# == Subfig C, second row left
	
	# a, p, r = auprc(binned_v1_results.predicted,
	#                 Int.(binned_v1_results.modified),
	#                 Set([1]))
	# lines!(ax_c, r, p,
	#        label = "RNA002 v1\nAUC=$(round(a; digits = 2))",
	#        color = COL_RNA002, linestyle = :dash)
	# p001 = precision(roc(binned_v1_results.modified,
	# 		     		 binned_v1_results.predicted .>= -log10(0.01)))
	# r001 = recall(roc(binned_v1_results.modified,
	#     		  	  binned_v1_results.predicted .>= -log10(0.01)))
	# scatter!(ax_c, [r001], [p001], color = COL_RNA002)

	ax_c = Axis(gc[1, 1],
			    title = "Precision/recall for v2 by chemistry",
			    aspect = 1,
		        xlabel = "Recall",
		        ylabel = "Precision")
	
	a, p, r = auprc(binned_rna002.predicted,
	                Int.(binned_rna002.modified),
	                Set([1]))
	lines!(ax_c, r[5:end], p[5:end],
	       label = "RNA002\nAUC=$(round(a; digits = 2))",
	       color = COL_RNA002)
	p001 = precision(roc(binned_rna002.modified,
			    		 binned_rna002.predicted .>= -log10(0.01)))
	r001 = recall(roc(binned_rna002.modified,
	    			  binned_rna002.predicted .>= -log10(0.01)))
	scatter!(ax_c, [r001], [p001], color = COL_RNA002)
	println("RNA002 precision = $p001, recall = $r001")
	
	a, p, r = auprc(binned_rna004.predicted,
	                Int.(binned_rna004.modified),
	                Set([1]))
	lines!(ax_c, r[5:end], p[5:end],
	       label = "RNA004\nAUC=$(round(a; digits = 2))",
	       color = COL_RNA004)
	p001 = precision(roc(binned_rna004.modified,
			    		 binned_rna004.predicted .>= -log10(0.01)))
	r001 = recall(roc(binned_rna004.modified,
	    		  binned_rna004.predicted .>= -log10(0.01)))
	scatter!(ax_c, [r001], [p001], color = COL_RNA004)
	println("RNA004 precision = $p001, recall = $r001")
	
	axislegend(ax_c, "Chemistry", position = :rt, labelsize = 12)
	
	# ==== Same coverage (c) =====
	ax_d1 = Axis(gd[1, 1:2],
	             title = "Precision at equal coverage",
	             xlabel = "Depth",
	             xticks = (1:3, ["100", "500", "1000"]),
	             ylabel = "Precision",
	             yticks = [0, 0.25, 0.5, 0.75, 1.0])
	ylims!(ax_d1, [0, 1.15])
	ax_d2 = Axis(gd[1, 3:4],
	             title = "Recall at equal coverage",
	             xlabel = "Depth",
	             xticks = (1:3, ["100", "500", "1000"]),
	             ylabel = "Recall",
	             yticks = [0, 0.25, 0.5, 0.75, 1.0])
	ylims!(ax_d2, [0, 1.15])
	
	# roc_rna002_100 = roc(binned_rna002_100.modified,
	#                      binned_rna002_100.predicted .>= -log10(0.01))
	# roc_rna002_500 = roc(binned_rna002_500.modified,
	#                      binned_rna002_500.predicted .>= -log10(0.01))
	# roc_rna002_1000 = roc(binned_rna002_1000.modified,
	#                       binned_rna002_1000.predicted .>= -log10(0.01))
	# roc_rna004_100 = roc(binned_rna004_100.modified,
	#                      binned_rna004_100.predicted .>= -log10(0.01))
	# roc_rna004_500 = roc(binned_rna004_500.modified,
	#                      binned_rna004_500.predicted .>= -log10(0.01))
	# roc_rna004_1000 = roc(binned_rna004_1000.modified,
	#                       binned_rna004_1000.predicted .>= -log10(0.01))
	roc_rna002_100 = roc(common_binned_100.modified,
	                     common_binned_100.predicted .>= -log10(0.01))
	roc_rna002_500 = roc(common_binned_500.modified,
	                     common_binned_500.predicted .>= -log10(0.01))
	roc_rna002_1000 = roc(common_binned_1000.modified,
	                      common_binned_1000.predicted .>= -log10(0.01))
	roc_rna004_100 = roc(common_binned_100.modified_1,
	                     common_binned_100.predicted_1 .>= -log10(0.01))
	roc_rna004_500 = roc(common_binned_500.modified_1,
	                     common_binned_500.predicted_1 .>= -log10(0.01))
	roc_rna004_1000 = roc(common_binned_1000.modified_1,
	                      common_binned_1000.predicted_1 .>= -log10(0.01))
	
	df = DataFrame(
	     # kit = ["RNA002", "RNA002", "RNA002", "RNA004", "RNA004", "RNA004"],
	     kit = [1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2],
	     # depth = [100, 500, 1000, 100, 500, 1000],
	     depth = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
	     metric = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
	     values = [precision(roc_rna002_100),
		       precision(roc_rna002_500),
		       precision(roc_rna002_1000),
		       precision(roc_rna004_100),
		       precision(roc_rna004_500),
		       precision(roc_rna004_1000),
		       recall(roc_rna002_100),
		       recall(roc_rna002_500),
		       recall(roc_rna002_1000),
		       recall(roc_rna004_100),
		       recall(roc_rna004_500),
		       recall(roc_rna004_1000)])
	
	colormap = Dict(1 => "#56B4E9",
		        	2 => "#009E73")
	
	barplot!(ax_d1,
	         df.depth[df.metric .== 1],
	         df.values[df.metric .== 1],
	         dodge = df.kit[df.metric .== 1],
	         color = map(k -> colormap[k], df.kit[df.metric .== 1]),
	         bar_labels = map(v -> (@sprintf("%.2f", v)), df.values[df.metric .== 1]),
			 label_size = 11)
	
	barplot!(ax_d2,
	         df.depth[df.metric .== 2],
	         df.values[df.metric .== 2],
	         dodge = df.kit[df.metric .== 2],
	         color = map(k -> colormap[k], df.kit[df.metric .== 2]),
	         bar_labels = map(v -> (@sprintf("%.2f", v)), df.values[df.metric .== 2]),
			 label_size = 11)

	
	Legend(ge[2, 1],
		   [LineElement(color = :black, linestyle = :solid),
		    LineElement(color = :black, linestyle = :dash)],
		   ["GLORI", "Nanocompore"],
		   valign = :top,
		   halign = :center,
		   orientation = :horizontal)
	
	Legend(gd[2, 4],
	       [PolyElement(color = COL_RNA002),
			PolyElement(color = COL_RNA004)],
	       ["RNA002", "RNA004"],
		   valign = :center,
		   halign = :right,
		   orientation = :horizontal)
	

	
	# ==== D2: RNA002/RNA002 barplot & venn ====
	ax_glori = Axis(gd[3, 1:4][2, 1],
		    yticks = (1:2, ["RNA002", "RNA004"]),
		    xlabel = "m⁶A sites in GLORI",
		    xticks = [0, 50000, 100000],
		    xtickformat = "{:n}",
		    height = 75)
	# xlims!(ax_glori, (0, 110000))
	xlims!(ax_glori, (0, nrow(glori)))
	
	ax_venn = Axis(ge[3, 1],
	               title = "Overlap of true positives",
   	               height = 76)
	
	# ax22 = Axis(f[1, 3:6][2, 2],
	#             xlabel = "GLORI modification stoichiometry",
	#             ylabel = "Density",
	#             height = 100)
	
	glori_not_covered_rna002 = size(
	      antijoin(glori, binned_rna002,
		       on = [:chrom => :chr,
			         :strand => :strand,
			         :bin => :bin]), 1)
	glori_not_covered_rna004 = size(
	      antijoin(glori, binned_rna004,
		       on = [:chrom => :chr,
			         :strand => :strand,
			         :bin => :bin]), 1)
	# local roc002 = roc(binned_rna002.modified .== 1, binned_rna002.predicted .>= -log10(0.01))
	# local roc004 = roc(binned_rna004.modified .== 1, binned_rna004.predicted .>= -log10(0.01))
	
	df = DataFrame(
	     #kit = ["RNA002", "RNA002", "RNA002", "RNA004", "RNA004", "RNA004"],
	     #class = ["Unmapped", "Missed", "Detected", "Unmapped", "Missed", "Detected"],
	     kit = [1, 1, 1, 2, 2, 2],
	     class = [0, 1, 2, 0, 1, 2],
	     number = [glori_not_covered_rna002,
		       false_negative(roc002),
		       true_positive(roc002),
		       glori_not_covered_rna004,
		       false_negative(roc004),
		       true_positive(roc004)])
	
	colormap = Dict(0 => :black,
	                1 => :grey,
	                2 => :lightgreen)
	
	barplot!(ax_glori, df.kit, df.number,
	         stack = df.class,
	         color = map(c -> colormap[c], df.class),
	         direction = :x,
	         bar_labels = [(@sprintf "%d" glori_not_covered_rna002),
	                       (@sprintf "%d" false_negative(roc002)),
	                       (@sprintf "%d" true_positive(roc002)),
	                       (@sprintf "%d" glori_not_covered_rna004),
	                       (@sprintf "%d" false_negative(roc004)),
	                       (@sprintf "%d" true_positive(roc004))],
	         #label_position = :end, #[:start, :start, :end, :start, :center, :end],
			 label_position = :center,
			 label_size = [12, 12, 12, 12, 12, 12],
	         # label_align = :lb,
	         #label_offset = [-160, -70, 5, -135, -80, 5],
	         label_color = [:white, :white, :black, :white, :white, :black])
	# bar_labels = map(n -> (@sprintf "%d" n), df.number),
	# label_position = :center,
	# label_color = :white)
	
	local rna002_tps = binned_rna002[binned_rna002.modified .== 1 .&& binned_rna002.predicted .>= -log10(0.01), :]
	local rna004_tps = binned_rna004[binned_rna004.modified .== 1 .&& binned_rna004.predicted .>= -log10(0.01), :]
	
	local common = size(innerjoin(rna002_tps, rna004_tps,
	                              on = [:chr, :strand, :bin],
	                              makeunique = true), 1)
	local rna002_only = size(antijoin(rna002_tps, rna004_tps,
	                                  on = [:chr, :strand, :bin],
	                                  makeunique = true), 1)
	local rna004_only = size(antijoin(rna004_tps, rna002_tps,
	                                  on = [:chr, :strand, :bin],
	                                  makeunique = true), 1)
	
	# axislegend(ax_d_bl, [PolyElement(color = :lightgreen)], ["TP"]; position=:rt)
	# ax_d_bl = Axis(gd[2, 1][2, 1],
	Legend(gd[2, 1:3][1, 1],
		   [PolyElement(color = :black),
			PolyElement(color = :grey),
			PolyElement(color = :lightgreen)],
	       ["no coverage", "false negatives", "true positives"],
		   orientation = :horizontal,
		   framevisible = false,
		   valign = :center,
		   halign = :left)
	
	# rowsize!(f.layout, 1, Relative(0.75))
	
	poly!(ax_venn, Circle(Point2f(-60, 0), sqrt(9851/pi)), color = (:white, 0), linestyle = :solid, strokewidth = 2, strokecolor = COL_RNA002)
	poly!(ax_venn, Circle(Point2f(0, 0), sqrt(21627/pi)), color = (:white, 0), linestyle = :solid, strokewidth = 2, strokecolor = COL_RNA004)
	
	# xlims!(ax21, (-1.6, 2))
	# ylims!(ax21, (-1.6, 2))
	
	text!(ax_venn, -95, 0, text = "$rna002_only", align = (:center, :center), fontsize = 12)
	text!(ax_venn, -50, 0, text = "$common", align = (:center, :center), fontsize = 12)
	text!(ax_venn, 35, 0, text = "$rna004_only", align = (:center, :center), fontsize = 12)
	
	# poly!(ax21, Circle(Point2f(-0.9, 0), 1.3*2.777/(2*pi)),
	#       color = (:white, 0),
	#       linestyle = :solid,
	#       strokewidth = 2,
	#       strokecolor = COL_RNA002)
	# poly!(ax21, Circle(Point2f(0.4, 0), 9.043/(2*pi)),
	#       color = (:white, 0),
	#       linestyle = :solid,
	#       strokewidth = 2,
	#       strokecolor = COL_RNA004)
	# 
	# xlims!(ax21, (-1.6, 2))
	# ylims!(ax21, (-1.6, 2))
	#
	# text!(ax21, -1.25, 0, text = "$rna002_only", align = (:center, :center), fontsize = 10)
	# text!(ax21, -0.74, 0, text = "$common", align = (:center, :center), fontsize = 10)
	# text!(ax21, 0.3, 0, text = "$rna004_only", align = (:center, :center), fontsize = 10)

	ax_stoich = Axis(ge[1, 1],
				title = "Modification stoichiometry density",
			   	xlabel = "Modification stoichiometry",
				ylabel = "Density",)
	
	density!(ax_stoich,
		     glori.ratio,
	         color = (COL_RNA004, 0),
	         strokecolor = :black,
	         strokewidth = 1,
	         # strokearound = true,
		     linestyle = :solid,
	         label = "Overall")
	density!(ax_stoich,
	         # leftjoin(rna002[rna002.modified .== 1, :],
	         #          binned_glori,
	         #          on = [:chr, :strand, :bin],
	         #          makeunique = true).ratio,
			 innerjoin(glori, ref002[:, [:chr, :strand, :genomicPos]],
					   on = [:chrom => :chr, :strand, :start => :genomicPos]).ratio,
	         color = (COL_RNA002, 0),
	         strokecolor = COL_RNA002,
	         strokewidth = 1,
	         # strokearound = true,
		     linestyle = :solid,
	         # boundary = (0, 1),
	         label = "RNA002")
	density!(ax_stoich,
	         innerjoin(rna002[rna002.modified .== 1 .&& rna002.predicted .>= -log10(0.01), :],
	                  binned_glori,
	                  on = [:chr, :strand, :bin],
	                  makeunique = true).ratio,
	         color = (COL_RNA002, 0),
	         strokecolor = COL_RNA002,
	         strokewidth = 2,
	         # strokearound = true,
		     linestyle = :dash,
	         # boundary = (0, 1),
	         label = "RNA002")
	density!(ax_stoich,
	         # leftjoin(rna004[rna004.modified .== 1, :],
	         #          glori,
	         #          on = [:chr, :strand, :bin],
	         #          makeunique = true).ratio,
			 innerjoin(glori, ref004[:, [:chr, :strand, :genomicPos]],
					   on = [:chrom => :chr, :strand, :start => :genomicPos]).ratio,
	         color = (COL_RNA004, 0),
	         strokecolor = COL_RNA004,
	         strokewidth = 1,
	         # strokearound = true,
		     linestyle = :solid,
	         label = "RNA004")
	density!(ax_stoich,
	         innerjoin(rna004[rna004.modified .== 1 .&& rna004.predicted .>= -log10(0.01), :],
	                  binned_glori,
	                  on = [:chr, :strand, :bin],
	                  makeunique = true).ratio,
	         color = (COL_RNA004, 0),
	         strokecolor = COL_RNA004,
	         strokewidth = 2,
	         # strokearound = true,
		     linestyle = :dash,
	         label = "RNA004")
	# density!(ax22,
	# glori.mean_m6A_level,
	# color = (:grey, 0),
	# linestyle = :dash,
	# strokecolor = :grey,
	# strokewidth = 2,
	# strokearound = true)
	
	hidedecorations!(ax_venn)
	hidespines!(ax_venn)
	f
end

# ╔═╡ af54cf26-5304-4250-a21c-61b08871defa
begin
	local f = Figure(size = (640, 340))
	# local ax_collapse = Axis(f[1, 1],
	# 						 aspect = 1/2.6,
	# 						 title = "Preprocessing\ntime",
	# 						 #titlesize = 13,
	# 						 ylabel = "Execution time (m) / 1,000 transcripts",
	# 						 xticks = ([1, 2], ["v1", "v2"]),
	# 					  	 yticks = 0:100:450)
	# ylims!(ax_collapse, (-1, 500))
	# local ax_comp = Axis(f[1, 2],
	# 					 aspect = 1/2.6,
	# 					 title = "Analysis\ntime",
	# 					 ylabel = "Execution time (m) / 1,000 transcripts",
	# 					 xticks = ([1, 2], ["v1", "v2"]),
	# 					 yticks = 0:100:450)
	# ylims!(ax_comp, (-1, 500))
	# local ax_space = Axis(f[1, 3],
	# 					  aspect = 1/2.6,
	# 					  title = "Preprocessing\ndisk space",
	# 					  ylabel = "Disk space (GB) / 1 million reads",
	# 					  xticks = ([1, 2], ["v1", "v2"]),
	# 					  yticks = 10:10:60)
	# ylims!(ax_space, (-1, 67))
	# local ax_space2 = Axis(f[1, 4],
	# 				 	   aspect = 1/2.6,
	# 					   title = "Analysis\ndisk space",
	# 					   ylabel = "Disk space (MB) / 1,000 tested sites",
	# 					   xticks = ([1, 2], ["v1", "v2"]),
	# 					   yticks = 0:5:20)
	# ylims!(ax_space2, (-0.05, 22.2))

	# barplot!(ax_collapse,
	# 		 [(35*60+3)/60, (7*60)/60],
	# 		 color = :black,
	# 		 bar_labels = ["35m", "7m"])
	
	# barplot!(ax_comp,
	# 		 # [59u"minute" + 57u"s", 3u"minute" + 21u"s"],
	# 		 # [(59*60+57)/60, (3*60+21)/60],
	# 		 [441, 22],
	# 		 color = :black,
	# 		 bar_labels = ["441m", "22m"])

	# barplot!(ax_space,
	# 		 [55.52, 19.43],
	# 		 color = :black,
	# 		 bar_labels = ["56GB", "19GB"])

	# barplot!(ax_space2,
	# 		 [17, 1],
	# 		 color = :black,
	# 		 bar_labels = ["17MB", "1MB"])

	local ax_collapse = Axis(f[1, 1],
							 aspect = 1/2.6,
							 title = "Preprocessing\ntime",
							 #titlesize = 13,
							 ylabel = "Execution time (s) / transcript",
							 xticks = ([1, 2], ["v1", "v2"]),
							  	 yticks = 0:0.5:2.5)
	ylims!(ax_collapse, (-0.1, 2.7))
	local ax_comp = Axis(f[1, 2],
						 aspect = 1/2.6,
						 title = "Analysis\ntime",
						 ylabel = "Execution time (s) / transcript",
						 xticks = ([1, 2], ["v1", "v2"]),
						 yticks = 0:5:30)
	ylims!(ax_comp, (-1, 31))
	local ax_space = Axis(f[1, 3],
						  aspect = 1/2.6,
						  title = "Preprocessing\ndisk space",
						  ylabel = "Disk space (GB) / 1 million reads",
						  xticks = ([1, 2], ["v1", "v2"]),
						  yticks = 10:10:60)
	ylims!(ax_space, (-1, 67))
	local ax_space2 = Axis(f[1, 4],
					 	   aspect = 1/2.6,
						   title = "Analysis\ndisk space",
						   ylabel = "Disk space (MB) / 1,000 tested sites",
						   xticks = ([1, 2], ["v1", "v2"]),
						   yticks = 0:5:20)
	ylims!(ax_space2, (-0.05, 22.2))

	local t1 = (35*60+3)/1000
	local t2 = (7*60)/1000
	barplot!(ax_collapse,
			 #[(35*60+3)/60, (7*60)/60],
			 [t1, t2],
			 color = :black,
			 bar_labels = ["$(@sprintf("%.2f", t1))s",
						   "$(@sprintf("%.2f", t2))s"])
	
	local t1 = (59*60+57)/136
	local t2 = (3*60+21)/136
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

	f	
end

# ╔═╡ eee69a22-1144-4215-ac6b-ce31937a60bd
begin
	local df002 = copy(binned_rna002[:, [:predicted, :modified, :ratio]])
	local df004 = copy(binned_rna004[:, [:predicted, :modified, :ratio]])
	df002[:, :ratio] = ifelse.(df002.ratio .=== missing, 0, df002.ratio)
	df002[:, :chemistry] .= "rna002"
	df004[:, :ratio] = ifelse.(df004.ratio .=== missing, 0, df004.ratio)
	df004[:, :chemistry] .= "rna004"

	local df = vcat(df002, df004)
	df[:, :ratio_bin] = map(r -> div(100*r, 10)*10, df.ratio)
	df = combine(groupby(df, [:chemistry, :ratio_bin]),
				 nrow => :count,
		   		 [:predicted, :modified] => ((p, m) -> recall(roc(m, p .>= 2))) => :recall)

	sort!(df, [:chemistry, :ratio_bin])

	# local f = Figure()
	# local ax = Axis(f[1, 1])

	# barplot!(ax,
	# 		 df.ratio_bin,
	# 		 df.recall,
	# 		 color = ifelse.(df.chemistry .== "rna002", 0, 1),
	# 	     dodge = ifelse.(df.chemistry .== "rna002", 0, 1))

	# f

	# data(df) * mapping(:recall, color = :chemistry, dodge = :chemistry) * frequency() |> draw
	data(df) * mapping(:ratio_bin => "Stoichiometry", :recall => "Recall", color = :chemistry) * (visual(Lines) + visual(Scatter)) |> draw(scales(Color = (; palette = [COL_RNA002, COL_RNA004])))
end

# ╔═╡ 7c45dd08-df7b-44aa-a8c4-4895d73e4bd9
binned_rna004[binned_rna004.modified .=== 1 .&& binned_rna004.ratio .> 0.9 .&& binned_rna004.predicted .< 2, :]

# ╔═╡ ec46df8d-4a25-4fe6-b82b-9292c8d345ef
minimum(skipmissing([missing]), init = 1)

# ╔═╡ 3c8398c3-9186-4cb3-89fd-95edd5f94d94
begin
	# local df = binned_rna004[binned_rna004.modified .=== 1 .&& binned_rna004.ratio .> 0.9, :]
	local df = binned_rna004[binned_rna004.modified .=== 1 .&& binned_rna004.ratio .> 0.9 .&& binned_rna004.ratio .< 1, :]
	sum(df.GMM_chi2_pvalue .=== missing) / nrow(df)
end

# ╔═╡ e25b7e46-9ea3-4210-a661-b55e69b8fc2c
# rna004[rna004.bin .== 102743, :]
rna004[rna004.bin .== 137185, :]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
MLBase = "f0e99cf1-93fa-52ec-9ecc-5026115318e0"
MultipleTesting = "f8716d33-7c4a-5097-896f-ce0ecbd3ef6b"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
AlgebraOfGraphics = "~0.10.2"
CSV = "~0.10.15"
CairoMakie = "~0.13.6"
DataFrames = "~1.7.0"
MLBase = "~0.9.2"
MultipleTesting = "~0.6.0"
StatsBase = "~0.34.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "2d72d72b0cd446ffe9f1ea0d8f6693826a99594b"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AlgebraOfGraphics]]
deps = ["Accessors", "Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "Isoband", "KernelDensity", "Loess", "Makie", "NaturalSort", "PlotUtils", "PolygonOps", "PooledArrays", "PrecompileTools", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "90983916b063f36f585f9f54b06e7fc50dfe49e8"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.10.2"

    [deps.AlgebraOfGraphics.extensions]
    AlgebraOfGraphicsDynamicQuantitiesExt = "DynamicQuantities"
    AlgebraOfGraphicsUnitfulExt = "Unitful"

    [deps.AlgebraOfGraphics.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "d116d9b54ff8a3b84b2ed0be32dff6304e9c7798"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.13.6"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "e771a63cc8b539eca78c85b0cabd9233d6c8f06f"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "5620ff4ee0084a6ab7097a27ba0c19290200b037"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.4"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "a86af9c4c4f33e16a2b2ff43c2113b2f390081fa"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.4.5"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3e6d038b77f22791b8e3472b7c633acea1ecac06"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.120"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "8cc47f299902e13f90405ddb5bf87e5d474c0d38"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "6.1.2+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "d52e255138ac21be31fa633200b65e4e71d26802"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.6"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "273bd1cd30768a2fddfa3fd63bbc746ed7249e5f"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.9.0"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "8e233d5167e63d708d41f87597433f59a0f213fe"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.4"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "294e99f19869d0b0cb71aef92f19d03649d028d5"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.4.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "2670cf32dcf0229c9893b895a9afe725edb23545"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.9"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "dc6bed05c15523624909b3953686c5f5ffa10adc"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "6a9fde685a7ac1eb3495f8e812c5a7c3711c2d5e"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.3"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Random", "RoundingEmulator"]
git-tree-sha1 = "694c52705f8b23dc5b39eeac629dc3059a168a40"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.35"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "f749e7351f120b3566e5923fefdf8e52ba5ec7f9"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.4"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MLBase]]
deps = ["IterTools", "Random", "Reexport", "StatsBase"]
git-tree-sha1 = "ac79beff4257e6e80004d5aee25ffeee79d91263"
uuid = "f0e99cf1-93fa-52ec-9ecc-5026115318e0"
version = "0.9.2"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "f79e47c8341376c283d3ff3b6eaeee2f73ce69a1"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.22.6"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "733d910c70805e7114c82508bae99c6cdf004466"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.9.3"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "f5a6805fb46c0285991009b526ec6fae43c6dec2"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultipleTesting]]
deps = ["Distributions", "SpecialFunctions", "StatsBase"]
git-tree-sha1 = "1e98f8f732e7035c4333135b75605b74f3462b9b"
uuid = "f8716d33-7c4a-5097-896f-ce0ecbd3ef6b"
version = "0.6.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+4"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9216a80ff3682833ac4b733caa8c00390620ba5d"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.0+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3b31172c032a1def20c98dae3f2cdc9d10e3b561"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "b81c5035922cc89c2d9523afc6c54be512411466"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.5"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "8e45cecc66f3b42633b8ce14d431e8e57a3e242e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "9022bcaa2fc1d484f1326eaa4db8db543ca8c66d"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.4"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "8ad2e38cbb812e29348719cc63580ec1dfeb9de4"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.1"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "f21231b166166bebc73b99cea236071eb047525b"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.3"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d2282232f8a4d71f79e85dc4dd45e5b12a6297fb"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.23.1"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "002748401f7b520273e2b506f61cab95d4701ccf"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.48+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "d2408cac540942921e7bd77272c32e58c33d8a77"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.5.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc541bb19ed5b0ede95581fb2e41ecf179527d2"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.6.0+0"
"""

# ╔═╡ Cell order:
# ╠═11b96b0e-33e1-11f0-23a0-f998c7e38e3f
# ╠═78c86cd4-1dce-4afb-bd82-633875018901
# ╠═b0e757dd-c8ce-45f6-b448-e8543168c7bf
# ╠═b7f89781-6332-4942-8d5b-7d95e0550b63
# ╠═6608e7a5-fde0-4f0c-9a33-a89c0fe0b858
# ╠═05ce3923-41da-4ab4-9ae5-860e0531f8f8
# ╠═672662f5-68a9-4c48-a3b8-45a9fb5683d0
# ╠═7d536c8b-d7c7-42cb-a086-2d76b3748380
# ╠═52c91271-3317-427e-9472-28ff7e9a46c9
# ╠═ac6bc167-4e00-44df-b32a-397eb8f2300f
# ╠═d5035c67-6f2e-4135-a9d1-bcd84a0f2c27
# ╠═f56b4ce6-47d5-48a4-8da8-4e068b2e4da9
# ╠═b6d5e7db-6e65-4120-901d-bcfce89b1537
# ╠═708e963c-f42c-415e-8d6b-f6836c878d87
# ╠═c87569b4-4ffd-44c2-9216-5193ddbcc736
# ╠═e2cb039e-8dd7-4983-965a-f52b7b718029
# ╠═749b06ac-a43b-4d1f-ba66-6a9b11ce5038
# ╠═8367cd69-1b11-4ac4-8b98-a3de8c5ad664
# ╠═a4751773-afb8-4518-a639-e5c9aa1e1765
# ╠═9ea0a845-e839-411c-99ad-051979516cb5
# ╠═725002c7-b937-4240-bbb0-22bd28972e90
# ╠═35650746-09b3-4529-ab73-28769d178437
# ╠═3ec16ac1-dbec-4ce9-af4d-a913bdc68222
# ╠═cc423cb8-7e32-4b76-b3e6-9fab041f2d1e
# ╠═a50423b1-d847-4429-9bb9-9d8a7523a459
# ╠═884239d8-5613-44bc-a14c-348b3e31d3c3
# ╠═866b5d33-c2a7-4400-845e-50575053be71
# ╠═9193fd14-d45c-4534-a029-75ab8171ea40
# ╠═7a74c048-770c-4782-9c19-ab800f66b7e1
# ╠═96a61659-2925-44ce-ae5a-a2e0c6f7b0f4
# ╠═9054c2f0-7050-4bb9-bd10-da57ec7fa69c
# ╠═9b8efe40-93b3-4dc9-8a01-ac06eae5024e
# ╠═99b9d173-618d-4fd2-9bda-89d6e991711b
# ╠═095e0c0f-2af2-4dc1-b852-bb554d1d8af1
# ╠═dc4ee9d0-7c6a-49fd-aa99-1fe3b5cb3f87
# ╠═d891c628-9785-400b-bdc4-bbc9831850f9
# ╠═0e079cb8-3014-4149-ab6d-38f4ed3c256d
# ╠═0697fd84-432e-401b-97c3-dc56b0cccbff
# ╠═e48d049b-4297-4c1e-9162-85f947d1dc20
# ╠═91c41c2e-5ae1-4cb2-86b9-d4fb953570d1
# ╠═2f1d642b-d646-416d-9550-bbe30f2044a1
# ╠═53df2546-878c-4027-bb5b-6f8a4af5bb5f
# ╠═b94447fc-762a-43b8-9e0c-987064a62618
# ╠═31b0f531-c921-46b2-aada-87098aa16d97
# ╠═af54cf26-5304-4250-a21c-61b08871defa
# ╠═71667c1d-cf1a-43a6-b3a6-ffafd87fc069
# ╠═eee69a22-1144-4215-ac6b-ce31937a60bd
# ╠═7c45dd08-df7b-44aa-a8c4-4895d73e4bd9
# ╠═ec46df8d-4a25-4fe6-b82b-9292c8d345ef
# ╠═3c8398c3-9186-4cb3-89fd-95edd5f94d94
# ╠═e25b7e46-9ea3-4210-a661-b55e69b8fc2c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
