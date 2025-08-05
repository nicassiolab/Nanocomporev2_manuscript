### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 4d706412-96c2-41f9-9c6b-be708fa6c1fd
using Pkg

# ╔═╡ 1152cca8-5b2f-11f0-0f50-ef51dcc2e765
Pkg.add(["DataFrames", "SQLite", "CairoMakie", "AlgebraOfGraphics"])

# ╔═╡ e71e0883-a160-4659-a16b-f2cbb55e93b1
begin
	using CSV
	using DataFrames
	using SQLite
	using CairoMakie
	using AlgebraOfGraphics
end

# ╔═╡ 7ff9edbd-d41b-44fb-a2fb-3309976646da
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

# ╔═╡ f72f9f01-38e0-4656-a139-d6cf2b714317
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
	

# ╔═╡ 3ed44c06-f642-446f-b263-67000654bb01
abl1_data = get_transcript_signals(dbs, "ENST00000372348.8|ENSG00000097007.19|OTTHUMG00000020813.5|OTTHUMT00000054685.3|ABL1-202|ABL1|6741|protein_coding|")

# ╔═╡ 6df7f0f7-1961-4ff1-8ea8-b3addde40d48
sample_ids = Dict(zip(keys(dbs), 1:length(dbs)))

# ╔═╡ 8e155c40-457d-43a1-a09d-083dedd05e16
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

# ╔═╡ 901b8243-6836-4a09-ad4b-8877d493025b
abl1_4838 = get_pos_data(abl1_data, 4839)

# ╔═╡ c96b3552-6b48-4cd0-a005-103f148a00dc
data(abl1_4838) * mapping(:dwell => "log10(dwell)", :intensity, color = :condition) * visual(Scatter, markersize = 6.5) |> draw(scales(
      			Color = (; palette = [:green, :orange, :purple]));
			 axis = (; title = "ABL1-202, position 4839 - STORM-insensitive m6A"))

# ╔═╡ b5262135-6641-4bc2-8ec3-efdbb830efb8
abl1_201_3670 = get_pos_data(get_transcript_signals(dbs, "ENST00000318560.6|ENSG00000097007.19|OTTHUMG00000020813.5|OTTHUMT00000054684.2|ABL1-201|ABL1|5578|protein_coding|"), 3671)

# ╔═╡ 7827e966-b151-4db4-b3b5-527de6293550
data(abl1_201_3670) * mapping(:dwell => "log10(dwell)", :intensity, color = :condition) * visual(Scatter, markersize = 6.5) |> draw(scales(
      			Color = (; palette = [:green, :orange, :purple]));
			 axis = (; title = "ABL1-201, position 3671"))

# ╔═╡ ef151ebd-f708-4472-b957-19face6a4848
ythdf1_data = get_transcript_signals(dbs, "ENST00000370339.8|ENSG00000149658.18|OTTHUMG00000032955.4|OTTHUMT00000080110.3|YTHDF1-202|YTHDF1|3198|protein_coding|")

# ╔═╡ c82f835d-7237-46b3-8607-b4c03e9c23f7
ythdf1_1043 = get_pos_data(ythdf1_data, 1044)

# ╔═╡ 7fc7b82a-59a0-4cf8-903c-92d2e7269f0c
data(ythdf1_1043) * mapping(:dwell, :intensity, color = :condition) * visual(Scatter, markersize = 5) |> draw(scales(
      			Color = (; palette = [:green, :orange, :purple])); axis = (; title = "YTHDF1-202, position 1044 - STORM-insensitive m6A"))

# ╔═╡ f15866d5-82db-4e63-8738-6dca6697ac68
begin
	slc5a6_2488 = get_pos_data(get_transcript_signals(dbs, "ENST00000310574.8|ENSG00000138074.15|OTTHUMG00000097075.10|OTTHUMT00000214194.2|SLC5A6-201|SLC5A6|3213|protein_coding|"), 2489)
	slc5a6_2488 = slc5a6_2488[slc5a6_2488.intensity .> 100, :]
end

# ╔═╡ 681a66d6-d63a-4bdd-8811-9f6d70ea797f
data(slc5a6_2488) * mapping(:dwell => "log10(dwell)", :intensity, color = :condition) * visual(Scatter, markersize = 4) |> draw(scales(
      			Color = (; palette = [:green, :orange, :purple]));
			 axis = (; title = "SLC5A6-201, position 2489"))

# ╔═╡ fe491dea-7146-4f89-a54a-71ef28164fcd
data(slc5a6_2488) * mapping(:dwell => "log10(dwell)", :intensity, color = :condition) * (AlgebraOfGraphics.density() * visual(Contour) + visual(Scatter, markersize = 2.5)) |> draw(scales(Color = (; palette = [:green, :orange, :purple]));
																																												   axis = (; title = "SLC5A6-201, position 2489"))
	
	
	# visual(Scatter, markersize = 5) |> draw(scales(
      			# Color = (; palette = [:green, :orange, :purple]));
			 # axis = (; title = "SLC5A6-201, position 2489"))

# ╔═╡ a630043b-7031-4eaa-b9ea-71fb437fb6d4
fuca2_data = get_transcript_signals(dbs, "ENST00000002165.11|ENSG00000001036.14|OTTHUMG00000015728.3|OTTHUMT00000042521.3|FUCA2-201|FUCA2|2385|protein_coding|")

# ╔═╡ e19b1b52-42a3-4fc5-8015-4323f568e230
fuca2_616 = get_pos_data(fuca2_data, 616)

# ╔═╡ 6f4d5f1e-515c-4322-946f-e387d292407f
data(fuca2_616[fuca2_616.condition .!= "IVT", :]) * mapping(:dwell => "log10(dwell)", :intensity, color = :condition) * (AlgebraOfGraphics.density() * visual(Contour) + visual(Scatter, markersize = 2.5)) |> draw(scales(Color = (; palette = [:green, :orange, :purple]));
																																												   axis = (; title = "FUCA2-201, position 616"))

# ╔═╡ 23f93df5-2387-40c6-8058-adcb87d821b9
esrra_data = get_transcript_signals(dbs, 
"ENST00000000442.11|ENSG00000173153.17|OTTHUMG00000150641.9|OTTHUMT00000319303.3|ESRRA-201|ESRRA|2274|protein_coding|")

# ╔═╡ b76ab0a9-eae8-47b6-a947-df797d45fe5b
esrra_1489 = get_pos_data(esrra_data, 1489)

# ╔═╡ 87236472-7532-4c35-bbc1-1179c0cb8a2a
data(esrra_1489[esrra_1489.condition .!= "IVT", :]) * mapping(:dwell => "log10(dwell)", :intensity, color = :condition) * (AlgebraOfGraphics.density() * visual(Contour) + visual(Scatter, markersize = 2.5)) |> draw(scales(Color = (; palette = [:green, :orange, :purple]));
																																												   axis = (; title = "ESRRA-202, position 1489"))

# ╔═╡ Cell order:
# ╠═4d706412-96c2-41f9-9c6b-be708fa6c1fd
# ╠═1152cca8-5b2f-11f0-0f50-ef51dcc2e765
# ╠═e71e0883-a160-4659-a16b-f2cbb55e93b1
# ╠═7ff9edbd-d41b-44fb-a2fb-3309976646da
# ╠═f72f9f01-38e0-4656-a139-d6cf2b714317
# ╠═3ed44c06-f642-446f-b263-67000654bb01
# ╠═6df7f0f7-1961-4ff1-8ea8-b3addde40d48
# ╠═8e155c40-457d-43a1-a09d-083dedd05e16
# ╠═901b8243-6836-4a09-ad4b-8877d493025b
# ╠═c96b3552-6b48-4cd0-a005-103f148a00dc
# ╠═b5262135-6641-4bc2-8ec3-efdbb830efb8
# ╠═7827e966-b151-4db4-b3b5-527de6293550
# ╠═ef151ebd-f708-4472-b957-19face6a4848
# ╠═c82f835d-7237-46b3-8607-b4c03e9c23f7
# ╠═7fc7b82a-59a0-4cf8-903c-92d2e7269f0c
# ╠═f15866d5-82db-4e63-8738-6dca6697ac68
# ╠═681a66d6-d63a-4bdd-8811-9f6d70ea797f
# ╠═fe491dea-7146-4f89-a54a-71ef28164fcd
# ╠═a630043b-7031-4eaa-b9ea-71fb437fb6d4
# ╠═e19b1b52-42a3-4fc5-8015-4323f568e230
# ╠═6f4d5f1e-515c-4322-946f-e387d292407f
# ╠═23f93df5-2387-40c6-8058-adcb87d821b9
# ╠═b76ab0a9-eae8-47b6-a947-df797d45fe5b
# ╠═87236472-7532-4c35-bbc1-1179c0cb8a2a
