include("lib.jl")


mod_ref = DataFrame(CSV.File("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/combined_mod_ref.tsv", delim = '\t'))

ivt = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_gmm_gpu_v0.2.8/out_nanocompore_results.tsv",
				   "GMM_chi2_qvalue",
				   shift=4,
				   genomic_collapse=false)


ivt_p12 = read_results("/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/nanocompore_output/WT_IVT_eventalign_p12_motor_gmm_gpu_v0.2.8_fixed_motor_assignent//out_nanocompore_results.tsv",
					   "GMM_chi2_qvalue",
					   shift=4,
					   genomic_collapse=false)


df = innerjoin(innerjoin(ivt, ivt_p12,
							   on = [:ref_id, :pos],
							   makeunique = true),
					 mod_ref,
					 on = [:chr, :strand, :genomicPos => :pos])

pseUs = df[df.mod .== "Y", :]

fig = Figure()
ax = Axis(fig[1, 1],
		  title = "Reference pseudouridine sites",
		   xlabel = "|LORmotor| - |LOR|",
		   ylabel = "-log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒 motor) + log₁₀(𝑃-𝑣𝑎𝑙𝑢𝑒)",
		   xticks = -3:1:3)

color = map(r -> if r.predicted >= 2 && r.predicted_1 >= 2
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

colorcounts = countmap(color)

Legend(fig[1, 2],
	   [PolyElement(color = :orange),
		PolyElement(color = :blue),
		PolyElement(color = :black),
		PolyElement(color = :grey)],
       ["Becomes significant ($(colorcounts[:orange]))", "Becomes insignificant ($(colorcounts[:blue]))", "Remains significant ($(colorcounts[:black]))", "Remains insignificant ($(colorcounts[:grey]))"],
	   orientation = :vertical,
	   framevisible = false)


save("supplemental_fig5_pseudouridine.png", fig)
save("supplemental_fig5_pseudouridine.pdf", fig)
