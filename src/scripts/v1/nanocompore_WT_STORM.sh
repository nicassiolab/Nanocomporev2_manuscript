#!/bin/bash
#PBS -l select=1:ncpus=24:mem=90gb
#PBS -l walltime=120:00:00
#PBS -j oe
#PBS -N nanocompore_WT_STORM
#PBS -q workq



source /home/mzdravkov/.bashrc


conda activate nanocomporev1
echo "Starting the analysis with v1"
date
nanocompore sampcomp \
	--file_list1 /projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/WT_1_eventalign_collapsed.tsv/out_eventalign_collapse.tsv,/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/WT_2_eventalign_collapsed.tsv/out_eventalign_collapse.tsv,/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/WT_3_eventalign_collapsed.tsv/out_eventalign_collapse.tsv \
	--file_list2 /projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/STORM_1_eventalign_collapsed.tsv/out_eventalign_collapse.tsv,/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/STORM_2_eventalign_collapsed.tsv/out_eventalign_collapse.tsv,/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/STORM_3_eventalign_collapsed.tsv/out_eventalign_collapse.tsv \
	--label1 WT \
	--label2 STORM \
	--fasta /work/mzdravkov/gencode.v41.transcripts.fa \
	--nthreads 24 \
	--outpath /projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1_part2
echo "The v1 run has completed"
date

