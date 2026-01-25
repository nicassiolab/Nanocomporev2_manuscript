#!/bin/bash
#PBS -l select=1:ncpus=16:mem=90gb:ngpus=2
#PBS -l walltime=120:00:00
#PBS -j oe
#PBS -N benchmark_v1_v2_nanocompore
#PBS -q workq



source /home/mzdravkov/.bashrc


conda activate nanocomporev1
echo "Starting the analysis with v1"
date
nanocompore sampcomp \
	--file_list1 /projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/WT_1_eventalign_collapsed.tsv/out_eventalign_collapse.tsv \
	--file_list2 /projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/STORM_1_eventalign_collapsed.tsv/out_eventalign_collapse.tsv \
	--label1 WT_1 \
	--label2 STORM_1 \
	--fasta /work/mzdravkov/gencode.v41.transcripts.fa \
	--nthreads 18 \
	--outpath /projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_1_STORM_1_benchmark_v1
echo "The v1 run has completed"
date


cd /home/mzdravkov/nanocompore
echo "Starting the analysis with v2"
date
uv run python /home/mzdravkov/nanocompore/nanocompore/__main__.py run /projects/CGS_shared/FN_shared_projects/nanocompore_v2/scripts/RNA002/WT_1_STORM_1_benchmark.yaml
echo "The v2 run has completed"
date

