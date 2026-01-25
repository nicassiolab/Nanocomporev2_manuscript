#!/bin/bash
#PBS -l select=1:ncpus=24:mem=64gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -N nanocomporev1_rna002_seed_42
#PBS -q workq

cd /work/mzdravkov                              # Change directory to the submission directory
source /home/mzdravkov/.bashrc
module load miniconda/23.1.0                    # Load the module my_module_name to set up the environment
conda init # bash
conda activate nanocomporev1

wt1_eventalign="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/WT_1_eventalign_collapsed.tsv/out_eventalign_collapse.tsv"
wt2_eventalign="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/WT_2_eventalign_collapsed.tsv/out_eventalign_collapse.tsv"
wt3_eventalign="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/WT_3_eventalign_collapsed.tsv/out_eventalign_collapse.tsv"

storm1_eventalign="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/STORM_1_eventalign_collapsed.tsv/out_eventalign_collapse.tsv"
storm2_eventalign="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/STORM_2_eventalign_collapsed.tsv/out_eventalign_collapse.tsv"
storm3_eventalign="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/STORM_3_eventalign_collapsed.tsv/out_eventalign_collapse.tsv"

output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/nanocompore_output/WT_STORM_v1_seed_42"

ref='/work/mzdravkov/gencode.v41.transcripts.fa'

nanocompore sampcomp --file_list1 $wt1_eventalign,$wt2_eventalign,$wt3_eventalign \
		     --file_list2 $storm1_eventalign,$storm2_eventalign,$storm3_eventalign \
		     --label1 "WT" \
		     --label2 "STORM" \
		     --fasta $ref \
		     --outpath $output \
		     --comparison_methods GMM,KS \
		     --logit \
		     --min_coverage 30 \
		     --pvalue_thr 0.01 \
		     -t 24 \
		     --overwrite

conda deactivate

exit 0
