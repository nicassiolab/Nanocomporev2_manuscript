#!/bin/bash
#PBS -l select=1:ncpus=16:mem=90gb:nodetype=gpu:scratch=450GB
#PBS -l walltime=120:00:00
#PBS -j oe
#PBS -N benchmark_v1_v2_nanocompore_collapse_WT_1
#PBS -q workq

source /home/mzdravkov/.bashrc


sample="WT_1"

blow5="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/tmp/${sample}_subset_1000.blow5"
aligned_bam="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/aligned/${sample}_aligned.bam"

ref='/work/mzdravkov/gencode.v41.transcripts.fa'

reads="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/${sample}_eventalign_collapsed.tsv/reads_1000_transcript_subset.txt"

eventalign_tsv="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/tmp/${sample}_subset_1000_eventalign.tsv"
summary="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/${sample}_eventalign_summary_subset_1000.txt"
v1_output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/${sample}_eventalign_collapsed_subset_1000_v1.tsv"
v2_output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/${sample}_eventalign_collapsed_subset_1000_v2.sqlite"

date

# echo 'Getting fasta'
full_fasta="${PBS_SCRATCHDIR}/${sample}.fa"
samtools fasta -@ 16 $aligned_bam > $full_fasta
echo 'Finished getting fasta'
ls -lh $full_fasta

fasta="${PBS_SCRATCHDIR}/${sample}_subset_1000.fa"
seqtk subseq $full_fasta $reads > $fasta

ls -lh $fasta
samtools faidx $fasta
echo "Indexed the fasta"


echo 'Indexing the blow5'
/home/mzdravkov/f5c-v1.5/f5c_x86_64_linux index -t 16 --slow5 $blow5 $fasta
echo 'Finished indexing the blow5'


echo "Running f5c eventalign"
/home/mzdravkov/f5c-v1.5/f5c_x86_64_linux eventalign --reads $fasta \
						     --bam $aligned_bam \
						     --genome $ref \
						     --slow5 $blow5 \
						     --scale-events \
						     --print-read-name \
						     --secondary=yes \
						     --min-mapq 0 \
						     --samples \
						     --summary $summary \
						     --threads 16 \
						     --rna \
						     --pore r9 \
						     -x desktop-mid \
						     -K 1024 \
	> $eventalign_tsv
echo "f5c eventalign has completed"
date

conda activate nanocomporev1


echo "Starting to collapse events with v1"
date

cat $eventalign_tsv | nanocompore eventalign_collapse -t 18 -o $v1_output

echo "Finished collapsing with v1"
date

conda activate py39 
cd /home/mzdravkov/nanocompore

echo "Starting to collapse events with v2"
date

cat $eventalign_tsv | uv run python /home/mzdravkov/nanocompore/nanocompore/__main__.py eventalign_collapse --ref $ref -o $v2_output --nthreads 17 --tmp $PBS_SCRATCHDIR


echo "Finished collapsing with v2"
date



exit 0
