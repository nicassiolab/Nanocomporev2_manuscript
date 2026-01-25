#!/bin/bash
#PBS -l select=1:ncpus=24:mem=96gb:scratch=300gb
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -N minimap2_wt_1_rna004
#PBS -q workq


conda activate py39 

source /home/mzdravkov/.bashrc

sample="WT_1"

transcriptome_fasta="/work/mzdravkov/gencode.v41.transcripts.fa"
basecalled_bam="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/basecalled/${sample}_unaligned.bam"
output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/aligned/${sample}_aligned_MM_ML.bam"


fasta="${PBS_SCRATCHDIR}/${sample}.fasta"
samtools fasta -@ 24 -T "mv,ts,pi,sp,ns,MM,ML" $basecalled_bam > $fasta

minimap2 -y -t 24 -ax map-ont -L $transcriptome_fasta $fasta --MD \
	| samtools view -bh -F 2324 \
 	| samtools sort -O bam \
 	> $output

samtools index -@ 24 $output

conda deactivate

exit 0
