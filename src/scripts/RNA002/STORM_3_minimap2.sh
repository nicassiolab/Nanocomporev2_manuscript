#!/bin/bash
#PBS -l select=1:ncpus=12:mem=32gb:scratch=300gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -N minimap2_STORM_3_rna002_MD
#PBS -q workq


conda activate py39 

source /home/mzdravkov/.bashrc

sample="STORM_3"

transcriptome_fasta="/work/mzdravkov/gencode.v41.transcripts.fa"
basecalled_bam="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/basecalled/${sample}_unaligned.bam"
output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/aligned/${sample}_aligned.bam"


fasta="${PBS_SCRATCHDIR}/${sample}.fasta"
samtools fasta -@ 12 -T "mv,ts,pi,sp,ns" $basecalled_bam > $fasta

minimap2 -y -t 12 -ax map-ont -L $transcriptome_fasta $fasta --MD \
	| samtools view -bh -F 2324 \
 	| samtools sort -O bam \
 	> $output

samtools index -@ 12 $output

conda deactivate

exit 0
