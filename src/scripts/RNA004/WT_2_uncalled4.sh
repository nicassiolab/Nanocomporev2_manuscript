#!/bin/bash
#PBS -l select=1:ncpus=14:mem=32gb:scratch=200gb:nodetype=cpu
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -N WT_2_uncalled4
#PBS -q workq


module load miniconda/23.1.0
source /home/mzdravkov/.bashrc
conda init
conda activate py39
cd /work/mzdravkov

sample="WT_2"

transcriptome_fasta="/work/mzdravkov/gencode.v41.transcripts.fa"
pod5="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/tmp/${sample}.pod5"
aligned_bam="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/aligned/${sample}_aligned.bam"
output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/${sample}_u4.bam"

date

echo 'Getting fasta'
fasta="${PBS_SCRATCHDIR}/${sample}.fa.gz"
samtools fasta -@ 14 -T "mv,ts,pi,sp,ns" $aligned_bam > $fasta
echo 'Finished getting fasta'
ls -lh $fasta

date

# tmp_output="${PBS_SCRATCHDIR}/${sample}_u4.bam"

uncalled4 align --ref $transcriptome_fasta --reads $pod5 --bam-in $aligned_bam --bam-out $output --kit SQK-RNA004 --flowcell FLO-PRO004RA -p 14


# date
# 
# echo "Finished resquiggling. Will copy the file from scratch to ${output}"
# cp $tmp_output $output

echo 'Done'

date


exit 0
