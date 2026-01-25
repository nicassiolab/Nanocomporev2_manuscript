#!/bin/bash
#PBS -l select=1:ncpus=14:mem=90gb:ngpus=0:scratch=250gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -N WT_SPK_remora
#PBS -q workq

conda activate py39 

source /home/mzdravkov/.bashrc
cd /home/mzdravkov/nanocompore


sample="WT_SPK"

pod5="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/tmp/${sample}.pod5"
aligned_bam="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/aligned/${sample}_aligned.bam"
output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA004/resquiggled/${sample}_remora.sqlite"

ref='/work/mzdravkov/gencode.v41.transcripts.fa'

tmp_output="${PBS_SCRATCHDIR}/${sample}_remora.sqlite"

uv run python /home/mzdravkov/nanocompore/nanocompore/__main__.py remora_resquiggle --ref $ref --pod5 $pod5 --bam $aligned_bam -o $tmp_output --nthreads 14 --kit RNA004

ls -lh $tmp_output

echo "Finished resquiggling. Will copy the file from scratch to ${output}"
cp $tmp_output $output
echo 'Done'


exit 0
