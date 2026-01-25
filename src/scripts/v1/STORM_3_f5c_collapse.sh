#!/bin/bash
#PBS -l select=1:ncpus=24:mem=150gb:ngpus=0:scratch=550gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -N v1_f5c_collapse_STORM_3_rna002
#PBS -q workq


source /home/mzdravkov/.bashrc
conda activate nanocomporev1


sample="STORM_3"

blow5="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/tmp/${sample}.blow5"
aligned_bam="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/aligned/${sample}_aligned.bam"
output="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/${sample}_eventalign_collapsed.tsv"
summary="/projects/CGS_shared/FN_shared_projects/nanocompore_v2/data/RNA002/resquiggled/${sample}_eventalign_summary_v1.txt"

ref='/work/mzdravkov/gencode.v41.transcripts.fa'


# We already have the fasta extracted for this sample
echo 'Getting fasta'
fasta="${PBS_SCRATCHDIR}/${sample}.fa"
samtools fasta -@ 24 -T "mv,ts,pi,sp,ns" $aligned_bam > $fasta
echo 'Finished getting fasta'
ls -lh $fasta

# tmp_summary="${PBS_SCRATCHDIR}/${sample}_eventalign_summary.txt"
# tmp_output="${PBS_SCRATCHDIR}/${sample}.sqlite"

echo 'Indexing the blow5'
/home/mzdravkov/f5c-v1.5/f5c_x86_64_linux index -t 24 --slow5 $blow5 $fasta
echo 'Finished indexing the blow5'


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
						     --threads 24 \
						     --rna \
						     --pore r9 \
						     -x desktop-mid \
						     -K 1024 \
	| nanocompore eventalign_collapse -t 24 -o $output

# ls -lh $tmp_output

# echo "Finished collapsing. Will copy the file from scratch to ${output}"
# cp $tmp_output $output
# cp $tmp_summary $summary
# echo 'Done'


exit 0
