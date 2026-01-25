#!/bin/bash
#PBS -l select=1:ncpus=18:mem=90gb:ngpus=2
#PBS -l walltime=120:00:00
#PBS -j oe
#PBS -N nanocompore_WT_STORM_RNA002_eventalign_p10_motor
#PBS -q workq

source /home/mzdravkov/.bashrc
cd /home/mzdravkov/nanocompore


uv run python /home/mzdravkov/nanocompore/nanocompore/__main__.py run /projects/CGS_shared/FN_shared_projects/nanocompore_v2/scripts/RNA002/WT_STORM_eventalign_p10_motor.yaml

