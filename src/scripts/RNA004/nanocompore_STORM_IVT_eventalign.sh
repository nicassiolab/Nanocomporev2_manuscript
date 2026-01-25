#!/bin/bash
#PBS -l select=1:ncpus=18:mem=90gb:ngpus=2
#PBS -l walltime=14:00:00
#PBS -j oe
#PBS -N nanocompore_STORM_IVT_RNA004_eventalign
#PBS -q workq

source /home/mzdravkov/.bashrc
cd /home/mzdravkov/nanocompore


uv run python /home/mzdravkov/nanocompore/nanocompore/__main__.py run /projects/CGS_shared/FN_shared_projects/nanocompore_v2/scripts/RNA004/STORM_IVT_eventalign.yaml

