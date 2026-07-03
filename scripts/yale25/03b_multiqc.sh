#!/bin/bash -e
#SBATCH --job-name=yale25_alignment_multiqc
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/clh162/henry/logs/yale25_hisat2_alignment_multiqc_%A.out
#SBATCH --error=/work/clh162/henry/logs/yale25_hisat2_alignment_multiqc_%A.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Activate conda environment with Multiqc loaded ##
source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate RNA-seq

## Set paths ##
HISAT2_SUMMARY=/work/clh162/henry/results/yale25/aligned
MULTIQC_OUT=/work/clh162/henry/results/yale25/aligned/multiqc

## Run MultiQC ##
echo "Running MultiQC on alignment summary files"

multiqc ${HISAT2_SUMMARY}/*_hisat2_summary.txt -o ${MULTIQC_OUT}

echo "MultiQC complete!"

conda deactivate