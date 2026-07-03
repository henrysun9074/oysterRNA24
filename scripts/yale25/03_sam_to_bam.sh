#!/bin/bash -e
#SBATCH --job-name=Yale25_sam_to_bam_array
#SBATCH --time=7-00:00:00
#SBATCH --array=1-18
#SBATCH --output=/work/clh162/henry/logs/Yale25_sam_to_bam_array_%A_%a.out
#SBATCH --error=/work/clh162/henry/logs/Yale25_sam_to_bam_array_%A_%a.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module - samtools already loaded on DCC ##
module load samtools
# Check version of samtools (should be samtools/1.21 -- newest version DCC has installed)
echo samtools

## Set Paths ##
SAM_DIR=/work/clh162/henry/results/yale25/aligned
BAM_DIR=/work/clh162/henry/results/yale25/aligned/aligned_bam
mkdir -p ${BAM_DIR}

## Set up direction/path to each sample ##
# Make list of trimmed sample names (without _R1/_R2 suffix)
SAMPLES=($(ls ${SAM_DIR}/*.sam | sed 's/.sam//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

## Convert from SAM to BAM file format ##
echo "Converting sample " ${SAMPLE} " from SAM to BAM..."

# First convert SAM --> BAM
samtools view -b ${SAM_DIR}/${SAMPLE}.sam -o ${BAM_DIR}/${SAMPLE}.bam
# Next sort the alignment (BAM) file
samtools sort -@ ${SLURM_CPUS_PER_TASK} ${BAM_DIR}/${SAMPLE}.bam -o ${BAM_DIR}/${SAMPLE}_sorted.bam
# Finally index on the sorted aligment file 
samtools index ${BAM_DIR}/${SAMPLE}_sorted.bam