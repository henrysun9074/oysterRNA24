#!/bin/bash -e
#SBATCH --job-name=featureCounts_array
#SBATCH --partition=schultzlab
#SBATCH --array=1-18
#SBATCH --output=/work/clh162/henry/logs/featureCounts_%a.out
#SBATCH --error=/work/clh162/henry/logs/featureCounts_%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module ##
module load Subread

## Set paths ##
GENOME=/work/clh162/OysterRNA24/hisat2_align/Cv_genome_RU_2025_shared
BAM_DIR=/work/clh162/OysterRNA24/hisat2_align/alignedreads_bam
COUNT_DIR=/work/clh162/henry/counts
mkdir -p ${COUNT_DIR}

## Set up direction/path to each sample ##
SAMPLES=($(ls ${BAM_DIR}/*_sorted.bam | sed 's/_sorted.bam//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SLURM_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$SLURM_INDEX]}


## Run featureCounts for gene quantification ##
echo "Running featureCounts for sample: ${SAMPLE}"
echo "Using file: ${BAM_DIR}/${SAMPLE}_sorted.bam"

# Hi Calista Callie Hundley
# Here are the changes I made:
# Removed -v so it actually runs instead of just printing the version name. That's why it was completing in 10s, it was just printing Subread v2.3 without running anything.
# Changed -o output name so each array task writes its own unique file (previous was all writing to the same file).
# Changed -a from gff to gtf so it reads the proper format. It was throwing an error with the gff. It looks like what you did to convert gff to gtf worked!
featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -p --countReadPairs -B \
    -a ${GENOME}/*.gtf \
    -o ${COUNT_DIR}/${SAMPLE}_gene_counts.txt \
    ${BAM_DIR}/${SAMPLE}_sorted.bam

echo "featureCounts completed for sample: ${SAMPLE}"