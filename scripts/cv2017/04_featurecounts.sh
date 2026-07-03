#!/bin/bash -e
#SBATCH --job-name=cv2017_featureCounts_array
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/clh162/henry/logs/featureCounts_gtf_%A.out
#SBATCH --error=/work/clh162/henry/logs/featureCounts_gtf_%A.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module ##
module load Subread

## Set paths ## 
GENOME=/work/clh162/henry/ref/cv2017
BAM_DIR=/work/clh162/henry/results/cv2017/aligned_bam
COUNT_DIR=/work/clh162/henry/results/cv2017/featureCounts_counts
mkdir -p ${COUNT_DIR}

## Set up direction/path to each sample ##
# Make list of trimmed sample names (without _R1/_R2 suffix)
SAMPLES=($(ls ${BAM_DIR}/*_sorted.bam | sed 's/_sorted.bam//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

## Run featureCounts for gene quantification ##
echo "Running featureCounts for sample" ${SAMPLE} 
echo "Using file: ${BAM_DIR}/${SAMPLE}_sorted.bam"

featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -p --countReadPairs -B \
    -a ${GENOME}/*.gtf \
    -o ${COUNT_DIR}/gene_counts.txt \
    ${BAM_DIR}/*_sorted.bam   

echo "featureCounts completed for samples in" ${BAM_DIR} "directory"