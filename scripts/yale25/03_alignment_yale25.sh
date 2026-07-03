#!/bin/bash -e
#SBATCH --job-name=Yale25_alignment_array
#SBATCH --time=7-00:00:00
#SBATCH --array=1-18
#SBATCH --output=/work/clh162/henry/logs/Yale25_hisat2_alignment_%A_%a.out
#SBATCH --error=/work/clh162/henry/logs/Yale25_hisat2_alignment_%A_%a.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module - HISAT2 already loaded on DCC ## 
module load HISAT2

## Set Paths ## 
RAW_DIR=/work/clh162/OysterRNA24/rawreads
TRIMMED_DIR=/work/clh162/OysterRNA24/trimmedreads
INDEX_DIR=/work/clh162/henry/ref/yale25/index
ALIGNED_DIR=/work/clh162/henry/results/yale25/aligned
mkdir -p ${ALIGNED_DIR} 

## Set up direction/path to each sample ##
# Make list of trimmed sample names (without _R1/_R2 suffix)
SAMPLES=($(ls ${RAW_DIR}/*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | xargs -n 1 basename))

# Index an individual sample from the list for this array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
# List the sample to see if naming the correct thing 
echo "Sample being processed here: " ${SAMPLE}

# Define R1 and R2 for the sample 
R1=${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz
R2=${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz
# List R1 and R2 to see if naming the correct thing 
echo "Path to R1 is " ${R1} "and path to R2 is " ${R2}

## Run Alignment ##
echo "Aligning sample:" ${SAMPLE}
hisat2 \
    -p ${SLURM_CPUS_PER_TASK} \
    -x ${INDEX_DIR}/c.virginica_yale25_HFM_index \
    -1 ${R1} \
    -2 ${R2} \
    -S ${ALIGNED_DIR}/${SAMPLE}.sam \
    --summary-file ${ALIGNED_DIR}/${SAMPLE}_hisat2_summary.txt

echo "Alignment of " ${SAMPLE} "complete!"
