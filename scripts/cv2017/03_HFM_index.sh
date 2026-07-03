#!/bin/bash -e
#SBATCH --job-name=CV2017_hisat2_HFM_index
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/clh162/henry/logs/CV2017_hisat2_HFM_index.out
#SBATCH --error=/work/clh162/henry/logs/CV2017_hisat2_HFM_index.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module - HISAT2 already loaded on DCC ## 
module load HISAT2

## Set Paths ## 
GENOME=/work/clh162/henry/ref/cv2017/GCF_002022765.2_C_virginica-3.0_genomic.fna
GTF=/work/clh162/henry/ref/cv2017/GCF_002022765.2_C_virginica-3.0_genomic.gtf
INDEX_DIR=/work/clh162/henry/ref/cv2017/index
mkdir -p ${INDEX_DIR}

## Create HFM (Hierarchical FM) index ## 
    # Aligns reads to a single reference genome 
hisat2-build \
    -p ${SLURM_CPUS_PER_TASK} \
    ${GENOME} ${INDEX_DIR}/c.virginica_2017_HFM_index 