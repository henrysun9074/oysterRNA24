#!/bin/bash -e
#SBATCH --job-name=yale25_featureCounts_array
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/clh162/henry/logs/yale25_featureCounts_gtf.out
#SBATCH --error=/work/clh162/henry/logs/yale25_featureCounts_gtf.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module ##
module load Subread

## Set paths ## 
GENOME=/work/clh162/henry/ref/yale25
BAM_DIR=/work/clh162/henry/results/yale25/aligned/aligned_bam
COUNT_DIR=/work/clh162/henry/results/yale25/featureCounts_counts
mkdir -p ${COUNT_DIR}

echo "Running featureCounts on all samples simultaneously..."

featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -p --countReadPairs -B \
    -t gene \
    -g gene_id \
    -a ${GENOME}/*.gtf \
    -o ${COUNT_DIR}/gene_counts_matrix.txt \
    ${BAM_DIR}/*_sorted.bam 