#!/bin/bash -e
#SBATCH --job-name=cv2017_featureCounts
#SBATCH --time=7-00:00:00
#SBATCH --output=/work/clh162/henry/logs/featureCounts_gtf.out
#SBATCH --error=/work/clh162/henry/logs/featureCounts_gtf.err
#SBATCH --partition=common
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

## Load module ##
module load Subread

GENOME=/work/clh162/henry/ref/cv2017
BAM_DIR=/work/clh162/henry/results/cv2017/aligned/aligned_bam
COUNT_DIR=/work/clh162/henry/results/cv2017/featureCounts_counts
mkdir -p ${COUNT_DIR}

echo "Running featureCounts on all samples simultaneously..."

## switch -t exon to -t gene for ncbi file (?)
featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -p --countReadPairs -B \
    -t exon \
    -g gene_id \
    -a ${GENOME}/fixed_ncbi_annotation.gtf \
    -o ${COUNT_DIR}/counts_matrix.txt \
    ${BAM_DIR}/*_sorted.bam
