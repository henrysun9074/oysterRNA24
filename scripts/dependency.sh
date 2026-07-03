################## for CV2017
cd /work/clh162/henry/scripts/cv2017

# 1. Submit sam_to_bam to run after alignment finishes
JOB_BAM=$(sbatch --dependency=afterok:49231290 03_sam_to_bam.sh | awk '{print $4}')

# 2. Submit multiqc to run after sam_to_bam finishes successfully
JOB_QC=$(sbatch --dependency=afterok:$JOB_BAM 03b_multiqc.sh | awk '{print $4}')

# 3. Submit featurecounts to run after multiqc finishes successfully
sbatch --dependency=afterok:$JOB_QC 04_featurecounts.sh

############### for Yale25
cd /work/clh162/henry/scripts/yale25

# 1. Submit sam_to_bam to run after featurecounts finishes for CV2017
JOB_BAM2=$(sbatch --dependency=afterok:49231937 03_sam_to_bam.sh | awk '{print $4}')

# 2. Submit multiqc to run after sam_to_bam finishes successfully
JOB_QC2=$(sbatch --dependency=afterok:$JOB_BAM2 03b_multiqc.sh | awk '{print $4}')

# 3. Submit featurecounts to run after multiqc finishes successfully
sbatch --dependency=afterok:$JOB_QC2 04_featurecounts.sh
