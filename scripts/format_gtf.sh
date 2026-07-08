# 2017
GENOME=/work/clh162/henry/ref/cv2017

awk 'BEGIN {FS="\t"; OFS="\t"} 
$3 == "exon" && $9 ~ /gene_id ""/ {
    # Extract the transcript_id to use as a backup
    match($9, /transcript_id "[^"]+"/, t);
    if (t[0] != "") {
        sub(/gene_id ""/, "gene_id " substr(t[0], 15), $9);
    } else {
        # If no transcript_id exists, use the product name
        match($9, /product "[^"]+"/, p);
        if (p[0] != "") {
            sub(/gene_id ""/, "gene_id " substr(p[0], 9), $9);
        }
    }
} {print}' ${GENOME}/GCF_002022765.2_C_virginica-3.0_genomic.gtf > ${GENOME}/fixed_ncbi_annotation.gtf

# 2025
GENOME2=/work/clh162/henry/ref/yale25

awk 'BEGIN {FS="\t"; OFS="\t"} 
$3 == "exon" && $9 ~ /gene_id ""/ {
    # Extract the transcript_id to use as a backup
    match($9, /transcript_id "[^"]+"/, t);
    if (t[0] != "") {
        sub(/gene_id ""/, "gene_id " substr(t[0], 15), $9);
    } else {
        # If no transcript_id exists, use the product name
        match($9, /product "[^"]+"/, p);
        if (p[0] != "") {
            sub(/gene_id ""/, "gene_id " substr(p[0], 9), $9);
        }
    }
} {print}' ${GENOME2}/GCF_053477285.1_ASM5347728v1_genomic.gtf > ${GENOME2}/fixed_ncbi_annotation_yale25.gtf

#