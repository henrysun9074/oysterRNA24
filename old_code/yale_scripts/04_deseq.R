if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("biomaRt")
# BiocManager::install("genefilter")
# BiocManager::install("pheatmap")
# BiocManager::install("WGCNA")
# BiocManager::install("flashClust")

# all of these packages sans ggplot2 install via Biocmanager
library(ggplot2)
library(DESeq2)
library(apeglm)
library(biomaRt)
library(goseq)
library(genefilter)
library(pheatmap)
library(flashClust)
library(WGCNA)

set.seed(123)

directory <- "/hpc/group/wonglab/OysterRNAseq2024/htseq_counts_newgenome"
setwd("/hpc/group/wonglab/OysterRNAseq2024")
sampleFiles <- list.files(directory, pattern = "counts.txt")
sampleCondition <- c(rep("CMAST", 9), rep("DAF", 9))
sampleTable <- data.frame(
  sampleName = gsub("_counts.txt", "", sampleFiles), 
  fileName = sampleFiles,
  condition = factor(sampleCondition)
)
sampleTable$condition <- relevel(sampleTable$condition, ref = "CMAST")
sampleTable$month <- sub(".*\\d{2}([A-Za-z]{3})\\d{2}.*", "\\1", sampleTable$fileName)
sampleTable$month <- factor(sampleTable$month)
head(sampleTable)

# Create the DESeq object from matrix and strip header of count file
l <- lapply(file.path(directory, sampleTable$fileName), function(x) {
  # change to read.table(x, skip = 1, header = FALSE, row.names = 1) if stripping header
  read.table(x, header = TRUE, row.names = 1)
})

countMatrix <- do.call(cbind, l)
colnames(countMatrix) <- sampleTable$sampleName

# Create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = sampleTable,
                              design = ~ condition)
# run the differential expression pipeline
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "CMAST", "DAF"))
resOrdered <- res[order(res$padj), ]
res$Symbol <- mcols(dds)$symbol

####change filesave for new runs
# write.csv(as.data.frame(resOrdered), file="/hpc/group/wonglab/henry/results/Jan19_deseq_newgenome.csv")

# plotMA 
plotMA(res, ylim=c(-5,5))
abline(h=0, lwd=1) # Add a line at zero
mtext("Higher in DAF", side=3, adj=1, line=-1.5, cex=1, col="darkred")
mtext("Higher in CMAST", side=1, adj=1, line=-1.5, cex=1, col="darkblue")

# plotPCA 
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition", "month"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = month)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() 

## heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

anno <- as.data.frame(colData(vsd)[, c("month", "condition")])
anno$month <- factor(anno$month, levels = c("Jun", "Jul", "Aug"))
anno$condition <- as.factor(anno$condition)
anno <- anno[order(anno$condition, anno$month), ]
mat_sorted <- mat[, rownames(anno)]
pheatmap(mat_sorted, 
         annotation_col = anno, 
         show_rownames = FALSE, 
         legend = FALSE,
         show_colnames = FALSE,
         cluster_cols = FALSE, # Required to keep our manual sorting
         scale = "none")

#plotMA with LFC shrinkage 
resultsNames(dds)
LFC <- lfcShrink(dds, coef = "condition_DAF_vs_CMAST", type = "apeglm")
plotMA(LFC, main = 'condition_DAF_vs_CMAST', cex = 0.5)
#histogram of p-values with and without FDR adjustment
hist(LFC$pvalue, breaks = 50, col = 'grey', main = 'condition_DAF_vs_CMAST', xlab = 'p-value')
hist(LFC$padj, breaks = 50, col = 'grey', main = 'condition_DAF_vs_CMAST', xlab = 'Adjusted p-value')

##volcano plot
par(mar = c(5, 4, 4, 4))
lfc = 2
pval = 0.05

tab = data.frame(logFC = LFC$log2FoldChange, negLogPval = -log10(LFC$padj)) #make a data frame with the log2 fold-changes and adjusted p-values

plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), main = 'condition_DAF_vs_CMAST') #replace main = with your title

#Genes with a fold-change greater than 2 and p-value<0.05:
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)

mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)



################################################################################
### GOSeq

mart <- useMart("metazoa_mart", 
                dataset = "cvgca002022765v4_eg_gene", 
                host = "https://metazoa.ensembl.org")
attrs <- getBM(attributes = c("ensembl_gene_id", "transcript_length", "go_id", "namespace_1003"), 
               mart = mart)

gene2cat <- attrs[attrs$go_id != "", c("ensembl_gene_id", "go_id")]
gene_lengths <- aggregate(transcript_length ~ ensembl_gene_id, data = attrs, max)
lengths_vec <- gene_lengths$transcript_length
names(lengths_vec) <- gene_lengths$ensembl_gene_id

# DE gene vector where 1 = significant and 0 = not
genes_vector <- as.integer(res$padj < 0.05 & !is.na(res$padj))
names(genes_vector) <- rownames(res)

# goseq requiresmtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.6, line = 0.5)

########### How many genes are differentially expressed?
attach(as.data.frame(LFC))
#The total number of DEGs with an adjusted p-value<0.05
summary(LFC, alpha=0.05)
#The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)

#Decreased expression:
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <(-2)) #fold-change greater than 2
#Increased expression
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) #fold-change greater than 2
# 61 overexpressed at DUML, 23 underexpressed that the DE vector and lengths vector have the same IDs
common_ids <- intersect(names(genes_vector), names(lengths_vec))
genes_vector <- genes_vector[common_ids]
lengths_vec <- lengths_vec[common_ids]

# account for the fact that longer genes are easier to sequence/detect
pwf <- nullp(genes_vector, bias.data = lengths_vec)
# run goseq
GO.wall <- goseq(pwf, gene2cat = gene2cat)
# multiple test correction with BH
GO.wall$padj <- p.adjust(GO.wall$over_represented_pvalue, method = "BH")

# Filter for significant terms (e.g., padj < 0.05)
enriched_GO <- GO.wall[GO.wall$padj < 0.05, ]

####################################
# testing revigo
common_ids <- intersect(rownames(res), names(lengths_vec))
res_sync <- res[common_ids, ]
lengths_sync <- lengths_vec[common_ids]

# downregulated genes
degs_down <- as.integer(res_sync$padj < 0.05 & !is.na(res_sync$padj) & res_sync$log2FoldChange < 0)
names(degs_down) <- rownames(res_sync)
pwf.dn <- nullp(degs_down, bias.data = lengths_sync)
go.results.dn <- goseq(pwf.dn, gene2cat = gene2cat)
go.results.dn$padj <- p.adjust(go.results.dn$over_represented_pvalue, method = "BH")

# upregulated genes
degs_up <- as.integer(res_sync$padj < 0.05 & !is.na(res_sync$padj) & res_sync$log2FoldChange > 0)
names(degs_up) <- rownames(res_sync)
pwf.up <- nullp(degs_up, bias.data = lengths_sync)
go.results.up <- goseq(pwf.up, gene2cat = gene2cat)
go.results.up$padj <- p.adjust(go.results.up$over_represented_pvalue, method = "BH")

# export for revigo http://revigo.irb.hr/ 
setwd("/hpc/group/wonglab/henry/results")
# write.table(go.results.up[go.results.up$padj < 0.05, 1:2], 
#             'Jan15-Oyster-GO-up-0.05.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
# write.table(go.results.dn[go.results.dn$padj < 0.05, 1:2], 
#             'Jan15-Oyster-GO-down-0.05.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
## same results - downregulated -- protein refolding, response to heat, unfolded protein binding 


################################################################################
## WGCNA

datExpr <- as.data.frame(t(assay(vsd)))

# gsg removes genes with too many missing values
gsg <- goodSamplesGenes(datExpr)
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}


####
# call network topology analysis function - this takes a while

# powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 3)
sft <- pickSoftThreshold(datExpr)
                         
# plot the results to find the "elbow"
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Model Fit,signed R^2",
     main = "Scale independence")
abline(h=0.90, col="red")

softPower <- 6 
net <- blockwiseModules(datExpr, power = softPower,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE, verbose = 3)

# plot dendrogram
moduleColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

traits <- colData(vsd)[, c("condition", "month")]
traits$condition <- as.numeric(as.factor(traits$condition))
traits$month <- as.numeric(as.factor(traits$month))

# ME = module eigengene
MEs <- net$MEs
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")

#plot
textMatrix_filtered <- paste(signif(moduleTraitCor, 2), "\n(",
                             signif(moduleTraitPvalue, 1), ")", sep = "")
# Replace text in non-significant cells with an empty string
textMatrix_filtered[moduleTraitPvalue > 0.05] <- ""

pdf_file <- "WGCNA_ModuleTrait_Heatmap.pdf"
pdf(file = pdf_file, width = 8, height = 12)

# set margins manually
# par(mar = c(6, 12, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_filtered,  # Using the significant-only version
  setStdMargins = FALSE,             # Crucial: allows our par(mar) to take effect
  cex.text = 0.5,                    # Size of correlation/p-value text
  cex.lab.y = 0.6,                   # Size of Module names
  cex.lab.x = 0.8,                   # Size of Trait names
  zlim = c(-1, 1),
  main = "Significant Module-Trait Relationships"
)
dev.off()
cat("Heatmap saved to:", getwd(), "/", pdf_file, sep="")


###########
# GOSeq with WGCNA module with lowest p value
cond_col <- which(colnames(moduleTraitPvalue) == "condition")
best_module_idx <- which.min(moduleTraitPvalue[, cond_col])
best_module_name <- rownames(moduleTraitPvalue)[best_module_idx]
target_color <- sub("ME", "", best_module_name)

#run goseq same as before
module_genes_bool <- as.integer(moduleColors == target_color)
names(module_genes_bool) <- colnames(datExpr)

common_ids <- intersect(names(module_genes_bool), names(lengths_vec))
module_genes_sync <- module_genes_bool[common_ids]
lengths_sync <- lengths_vec[common_ids]

pwf.module <- nullp(module_genes_sync, bias.data = lengths_sync)
go.module <- goseq(pwf.module, gene2cat = gene2cat)

go.module$padj <- p.adjust(go.module$over_represented_pvalue, method = "BH")
sig_go_module <- go.module[go.module$padj < 0.05, ]
sig_go_module ## NOTHING!!