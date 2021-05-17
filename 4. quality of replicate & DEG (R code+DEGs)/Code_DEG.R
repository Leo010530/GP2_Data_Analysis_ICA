library("DESeq2")
library("dplyr")
library("ggplot2")
library("hexbin")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("car")

# Import data
readcount=read.delim("matrix_paired.txt")
# The last column is the gene length
names(readcount)[1:12] <- c("K562_CTRL_replicate1",
                            "K562_CTRL_replicate2",
                            "K562_CTRL_replicate3",
                            "K562_CTRL_replicate4",
                            "K562_AKB_replicate1",
                            "K562_AKB_replicate2",
                            "K562_AKB_replicate3",
                            "K562_AKB_replicate4",
                            "K562_Pyruvate_replicate1",
                            "K562_Pyruvate_replicate2",
                            "K562_Pyruvate_replicate3",
                            "K562_Pyruvate_replicate4")
# Export all ensembl gene ids and convert them to ENTREZ_GENE_ID
Export(data.frame(row.names(readcount)),"gene_id.csv")

# Part 1: Data quality control for read counts

# Sum of the read counts for each sample
total.cov=apply(readcount[,-13],2,sum) # Sums read counts column-wise
coverage <- data.frame(names(readcount)[1:12],total.cov)
names(coverage)[1] <- "Sample"
row.names(coverage) <- seq(1,12,1)
g <- ggplot(coverage,aes(x=Sample,y=total.cov)) 
g <- g + geom_col()
g <- g + labs(y="log10 (total counts over genes)", 
              title = "Sum of the read counts for each sample")
g <- g + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
g

# Number of genes covered
genecovered=apply(readcount[,-13],2,function(x) sum(x > 0)) # Count all genes with a readcount >0
coverage <- data.frame(names(readcount)[1:12],genecovered)
names(coverage)[1] <- "Sample"
row.names(coverage) <- seq(1,12,1)
g <- ggplot(coverage,aes(x=Sample,y=genecovered)) 
g <- g + geom_col()
g <- g + labs(y="genes covered", 
              title = "Genes covered in each sample")
g <- g + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
g

# Part 2: replicate quality assessment 

# Load data into DESeq2

# Import control and AKB group
configure = 
  data.frame(condition=factor(c("ctrl","ctrl","ctrl","ctrl","akb","akb","akb","akb")),
            type=c("r1","r2","r3","r4","r1","r2","r3","r4"))
dds_akb = DESeqDataSetFromMatrix(countData = readcount[,1:8], 
                              colData = configure,
                              design = ~ condition)
dds_akb

# Import control and pyruvate group
configure =
 data.frame(condition=factor(c("ctrl","ctrl","ctrl","ctrl","pyru","pyru","pyru","pyru")),
            type=c("r1","r2","r3","r4","r1","r2","r3","r4"))
dds_pyru = DESeqDataSetFromMatrix(countData = readcount[,c(1:4,9:12)], 
                            colData = configure,
                            design = ~ condition)
dds_pyru

# Keep the genes that have total read count >=15
# pyru
keep = rowSums(counts(dds_pyru)) >= 15 
dds_pyru = dds_pyru[keep,]
dim(counts(dds_pyru))
# akb
keep = rowSums(counts(dds_akb)) >= 15
dds_akb = dds_akb[keep,]
dim(counts(dds_akb))

# The variance stabilizing transformation and the rlog

# vst
vsd_akb = vst(dds_akb, blind = FALSE)
head(assay(vsd_akb), 3)
vsd_pyru = vst(dds_pyru, blind = FALSE)
head(assay(vsd_pyru), 3)

# rlog
rld_akb = rlog(dds_akb, blind = FALSE)
head(assay(rld_akb), 3)
rld_pyru = rlog(dds_pyru, blind = FALSE)
head(assay(rld_pyru), 3)

dds_akb = estimateSizeFactors(dds_akb)
dds_pyru = estimateSizeFactors(dds_pyru)

df_akb = bind_rows(
  as_tibble(log2(counts(dds_akb, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd_akb)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld_akb)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df_akb)[1:2] = c("x", "y") 
ggplot(df_akb, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)+ggtitle("Scatterplot of transformed counts from two samples")

df_pyru = bind_rows(
  as_tibble(log2(counts(dds_pyru, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd_pyru)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld_pyru)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df_pyru)[1:2] = c("x", "y") 
ggplot(df_pyru, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

# Sample distance (using vst)
sampleDists_akb = dist(t(assay(vsd_akb)))
sampleDists_akb
sampleDists_pyru = dist(t(assay(vsd_pyru)))
sampleDists_pyru

sampleDistMatrix_akb = as.matrix(sampleDists_akb)
sampleDistMatrix_pyru = as.matrix(sampleDists_pyru)
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_akb,
         clustering_distance_rows = sampleDists_akb,
         clustering_distance_cols = sampleDists_akb,
         col = colors)
pheatmap(sampleDistMatrix_pyru,
         clustering_distance_rows = sampleDists_pyru,
         clustering_distance_cols = sampleDists_pyru,
         col = colors)

# PCA (using vst)
plotPCA(vsd_akb)
plotPCA(vsd_pyru)

dim(counts(dds))

# Use heatmaps to plot highly variable genes (using vst result)
# akb
topVarGenes = head(order(rowVars(assay(vsd_akb)), decreasing = TRUE), 20)
mat = assay(vsd_akb)[ topVarGenes, ]
pheatmap(mat)
# pyru
topVarGenes = head(order(rowVars(assay(vsd_pyru)), decreasing = TRUE), 20)
mat = assay(vsd_pyru)[ topVarGenes, ]
pheatmap(mat)


# Differential expression analysis with DESeq2

# Call DEGs with DESeq2
dds_pyru = DESeq(dds_pyru)
dds_akb = DESeq(dds_akb)

# genes with padj <0.05 were outputted
res_akb = results(dds_akb,alpha = 0.05)
res_akb_sig = subset(res_akb, padj < 0.05)
summary(res_akb_sig)
write.table(res_akb_sig,file = "DEseq2_output_ctrl_akb.txt",sep = "\t")
res_pyru = results(dds_pyru,alpha = 0.05)
res_pyru_sig = subset(res_pyru, padj < 0.05)
summary(res_pyru_sig)
write.table(res_pyru_sig,file = "DEseq2_output_ctrl_pyru.txt",sep = "\t")

# The following analysis is not included in the report

# MA plot
plotMA(res_akb)
plotMA(res_pyru)
 
# Plot the p-value distribution (very low expressed genes (basemean<1) are excluded)
hist(res_akb$pvalue[res_akb$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white",main = "p value distribution",xlab = "p value")

# volcano <- data.frame(res$pvalue,res$log2FoldChange)
# names(volcano) <- c("pvalue","log2FoldChange")
# threshold <- as.factor(((volcano$log2FoldChange>5)|(volcano$log2Foldchange < -5)) & (volcano$pvalue<0.05))
# g <- ggplot(volcano,aes(log2FoldChange,pvalue,colour=threshold)) + geom_point()
# g <- g + labs(title="Volcano plot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-10,10)
# g <- g + geom_vline(xintercept=c(-1.5,1.5),linetype="dotted",size=1) + geom_hline(yintercept=-log2(0.05),col="blue")
# g