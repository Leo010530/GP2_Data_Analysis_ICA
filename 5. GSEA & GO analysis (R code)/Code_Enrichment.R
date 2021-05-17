library("ggplot2")
library("clusterProfiler")
library("reactome.db")
library("org.Hs.eg.db")
library("enrichplot")
library("DOSE")

# Prepare Rank

# ENTREZ ID ranked genes
ctrl_akb = read.table("DEseq2_output_ctrl_akb.txt",header=T,row.names=1)
ctrl_pyru = read.table("DEseq2_output_ctrl_pyru.txt",header=T,row.names=1)
convert = read.csv("ensembl2enterz.csv")
convert <- convert[!is.na(convert$To),]
names(convert)[1] <- "From"
# akb
row_ensembl <- row.names(ctrl_akb) 
row_enterz <- c()
row_ensembl_new <- c()
for (i in row_ensembl) {
   if (sum(convert[which(convert$From==i),"To"])>0){
     row_ensembl_new <- c(i,row_ensembl_new)
     row_enterz <- c(convert[which(convert$From==i),"To"][1],row_enterz)
   }
}
ctrl_akb <- ctrl_akb[row_ensembl_new,]
ctrl_akb$enterzid <- row_enterz
ctrl_akb <- ctrl_akb[sort(ctrl_akb$log2FoldChange,index.return=TRUE,decreasing = T)$ix,]
Rank_entrez_akb <- ctrl_akb$log2FoldChange
names(Rank_entrez_akb) <- ctrl_akb$enterzid
Rank_entrez_akb <- Rank_entrez_akb[!duplicated(ctrl_akb$enterzid)]

# pyru
row_ensembl <- row.names(ctrl_pyru) 
row_enterz <- c()
row_ensembl_new <- c()
for (i in row_ensembl) {
   if (sum(convert[which(convert$From==i),"To"])>0){
      row_ensembl_new <- c(i,row_ensembl_new)
      row_enterz <- c(convert[which(convert$From==i),"To"][1],row_enterz)
   }
}
ctrl_pyru <- ctrl_pyru[row_ensembl_new,]
ctrl_pyru$enterzid <- row_enterz
ctrl_pyru <- ctrl_pyru[sort(ctrl_pyru$log2FoldChange,index.return=TRUE,decreasing = T)$ix,]
Rank_entrez_pyru <- ctrl_pyru$log2FoldChange
names(Rank_entrez_pyru) <- ctrl_pyru$enterzid
Rank_entrez_pyru <- Rank_entrez_pyru[!duplicated(ctrl_pyru$enterzid)]

# ENSEMBL ID ranked genes

# akb
ctrl_akb = read.table("DEseq2_output_ctrl_akb.txt",header=T,row.names=1)
ctrl_akb <- ctrl_akb[sort(ctrl_akb$log2FoldChange,index.return=TRUE,decreasing = T)$ix,]
Rank_ensembl_akb <- ctrl_akb$log2FoldChange
names(Rank_ensembl_akb) <- row.names(ctrl_akb)

# pyru
ctrl_pyru = read.table("DEseq2_output_ctrl_pyru.txt",header=T,row.names=1)
ctrl_pyru <- ctrl_pyru[sort(ctrl_pyru$log2FoldChange,index.return=TRUE,decreasing = T)$ix,]
Rank_ensembl_pyru <- ctrl_pyru$log2FoldChange
names(Rank_ensembl_pyru) <- row.names(ctrl_pyru)

# Enrichment analysis (Ranked gene list)

# gseKEGG (not included in the report)

# akb
KEGG_akb <- gseKEGG(Rank_entrez_akb,organism = "hsa")
sortKEGG_akb <- KEGG_akb[order(KEGG_akb$enrichmentScore,decreasing = T),]
gseaplot2(KEGG_akb,"hsa03010")

# pyru
KEGG_pyru <- gseKEGG(Rank_entrez_pyru,organism = "hsa")
sortKEGG_pyru <- KEGG_pyru[order(KEGG_pyru$enrichmentScore,decreasing = T),]

# gseGO ("https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html")

go_akb <- gseGO(geneList = Rank_ensembl_akb,
            OrgDb = org.Hs.eg.db,
            keyType = "ENSEMBL")
dotplot(go_akb)

go_akb[grep("adhesion",go_akb@result$Description),]
go_pyru <- gseGO(geneList = Rank_ensembl_pyru,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL")
dotplot(go_pyru)

# gene ontology analysis (log2 fold change>1)
ego_akb <- enrichGO(
   gene          = names(Rank_entrez_akb[Rank_entrez_akb>1]),
   keyType = "ENTREZID",
   OrgDb         = org.Hs.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE)
# select apoptotic process
ego_akb@result[grep("apoptotic",ego_akb@result$Description),]

ego_pyru <- enrichGO(
   gene          = names(Rank_entrez_pyru[Rank_entrez_pyru>1]),
   keyType = "ENTREZID",
   OrgDb         = org.Hs.eg.db,
   ont           = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.01,
   qvalueCutoff  = 0.05,
   readable      = TRUE)
plotGOgraph(ego_pyru)
# select apoptotic process
ego_pyru@result[grep("apoptotic",ego_pyru@result$Description),]

# GSEA analysis (only some of the results are included in the report)

gmt <- read.gmt("msigdb.v7.4.entrez.gmt")
GSEAres_akb <- GSEA(Rank_entrez_akb,TERM2GENE=gmt)
dotplot(GSEAres_akb)
gseaplot2(GSEAres_akb,909,color="red",pvalue_table = T) 

GSEAres_pyru <- GSEA(Rank_entrez_pyru,TERM2GENE=gmt)
dotplot(GSEAres_pyru)
gseaplot2(GSEAres_pyru,909,color="red",pvalue_table = T) 

oncogenicgmt <- read.gmt("c6.all.v7.4.entrez.gmt")
Oncogenic_akb <- GSEA(Rank_entrez_akb,TERM2GENE=oncogenicgmt)
gseaplot2(Oncogenic_akb,2,color="red",pvalue_table = T) 
dotplot(Oncogenic_akb)

Oncogenic_pyru <- GSEA(Rank_entrez_pyru,TERM2GENE=oncogenicgmt)
gseaplot2(Oncogenic_pyru,2,color="red",pvalue_table = T) 


hallmarkgmt <- read.gmt("h.all.v7.4.entrez.gmt")
hallmark_akb <- GSEA(Rank_entrez_akb,TERM2GENE = hallmarkgmt)
gseaplot2(hallmark_akb,1,color="red",pvalue_table = T) # UNFOLDED_PROTEIN_RESPONSE
gseaplot2(hallmark_akb,9,color="red",pvalue_table = T) # apical junction 
gseaplot2(hallmark_akb,2,color="red",pvalue_table = T) # HALLMARK_MYC_TARGETS_V1
dotplot(hallmark_akb)
sorted_hallmark_akb <- hallmark_akb@result[sort(hallmark_akb@result$NES,index.return=TRUE,decreasing = T)$ix,]
ggplot(sorted_hallmark_akb,aes(x=ID,y=NES))+geom_col()+coord_flip()

hallmark_pyru <- GSEA(Rank_entrez_pyru,TERM2GENE = hallmarkgmt)
gseaplot2(hallmark_pyru,1,color="red",pvalue_table = T) # UNFOLDED_PROTEIN_RESPONSE
dotplot(hallmark_pyru)
sorted_hallmark_pyru <- hallmark_pyru@result[sort(hallmark_pyru@result$NES,index.return=TRUE,decreasing = T)$ix,]
ggplot(sorted_hallmark_pyru,aes(x=ID,y=NES))+geom_col()+coord_flip()

micrornagmt <- read.gmt("c3.mir.v7.4.entrez.gmt")
microrna <- GSEA(Rank,TERM2GENE = micrornagmt)
gseaplot2(microrna,3,color="red",pvalue_table = T)
