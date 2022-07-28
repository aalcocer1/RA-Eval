####Using DESeq for differential gene expression analysis####


#libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(tidyverse)
library(enrichR)
library(msigdbr)

#read Data
ews <- readRDS("EWS.rds")

#Explore data
head(assay(ews))
rowData(ews)
coldata <-colData(ews)

###Create DESeq2 object###
counts.data <- assay(ews)

all(colnames(counts.data)%in% rownames(coldata))#Checking to see if rownames in coldata match with colnames in counts data
all(colnames(counts.data) == rownames(coldata)) #Cheking to see if they are in the same order


assay(ews, "counts")
sampleInfo <- as.data.frame(colData(ews))
colnames(counts.data)
dds <- DESeqDataSetFromMatrix(countData = counts.data,
                                 colData = coldata,
                                 design = ~ condition)

#Filtering: removing rows with low gene counts
keep <- rowSums(counts(dds)) >= 10
keep
ddsMat <- dds[keep,]
head(ddsMat)
ddsMat$condition

#setting the factor level
ddsMat$condition <- relevel(ddsMat$condition, ref = "shCTR")
ddsMat$condition

###Run DESeq###
dds<- DESeq(ddsMat)

#Table of differently expressed genes
res <- results(dds)
res
sigs <- na.omit(res) #omiting the NA values
sigs <- sigs[sigs$padj <0.05,] #Showing genes with p-value <0.05
sigs
write.csv(sigs, "DEGs.csv") #writing the table as a csv
vsdata <- vst(dds, blind=FALSE)
head(vsdata)
#Explore Results
summary(sigs)

#MA plot
plotMA(sigs)

###PCA Plot###
plotPCA(vsdata, intgroup = "condition")

#Dispersion plot
plotDispEsts(dds)

#converting significant results to data frame
df <- as.data.frame(sigs)
df
#Get Gene ID's and add to data frame
tmp=gsub("\\..*","",row.names(df))
tmp
df$symbol <- mapIds(org.Hs.eg.db, keys = tmp, keytype = "ENSEMBL", column = "SYMBOL") #Adding GENE ID to data frame
#Getting normalized counts from dds
mat <- counts(dds, normalized = T)[rownames(df),]

#Get z-score for each row
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- colnames(mat)
mat.z.10 <- head(mat.z, n = 10L)
###Heatmap Plot###
tmp.1=gsub("\\..*","",row.names(mat.z.10))
tmp.1
tmp.1$symbol <- mapIds(org.Hs.eg.db, keys = tmp.1, keytype = "ENSEMBL", column = "SYMBOL")
tmp.1$symbol
rownames(mat.z.10) <- tmp.1$symbol

heatmap(mat.z.10, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z.10),
        name = "Z-score", main = "Heatmap of Expressed DEGs")

###Volcano Plot###
EnhancedVolcano(df, x = "log2FoldChange", y = "padj", lab = df$symbol) 

###Enrichment Analysis###
H <- msigdbr(species = "Homo sapiens", subcategory = "CP:KEGG") #KEGG pathway database
class(H)
?msigdbr
#Define significant genes based on P-value by looking at distribution
ggplot(df, aes(x = padj))+
  geom_histogram()

#Set Cutoff
signif <- df %>%
  filter(pvalue <= 0.01)
head(signif)
 
#Get Gene ID
signif.ensembl <- unique(signif$symbol)
H.ensembl <- dplyr::select(H, gs_name, gene_symbol)

#Run Enrichment
enrich.signif.df <- enricher(gene = signif.ensembl, TERM2GENE = H.ensembl)
head(enrich.signif.df@result) #View Results

#Writing results to CSV
write.csv(enrich.signif.df, "KEGGEnrichmentAnalysis.csv")
#Figure to summarize Results (Figure to Summarize Enrichment Analysis?)
barplot(enrich.signif.df, showCategory = 10)




