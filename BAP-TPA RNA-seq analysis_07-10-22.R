# 
# Author: Hsiao-Chen Kuo 
# Email: dine820410@gmail.com
# 

# 1. CRAN packages

#install.packages("ggplot2")
#install.packages("pheatmap")
#install.packages("PoiClaClu")
#install.packages("dplyr")

# 2. Install Bioconductor packages
#BiocManager::install("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("apeglm")
#BiocManager::install("ashr")


library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(apeglm)
library(ashr)


# ------- Prepare feature counts matrix ------- 

# Load data from a csv file
setwd("C:/Rutgers/Github/BAP-TPA_UA_skin-in vivo")
data <- read.csv("featurecounts.results_20200508-150414_transformed.csv", sep = ",", header = TRUE, row.names=1)
str(data)
head(data)

# only keep the columns with feature counts of each sample
colnames(data)
data <- data[,c(6:35)] 
head(data)

# Match column names from count table with the row names from design table
colnames(data) <- c(1:30)
colnames(data)

# Transform data into a matrix
cts <-as.matrix(data)
head(cts)
str(cts)


## ------- Prepare design table -------
# design table includes the treatment group, treatment duration, or sample replicates by which you want to compare
# design table can include more than one condition by which you want to compare

coldata <- read.csv("design table_RNA-seq-07-05-22.csv", sep = ",", header = TRUE, row.names = 1)
coldata
coldata <- coldata[,c("trt","wk")] # only keep the columns of conditions by which you want to compare
coldata$trt <- factor(coldata$trt)
coldata$wk <- factor(coldata$wk)

## double check if row names of the design matches the column names of the feature counts matrix ##
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

# if not, coerce the row names of the design to match the column names of the feature counts matrix
# cts <- cts[, rownames(coldata)]
# all(rownames(coldata) == colnames(cts))

## Transform data into a matrix
cts <-as.matrix(data)
head(cts)
str(cts)

## ------- Prepare DESeq DataSet by DESeqDataSetFromMatrix -------
# 1. pre-filter the row, to keep only rows that: 
# (1) do not have 0 ##
# (2) have at least 50 reads in total ##

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ wk + trt + wk:trt)
#(1) 
dds <- dds[apply(counts(dds), 1, function(row) all(row !=0)),]
#(2)
dds <- dds[rowSums(counts(dds)) >= 50,]

dds2 <- dds # copy dds to dds2 to set the ref to BAP+TPA

# 2. explicitly set the factors level and specify the reference level 

#(1) set ref to Control
dds$trt <- factor(dds$trt, levels = c("Control", "BAP+TPA", "BAP+TPA+UA"))
dds$trt <- relevel(dds$trt, ref = "Control")

#(2) set ref to BAP+TPA

dds2$trt <- factor(dds2$trt, levels = c("Control", "BAP+TPA", "BAP+TPA+UA"))
dds2$trt <- relevel(dds2$trt, ref = "BAP+TPA")


## ------- Plot the heatmap of Poisson distances between samples ------- 

# *Background info:
# Poisson distance is used to examine overall similarity between samples
# Poisson distance consider variation between sample counts when calculating distances
# Euclidian distance do NOT onsider variation between sample counts when calculating distances
# use Poisson distance for raw (not normalized) count matrix
# use Euclidean distance for data normalized by regularized-logarithm transformation (rlog) or variance stablization transfromation (vst)
# Reference: https://www.biostars.org/p/278391/

# 1. transpose the counts and transform the poisson distance into a matrix
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)

# 2. set the row names and remove the column names
rownames(samplePoisDistMatrix) <- paste(dds$trt, dds$wk, c("wk"), sep = "." )
colnames(samplePoisDistMatrix) <- NULL

colors = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(255)



# 3. Plot heatmap

# without Clustering

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         cluster_rows=F, 
         cluster_cols=F,
         col = colors)


# with Clustering

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         cluster_cols=F,
         col = colors)

## ------- Perform DESeq data analysis ------- 
# Reference: 
# Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

dds <- DESeq(dds) # ref is control
dds2 <- DESeq(dds2) # ref is BAP+TPA

res <- results(dds, alpha = 0.05)
summary(res)

res_TPA.vs.C <- results(dds, contrast=c("trt", "BAP+TPA", "Control"), alpha = 0.05)
summary(res_TPA.vs.C)
mcols(res_TPA.vs.C, use.names=TRUE)

res_UA.vs.TPA <- results(dds, contrast=c("trt", "BAP+TPA+UA", "BAP+TPA"), alpha = 0.05)
summary(res_UA.vs.TPA)

write.csv(res_TPA.vs.C, file = "deseq2_results_TPAvsC.csv")
write.csv(res_UA.vs.TPA, file = "deseq2_results_UAvsTPA.csv")



## ------- MA plot and Data Shrinkage------- 
# M stands for magnitude and A stands for average
# *Background info:
# MA plot visualizes the log2 fold change vs. average of normalized count  
# MA plot points are colored blue if Padj < 0.1 (default) -> can change to 0.05 by alpha = 0.05
# the triangle symbols are the points that outside of the y limit

# 1. MA plot of BAP+TPA vs. Control
plotMA(res_TPA.vs.C, ylim = c(-4, 4), alpha = 0.05)

# 2. MA plot of BAP+TPA+UA vs. BAP+TPA
plotMA(res_UA.vs.TPA, ylim = c(-4, 4), alpha = 0.05)

# 3. Log fold change/effect size shrinkage (for visualization and ranking)
# The shrunken log fold changes are useful for ranking and visualization, without the need for arbitrary filters on low count genes.
# The normal prior can sometimes produce too strong of shrinkage for certain datasets. 
# In DESeq2 (ver1.18): two additional adaptive shrinkage estimators, apeglm and ashr (available via the type argument of lfcShrink). 

# (1) apeglm shrinkage
# Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 10.1093/bioinformatics/bty895

# BAP+TPA vs Control
resultsNames(dds)
resLFC_TPA.vs.C_ape <- lfcShrink(dds, coef="trt_BAP.TPA_vs_Control", type="apeglm")
plotMA(resLFC_TPA.vs.C_ape, ylim = c(-4, 4), alpha = 0.05)

# BAP+TPA+UA vs BAP+TPA
resultsNames(dds2)
resLFC_UA.vs.TPA_ape <- lfcShrink(dds2, coef="trt_BAP.TPA.UA_vs_BAP.TPA", type="apeglm")
plotMA(resLFC_UA.vs.TPA_ape, ylim = c(-4, 4), alpha = 0.05)


# (2) ashr shrinkage
# Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2. 10.1093/biostatistics/kxw041

# BAP+TPA vs Control
resultsNames(dds)
resLFC_TPA.vs.C_ash <- lfcShrink(dds, coef="trt_BAP.TPA_vs_Control", type="ashr")
plotMA(resLFC_TPA.vs.C_ash, ylim = c(-4, 4), alpha = 0.05)

# using contrast gives the same results as above
# resultsNames(dds)
# resLFC_TPA.vs.C_ash2 <- lfcShrink(dds, contrast = c("trt", "BAP+TPA", "Control"), type="ashr")
# plotMA(resLFC_TPA.vs.C_ash2, ylim = c(-4, 4), alpha = 0.05)

# BAP+TPA+UA vs BAP+TPA
resultsNames(dds)

resLFC_UA.vs.TPA_ash <- lfcShrink(dds2, coef="trt_BAP.TPA.UA_vs_BAP.TPA", type="ashr")
plotMA(resLFC_UA.vs.TPA_ash, ylim = c(-4, 4), alpha = 0.05)

# this gives the same MA plot
# resLFC_UA.vs.TPA_ash <- lfcShrink(dds, contrast = c("trt", "BAP+TPA+UA", "BAP+TPA"), type="ashr")
# plotMA(resLFC_UA.vs.TPA_ash2, ylim = c(-4, 4), alpha = 0.05)
resultsNames(dds2)


##  ------- Determine Differential Expressed Genes (DEGs) ------- 

# 1. filter by 
# (1) adjust p value of 0.05
# (2) sort by log2 fold change in ascending/descending order

# sort in ascending order: the DEGs with the strongest down-regulation:

    # BAP/TPA vs Control
resSig.TPAvsC.down <- subset(resLFC_TPA.vs.C_ash, padj < 0.05)
head(resSig[ order( resSig$log2FoldChange ), ])
    # BAP/TPA+UA vs BAP/TPA
resSig.UAvsTPA.down <- subset(resLFC_UA.vs.TPA_ash, padj < 0.05)
head(resSig[ order( resSig$log2FoldChange ), ])

# sort in descending order: the DEGs with the strongest up-regulation:

    # BAP/TPA vs Control
resSig.TPAvsC.up <- subset(resLFC_TPA.vs.C_ash, padj < 0.05)
head(resSig[ order( -resSig$log2FoldChange ), ])
    # BAP/TPA+UA vs BAP/TPA
resSig.UAvsTPA.up <- subset(resLFC_UA.vs.TPA_ash, padj < 0.05)
head(resSig[ order( -resSig$log2FoldChange ), ])

# 2. filter by 
# (1) adjust p value of 0.05
# (2) log2 fold change of 1 
# (3) baseMean >= 20
# (4) omit NA

resSig.TPAvsC.method2 <- resLFC_TPA.vs.C_ash[which(resLFC_TPA.vs.C_ash$padj < 0.05 & abs(resLFC_TPA.vs.C_ash$log2FoldChange) >= 1 & resLFC_TPA.vs.C_ash$baseMean >= 20), ]
resSig.TPAvsC.method2 <- na.omit(resultsSig)
resSig.TPAvsC.method2

resSig.UAvsTPA.method2 <- resLFC_UA.vs.TPA_ash[which(resLFC_UA.vs.TPA_ash$padj < 0.05 & abs(resLFC_UA.vs.TPA_ash$log2FoldChange) >= 1 & resLFC_UA.vs.TPA_ash$baseMean >= 20), ]
resSig.UAvsTPA.method2 <- na.omit(resultsSig)
resSig.UAvsTPA.method2

# mini-conclusion:
# gene of interest: 
# genes which show high log2 fold change that may be the first priority to look into
# method: filter by padj first and sort by DESeq2
# 
# can also consider some genes which have high p value (sort by padj)
# the upregulation and downregulation depends on the threshold you set


# ------- Plot Volcano Plot -------

library(EnhancedVolcano)
# x axis is the how much the genes differ from the control
# y axis means how statistically significant the change as compared to control

   # BAP/TPA vs Control
EnhancedVolcano(resLFC_TPA.vs.C_ash, lab = rownames(resLFC_TPA.vs.C_ash), 
                x = "log2FoldChange", y="padj",
                # add custom log2FC and adjusted P-value cutoffs
                pCutoff = 0.05, FCcutoff = 2, #these are Default thresholds
                # Adjust axis limits
                xlim = c(-8, 8),
                # Modify border
                border = "full", borderWidth = 1.5, borderColour = "black", 
                # remove gridlines
                gridlines.major = FALSE, gridlines.minor = FALSE, 
                # Add title
                title = "BAP+TPA versus Control", subtitle = "Differential Expressed Genes",
                # Add connector between labels and points
                drawConnectors = TRUE, widthConnectors = 0.75,
                # Customize legend position and point shape
                legendPosition = 'bottom', 
                shape = 19,
                # Select the genes to show labels
                selectLab = c('Nfe2l2','Nqo1','Hmox1','Tnfrsf1b')
                )
  # BAP/TPA+UA vs BAP/TPA
EnhancedVolcano(resLFC_UA.vs.TPA_ash, lab = rownames(resLFC_UA.vs.TPA_ash), 
                x = "log2FoldChange", y="padj",
                pCutoff = 0.05, FCcutoff = 2,
                xlim = c(-8, 8),
                border = "full", borderWidth = 1.5, borderColour = "black", 
                gridlines.major = FALSE, gridlines.minor = FALSE, 
                title = "BAP+TPA versus Control", subtitle = "Differential Expressed Genes",
                drawConnectors = TRUE, widthConnectors = 0.75,
                legendPosition = 'bottom', 
                shape = 19,
                selectLab = c('Nfe2l2','Nqo1','Hmox1','Tnfrsf1b')
                )

# ------- Transform read counts -------
# transformed values can be used for visualization and clustering, not for differential expression analysis

# sample > 30: perform vst (variance stabilizing transformation) transformation
# vsdata <- vst(dds, blind = FALSE)

# smaller samples: perform regularized-logarithm transformation(rlog) transformation

# vst considers (1) size factor(depth normalization) and (2) mean dispersion across the sample,
# while rld only accounts for size factor.
# Reference:
# Michael I Love, Wolfgang Huber, Simon Anders: Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 2014, 15:550. http://dx.doi.org/10.1186/s13059-014-0550-8

rld <- rlog(dds)
colnames(dds)
rld_TPA.vs.C <- rlog(dds[,1:18])
rld_UA.vs.TPA <- rlog(dds[,7:30])

rld_w5 <- rlog(dds[,dds$wk == 5])

rld_w15 <- rlog(dds[,dds$wk == 15])

rld_w20 <- rlog(dds[,dds$wk == 20])

rld_w26 <- rlog(dds[,dds$wk == 26])

plotPCA(rld_w5, intgroup = c("trt"))
plotPCA(rld_w15, intgroup = c("trt"))
plotPCA(rld_w20, intgroup = c("trt"))
plotPCA(rld_w26, intgroup = c("trt"))


# ------- Plot Principal Component Analysis ------- 
# project variation of the data into 2D format

plotPCA(rld_TPA.vs.C, intgroup = c("trt"))
plotPCA(rld_UA.vs.TPA, intgroup = c("trt"))


plotPCA(rld, intgroup = c("trt"))

plotPCA(rld, intgroup = c("wk"))
plotPCA(rld, intgroup = c("repl"))

plotPCA(rld, intgroup = c("trt","repl"))
plotPCA(rld, intgroup = c("trt","wk"))
?plotPCA

##
pcaData <- plotPCA(rld, intgroup=c("trt", "wk"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$trt <- factor(pcaData$trt, levels = c("Control","BAP+TPA", "BAP+TPA+UA"))
ggplot(pcaData, aes(PC1, PC2, color=trt, shape=wk)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size =1.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size =10, face = "bold")) +
  scale_shape_manual(values = c(16, 17, 15, 10))


# ------- Plot Dispersion of dds -------
# plot variability between replicates as function of normalized read count
# as read counts gets higher, there is less variation
# red line below 1 is preferred, and the lower the better

plotDispEsts(dds)

# ------- Plot DEGs -------
# cannot tell the reference group from here

head(resultsSig[order(resultsSig$log2FoldChange ), ] ) # DEGs downregulated compared to reference group
tail(resultsSig[order(resultsSig$log2FoldChange ), ] ) # DEGs upregulated compared to reference group

# ------- plot Expression of Individual Genes------- 
# extract gene with largest positive log2 Fold Change 
plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup="trt", normalized = TRUE)
plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup="wk", normalized = TRUE)


# extract gene counts for a specific gene of interest
geneCounts <- plotCounts(dds, gene = "Nfe2l2", intgroup = c("trt","wk"), returnData = TRUE)
geneCounts

text <- element_text(size = 14)

ggplot(geneCounts, aes(x = trt, y = count, color = wk, group = wk)) + 
  geom_point(size = 4) + 
  geom_line(size = 1,linetype = 2)  + 
  theme_classic() +
  scale_color_manual(name = "Treatment Weeks", values=c("orangered", "limegreen", "deepskyblue", "deeppink")) + 
  theme(axis.title = text) + 
  theme(axis.text.x = element_text(color="black", size = 12), axis.text.y = element_text(color="black", size = 12)) + 
  scale_x_discrete(name="Treatment", limits=c("Control", "BAP+TPA", "BAP+TPA+UA"), labels=c("Control", "BAP+TPA", "BAP+TPA+UA")) + 
  scale_y_log10(name="Normalized count") + 
  ggtitle("Nrf2") + 
  facet_wrap(vars(wk))


ggsave(file="paired_dots.png", height=5, width=5)


#
dev.off()

# boxplot for unpaired observations
p <- ggplot(geneCounts, aes(x = dex, y = count, fill = dex)) + geom_boxplot(outlier.shape = NA) + theme_classic() + geom_point(size = 3, position=position_jitterdodge(),alpha=0.5)
p <- p + scale_fill_manual(name = "Treatment", limits=c("untrt", "trt"), labels = c("None", "Dex"), values=c("orangered", "deepskyblue"))
p <- p + theme(axis.title = text) + theme(axis.text.x = element_text(color="black", size = 12), axis.text.y = element_text(color="black", size = 12))
p <- p + scale_x_discrete(name="Treatment", limits=c("untrt", "trt"), labels=c("None", "Dex")) + scale_y_log10(name="Normalized count")  + ggtitle("ENSG00000127954")
p
ggsave(file="boxplot_with_dots.png", height=5, width=5)

geneCounts
geneCounts %>% 
  ggplot(mapping = aes(x = trt, y = count, group = trt, fill = trt)) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(), color = "black", size = 4, alpha = 0.5) +
  xlab("Treatment") +
  ylab("Normalized count") +
  labs(title = "Nrf2")+
  scale_x_discrete(name="Treatment", limits=c("Control", "BAP+TPA", "BAP+TPA+UA"), labels=c("Control", "BAP+TPA", "BAP+TPA+UA")) + 
  scale_y_continuous(trans = "log10")+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
  ) +
  theme_classic() + 
  facet_wrap(vars(wk))

# ------- plot heatmap of all DE genes ------- 
mat <- assay(rld)
mat_w5 <- assay(rld_w5)

head(mat)
colnames(resSig)
rownames(resSig)
#idx <- resultsSig$ENSEMBL
DEgenes <- mat[rownames(resSig),]

colData(rld)
annotation <- as.data.frame(colData(rld)[, c("trt","wk")])
annotation


DEgenes_w5 <- mat_w5[rownames(resultsSig), ]
annotation_w5 <- as.data.frame(colData(rld_w5)[, c("trt")])


# scale by "row" is to help visualize the relative difference between groups
# otherwise (if set scale to "none"), you will only see absolute individual gene expression levels 
pheatmap(DEgenes_w5, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation", annotation_col = annotation_w5, main="Differentially Expressed genes")

pheatmap(DEgenes, scale = "row", show_rownames = FALSE, 
         clustering_distance_rows = "correlation", 
         annotation_col = annotation, 
         main="Differentially Expressed genes", 
         cluster_cols = F,
         show_colnames=F
)

# ------- Plot Heatmap of Top 20 DE genes ------- 
orderedSig <- resSig[order(resSig$padj), ]
orderedSig
head(orderedSig, n = 20)
id1 <- mat[rownames(orderedSig),]
topDE <- mat[id1,]

orderedSig_LFC <- resSig[order(resSig$), ]
orderedSig_desc
head(orderedSig_desc, n = 20)


id1 <- mat[rownames(orderedSig),]
topDE <- mat[id1,]

id2 <- mat[rownames(orderedSig_desc),]
topDE2 <- mat[id2,]

head(topDE)

head(topDE2)
# id1 <- orderedSig$ENSEMBL
# 
# id2 <- orderedSig$SYMBOL

#topDE <- mat[id1,]

#rownames(topDE) <- id2
top20DE <- head(id1, n=20)
top20DE2 <- head(id2, n=20)

pheatmap(top20DE, scale = "row", clustering_distance_rows = "correlation", annotation_col = annotation, main="Top 20 Differentially Expressed genes")
pheatmap(top20DE, scale = "row", clustering_distance_rows = "correlation", annotation_col = annotation, main="Top 20 Differentially Expressed genes", cluster_cols = F)

pheatmap(top20DE2, scale = "row", clustering_distance_rows = "correlation", annotation_col = annotation, main="Top 20 Differentially Expressed genes", cluster_cols = F)


##
orderedSig <- res_TPA.vs.C[order(res_TPA.vs.C$padj),]

DGEgenes <- rownames(subset(orderedSig, padj <0.05))

rlog.dge <- rld[DGEgenes,] %>% assay

pheatmap(rlog.dge, scale = "row", annotation_col = annotation, show_rownames = F, show_colnames = F, cluster_cols = F, treeheight_row = 0)

##
orderedSig2 <- res_TPA.vs.C[order(res_UA.vs.TPA$padj),]

DGEgenes2 <- rownames(subset(orderedSig2, padj <0.05))

rlog.dge2 <- rld[DGEgenes2,] %>% assay

pheatmap(rlog.dge2, scale = "row", annotation_col = annotation, show_rownames = F, show_colnames = F, cluster_cols = F, treeheight_row = 0)

## ashr
orderedSig3 <- resLFC_TPA.vs.C_ash[order(resLFC_TPA.vs.C_ash$padj),]

DGEgenes3 <- rownames(subset(orderedSig3, padj <0.05))

rlog.dge3 <- rld[DGEgenes3,] %>% assay

pheatmap(rlog.dge3, scale = "row", annotation_col = annotation, show_rownames = F, show_colnames = F, cluster_cols = F, treeheight_row = 0)
