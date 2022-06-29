#######################################
#INSTALL AND LOAD LIBRARIES/PACKAGES#
#######################################
BiocManager::install("DESeq2")
BiocManager::install("lazyeval") 
BiocManager::install("ggplot2")
BiocManager::install("BiocGenerics")
BiocManager::install("parallel")
BiocManager::install("affy")
BiocManager::install("vegan")
BiocManager::install("pheatmap")
BiocManager::install("IHW")
BiocManager::install("gplots")
BiocManager::install("apeglm")
BiocManager::install("ashr")
BiocManager::install("ggrepel")
BiocManager::install("EnhancedVolcano")
BiocManager::install("cowplot")
BiocManager::install("plyr")
BiocManager::install("dplyr")
install.packages("rJava",type='source')
BiocManager::install("xlsxjars")
BiocManager::install("XLConnect")
BiocManager::install("xlsx")
BiocManager::install("matrixStats")
BiocManager::install("ggfortify")
BiocManager::install("factoextra")
BiocManager::install("tidyverse")
BiocManager::install("ggrepel")


library(BiocGenerics)
library(matrixStats)
library(parallel)
library(affy)
library(DESeq2)
library(lazyeval)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library(ggrepel)
library(IHW)
library(apeglm)
library(ashr)
library(stringr)
library(EnhancedVolcano)
library(magrittr)
library(gridExtra)
library(grid)
library(cowplot)
library(gridExtra)
library(gtable)
library(reshape2)
library(plyr)
library(dplyr)
library(rJava)
.jcall()
.jinit(parameters="-Xmx6g")
.jinit()
.jcheck()
library(xlsxjars)
library(XLConnect)
library(xlsx)

####################
#UPLOAD COUNTS DATA
####################
#setting working directory
#universitate
setwd('/Users/nadezdajanina/Desktop/LAB_02.09.20/1_PROJECTS/DPH1/ReadsCounts') 
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/ReadsCounts')

#UPLOAD COUNTS DATA
filelist = list.files(pattern = "*.txt") # Make a file list from all count text files
print(filelist)

datafr = do.call(cbind,lapply(filelist,function(fn)read.table(fn,header=FALSE, sep="\t")[,2])) #merge all secoond columns of count files

genes <- read.delim(filelist[1], head = FALSE)[,1] # for the first column
genes2 <- as.character(genes)

count0 = cbind(genes2, datafr) #merge with gene name

countdata <- count0[,-1]
rownames(countdata) <- count0[,1]

colnames(countdata)=c("M03_R1",	"M03_R2", "M03_R3", "M52_R1",	"M52_R2", "M52_R3", "WT_R1", "WT_R2", "WT_R3") # add column names

D=as.data.frame(countdata)
str(D)
is.recursive(D)
hist(as.numeric(D$WT_R1), breaks=50, col="red")
hist(as.numeric(D$M03_R1), breaks=50, col="red")

indx <- sapply(D, is.factor)
D[indx] <- lapply(D[indx], function(x) as.integer(as.character(x)))
str(D)

COUNTS=D[c(7, 8, 9, 1, 2, 3, 4, 5, 6)]


hist(D$WT_R1, breaks=20000, col="blue", border = "white", xlim=c(0,1500), main="Counts per gene", xlab="Counts (truncated axis)", ylab="Number of genes", las=1, cex.axis=0.7)

epsilon=1 # pseudo-count to avoid problems with log(0)
hist(log2((D$WT_R1) + epsilon), breaks=100, col="blue", border = "white", main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", las=1, cex.axis=0.7)
hist(log2((D$M03_R1) + epsilon), breaks=100, col="blue", border = "white", main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", las=1, cex.axis=0.7)

###############
#COUNTS PLOTS
###############
cvet=c(rep("orange", 3), rep("steelblue3", 3), rep("olivedrab", 3))

boxplot(log2(COUNTS + epsilon), pch=".", col=cvet,
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(counts +1)", main="Box plots of non-normalized log2(counts+1) per sample")


plotDensity(log2(COUNTS + epsilon), lty=1, col=cvet, lwd=2, las=1, cex.axis=0.5, main="Density plot for non-normalized counts")
grid()
legend("topright", legend=names(COUNTS), col=cvet, lwd=2, cex=0.8)

#SCATTER PLOTS
# Define a function to draw a scatter plot for a pair of variables (samples) with density colors
plotFun <- function(x,y){ 
  dns <- densCols(x,y); 
  points(x,y, col=dns, pch=".", panel.first=grid());  
  #  abline(a=0, b=1, col="brown")
}

# Plot the scatter plot for a few pairs of variables selected at random
pairs(log2(COUNTS[1:9] + epsilon), 
      panel=plotFun, lower.panel = NULL, main="Scatter plot of non-normalized log2-counts for all samples", cex.main=0.9) 



#############################
# Geting info about genes #
#############################

#GENES LENGTHS
#uni
Ath_Anno = read.table("/Users/nadezdajanina/Desktop/LAB_02.09.20/1_PROJECTS/DPH1/TAIR10_representative_exon_formated_forR.gtf", sep = "\t", header = FALSE)      
#my laptop
Ath_Anno = read.table("/Users/nadezhdayanina/Desktop/DPH1/TAIR10_representative_exon_formated_forR.gtf", sep = "\t", header = FALSE)

nrow(Ath_Anno) # 153621

head(Ath_Anno)
names(Ath_Anno) <- c("Chr","meta","type","start","end","F","G","H","geneID")                                 
Ath_Anno$size = (Ath_Anno$end - Ath_Anno$start)+1                                                         
GeneSizes = aggregate(Ath_Anno$size, by=list(Category=Ath_Anno$geneID), FUN=sum) 
names(GeneSizes) <- c("GeneID","GeneSize")
nrow(GeneSizes) #32678
length(GeneSizes$GeneSize) #32678
head(GeneSizes)

#GENES DESCRIPTION(TAIR10)
#uni
Ath_Anno_TAIR10<-read.delim("/Users/nadezdajanina/Desktop/LAB_02.09.20/1_PROJECTS/DPH1/TAIR10_functional_descriptions.txt", header=TRUE, sep = "\t", fill=TRUE, quote="")
#my laptop
Ath_Anno_TAIR10<-read.table("/Users/nadezhdayanina/Desktop/DPH1/TAIR10_functional_descriptions.txt", header=TRUE, sep = "\t", fill=TRUE, quote="")

nrow(Ath_Anno_TAIR10) #41671
Ath_Anno_TAIR10$Ath_geneID = substr(Ath_Anno_TAIR10$Model_name, 1,9)


#uni
Ath_Anno_TAIR10<-read.delim("/Users/nadezdajanina/Desktop/DPH1/Annotation.txt", header=TRUE, sep = "\t", fill=TRUE, quote="")
#my laptop
Ath_Anno_TAIR10<-read.table("/Users/nadezhdayanina/Desktop/DPH1/Annotation.txt", header=TRUE, sep = "\t", fill=TRUE, quote="")

nrow(Ath_Anno_TAIR10) #33602
names(Ath_Anno_TAIR10)[1]<-paste("Ath_geneID")

######################
# NORMALIZATION TPM #
######################
#Convert counts to transcripts per million (TPM).
#Convert a numeric matrix of features (rows) and conditions (columns) with raw feature counts to transcripts per million.
#tpm A numeric matrix normalized by library size and feature length.
#TPM: 1)normalize for gene length(RPK) 2)normalize for sequencing depth(transcript per million)

sCOUNTS <- COUNTS[ order(row.names(COUNTS)), ]
head(sCOUNTS)
pairs(log2(sCOUNTS[1:9] + epsilon), 
      panel=plotFun, lower.panel = NULL, main="Scatter plot of Non-norm. log2-sorted counts for all samples") 


LEN=as.data.frame.vector((GeneSizes$GeneSize)/1000)

# Average gene size:
mean(LEN$`(GeneSizes$GeneSize)/1000`)*1000 #1536.057 bp


# Process one column at a time
TPM <- do.call(cbind, lapply(1:ncol(sCOUNTS), function(i) {
  rate = log(sCOUNTS[,i]) - log(LEN[,1])
  denom = logSumExp(rate)
  exp(rate - denom + log(1e6))
}))
# Copy the row and column names from the original matrix.
colnames(TPM) <- colnames(sCOUNTS)
rownames(TPM) <- rownames(sCOUNTS)
TPM=as.data.frame(TPM)
head(TPM)

#Plots for normalited data
pairs(log2(TPM[1:9] + epsilon), 
      panel=plotFun, lower.panel = NULL, main="Scatter plot of normalized (TPM) log2-transfermed counts", cex.main=0.9) 

boxplot(log2(TPM + epsilon), pch=".", col=cvet,
        horizontal=TRUE, cex.axis=0.6,
        las=1, ylab="Samples", xlab="log2(counts +1)", main="Box plots of normalized (TPM) log2-transformed counts", cex.main=1)

plotDensity(log2(TPM + epsilon), lty=1, col=cvet, lwd=2, las=1, cex.axis=0.8, main="Density plot for normalized (TPM) counts", xlab="log2(counts +1)")
grid()
legend("topright", legend=names(TPM), col=cvet, lwd=2, cex=0.8)


#Additional data for TPM
colnames(TPM)
TPM_plus = as.data.frame(TPM)
TPM_plus[, c(1:9)] <- sapply(TPM_plus[, c(1:9)], as.numeric)
is.numeric(TPM_plus$WT_R3)

TPM_plus$mean_WT = rowMeans(TPM_plus[, c(1:3)], na.rm = FALSE, dims = 1)
TPM_plus$mean_M03 = rowMeans(TPM_plus[, c(4:6)], na.rm = FALSE, dims = 1)
TPM_plus$mean_M52 = rowMeans(TPM_plus[, c(7:9)], na.rm = FALSE, dims = 1)

TPM_plus$SD_WT = rowSds(as.matrix(TPM_plus[, c(1:3)]), na.rm = FALSE)
TPM_plus$SD_M03 = rowSds(as.matrix(TPM_plus[, c(4:6)]), na.rm = FALSE)
TPM_plus$SD_M52 = rowSds(as.matrix(TPM_plus[, c(7:9)]), na.rm = FALSE)

colnames(TPM_plus)
TPM_plus$Ath_geneID = substr(rownames(TPM_plus), 1,9)
TPM_plus = TPM_plus[c(16,10:12,13:15, 1:9)]
nrow(TPM_plus) #32678

setwd('/Users/nadezdajanina/Desktop/DPH1/') 
write.xlsx2(TPM_plus, "TPM_plus.xlsx")

################
# TPM threshold
################
plotDensity(log2(TPM + epsilon), lty=1, col=cvet, lwd=2, las=1, cex.axis=0.85, cex.main=1, main="Expression threshold on normalized (TPM) counts density plot", xlab="log2(counts +1)")
grid()
legend("right", legend=names(TPM), col=cvet, lwd=2, cex=0.9)
abline(v=1.3,col="red")
legend("topright", legend="Threshold (TPM=2.46)", lwd=2, cex=0.9, col="red")
threshold = 2**1.3 - 1 #1.462289 (1.46) 
colnames(TPM)

nrow(TPM) #32678
TPM_th = TPM[apply(TPM[,1:9],1,function(x) any(x>1.46)),]
nrow(TPM_th) #17413

TPM_WT_th = TPM[apply(TPM[,1:3],1,function(x) any(x>1.46)),]
nrow(TPM_WT_th) #17030
TPM_M_th = TPM[apply(TPM[,4:9],1,function(x) any(x>1.46)),]
nrow(TPM_th) #17413

TPM_th$Ath_geneID = substr(rownames(TPM_th), 1,9)
TPM_th_MM = as.data.frame(TPM_th$Ath_geneID)
names(TPM_th_MM) = "Ath_geneID"

#save
#universitate
setwd('/Users/nadezdajanina/Desktop/DPH1') 
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1')

write.table(TPM_th,"TPM_th.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.xlsx2(TPM_th, "TPM_th.xlsx")
write.table(TPM_th_MM,"TPM_th_MM.txt",sep="\t",row.names=F,col.names=T,quote=F)


#########################################
# HIERARCHICAL CLUSTERING before DESeq2
#########################################
#To check how cluster normalized reads

#1_NonStand
d <- dist(t(as.matrix(TPM)))   # find distance matrix
hc <- hclust(d)             # apply hirarchical clustering 
plot(hc, cex = 0.8, hang = -1, main="Hierarchical clustering from non-standardized TMP")                    # plot the dendrogram 

#2_Standardize TPM
TPM_stand<-decostand(log2(na.omit(TPM)+1),"standardize") #Standardize. Each value in a column is standardized to a mean of 0 and standard deviation of 1. 
TPM_stand_dist<-dist(t(TPM_stand))
hcluster<-hclust(TPM_stand_dist,method="ward.D")
plot(hcluster, cex = 0.8, hang = -1, main="Hierarchical clustering from standardized TMP")

##
TPM_stand_dist_M=as.matrix(TPM_stand_dist)
rownames(TPM_stand_dist_M) <- colnames(TPM_stand)
colnames(TPM_stand_dist_M) <- colnames(TPM_stand)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(TPM_stand_dist_M,
         clustering_distance_rows=TPM_stand_dist,
         clustering_distance_cols=TPM_stand_dist,
         col=colors, cex = 1, main = "Heatmap from standardized TPM", cex.main=1)

TPM_stand2<-decostand(log2(na.omit(TPM)+1),"max") #Maximum per column. Each value in a column is divided by the maximum of the column. This is a simple standardization that can even up the influence of each of the columns. This is effective for variables that are measured on the same scale and that have similar spread (variance). Standardized values will range from 0 to 1. Note this standardization will not address changes in the spread or variability of variables measured on different scales.  
TPM_stand_dist2<-dist(t(TPM_stand2))
hcluster2<-hclust(TPM_stand_dist2,method="ward.D")
plot(hcluster2, cex = 0.8, hang = -1, main="Hierarchical clustering from standardized (max metod) TMP", cex.main=1)



##########################
# CREATING DATASET #
##########################
sCOUNTS$geneID = substr(rownames(sCOUNTS), 1,9)
rownames(sCOUNTS) = sCOUNTS$geneID
sCOUNTS = sCOUNTS[1:9]
  
genotype = factor(c(rep("WT", 3), rep("dph1_1", 3), rep("dph1_2", 3)))
phenotype = factor(c(rep("control", 3), rep("mutant", 6)))
coldata_1 = data.frame(row.names=colnames(sCOUNTS[1:9]), genotype)
coldata_2 = data.frame(row.names=colnames(sCOUNTS[1:9]), phenotype)
#coldata_2 = data.frame(row.names=colnames(sCOUNTS[1:9]), genotype, phenotype)
all(rownames(coldata_1)==colnames(sCOUNTS))

#Make a data frame form a matrix
dds = DESeqDataSetFromMatrix(countData=sCOUNTS, colData=coldata_1, design = ~genotype)
dds2 = DESeqDataSetFromMatrix(countData=sCOUNTS, colData=coldata_2, design = ~phenotype)
#dds = DESeqDataSetFromMatrix(countData=sCOUNTS, colData=coldata_1, design = ~genotype,  tidy = TRUE, ignoreRank = FALSE)
print(dds)
is(dds) # What kind of object is it?
isS4(dds)
slotNames(dds)  # What does it contain? The list of slot names

print(dds2)
is(dds2) # What kind of object is it?
isS4(dds2)
slotNames(dds2)  # What does it contain? The list of slot names



# DISPERTION ESTIMATION by DESeq2
# Estimates genewise dispersion using maximum likelihood
# Fits a curve to capture the dependence of these estimates on the average expression strength
# Shrinks genewise values towards the curve using an empirical Bayes approach (The amount of shrinkage depends on several things including sample size. Genes with high gene-wise dispersion estimates are dispersion outliers (blue circles abouve the cloud) and they are not shrunk)
#Plot dispersion plots

ddsDefault <- DESeq(dds)
plotDispEsts(ddsDefault, main= "Dispersion estimation_default")

ddsloc <- DESeq(dds ,fitType = "local") #Test another option for fitType as parametric
plotDispEsts(ddsloc, main= "Dispersion estimation_local") #Dispersion plot

ddsPra <- DESeq(dds ,fitType = "parametric") #Test another option for fitType as parametric
plotDispEsts(ddsPra, main= "Dispersion estimation_parametric") #Dispersion plot for parametric

ddsmean <- DESeq(dds ,fitType = "mean") #Test another option for fitType as mean
plotDispEsts(ddsmean, main= "All samples_mean") #Dispersion plot for mean


dds.SF <- estimateSizeFactors(dds) #Check the size factors calculated by geometric mean
print(sizeFactors(dds.SF))

plot(sizeFactors(dds.SF), colSums(counts(dds.SF)),
     xlab = "Size factors", xlim = c(0.73,1.35), cex.lab=1.2,
     ylab = "Library sizes", main = "All samples_Correlation between size factor and library size", cex.main = 1) #Correlation between size factor and library size
abline(lm(colSums(counts(dds.SF)) ~ sizeFactors(dds.SF) + 0), col = "red")
legend("bottomright", legend="R=0.991186", lwd=2, cex=0.9, col="red")


#Correlation between size factor and library size with tissue info
df.SF <- as.data.frame(sizeFactors(dds.SF))
df.colsum <- as.data.frame(colSums(counts(dds.SF)))

plot(df.SF[c(1:3),], df.colsum[c(1:3),],
     xlab = "Size factors",  xlim = c(0.73,1.35), cex.lab=1.2,
     ylab = "Library sizes", ylim=c(16500000,31000000),col="orange") 
par(new=TRUE)
plot(df.SF[c(4:6),], df.colsum[c(4:6),],
     xlab = "Size factors", xlim = c(0.73,1.35),cex.lab=1.2,
     ylab = "Library sizes", ylim=c(16500000,31000000), col="blue")
par(new=TRUE)
plot(df.SF[c(7:9),], df.colsum[c(7:9),],
     xlab = "Size factors", xlim = c(0.73,1.35),cex.lab=1.2,
     ylab = "Library sizes", ylim=c(16500000,31000000), col="darkgreen")
abline(lm(colSums(counts(dds.SF)) ~ sizeFactors(dds.SF) + 0), col = "red")
legend(x=1.23, y=19500000, legend=c('WT', 'dph1_1', 'dph1_2'), col=c('orange', 'blue', 'darkgreen'), pch = c(1,1), box.lty=0, cex= 1)
legend("topleft", legend="R=0.991186", lwd=2, cex=1, col="red")

#Correlaton test for all samples and subsets
cor.test(sizeFactors(dds.SF), colSums(counts(dds.SF)))
#Result for all: 0.991186


##########################
# Normalization by DESeq
##########################
# Aim to make normalized counts for non-differentially expressed genes similas between samples. Do not aim to adjust count distributions between samples.
# Assume that: 1)Most genes are not differentially expressed 2)Differentially expressed genes are divided equally between up- and down-regulation.
# Do not transform data, but use normalization factors within statistical testing.
# Procedure:
# 1) Take geometric mean of gene#s counts across all samples
# 2) Divide gene's counts in a sample by the geometric mean
# 3) Take median of these ratios - sample's normalization factor

Log <- rlog(dds, blind = FALSE) #regularized logarithm transformation
sampleDists <- dist(t(assay(Log)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(Log)
colnames(sampleDistMatrix) <- colnames(Log)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Hierarchical_clustering_Cu_rlog.pdf",width=4.72441,height=4.72441,paper="special")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, cex = 1, main = "Heatmap from regularized logarithm transformation (rlog)")


dev.off()


vsd = vst(dds, blind = TRUE) # Variance stabilizing transformation

sampleDists_vsd <- dist( t( assay(vsd) ) ) #Calculate sample distance
sampleDistMatrix_vsd <- as.matrix( sampleDists_vsd )
colours_vsd = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2(sampleDistMatrix_vsd, trace="none", col=colours_vsd, main = "Heatmap from Variance stabilizing transformation (vsd)", cexCol = 1.2, cexRow = 1.2, cex.main = 0.5)

########################
##    PCA    ##
########################

#Test PCA
plotPCA(vsd, intgroup=c("genotype"))
plotPCA(Log, intgroup=c("genotype"))


####
TPM_stand_T <- data.frame(TPM_stand)
TPM_stand_T <- data.frame(t(TPM_stand_T))
TPM_stand_T <- data.frame(t(TPM_stand[-1]))
rownames(TPM_stand_T) <- colnames(TPM_stand)
colnames(TPM_stand_T) <- rownames(TPM_stand)

pca <- prcomp(TPM_stand_T)

pca.proportionvariances <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100
plot(pca$x, type="n", main="Principal components analysis_TPM-stand", xlab=paste("PC1, ", round(pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(pca.proportionvariances[2], 2), "%"))
points(pca$x, col="black", pch=16, cex=1)

autoplot(pca)

percentVar <- round(pca.proportionvariances)

ggplot(pca, aes(PC1, PC2, color=genotype)) +
  geom_point(size=1.5) + geom_text_repel(aes(label=rownames(coldata_1)), size=4, segment.color = "grey50") +
  xlab(paste0("PC1 (",percentVar[1],"%) ")) +
  ylab(paste0("PC2 (",percentVar[2],"%) ")) + 
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA plot_TPM standardize") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = pcaData$genotype,
                     values=c( "#0000FF", "#66CC00", "#FF0000")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12)) +
  theme(panel.grid.minor = element_blank())



#Visualization of PCA
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  geom_point(size=1.5) + geom_text_repel(aes(label=rownames(coldata_1)), size=4, segment.color = "grey50") +
  xlab(paste0("PC1 (",percentVar[1],"%) ")) +
  ylab(paste0("PC2 (",percentVar[2],"%) ")) + 
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA plot") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = pcaData$genotype,
                     values=c( "#0000FF", "#66CC00", "#FF0000")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12)) +
  theme(panel.grid.minor = element_blank())



#Visualization of PCA (FOR ARTICLE!)
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, shape=genotype)) +
  geom_point(size=4, color="Black", fill="white") +
  scale_shape_manual(values=c(0, 2, 19)) +
  xlab(paste0("PC1 (",percentVar[1],"%) ")) +
  ylab(paste0("PC2 (",percentVar[2],"%) ")) + 
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA plot") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = pcaData$genotype,
                     values=c( "#999999", "#66CC00", "#FF0000")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12)) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank())


##################################
# DESeq2#  ---  VARIANT 1
##################################
##############
dds1 = DESeq(dds)
result_dph1_1_WT <- results(dds1,contrast=c("genotype","dph1_1","WT"))
result_dph1_2_WT <- results(dds1,contrast=c("genotype","dph1_2","WT"))
result_dph1_1_dph1_2 <- results(dds1,contrast=c("genotype","dph1_1","dph1_2"))

#Choose significantvariants
res_dph1_1_WT_FC_2  = subset(result_dph1_1_WT, (padj<=0.05 & log2FoldChange > 1) | (padj<=0.05 & log2FoldChange < -1))
res_dph1_2_WT_FC_2  = subset(result_dph1_2_WT, (padj<=0.05 & log2FoldChange > 1) | (padj<=0.05 & log2FoldChange < -1))
res_dph1_1_dph1_2_FC_2  = subset(result_dph1_1_dph1_2, (padj<=0.05 & log2FoldChange > 1) | (padj<=0.05 & log2FoldChange < -1))
res_dph1_1_WT_FC_1.5  = subset(result_dph1_1_WT, (padj<=0.05 & log2FoldChange > 0.5849625) | (padj<=0.05 & log2FoldChange < -0.5849625))
res_dph1_2_WT_FC_1.5  = subset(result_dph1_2_WT, (padj<=0.05 & log2FoldChange > 0.5849625) | (padj<=0.05 & log2FoldChange < -0.5849625))
res_dph1_1_dph1_2_FC_1.5  = subset(result_dph1_1_dph1_2, (padj<=0.05 & log2FoldChange > 0.5849625) | (padj<=0.05 & log2FoldChange < -0.5849625))

summary(res_dph1_1_WT_FC_2) #353 genes (204/149)
summary(res_dph1_2_WT_FC_2) #336 genes (191/145)
summary(res_dph1_1_dph1_2_FC_2) #0 genes (0/0)
summary(res_dph1_1_WT_FC_1.5) #814 genes (459/355)
summary(res_dph1_2_WT_FC_1.5) #808 genes (424/384)
summary(res_dph1_1_dph1_2_FC_1.5) #4 genes (0/4)

#setting working directory
#universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results') 
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results')
write.table(res_dph1_1_WT_FC_2, "res_dph1_1_WT_FC_2.txt", sep = "\t", row.names = TRUE)
write.table(res_dph1_2_WT_FC_2, "res_dph1_2_WT_FC_2.txt", sep = "\t", row.names = TRUE)
write.table(res_dph1_1_dph1_2_FC_2, "res_dph1_1_dph1_2_FC_2.txt", sep = "\t", row.names = TRUE)
write.table(res_dph1_1_WT_FC_1.5, "res_dph1_1_WT_FC_1.5.txt", sep = "\t", row.names = TRUE)
write.table(res_dph1_2_WT_FC_1.5, "res_dph1_2_WT_FC_1.5.txt", sep = "\t", row.names = TRUE)
write.table(res_dph1_1_dph1_2_FC_1.5, "res_dph1_1_dph1_2_FC_1.5.txt", sep = "\t", row.names = TRUE)


#Adding genes info to the result tables
res_dph1_1_WT_FC_2_Ano = as.data.frame(res_dph1_1_WT_FC_2)
res_dph1_2_WT_FC_2_Ano  = as.data.frame(res_dph1_2_WT_FC_2)
res_dph1_1_WT_FC_1.5_Ano  = as.data.frame(res_dph1_1_WT_FC_1.5)
res_dph1_2_WT_FC_1.5_Ano  = as.data.frame(res_dph1_2_WT_FC_1.5)
res_dph1_1_dph1_2_FC_1.5_Ano  = as.data.frame(res_dph1_1_dph1_2_FC_1.5)

res_dph1_1_WT_FC_2_Ano$Ath_geneID = substr(rownames(res_dph1_1_WT_FC_2), 1,9)
res_dph1_2_WT_FC_2_Ano$Ath_geneID  = substr(rownames(res_dph1_2_WT_FC_2), 1,9)
res_dph1_1_WT_FC_1.5_Ano$Ath_geneID  = substr(rownames(res_dph1_1_WT_FC_1.5), 1,9)
res_dph1_2_WT_FC_1.5_Ano$Ath_geneID  = substr(rownames(res_dph1_2_WT_FC_1.5), 1,9)
res_dph1_1_dph1_2_FC_1.5_Ano$Ath_geneID  = substr(rownames(res_dph1_1_dph1_2_FC_1.5), 1,9)

res_dph1_1_WT_FC_2_Ano = merge(res_dph1_1_WT_FC_2_Ano, Ath_Anno_TAIR10, by ="Ath_geneID", all.x = TRUE)
res_dph1_2_WT_FC_2_Ano  = merge(res_dph1_2_WT_FC_2_Ano, Ath_Anno_TAIR10, by ="Ath_geneID", all.x = TRUE)
res_dph1_1_WT_FC_1.5_Ano  = merge(res_dph1_1_WT_FC_1.5_Ano, Ath_Anno_TAIR10, by ="Ath_geneID", all.x = TRUE)
res_dph1_2_WT_FC_1.5_Ano  = merge(res_dph1_2_WT_FC_1.5_Ano, Ath_Anno_TAIR10, by ="Ath_geneID", all.x = TRUE)
res_dph1_1_dph1_2_FC_1.5_Ano  = merge(res_dph1_1_dph1_2_FC_1.5_Ano, Ath_Anno_TAIR10, by ="Ath_geneID", all.x = TRUE)

#FC > |2.0|
nrow(res_dph1_1_WT_FC_2_Ano) #353
length(unique(res_dph1_1_WT_FC_2_Ano$Ath_geneID)) #353

nrow(res_dph1_2_WT_FC_2_Ano) #336
length(unique(res_dph1_2_WT_FC_2_Ano$Ath_geneID)) #336

#FC > |1.5|
nrow(res_dph1_1_WT_FC_1.5_Ano) #814
length(unique(res_dph1_1_WT_FC_1.5_Ano$Ath_geneID)) #814

nrow(res_dph1_2_WT_FC_1.5_Ano) #808
length(unique(res_dph1_2_WT_FC_1.5_Ano$Ath_geneID)) #808

nrow(res_dph1_1_dph1_2_FC_1.5_Ano) #4
length(unique(res_dph1_1_dph1_2_FC_1.5_Ano$Ath_geneID)) #4



#full results after thresholds
res_dph1_1_WT_FC_2_Ano_th = merge(res_dph1_1_WT_FC_2_Ano, TPM_th_MM, by="Ath_geneID")
res_dph1_2_WT_FC_2_Ano_th  = merge(res_dph1_2_WT_FC_2_Ano, TPM_th_MM, by="Ath_geneID")
res_dph1_1_WT_FC_1.5_Ano_th  = merge(res_dph1_1_WT_FC_1.5_Ano, TPM_th_MM, by="Ath_geneID")
res_dph1_2_WT_FC_1.5_Ano_th  = merge(res_dph1_2_WT_FC_1.5_Ano, TPM_th_MM, by="Ath_geneID")
res_dph1_1_dph1_2_FC_1.5_Ano_th  = merge(res_dph1_1_dph1_2_FC_1.5_Ano, TPM_th_MM, by="Ath_geneID")

nrow(res_dph1_1_WT_FC_2_Ano_th) #273 (from 353 genes)
nrow(res_dph1_2_WT_FC_2_Ano_th) #263 (from 336 genes)
nrow(res_dph1_1_WT_FC_1.5_Ano_th) #732 (from 814 genes)
nrow(res_dph1_2_WT_FC_1.5_Ano_th) #733 (from 808 genes)
nrow(res_dph1_1_dph1_2_FC_1.5_Ano_th) #4 (from 4 genes)



#PREPARE FILES FOR MAPMAN INPUT #Extracting crusial rows for MM
res_dph1_1_WT_FC_2_Ano_th_MM = res_dph1_1_WT_FC_2_Ano_th[c(1,3)]
res_dph1_2_WT_FC_2_Ano_th_MM  = res_dph1_2_WT_FC_2_Ano_th[c(1,3)]
res_dph1_1_WT_FC_1.5_Ano_th_MM  = res_dph1_1_WT_FC_1.5_Ano_th[c(1,3)]
res_dph1_2_WT_FC_1.5_Ano_th_MM  = res_dph1_2_WT_FC_1.5_Ano_th[c(1,3)]
res_dph1_1_dph1_2_FC_1.5_Ano_th_MM  = res_dph1_1_dph1_2_FC_1.5_Ano_th[c(1,3)]


#PREPARE FILES FOR VENNY2.1 INPUT
#Genes with LFC>0 
res_dph1_1_WT_FC_2_Ano_VennyUP = res_dph1_1_WT_FC_2_Ano[res_dph1_1_WT_FC_2_Ano$log2FoldChange > 0,]
res_dph1_2_WT_FC_2_Ano_VennyUP = res_dph1_2_WT_FC_2_Ano[res_dph1_2_WT_FC_2_Ano$log2FoldChange > 0,]
res_dph1_1_WT_FC_1.5_Ano_VennyUP = res_dph1_1_WT_FC_1.5_Ano[res_dph1_1_WT_FC_1.5_Ano$log2FoldChange > 0,]
res_dph1_2_WT_FC_1.5_Ano_VennyUP = res_dph1_2_WT_FC_1.5_Ano[res_dph1_2_WT_FC_1.5_Ano$log2FoldChange > 0,]
res_dph1_1_dph1_2_FC_1.5_Ano_VennyUP = res_dph1_1_dph1_2_FC_1.5_Ano[res_dph1_1_dph1_2_FC_1.5_Ano$log2FoldChange > 0,]
# after threshhold
res_dph1_1_WT_FC_2_Ano_th_VennyUP = res_dph1_1_WT_FC_2_Ano_th[res_dph1_1_WT_FC_2_Ano_th$log2FoldChange > 0,]
res_dph1_2_WT_FC_2_Ano_th_VennyUP = res_dph1_2_WT_FC_2_Ano_th[res_dph1_2_WT_FC_2_Ano_th$log2FoldChange > 0,]
res_dph1_1_WT_FC_1.5_Ano_th_VennyUP = res_dph1_1_WT_FC_1.5_Ano_th[res_dph1_1_WT_FC_1.5_Ano_th$log2FoldChange > 0,]
res_dph1_2_WT_FC_1.5_Ano_th_VennyUP = res_dph1_2_WT_FC_1.5_Ano_th[res_dph1_2_WT_FC_1.5_Ano_th$log2FoldChange > 0,]
res_dph1_1_dph1_2_FC_1.5_Ano_th_VennyUP = res_dph1_1_dph1_2_FC_1.5_Ano_th[res_dph1_1_dph1_2_FC_1.5_Ano_th$log2FoldChange > 0,]


#Genes with LFC<0 
res_dph1_1_WT_FC_2_Ano_VennyDOWN = res_dph1_1_WT_FC_2_Ano[res_dph1_1_WT_FC_2_Ano$log2FoldChange < 0,]
res_dph1_2_WT_FC_2_Ano_VennyDOWN = res_dph1_2_WT_FC_2_Ano[res_dph1_2_WT_FC_2_Ano$log2FoldChange < 0,]
res_dph1_1_WT_FC_1.5_Ano_VennyDOWN = res_dph1_1_WT_FC_1.5_Ano[res_dph1_1_WT_FC_1.5_Ano$log2FoldChange < 0,]
res_dph1_2_WT_FC_1.5_Ano_VennyDOWN = res_dph1_2_WT_FC_1.5_Ano[res_dph1_2_WT_FC_1.5_Ano$log2FoldChange < 0,]
res_dph1_1_dph1_2_FC_1.5_Ano_VennyDOWN = res_dph1_1_dph1_2_FC_1.5_Ano[res_dph1_1_dph1_2_FC_1.5_Ano$log2FoldChange < 0,]
# after threshhold
res_dph1_1_WT_FC_2_Ano_th_VennyDOWN = res_dph1_1_WT_FC_2_Ano_th[res_dph1_1_WT_FC_2_Ano_th$log2FoldChange < 0,]
res_dph1_2_WT_FC_2_Ano_th_VennyDOWN = res_dph1_2_WT_FC_2_Ano_th[res_dph1_2_WT_FC_2_Ano_th$log2FoldChange < 0,]
res_dph1_1_WT_FC_1.5_Ano_th_VennyDOWN = res_dph1_1_WT_FC_1.5_Ano_th[res_dph1_1_WT_FC_1.5_Ano_th$log2FoldChange < 0,]
res_dph1_2_WT_FC_1.5_Ano_th_VennyDOWN = res_dph1_2_WT_FC_1.5_Ano_th[res_dph1_2_WT_FC_1.5_Ano_th$log2FoldChange < 0,]
res_dph1_1_dph1_2_FC_1.5_Ano_th_VennyDOWN = res_dph1_1_dph1_2_FC_1.5_Ano_th[res_dph1_1_dph1_2_FC_1.5_Ano_th$log2FoldChange < 0,]


#### EXTRACTING RESULTS

#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/DEGs_full/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/DEGs_full/')

#full results
write.table(as.data.frame(res_dph1_1_WT_FC_2_Ano), file = "res_dph1_1_WT_FC_2_Ano.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_2_Ano), file = "res_dph1_2_WT_FC_2_Ano.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_WT_FC_1.5_Ano), file = "res_dph1_1_WT_FC_1.5_Ano.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_1.5_Ano), file = "res_dph1_2_WT_FC_1.5_Ano.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_dph1_2_FC_1.5_Ano), file = "res_dph1_1_dph1_2_FC_1.5_Ano.txt", sep = "\t", row.names = FALSE)

#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/DEGs_full_th/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/DEGs_full_th/')

write.table(as.data.frame(res_dph1_1_WT_FC_2_Ano_th), file = "res_dph1_1_WT_FC_2_Ano_th.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_2_Ano_th), file = "res_dph1_2_WT_FC_2_Ano_th.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_WT_FC_1.5_Ano_th), file = "res_dph1_1_WT_FC_1.5_Ano_th.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_1.5_Ano_th), file = "res_dph1_2_WT_FC_1.5_Ano_th.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_dph1_2_FC_1.5_Ano_th), file = "res_dph1_1_dph1_2_FC_1.5_Ano_th.txt", sep = "\t", row.names = FALSE)


#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/forMM_th/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/forMM_th/')

#results for MapMan input(only crusial columns and no NA)

#UP
write.table(as.data.frame(res_dph1_1_WT_FC_2_Ano_th_MM[res_dph1_1_WT_FC_2_Ano_th_MM$log2FoldChange > 0,]), file = "res_dph1_1_WT_FC_2_Ano_th_MM_UP.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_2_Ano_th_MM[res_dph1_2_WT_FC_2_Ano_th_MM$log2FoldChange > 0,]), file = "res_dph1_2_WT_FC_2_Ano_th_MM_UP.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_WT_FC_1.5_Ano_th_MM[res_dph1_1_WT_FC_1.5_Ano_th_MM$log2FoldChange > 0,]), file = "res_dph1_1_WT_FC_1.5_Ano_th_MM_UP.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_1.5_Ano_th_MM[res_dph1_2_WT_FC_1.5_Ano_th_MM$log2FoldChange > 0,]), file = "res_dph1_2_WT_FC_1.5_Ano_th_MM_UP.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_dph1_2_FC_1.5_Ano_th_MM[res_dph1_1_dph1_2_FC_1.5_Ano_th_MM$log2FoldChange > 0,]), file = "res_dph1_1_dph1_2_FC_1.5_Ano_th_MM_UP.txt", sep = "\t", row.names = FALSE)

#DOWN
write.table(as.data.frame(res_dph1_1_WT_FC_2_Ano_th_MM[res_dph1_1_WT_FC_2_Ano_th_MM$log2FoldChange < 0,]), file = "res_dph1_1_WT_FC_2_Ano_th_MM_DOWN.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_2_Ano_th_MM[res_dph1_2_WT_FC_2_Ano_th_MM$log2FoldChange < 0,]), file = "res_dph1_2_WT_FC_2_Ano_th_MM_DOWN.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_WT_FC_1.5_Ano_th_MM[res_dph1_1_WT_FC_1.5_Ano_th_MM$log2FoldChange < 0,]), file = "res_dph1_1_WT_FC_1.5_Ano_th_MM_DOWN.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_2_WT_FC_1.5_Ano_th_MM[res_dph1_2_WT_FC_1.5_Ano_th_MM$log2FoldChange < 0,]), file = "res_dph1_2_WT_FC_1.5_Ano_th_MM_DOWN.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(res_dph1_1_dph1_2_FC_1.5_Ano_th_MM[res_dph1_1_dph1_2_FC_1.5_Ano_th_MM$log2FoldChange < 0,]), file = "res_dph1_1_dph1_2_FC_1.5_Ano_th_MM_DOWN.txt", sep = "\t", row.names = FALSE)


#setting working directory
#universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/Venny/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/Venny/')

#results for Venny input(separate by LFC)
#Genes with LFC>0 
write.xlsx2(unique(res_dph1_1_WT_FC_2_Ano_VennyUP$Ath_geneID), "VennyUP_dph1_WT_FC_2.xlsx", sheetName = "dph1_1_WT_FC2", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_2_Ano_VennyUP$Ath_geneID), "VennyUP_dph1_WT_FC_2.xlsx", sheetName = "dph1_2_WT_FC2", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_2_Ano_VennyUP$Ath_geneID), unique(res_dph1_2_WT_FC_2_Ano_VennyUP$Ath_geneID)), "VennyUP_dph1_WT_FC_2.xlsx", sheetName = "overlap_dph1_1&2_WT_FC2", append=TRUE)
write.xlsx2(unique(res_dph1_1_WT_FC_1.5_Ano_VennyUP$Ath_geneID), "VennyUP_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_1_WT_FC_1.5", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_1.5_Ano_VennyUP$Ath_geneID), "VennyUP_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_2_WT_FC_1.5", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_1.5_Ano_VennyUP$Ath_geneID), unique(res_dph1_2_WT_FC_1.5_Ano_VennyUP$Ath_geneID)), "VennyUP_dph1_WT_FC_1.5.xlsx", sheetName = "overlap_dph1_1&2_WT_FC_1.5", append=TRUE)

#Genes with LFC<0 
write.xlsx2(unique(res_dph1_1_WT_FC_2_Ano_VennyDOWN$Ath_geneID), "VennyDOWN_dph1_WT_FC_2.xlsx", sheetName = "dph1_1_WT_FC2", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_2_Ano_VennyDOWN$Ath_geneID), "VennyDOWN_dph1_WT_FC_2.xlsx", sheetName = "dph1_2_WT_FC2", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_2_Ano_VennyDOWN$Ath_geneID), unique(res_dph1_2_WT_FC_2_Ano_VennyDOWN$Ath_geneID)), "VennyDOWN_dph1_WT_FC_2.xlsx", sheetName = "overlap_dph1_1&2_WT_FC2", append=TRUE)
write.xlsx2(unique(res_dph1_1_WT_FC_1.5_Ano_VennyDOWN$Ath_geneID), "VennyDOWN_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_1_WT_FC_1.5", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_1.5_Ano_VennyDOWN$Ath_geneID), "VennyDOWN_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_2_WT_FC_1.5", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_1.5_Ano_VennyDOWN$Ath_geneID), unique(res_dph1_2_WT_FC_1.5_Ano_VennyDOWN$Ath_geneID)), "VennyDOWN_dph1_WT_FC_1.5.xlsx", sheetName = "overlap_dph1_1&2_WT_FC_1.5", append=TRUE)


#setting working directory
#universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/Venny_th/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/Venny_th/')

#results for Venny input(separate by LFC, after threshfold)
#Genes with LFC>0 
write.xlsx2(unique(res_dph1_1_WT_FC_2_Ano_th_VennyUP$Ath_geneID), "VennyUP_th_dph1_WT_FC_2.xlsx", sheetName = "dph1_1_WT_FC2_th", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_2_Ano_th_VennyUP$Ath_geneID), "VennyUP_th_dph1_WT_FC_2.xlsx", sheetName = "dph1_2_WT_FC2_th", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_2_Ano_th_VennyUP$Ath_geneID), unique(res_dph1_2_WT_FC_2_Ano_th_VennyUP$Ath_geneID)), "VennyUP_th_dph1_WT_FC_2.xlsx", sheetName = "overlap_dph1_1&2_WT_FC2_th", append=TRUE)
write.xlsx2(unique(res_dph1_1_WT_FC_1.5_Ano_th_VennyUP$Ath_geneID), "VennyUP_th_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_1_WT_FC_1.5_th", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_1.5_Ano_th_VennyUP$Ath_geneID), "VennyUP_th_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_2_WT_FC_1.5_th", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_1.5_Ano_th_VennyUP$Ath_geneID), unique(res_dph1_2_WT_FC_1.5_Ano_th_VennyUP$Ath_geneID)), "VennyUP_th_dph1_WT_FC_1.5.xlsx", sheetName = "overlap_dph1_1&2_WT_FC_1.5_th", append=TRUE)

#Genes with LFC<0 
write.xlsx2(unique(res_dph1_1_WT_FC_2_Ano_th_VennyDOWN$Ath_geneID), "VennyDOWN_th_dph1_WT_FC_2.xlsx", sheetName = "dph1_1_WT_FC2_th", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_2_Ano_th_VennyDOWN$Ath_geneID), "VennyDOWN_th_dph1_WT_FC_2.xlsx", sheetName = "dph1_2_WT_FC2_th", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_2_Ano_th_VennyDOWN$Ath_geneID), unique(res_dph1_2_WT_FC_2_Ano_th_VennyDOWN$Ath_geneID)), "VennyDOWN_th_dph1_WT_FC_2.xlsx", sheetName = "overlap_dph1_1&2_WT_FC2_th", append=TRUE)
write.xlsx2(unique(res_dph1_1_WT_FC_1.5_Ano_th_VennyDOWN$Ath_geneID), "VennyDOWN_th_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_1_WT_FC_1.5_th", append=TRUE)
write.xlsx2(unique(res_dph1_2_WT_FC_1.5_Ano_th_VennyDOWN$Ath_geneID), "VennyDOWN_th_dph1_WT_FC_1.5.xlsx", sheetName = "dph1_2_WT_FC_1.5_th", append=TRUE)
write.xlsx2(intersect(unique(res_dph1_1_WT_FC_1.5_Ano_th_VennyDOWN$Ath_geneID), unique(res_dph1_2_WT_FC_1.5_Ano_th_VennyDOWN$Ath_geneID)), "VennyDOWN_th_dph1_WT_FC_1.5.xlsx", sheetName = "overlap_dph1_1&2_WT_FC_1.5_th", append=TRUE)


##################
## FULL TABLES ##
##################

#UNIVERSAL TABLE
################
dph1_1.vs.WT_FC_2 = as.data.frame(res_dph1_1_WT_FC_2)
dph1_2.vs.WT_FC_2 = as.data.frame(res_dph1_2_WT_FC_2)
dph1_1.vs.WT_FC_1.5 = as.data.frame(res_dph1_1_WT_FC_1.5)
dph1_2.vs.WT_FC_1.5 = as.data.frame(res_dph1_2_WT_FC_1.5)
dph1_1.vs.dph1_2_FC_1.5 = as.data.frame(res_dph1_1_dph1_2_FC_1.5)

dph1_1.vs.WT_FC_2$Ath_geneID = rownames(dph1_1.vs.WT_FC_2)
dph1_2.vs.WT_FC_2$Ath_geneID = rownames(dph1_2.vs.WT_FC_2)
dph1_1.vs.WT_FC_1.5$Ath_geneID = rownames(dph1_1.vs.WT_FC_1.5)
dph1_2.vs.WT_FC_1.5$Ath_geneID = rownames(dph1_2.vs.WT_FC_1.5)
dph1_1.vs.dph1_2_FC_1.5$Ath_geneID = rownames(dph1_1.vs.dph1_2_FC_1.5)


names(dph1_1.vs.WT_FC_2) = c("dph1_1.vs.WT_FC_2_baseMean", "dph1_1.vs.WT_FC_2_log2FoldChange", "dph1_1.vs.WT_FC_2_lfcSE", "dph1_1.vs.WT_FC_2_stat", "dph1_1.vs.WT_FC_2_pvalue", "dph1_1.vs.WT_FC_2_padj", "Ath_geneID")
names(dph1_2.vs.WT_FC_2) = c("dph1_2.vs.WT_FC_2_baseMean", "dph1_2.vs.WT_FC_2_log2FoldChange", "dph1_2.vs.WT_FC_2_lfcSE", "dph1_2.vs.WT_FC_2_stat", "dph1_2.vs.WT_FC_2_pvalue", "dph1_2.vs.WT_FC_2_padj", "Ath_geneID")
names(dph1_1.vs.WT_FC_1.5) = c("dph1_1.vs.WT_FC_1.5_baseMean", "dph1_1.vs.WT_FC_1.5_log2FoldChange", "dph1_1.vs.WT_FC_1.5_lfcSE", "dph1_1.vs.WT_FC_1.5_stat", "dph1_1.vs.WT_FC_1.5_pvalue", "dph1_1.vs.WT_FC_1.5_padj", "Ath_geneID")
names(dph1_2.vs.WT_FC_1.5) = c("dph1_2.vs.WT_FC_1.5_baseMean", "dph1_2.vs.WT_FC_1.5_log2FoldChange", "dph1_2.vs.WT_FC_1.5_lfcSE", "dph1_2.vs.WT_FC_1.5_stat", "dph1_2.vs.WT_FC_1.5_pvalue", "dph1_2.vs.WT_FC_1.5_padj", "Ath_geneID")
names(dph1_1.vs.dph1_2_FC_1.5) = c("dph1_1.vs.dph1_2_FC_1.5_baseMean", "dph1_1.vs.dph1_2_FC_1.5_log2FoldChange", "dph1_1.vs.dph1_2_FC_1.5_lfcSE", "dph1_1.vs.dph1_2_FC_1.5_stat", "dph1_1.vs.dph1_2_FC_1.5_pvalue", "dph1_1.vs.dph1_2_FC_1.5_padj", "Ath_geneID")


UNIVERSAL = TPM_plus
nrow(UNIVERSAL) #32678
UNIVERSAL = merge(UNIVERSAL, dph1_1.vs.WT_FC_2, by = "Ath_geneID", all.x = TRUE)
UNIVERSAL = merge(UNIVERSAL, dph1_2.vs.WT_FC_2, by = "Ath_geneID", all.x = TRUE)
UNIVERSAL = merge(UNIVERSAL, dph1_1.vs.WT_FC_1.5, by = "Ath_geneID", all.x = TRUE)
UNIVERSAL = merge(UNIVERSAL, dph1_2.vs.WT_FC_1.5, by = "Ath_geneID", all.x = TRUE)
UNIVERSAL = merge(UNIVERSAL, dph1_1.vs.dph1_2_FC_1.5, by = "Ath_geneID", all.x = TRUE)

# EXTRACTING GENES TABLES
#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results')

options(java.parameters = "-Xmx6g")  ## memory set to 6 GB
write.xlsx2(UNIVERSAL, "UNIVERSAL.xlsx", sheetName = "Variant1", append=TRUE)

##########
##########

#Genes with LFC>0
res_dph1_1_WT_FC_2_UP_table = res_dph1_1_WT_FC_2_Ano[res_dph1_1_WT_FC_2_Ano$log2FoldChange > 0,] #dph1_1.vs.WT_FC_2
res_dph1_2_WT_FC_2_UP_table = res_dph1_2_WT_FC_2_Ano[res_dph1_2_WT_FC_2_Ano$log2FoldChange > 0,] #dph1_2.vs.WT_FC_2
res_dph1_1_WT_FC_1.5_UP_table = res_dph1_1_WT_FC_1.5_Ano[res_dph1_1_WT_FC_1.5_Ano$log2FoldChange > 0,] #dph1_1.vs.WT_FC_1.5
res_dph1_2_WT_FC_1.5_UP_table = res_dph1_2_WT_FC_1.5_Ano[res_dph1_2_WT_FC_1.5_Ano$log2FoldChange > 0,] #dph1_2.vs.WT_FC_1.5
res_dph1_1_dph1_2_FC_1.5_UP_table = res_dph1_1_dph1_2_FC_1.5_Ano[res_dph1_1_dph1_2_FC_1.5_Ano$log2FoldChange > 0,] #dph1_1.vs.dph1_2_FC_1.5

#Intersection
#res_dph1_1and2_WT_FC_2
res_dph1_1and2_WT_FC_2_UP_table = intersect(res_dph1_1_WT_FC_2_UP_table$Ath_geneID, res_dph1_2_WT_FC_2_UP_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_2_UP_table) #129
res_dph1_1and2_WT_FC_2_UP_table = as.data.frame(res_dph1_1and2_WT_FC_2_UP_table)
names(res_dph1_1and2_WT_FC_2_UP_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_2_UP_table = merge(res_dph1_1and2_WT_FC_2_UP_table, res_dph1_1_WT_FC_2_Ano, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_2_UP_table = merge(res_dph1_1and2_WT_FC_2_UP_table, res_dph1_2_WT_FC_2_UP_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_2_UP_table)
res_dph1_1and2_WT_FC_2_UP_table = res_dph1_1and2_WT_FC_2_UP_table[, c(1:7,12:17,8:11)]

#res_dph1_1and2_WT_FC_1.5
res_dph1_1and2_WT_FC_1.5_UP_table = intersect(res_dph1_1_WT_FC_1.5_UP_table$Ath_geneID, res_dph1_2_WT_FC_1.5_UP_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_1.5_UP_table) #300
res_dph1_1and2_WT_FC_1.5_UP_table = as.data.frame(res_dph1_1and2_WT_FC_1.5_UP_table)
names(res_dph1_1and2_WT_FC_1.5_UP_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_1.5_UP_table = merge(res_dph1_1and2_WT_FC_1.5_UP_table, res_dph1_1_WT_FC_1.5_Ano, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_1.5_UP_table = merge(res_dph1_1and2_WT_FC_1.5_UP_table, res_dph1_2_WT_FC_1.5_UP_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_1.5_UP_table)
res_dph1_1and2_WT_FC_1.5_UP_table = res_dph1_1and2_WT_FC_1.5_UP_table[, c(1:7,12:17,8:11)]


#Genes with LFC<0
res_dph1_1_WT_FC_2_DOWN_table = res_dph1_1_WT_FC_2_Ano[res_dph1_1_WT_FC_2_Ano$log2FoldChange < 0,] #dph1_1.vs.WT_FC_2
res_dph1_2_WT_FC_2_DOWN_table = res_dph1_2_WT_FC_2_Ano[res_dph1_2_WT_FC_2_Ano$log2FoldChange < 0,] #dph1_2.vs.WT_FC_2
res_dph1_1_WT_FC_1.5_DOWN_table = res_dph1_1_WT_FC_1.5_Ano[res_dph1_1_WT_FC_1.5_Ano$log2FoldChange < 0,] #dph1_1.vs.WT_FC_1.5
res_dph1_2_WT_FC_1.5_DOWN_table = res_dph1_2_WT_FC_1.5_Ano[res_dph1_2_WT_FC_1.5_Ano$log2FoldChange < 0,] #dph1_2.vs.WT_FC_1.5
res_dph1_1_dph1_2_FC_1.5_DOWN_table = res_dph1_1_dph1_2_FC_1.5_Ano[res_dph1_1_dph1_2_FC_1.5_Ano$log2FoldChange < 0,] #dph1_1.vs.dph1_2_FC_1.5

#Intersection
#res_dph1_1and2_WT_FC_2
res_dph1_1and2_WT_FC_2_DOWN_table = intersect(res_dph1_1_WT_FC_2_DOWN_table$Ath_geneID, res_dph1_2_WT_FC_2_DOWN_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_2_DOWN_table) #94
res_dph1_1and2_WT_FC_2_DOWN_table = as.data.frame(res_dph1_1and2_WT_FC_2_DOWN_table)
names(res_dph1_1and2_WT_FC_2_DOWN_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_2_DOWN_table = merge(res_dph1_1and2_WT_FC_2_DOWN_table, res_dph1_1_WT_FC_2_Ano, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_2_DOWN_table = merge(res_dph1_1and2_WT_FC_2_DOWN_table, res_dph1_2_WT_FC_2_DOWN_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_2_DOWN_table)
res_dph1_1and2_WT_FC_2_DOWN_table = res_dph1_1and2_WT_FC_2_DOWN_table[, c(1:7,12:17,8:11)]

#res_dph1_1and2_WT_FC_1.5
res_dph1_1and2_WT_FC_1.5_DOWN_table = intersect(res_dph1_1_WT_FC_1.5_DOWN_table$Ath_geneID, res_dph1_2_WT_FC_1.5_DOWN_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_1.5_DOWN_table) #249
res_dph1_1and2_WT_FC_1.5_DOWN_table = as.data.frame(res_dph1_1and2_WT_FC_1.5_DOWN_table)
names(res_dph1_1and2_WT_FC_1.5_DOWN_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_1.5_DOWN_table = merge(res_dph1_1and2_WT_FC_1.5_DOWN_table, res_dph1_1_WT_FC_1.5_Ano, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_1.5_DOWN_table = merge(res_dph1_1and2_WT_FC_1.5_DOWN_table, res_dph1_2_WT_FC_1.5_DOWN_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_1.5_DOWN_table)
res_dph1_1and2_WT_FC_1.5_DOWN_table = res_dph1_1and2_WT_FC_1.5_DOWN_table[, c(1:7,12:17,8:11)]



#Add TPM to tables
res_dph1_1_WT_FC_2_Ano_table = merge(res_dph1_1_WT_FC_2_Ano, UNIVERSAL[,c(1,2,3,5,6,8:13)], by = "Ath_geneID", all.x = TRUE)
res_dph1_2_WT_FC_2_Ano_table = merge(res_dph1_2_WT_FC_2_Ano, UNIVERSAL[,c(1,2,4,5,7,8:10,14:16)], by = "Ath_geneID", all.x = TRUE)
res_dph1_1_WT_FC_1.5_Ano_table = merge(res_dph1_1_WT_FC_1.5_Ano, UNIVERSAL[,c(1,2,3,5,6,8:13)], by = "Ath_geneID", all.x = TRUE)
res_dph1_2_WT_FC_1.5_Ano_table = merge(res_dph1_2_WT_FC_1.5_Ano, UNIVERSAL[,c(1,2,4,5,7,8:10,14:16)], by = "Ath_geneID", all.x = TRUE)
res_dph1_1_dph1_2_FC_1.5_Ano_table = merge(res_dph1_1_dph1_2_FC_1.5_Ano, UNIVERSAL[,c(1,3,4,6,7,11:16)], by = "Ath_geneID", all.x = TRUE)


#Concatenating Intersection
res_dph1_1and2_WT_FC_2_Ano_table = rbind(res_dph1_1and2_WT_FC_2_UP_table, res_dph1_1and2_WT_FC_2_DOWN_table)
res_dph1_1and2_WT_FC_1.5_Ano_table = rbind(res_dph1_1and2_WT_FC_1.5_UP_table, res_dph1_1and2_WT_FC_1.5_DOWN_table)

res_dph1_1and2_WT_FC_2_Ano_table = merge(res_dph1_1and2_WT_FC_2_Ano_table, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)
res_dph1_1and2_WT_FC_1.5_Ano_table = merge(res_dph1_1and2_WT_FC_1.5_Ano_table, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)


# EXTRACTING GENES TABLES
#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/res_xlsx_TABLES')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/res_xlsx_TABLES')

options(java.parameters = "-Xmx6g")  ## memory set to 6 GB

write.xlsx2(res_dph1_1_WT_FC_2_Ano_table, "resTable_1.xlsx", sheetName = "res_dph1_1.vs.WT_FC=2", append=TRUE)
gc()
write.xlsx2(res_dph1_2_WT_FC_2_Ano_table, "resTable_1.xlsx", sheetName = "res_dph2_1.vs.WT_FC=2", append=TRUE)
gc()
write.xlsx2(res_dph1_1_WT_FC_1.5_Ano_table, "resTable_1.xlsx", sheetName = "res_dph1_1.vs.WT_FC=1.5", append=TRUE)
gc()
write.xlsx2(res_dph1_2_WT_FC_1.5_Ano_table, "resTable_1.xlsx", sheetName = "res_dph1_2.vs.WT_FC=1.5", append=TRUE)
gc()
write.xlsx2(res_dph1_1_dph1_2_FC_1.5_Ano_table, "resTable_1.xlsx", sheetName = "res_dph1_1.vs.dph1_2_FC=1.5", append=TRUE)
gc()
write.xlsx2(res_dph1_1and2_WT_FC_2_Ano_table, "resTable_1.xlsx", sheetName = "overlap_dph1_1&2_WT_FC=2", append=TRUE)
gc()
write.xlsx2(res_dph1_1and2_WT_FC_1.5_Ano_table, "resTable_1.xlsx", sheetName = "overlap_dph1_1&2_WT_FC=1.5", append=TRUE)
gc()



# FULL TABLES AFTER THRESHOLD
###############################
#Genes with LFC>0
res_dph1_1_WT_FC_2_UP_table = res_dph1_1_WT_FC_2_Ano_th[res_dph1_1_WT_FC_2_Ano_th$log2FoldChange > 0,] #dph1_1.vs.WT_FC_2
res_dph1_2_WT_FC_2_UP_table = res_dph1_2_WT_FC_2_Ano_th[res_dph1_2_WT_FC_2_Ano_th$log2FoldChange > 0,] #dph1_2.vs.WT_FC_2
res_dph1_1_WT_FC_1.5_UP_table = res_dph1_1_WT_FC_1.5_Ano_th[res_dph1_1_WT_FC_1.5_Ano_th$log2FoldChange > 0,] #dph1_1.vs.WT_FC_1.5
res_dph1_2_WT_FC_1.5_UP_table = res_dph1_2_WT_FC_1.5_Ano_th[res_dph1_2_WT_FC_1.5_Ano_th$log2FoldChange > 0,] #dph1_2.vs.WT_FC_1.5
res_dph1_1_dph1_2_FC_1.5_UP_table = res_dph1_1_dph1_2_FC_1.5_Ano_th[res_dph1_1_dph1_2_FC_1.5_Ano_th$log2FoldChange > 0,] #dph1_1.vs.dph1_2_FC_1.5

#Intersection
#res_dph1_1and2_WT_FC_2
res_dph1_1and2_WT_FC_2_UP_table = intersect(res_dph1_1_WT_FC_2_UP_table$Ath_geneID, res_dph1_2_WT_FC_2_UP_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_2_UP_table) #97
res_dph1_1and2_WT_FC_2_UP_table = as.data.frame(res_dph1_1and2_WT_FC_2_UP_table)
names(res_dph1_1and2_WT_FC_2_UP_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_2_UP_table = merge(res_dph1_1and2_WT_FC_2_UP_table, res_dph1_1_WT_FC_2_Ano_th, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_2_UP_table = merge(res_dph1_1and2_WT_FC_2_UP_table, res_dph1_2_WT_FC_2_UP_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_2_UP_table)
res_dph1_1and2_WT_FC_2_UP_table = res_dph1_1and2_WT_FC_2_UP_table[, c(1:7,12:17,8:11)]

#res_dph1_1and2_WT_FC_1.5
res_dph1_1and2_WT_FC_1.5_UP_table = intersect(res_dph1_1_WT_FC_1.5_UP_table$Ath_geneID, res_dph1_2_WT_FC_1.5_UP_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_1.5_UP_table) #267
res_dph1_1and2_WT_FC_1.5_UP_table = as.data.frame(res_dph1_1and2_WT_FC_1.5_UP_table)
names(res_dph1_1and2_WT_FC_1.5_UP_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_1.5_UP_table = merge(res_dph1_1and2_WT_FC_1.5_UP_table, res_dph1_1_WT_FC_1.5_Ano_th, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_1.5_UP_table = merge(res_dph1_1and2_WT_FC_1.5_UP_table, res_dph1_2_WT_FC_1.5_UP_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_1.5_UP_table)
res_dph1_1and2_WT_FC_1.5_UP_table = res_dph1_1and2_WT_FC_1.5_UP_table[, c(1:7,12:17,8:11)]


#Genes with LFC<0
res_dph1_1_WT_FC_2_DOWN_table = res_dph1_1_WT_FC_2_Ano_th[res_dph1_1_WT_FC_2_Ano_th$log2FoldChange < 0,] #dph1_1.vs.WT_FC_2
res_dph1_2_WT_FC_2_DOWN_table = res_dph1_2_WT_FC_2_Ano_th[res_dph1_2_WT_FC_2_Ano_th$log2FoldChange < 0,] #dph1_2.vs.WT_FC_2
res_dph1_1_WT_FC_1.5_DOWN_table = res_dph1_1_WT_FC_1.5_Ano_th[res_dph1_1_WT_FC_1.5_Ano_th$log2FoldChange < 0,] #dph1_1.vs.WT_FC_1.5
res_dph1_2_WT_FC_1.5_DOWN_table = res_dph1_2_WT_FC_1.5_Ano_th[res_dph1_2_WT_FC_1.5_Ano_th$log2FoldChange < 0,] #dph1_2.vs.WT_FC_1.5
res_dph1_1_dph1_2_FC_1.5_DOWN_table = res_dph1_1_dph1_2_FC_1.5_Ano_th[res_dph1_1_dph1_2_FC_1.5_Ano_th$log2FoldChange < 0,] #dph1_1.vs.dph1_2_FC_1.5

#Intersection
#res_dph1_1and2_WT_FC_2
res_dph1_1and2_WT_FC_2_DOWN_table = intersect(res_dph1_1_WT_FC_2_DOWN_table$Ath_geneID, res_dph1_2_WT_FC_2_DOWN_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_2_DOWN_table) #87
res_dph1_1and2_WT_FC_2_DOWN_table = as.data.frame(res_dph1_1and2_WT_FC_2_DOWN_table)
names(res_dph1_1and2_WT_FC_2_DOWN_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_2_DOWN_table = merge(res_dph1_1and2_WT_FC_2_DOWN_table, res_dph1_1_WT_FC_2_Ano_th, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_2_DOWN_table = merge(res_dph1_1and2_WT_FC_2_DOWN_table, res_dph1_2_WT_FC_2_DOWN_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_2_DOWN_table)
res_dph1_1and2_WT_FC_2_DOWN_table = res_dph1_1and2_WT_FC_2_DOWN_table[, c(1:7,12:17,8:11)]

#res_dph1_1and2_WT_FC_1.5
res_dph1_1and2_WT_FC_1.5_DOWN_table = intersect(res_dph1_1_WT_FC_1.5_DOWN_table$Ath_geneID, res_dph1_2_WT_FC_1.5_DOWN_table$Ath_geneID)
length(res_dph1_1and2_WT_FC_1.5_DOWN_table) #242
res_dph1_1and2_WT_FC_1.5_DOWN_table = as.data.frame(res_dph1_1and2_WT_FC_1.5_DOWN_table)
names(res_dph1_1and2_WT_FC_1.5_DOWN_table)[1]<-"Ath_geneID"
res_dph1_1and2_WT_FC_1.5_DOWN_table = merge(res_dph1_1and2_WT_FC_1.5_DOWN_table, res_dph1_1_WT_FC_1.5_Ano_th, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
res_dph1_1and2_WT_FC_1.5_DOWN_table = merge(res_dph1_1and2_WT_FC_1.5_DOWN_table, res_dph1_2_WT_FC_1.5_DOWN_table, by.x="Ath_geneID", all.x = TRUE,by.y="Ath_geneID")
colnames(res_dph1_1and2_WT_FC_1.5_DOWN_table)
res_dph1_1and2_WT_FC_1.5_DOWN_table = res_dph1_1and2_WT_FC_1.5_DOWN_table[, c(1:7,12:17,8:11)]



#Add TPM to tables
res_dph1_1_WT_FC_2_Ano_th_table = merge(res_dph1_1_WT_FC_2_Ano_th, UNIVERSAL[,c(1,2,3,5,6,8:13)], by = "Ath_geneID", all.x = TRUE)
res_dph1_2_WT_FC_2_Ano_th_table = merge(res_dph1_2_WT_FC_2_Ano_th, UNIVERSAL[,c(1,2,4,5,7,8:10,14:16)], by = "Ath_geneID", all.x = TRUE)
res_dph1_1_WT_FC_1.5_Ano_th_table = merge(res_dph1_1_WT_FC_1.5_Ano_th, UNIVERSAL[,c(1,2,3,5,6,8:13)], by = "Ath_geneID", all.x = TRUE)
res_dph1_2_WT_FC_1.5_Ano_th_table = merge(res_dph1_2_WT_FC_1.5_Ano_th, UNIVERSAL[,c(1,2,4,5,7,8:10,14:16)], by = "Ath_geneID", all.x = TRUE)
res_dph1_1_dph1_2_FC_1.5_Ano_th_table = merge(res_dph1_1_dph1_2_FC_1.5_Ano_th, UNIVERSAL[,c(1,3,4,6,7,11:16)], by = "Ath_geneID", all.x = TRUE)


#Concatenating Intersection
res_dph1_1and2_WT_FC_2_Ano_th_table = rbind(res_dph1_1and2_WT_FC_2_UP_table, res_dph1_1and2_WT_FC_2_DOWN_table)
res_dph1_1and2_WT_FC_1.5_Ano_th_table = rbind(res_dph1_1and2_WT_FC_1.5_UP_table, res_dph1_1and2_WT_FC_1.5_DOWN_table)

res_dph1_1and2_WT_FC_2_Ano_th_table = merge(res_dph1_1and2_WT_FC_2_Ano_th_table, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)
res_dph1_1and2_WT_FC_1.5_Ano_th_table = merge(res_dph1_1and2_WT_FC_1.5_Ano_th_table, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)


# EXTRACTING GENES TABLES
#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/res_xlsx_TABLES')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/res_xlsx_TABLES')

options(java.parameters = "-Xmx6g")  ## memory set to 6 GB

write.xlsx2(res_dph1_1_WT_FC_2_Ano_th_table, "resTable_1_th.xlsx", sheetName = "res_dph1_1.vs.WT_FC=2", append=TRUE)
gc()
write.xlsx2(res_dph1_2_WT_FC_2_Ano_th_table, "resTable_1_th.xlsx", sheetName = "res_dph2_1.vs.WT_FC=2", append=TRUE)
gc()
write.xlsx2(res_dph1_1_WT_FC_1.5_Ano_th_table, "resTable_1_th.xlsx", sheetName = "res_dph1_1.vs.WT_FC=1.5", append=TRUE)
gc()
write.xlsx2(res_dph1_2_WT_FC_1.5_Ano_th_table, "resTable_1_th.xlsx", sheetName = "res_dph1_2.vs.WT_FC=1.5", append=TRUE)
gc()
write.xlsx2(res_dph1_1_dph1_2_FC_1.5_Ano_th_table, "resTable_1_th.xlsx", sheetName = "res_dph1_1.vs.dph1_2_FC=1.5", append=TRUE)
gc()
write.xlsx2(res_dph1_1and2_WT_FC_2_Ano_th_table, "resTable_1_th.xlsx", sheetName = "overlap_dph1_1&2_WT_FC=2", append=TRUE)
gc()
write.xlsx2(res_dph1_1and2_WT_FC_1.5_Ano_th_table, "resTable_1_th.xlsx", sheetName = "overlap_dph1_1&2_WT_FC=1.5", append=TRUE)
gc()




##################################
# DESeq2#  ---  VARIANT 2
##################################
##############
dds2 = DESeq(dds2)
#phenotype = factor(c(rep("control", 3), rep("mutant", 6)))
#dds2 = DESeqDataSetFromMatrix(countData=sCOUNTS, colData=coldata_2, design = ~phenotype)
result_mutant_WT <- results(dds2,contrast=c("phenotype","mutant","control"))


#Choose significantvariants
result_mutant_WT_FC_2  = subset(result_mutant_WT, (padj<=0.05 & log2FoldChange > 1) | (padj<=0.05 & log2FoldChange < -1))
result_mutant_WT_FC_1.5  = subset(result_mutant_WT, (padj<=0.05 & log2FoldChange > 0.5849625) | (padj<=0.05 & log2FoldChange < -0.5849625))


summary(result_mutant_WT_FC_2) #435 genes (264/171)
summary(result_mutant_WT_FC_1.5) #969 genes (549/420)


#setting working directory
#universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results') 
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results')
write.table(result_mutant_WT_FC_2, "result_mutant_WT_FC_2.txt", sep = "\t", row.names = TRUE)
write.table(result_mutant_WT_FC_1.5, "result_mutant_WT_FC_1.5.txt", sep = "\t", row.names = TRUE)



#Adding genes info to the result tables
result_mutant_WT_FC_2_Ano = as.data.frame(result_mutant_WT_FC_2)
result_mutant_WT_FC_1.5_Ano  = as.data.frame(result_mutant_WT_FC_1.5)

result_mutant_WT_FC_2_Ano$Ath_geneID = substr(rownames(result_mutant_WT_FC_2), 1,9)
result_mutant_WT_FC_1.5_Ano$Ath_geneID  = substr(rownames(result_mutant_WT_FC_1.5), 1,9)

result_mutant_WT_FC_2_Ano = merge(result_mutant_WT_FC_2_Ano, Ath_Anno_TAIR10, by ="Ath_geneID", all.x = TRUE)
result_mutant_WT_FC_1.5_Ano  = merge(result_mutant_WT_FC_1.5_Ano, Ath_Anno_TAIR10, by ="Ath_geneID", all.x = TRUE)


#FC > |2.0|
nrow(result_mutant_WT_FC_2_Ano) #435
length(unique(result_mutant_WT_FC_2_Ano$Ath_geneID)) #435
#FC > |1.5|
nrow(result_mutant_WT_FC_1.5_Ano) #969
length(unique(result_mutant_WT_FC_1.5_Ano$Ath_geneID)) #969


#full results after thresholds
result_mutant_WT_FC_2_Ano_th = merge(result_mutant_WT_FC_2_Ano, TPM_th_MM, by="Ath_geneID")
result_mutant_WT_FC_1.5_Ano_th  = merge(result_mutant_WT_FC_1.5_Ano, TPM_th_MM, by="Ath_geneID")
nrow(result_mutant_WT_FC_2_Ano_th) #291 (from 435 genes)
nrow(result_mutant_WT_FC_1.5_Ano_th) #805 (from 969 genes)


#PREPARE FILES FOR MAPMAN INPUT #Extracting crusial rows for MM
result_mutant_WT_FC_2_Ano_th_MM = result_mutant_WT_FC_2_Ano_th[c(1,3)]
result_mutant_WT_FC_1.5_Ano_th_MM  = result_mutant_WT_FC_1.5_Ano_th[c(1,3)]


#PREPARE FILES FOR VENNY2.1 INPUT
#Genes with LFC>0 
result_mutant_WT_FC_2_Ano_VennyUP = result_mutant_WT_FC_2_Ano[result_mutant_WT_FC_2_Ano$log2FoldChange > 0,]
result_mutant_WT_FC_1.5_Ano_VennyUP = result_mutant_WT_FC_1.5_Ano[result_mutant_WT_FC_1.5_Ano$log2FoldChange > 0,]
# after threshhold
result_mutant_WT_FC_2_Ano_th_VennyUP = result_mutant_WT_FC_2_Ano_th[result_mutant_WT_FC_2_Ano_th$log2FoldChange > 0,]
result_mutant_WT_FC_1.5_Ano_th_VennyUP = result_mutant_WT_FC_1.5_Ano_th[result_mutant_WT_FC_1.5_Ano_th$log2FoldChange > 0,]

#Genes with LFC<0 
result_mutant_WT_FC_2_Ano_VennyDOWN = result_mutant_WT_FC_2_Ano[result_mutant_WT_FC_2_Ano$log2FoldChange < 0,]
result_mutant_WT_FC_1.5_Ano_VennyDOWN = result_mutant_WT_FC_1.5_Ano[result_mutant_WT_FC_1.5_Ano$log2FoldChange < 0,]
# after threshhold
result_mutant_WT_FC_2_Ano_th_VennyDOWN = result_mutant_WT_FC_2_Ano_th[result_mutant_WT_FC_2_Ano_th$log2FoldChange < 0,]
result_mutant_WT_FC_1.5_Ano_th_VennyDOWN = result_mutant_WT_FC_1.5_Ano_th[result_mutant_WT_FC_1.5_Ano_th$log2FoldChange < 0,]


#### EXTRACTING RESULTS

#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/DEGs_full/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/DEGs_full/')

#full results
write.table(as.data.frame(result_mutant_WT_FC_2_Ano), file = "result_mutant_WT_FC_2_Ano.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(result_mutant_WT_FC_1.5_Ano), file = "result_mutant_WT_FC_1.5_Ano.txt", sep = "\t", row.names = FALSE)

#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/DEGs_full_th/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/DEGs_full_th/')

write.table(as.data.frame(result_mutant_WT_FC_2_Ano_th), file = "result_mutant_WT_FC_2_Ano_th.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(result_mutant_WT_FC_1.5_Ano_th), file = "result_mutant_WT_FC_1.5_Ano_th.txt", sep = "\t", row.names = FALSE)


#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/forMM_th/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/forMM_th/')

#results for MapMan input(only crusial columns and no NA)

#UP
write.table(as.data.frame(result_mutant_WT_FC_2_Ano_th_MM[result_mutant_WT_FC_2_Ano_th_MM$log2FoldChange > 0,]), file = "result_mutant_WT_FC_2_Ano_th_MM_UP.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(result_mutant_WT_FC_1.5_Ano_th_MM[result_mutant_WT_FC_1.5_Ano_th_MM$log2FoldChange > 0,]), file = "result_mutant_WT_FC_1.5_Ano_th_MM_UP.txt", sep = "\t", row.names = FALSE)
#DOWN
write.table(as.data.frame(result_mutant_WT_FC_2_Ano_th_MM[result_mutant_WT_FC_2_Ano_th_MM$log2FoldChange < 0,]), file = "result_mutant_WT_FC_2_Ano_th_MM_DOWN.txt", sep = "\t", row.names = FALSE)
write.table(as.data.frame(result_mutant_WT_FC_1.5_Ano_th_MM[result_mutant_WT_FC_1.5_Ano_th_MM$log2FoldChange < 0,]), file = "result_mutant_WT_FC_1.5_Ano_th_MM_DOWN.txt", sep = "\t", row.names = FALSE)



##################
## FULL TABLES ##
##################

#UNIVERSAL TABLE
################
mutant.vs.WT_FC_2 = as.data.frame(result_mutant_WT_FC_2)
mutant.vs.WT_FC_1.5 = as.data.frame(result_mutant_WT_FC_1.5)

mutant.vs.WT_FC_2$Ath_geneID = rownames(mutant.vs.WT_FC_2)
mutant.vs.WT_FC_1.5$Ath_geneID = rownames(mutant.vs.WT_FC_1.5)

names(mutant.vs.WT_FC_2) = c("mutant.vs.WT_FC_2_baseMean", "mutant.vs.WT_FC_2_log2FoldChange", "mutant.vs.WT_FC_2_lfcSE", "mutant.vs.WT_FC_2_stat", "mutant.vs.WT_FC_2_pvalue", "mutant.vs.WT_FC_2_padj", "Ath_geneID")
names(mutant.vs.WT_FC_1.5) = c("mutant.vs.WT_FC_1.5_baseMean", "mutant.vs.WT_FC_1.5_log2FoldChange", "mutant.vs.WT_FC_1.5_lfcSE", "mutant.vs.WT_FC_1.5_stat", "mutant.vs.WT_FC_1.5_pvalue", "mutant.vs.WT_FC_1.5_padj", "Ath_geneID")

nrow(UNIVERSAL) #32678
UNIVERSAL = merge(UNIVERSAL, mutant.vs.WT_FC_2, by = "Ath_geneID", all.x = TRUE)
UNIVERSAL = merge(UNIVERSAL, mutant.vs.WT_FC_1.5, by = "Ath_geneID", all.x = TRUE)


# EXTRACTING GENES TABLES
#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results')

options(java.parameters = "-Xmx6g")  ## memory set to 6 GB
write.xlsx2(UNIVERSAL, "UNIVERSAL_2.xlsx", sheetName = "Variant1", append=TRUE)

##########
##########

#Add TPM to tables
result_mutant_WT_FC_2_Ano_table = merge(result_mutant_WT_FC_2_Ano, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)
result_mutant_WT_FC_1.5_Ano_table = merge(result_mutant_WT_FC_1.5_Ano, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)

# EXTRACTING GENES TABLES
#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/res_xlsx_TABLES')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/res_xlsx_TABLES')

options(java.parameters = "-Xmx6g")  ## memory set to 6 GB

write.xlsx2(result_mutant_WT_FC_2_Ano_table, "resTable_2.xlsx", sheetName = "res_mutant.vs.WT_FC=2", append=TRUE)
gc()
write.xlsx2(result_mutant_WT_FC_1.5_Ano_table, "resTable_2.xlsx", sheetName = "res_mutant.vs.WT_FC=1.5", append=TRUE)
gc()


# FULL TABLES AFTER THRESHOLD
###############################
#Add TPM to tables
result_mutant_WT_FC_2_Ano_th_table = merge(result_mutant_WT_FC_2_Ano_th, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)
result_mutant_WT_FC_1.5_Ano_th_table = merge(result_mutant_WT_FC_1.5_Ano_th, UNIVERSAL[,c(1:16)], by = "Ath_geneID", all.x = TRUE)

# EXTRACTING GENES TABLES
#setting working directory _ universitate
setwd('/Users/nadezdajanina/Desktop/DPH1/Results/res_xlsx_TABLES')
#my laptop
setwd('/Users/nadezhdayanina/Desktop/DPH1/Results/res_xlsx_TABLES')

options(java.parameters = "-Xmx6g")  ## memory set to 6 GB

write.xlsx2(result_mutant_WT_FC_2_Ano_th_table, "resTable_2_th.xlsx", sheetName = "res_mutant.vs.WT_FC=2", append=TRUE)
gc()
write.xlsx2(result_mutant_WT_FC_1.5_Ano_th_table, "resTable_2_th.xlsx", sheetName = "res_mutant.vs.WT_FC=1.5", append=TRUE)
gc()


#RETURN TO DEFAULT WORKING FOLDER
setwd('/Users/nadezdajanina/Desktop/DPH1') #setting working directory _ universitate
setwd('/Users/nadezhdayanina/Desktop/DPH1') #my laptop






