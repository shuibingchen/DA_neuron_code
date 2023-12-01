# load the counts of each sample
DMSO<-read.table(file = "DMSO_gene_count.tsv",header = TRUE,row.names = 1)

Metformin<-read.table(file = "Metformin_gene_count.tsv",header = TRUE,row.names = 1)

DMSO_Metformin<-cbind(DMSO,Metformin)

DMSO_Metformin<-DMSO_Metformin[-(60672:60682),]

# rename the colname to the group name
group<-rep(c("DMSO","mf"),each=3)

group<-as.factor(group)

colnames(DMSO_Metformin)<-paste(group,1:3,sep = "-")

# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(DMSO_Metformin)

y<-DGEList(counts=DMSO_Metformin,genes=genelist)

# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]


# Normalization
y <- calcNormFactors(y, method="TMM")

y$samples$group <- group

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# exploring differences between libraries
# MDS plot
pch<-c(15,16)
colors<-rep(c("red","green"),3)

plotMDS(y,col=colors[group],pch=pch[group])

legend("topright",legend=levels(group),pch=pch,col=colors)


# samples distance plot
library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(DMSO_Metformin))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(DMSO_Metformin))
colnames(sampleDistMatrix) <- paste(colnames(DMSO_Metformin))
colors2 <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sample_distance_plot<-
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors2)

# PCA plot
library(ggfortify)
pcDat <- prcomp(t(y$counts))

groupCol <- as.numeric(factor(y$samples$group)) + 1

PCA_plot<-autoplot(pcDat,
                   data = y$samples, 
                   colour="group", 
                   size=3)+ theme_bw()

PCA_plot


tiff(filename = "DMSO_Metformin_sample_distance_plot.tif",width = 12,height = 8,units ="cm",compression="lzw",bg="white",res=1024)
sample_distance_plot
dev.off()

tiff(filename = "DMSO_Metformin_PCA_plot.tif",width = 12,height = 8,units ="cm",compression="lzw",bg="white",res=1024)
PCA_plot
dev.off()
