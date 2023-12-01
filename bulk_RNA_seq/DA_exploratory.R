library(ggplot2)
# load the data
#load the phenotype table
pdata=read.table("phenotype_dopamine.txt",sep = "",row.names = 1)
# load the expression data
edata=read.table("expression_dopamine.txt",sep = "",row.names = 1)

# table for factor/character variables
pdata$group=as.factor(pdata$group)

table(pdata$group)


# remove low expression data
edata = edata[rowMeans(edata) > 10, ]


# look at overall distributions
boxplot(edata[,1])

boxplot(log2(edata[,1]+1))

boxplot(log2(edata+1),col=2,range=0)

# transform the edata
edata=log2(edata+1)


# get the exact principal components use prcomp
pc=prcomp(t(edata))

pcaData<-as.data.frame(pc$x)
pcaData$group<-pdata$group
percentVar<-round(100*(pc$sdev^2/sum(pc$sdev^2)))
p<-ggplot(data=pcaData,aes(x=PC1,y=PC2,color=group))+
  geom_point(size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()
p



# clustering
# sample distance plot
library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(edata))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(edata))
colnames(sampleDistMatrix) <- paste(colnames(edata))
colors2 <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sample_distance_plot<-
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors2)


# tree clustering
hclust1=hclust(sampleDists)
plot(hclust1)
plot(hclust1,hang=0.5)


