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

# the DEA result for all the genes
# dea <- lrt$table

y <- estimateDisp(y, design, robust = TRUE)

fit<-glmQLFit(y,design,robust = TRUE)

mf_vs_dm<-makeContrasts(mf-DMSO,levels = design)

res<-glmQLFTest(fit,contrast = mf_vs_dm)

toptag <- topTags(res, n = nrow(y$genes), p.value = 1)

dea <- toptag$table 

dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

write.table(dea,file = "Metformin_dmso_DEA.tsv",quote=FALSE,sep = "\t")

