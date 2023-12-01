# load the counts of each sample
DMSO<-read.table(file = "DMSO_gene_count.tsv",header = TRUE,row.names = 1)

Riluzole<-read.table(file = "Riluzole_gene_count.tsv",header = TRUE,row.names = 1)

DMSO_Riluzole<-cbind(DMSO,Riluzole)

DMSO_Riluzole<-DMSO_Riluzole[-(60672:60682),]

# rename the colname to the group name
group<-rep(c("DMSO","rz"),each=3)

group<-as.factor(group)

colnames(DMSO_Riluzole)<-paste(group,1:3,sep = "-")

# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(DMSO_Riluzole)

y<-DGEList(counts=DMSO_Riluzole,genes=genelist)

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

rz_vs_dm<-makeContrasts(rz-DMSO,levels = design)

res<-glmQLFTest(fit,contrast = rz_vs_dm)

toptag <- topTags(res, n = nrow(y$genes), p.value = 1)

dea <- toptag$table 

dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

write.table(dea,file = "Riluzole_dmso_DEA.tsv",quote=FALSE,sep = "\t")
