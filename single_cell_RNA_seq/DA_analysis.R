library(dplyr)
library(Seurat)
library(patchwork)

# load the mice scRNAseq data
DA_mock1.data <- Read10X(data.dir = "mock_1/filtered_feature_bc_matrix/")

DA_mock1 <- CreateSeuratObject(counts = DA_mock1.data, project = "DA_mock1", min.cells = 3, min.features = 100)

DA_mock1


# load the mice scRNAseq data
DA_infect.data <- Read10X(data.dir = "infect/filtered_feature_bc_matrix/")

DA_infect <- CreateSeuratObject(counts = DA_infect.data, project = "DA_infect", min.cells = 3, min.features = 100)

DA_infect




# merge two dataset
DA <- merge(DA_mock1, y = DA_infect, add.cell.ids = c("mock","infect"), project = "DA")


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
DA[["percent.mt"]] <- PercentageFeatureSet(DA, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(DA@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(DA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(DA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(DA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

# filter the cells based on nFeature_RNA and percent.mt
DA <- subset(DA, subset = nCount_RNA < 20000 & nFeature_RNA < 6000 & percent.mt < 10)


# normalization of the data
DA <- NormalizeData(DA, normalization.method = "LogNormalize", scale.factor = 10000)

DA <- FindVariableFeatures(DA, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(DA), 10)

plot3 <- VariableFeaturePlot(DA)
plot3
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot4


#Scaling the data including all the genes
all.genes <- rownames(DA)
DA <- ScaleData(DA, features = all.genes)

#'regress out' heterogeneity associated with mitochondrial contamination
DA <- ScaleData(DA, vars.to.regress = "percent.mt")


#Perform linear dimensional reduction
selected_features<-setdiff(VariableFeatures(object = DA),rownames(DA_infect.data)[36602:36613])
DA <- RunPCA(DA, features = selected_features)


# Examine and visualize PCA results a few different ways
print(DA[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(DA, dims = 1:2, reduction = "pca")

DimPlot(DA, reduction = "pca")

DimPlot(DA, reduction = "pca",split.by  = "orig.ident")


DimHeatmap(DA, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(DA, dims = 1:15, cells = 500, balanced = TRUE)


# Determine the 'dimensionality' of the dataset
DA <- JackStraw(DA, num.replicate = 100)

DA <- ScoreJackStraw(DA, dims = 1:20)

JackStrawPlot(DA, dims = 1:20,ymax = 0.8)

ElbowPlot(DA)


#Cluster the cells
DA <- FindNeighbors(DA, dims = 1:10)
DA <- FindClusters(DA, resolution = 0.2)
# resolution can adjust from 0.2 to 1.5

head(Idents(DA), 5)


#Run non-linear dimensional reduction (UMAP/tSNE)
DA <- RunUMAP(DA, dims = 1:10)

DimPlot(DA, reduction = "umap",label = TRUE)

DimPlot(DA, reduction = "umap",group.by = "orig.ident")
DimPlot(DA, reduction = "umap",split.by  = "orig.ident")

DA <- RunTSNE(DA, dims = 1:10)

DimPlot(DA, reduction = "tsne")
DimPlot(DA, reduction = "tsne",split.by  = "orig.ident")

# rename the idents of each cluster
DA<-RenameIdents(DA,'0'='0','1'='1','2'='1','3'='1','4'='2','5'='3')


DA$orig.ident<-factor(x=DA$orig.ident,levels = c("DA_mock1","DA_infect"))
# find all markers of cluster 0
cluster0.markers <- FindMarkers(DA, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 20)

# find the ratio of NR4A2 and TH positive cell in each group
DA_infect<-subset(DA,subset = orig.ident=="DA_infect")
DA_mock<-subset(DA,subset = orig.ident=="DA_mock1")

NR4A2_mock_positive_cells<-sum(GetAssayData(object = DA_mock, slot = "data")["NR4A2",]>0)
NR4A2_mock_negative_cells<-sum(GetAssayData(object = DA_mock, slot = "data")["NR4A2",]==0)

TH_mock_positive_cells<-sum(GetAssayData(object = DA_mock, slot = "data")["TH",]>0)
TH_mock_negative_cells<-sum(GetAssayData(object = DA_mock, slot = "data")["TH",]==0)

mock_NR4A2<-NR4A2_mock_positive_cells/(NR4A2_mock_positive_cells+NR4A2_mock_negative_cells)

mock_TH<-TH_mock_positive_cells/(TH_mock_positive_cells+TH_mock_negative_cells)


NR4A2_infect_positive_cells<-sum(GetAssayData(object = DA_infect, slot = "data")["NR4A2",]>0)
NR4A2_infect_negative_cells<-sum(GetAssayData(object = DA_infect, slot = "data")["NR4A2",]==0)

TH_infect_positive_cells<-sum(GetAssayData(object = DA_infect, slot = "data")["TH",]>0)
TH_infect_negative_cells<-sum(GetAssayData(object = DA_infect, slot = "data")["TH",]==0)

infect_NR4A2<-NR4A2_infect_positive_cells/(NR4A2_infect_positive_cells+NR4A2_infect_negative_cells)

infect_TH<-TH_infect_positive_cells/(TH_infect_positive_cells+TH_infect_negative_cells)

DA_ratio<-data.frame(NR4A2=c(mock_NR4A2,infect_NR4A2),
                     TH=c(mock_TH,infect_TH))

rownames(DA_ratio)<-c("mock","infect")

write.csv(DA_ratio,file = "DA_mock_infect_ratio.csv")

#  find the ratio of cells express SARSCOV2 gene in NR4A2 and TH positive cells
NR4A2_positive<-subset(DA_infect,subset = NR4A2>0)

NR4A2_negative<-subset(DA_infect,subset = NR4A2==0)

NR4A2_LMO3_positive<-subset(DA_infect,subset = NR4A2>0&LMO3>0)

NR4A2_CALB1_positive<-subset(DA_infect,subset = NR4A2>0&CALB1>0)

NR4A2_positive_number<-sum(GetAssayData(object = DA_infect_NR4A2_positive, slot = "data")["CoV2-N",]>=0)

NR4A2_negative_number<-sum(GetAssayData(object = DA_infect_NR4A2_negative, slot = "data")["CoV2-N",]>=0)

NR4A2_LMO3_number<-sum(GetAssayData(object = NR4A2_LMO3_positive, slot = "data")["CoV2-N",]>=0)

NR4A2_CALB1_number<-sum(GetAssayData(object = NR4A2_CALB1_positive, slot = "data")["CoV2-N",]>=0)

NR4A2_positive_N<-sum(GetAssayData(object = DA_infect_NR4A2_positive, slot = "data")["CoV2-N",]>0)
NR4A2_negetive_N<-sum(GetAssayData(object = DA_infect_NR4A2_negative, slot = "data")["CoV2-N",]>0)

NR4A2_LMO3_N<-sum(GetAssayData(object = NR4A2_LMO3_positive, slot = "data")["CoV2-N",]>0)

NR4A2_CALB1_N<-sum(GetAssayData(object = NR4A2_CALB1_positive, slot = "data")["CoV2-N",]>0)

cell_number<-data.frame(cell_number=c(NR4A2_positive_number,
                                 NR4A2_negative_number,
                                 NR4A2_LMO3_number,
                                 NR4A2_CALB1_number
                                 ),
                        CoV2_N_cell_number=c(NR4A2_positive_N,
                                             NR4A2_negetive_N,
                                             NR4A2_LMO3_N,
                                             NR4A2_CALB1_N))

rownames(cell_number)<-c("NR4A2+","NR4A2-","NR4A2+LMO3+","NR4A2+CALB1+")

write.csv(cell_number,file = "cell_number_in_subpopulation.csv")

# find markers for every cluster compared to all remaining cells, report only the positive ones
DA.markers <- FindAllMarkers(DA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DA.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top15 <- DA.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(DA, features = top10$gene) + NoLegend()

# find the cell number of each cluster
table(Idents(DA))


# features of A9/A10/senescence 
features_A9<-c("TUJ1","TH","LMO3","LMX1A","FOXA2","NR4A2","KCNJ6","MAP2")
features_A10<-c("TUJ1","TH","CALB1","NR4A2","LMX1A","FOXA2","MAP2")
features_DA<-c("LMO3","LMX1A","FOXA2","NR4A2","KCNJ6","MAP2")
features_senescence<-c("CDKN1A","LMNB1","IGFBP7","SNCA")
features_COV<-c("CoV2-M","CoV2-S","CoV2-N","CoV2-E","CoV2_orf1ab","CoV2_ORF3a","CoV2_ORF6","CoV2_ORF7a","CoV2_ORF7b","CoV2_ORF8","CoV2_ORF10")


RidgePlot(DA, features = features_A9, ncol = 2)

DotPlot(DA, features = features_A9, split.by = "orig.ident") + RotatedAxis()

DotPlot(DA, features = features_A9) + RotatedAxis()

DotPlot(DA, features = features_DA, group.by = "orig.ident") + RotatedAxis()

DotPlot(DA, features = features_A10, group.by = "orig.ident") + RotatedAxis()

DotPlot(DA, features = features_COV, group.by = "orig.ident") + RotatedAxis()

DotPlot(DA, features = features_senescence, group.by = "orig.ident") + RotatedAxis()

#find each cluster characters
VlnPlot(DA, features = c("TH"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("TH"),split.by = "orig.ident")


VlnPlot(DA, features = c("NR4A2"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("NR4A2"),split.by = "orig.ident")

VlnPlot(DA, features = c("NR4A2"))
FeaturePlot(DA,features = c("NR4A2"))

VlnPlot(DA, features = c("SLC18A2"))
FeaturePlot(DA,features = c("SLC18A2"))

VlnPlot(DA, features = c("SLC6A3"))
FeaturePlot(DA,features = c("SLC6A3"))

VlnPlot(DA, features = c("LMO3"))
FeaturePlot(DA,features = c("LMO3"))
VlnPlot(DA, features = c("LMO3"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("LMO3"),split.by = "orig.ident")

VlnPlot(DA, features = c("CALB1"))
FeaturePlot(DA,features = c("CALB1"))
VlnPlot(DA, features = c("CALB1"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("CALB1"),split.by = "orig.ident")

VlnPlot(DA, features = c("SNCA"))
FeaturePlot(DA,features = c("SNCA"))
VlnPlot(DA, features = c("SNCA"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("SNCA"),split.by = "orig.ident")

VlnPlot(DA, features = c("SOX6"))
FeaturePlot(DA,features = c("SOX6"))
VlnPlot(DA, features = c("SOX6"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("SOX6"),split.by = "orig.ident")



VlnPlot(DA, features = c("IGFBP7"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("IGFBP7"),split.by = "orig.ident")

VlnPlot(DA, features = c("LMNB1"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("LMNB1"),split.by = "orig.ident")

VlnPlot(DA, features = c("CDKN1A"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("CDKN1A"),split.by = "orig.ident")

VlnPlot(DA, features = c("CDKN2A"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("CDKN2A"),split.by = "orig.ident")


VlnPlot(DA, features = c("CoV2-E"))
FeaturePlot(DA,features = c("CoV2-E"))

VlnPlot(DA, features = c("CoV2-S"))
FeaturePlot(DA,features = c("CoV2-S"))

VlnPlot(DA, features = c("CoV2-N"))
FeaturePlot(DA,features = c("CoV2-N"))

VlnPlot(DA, features = c("CoV2-M"))
FeaturePlot(DA,features = c("CoV2-M"))

VlnPlot(DA, features = c("CoV2-M"),split.by  = "orig.ident")
FeaturePlot(DA,features = c("CoV2-M"),split.by = "orig.ident")

saveRDS(DA, file = "./DA_infect_mock1.rds")

gene_list <- read.table("gene_list.txt",sep = "",header = TRUE)
genes <- gene_list$Genes
for (gene in genes) {
  
  tiff(filename = paste0(gene,"_violin",".tiff"),width = 480,height = 480,units = "px",bg="white")
  
  p<-VlnPlot(DA, features = gene)
  
  print(p)

  dev.off()
}


for (gene in genes) {
  tiff(filename = paste0(gene,"_UMAP",".tiff"),width = 480,height = 480,units = "px",bg="white")
  p<-FeaturePlot(DA, features = gene)
  print(p)
  dev.off()
}


# select only the cells express TH or NR4A2
DA_TH_NR4A2<-subset(DA,subset = TH > 0 | NR4A2 > 0)

DimPlot(DA_TH_NR4A2, reduction = "umap",label = TRUE)

DimPlot(DA_TH_NR4A2, reduction = "umap",group.by  = "orig.ident")

DimPlot(DA_TH_NR4A2, reduction = "umap",split.by  = "orig.ident")


FeaturePlot(DA_TH_NR4A2,features = c("MAP2","NR4A2","FOXA2","LMX1A","LMO3","CALB1"),ncol = 3) + RotatedAxis()

DotPlot(DA_TH_NR4A2,features = c("CoV2-M","CoV2-S","CoV2-N","CoV2-E"),group.by = "orig.ident") + RotatedAxis()

DotPlot(DA_TH_NR4A2,features = c("CoV2-M","CoV2-S","CoV2-N","CoV2-E"),group.by = "Condition") + RotatedAxis()

DotPlot(DA_TH_NR4A2,features = c("LMO3","LMX1A","FOXA2","NR4A2","KCNJ6","MAP2"),group.by = "orig.ident") + RotatedAxis()

VlnPlot(DA_TH_NR4A2,features = c("CoV2-M","CoV2-S","CoV2-N","CoV2-E"),split.by = "orig.ident",ncol = 4) + RotatedAxis()

VlnPlot(DA_TH_NR4A2,features = c("CoV2-N"),group.by = "orig.ident")

DotPlot(DA_TH_NR4A2,features = c("CDKN1A","LMNB1","IGFBP7"),group.by = "orig.ident") + RotatedAxis()

VlnPlot(DA_TH_NR4A2,features = c("CDKN1A","LMNB1","IGFBP7"),group.by = "orig.ident") + RotatedAxis()

DotPlot(DA_TH_NR4A2,features = c("TP53","BAX","BCL2"),group.by = "orig.ident") + RotatedAxis()

VlnPlot(DA_TH_NR4A2,features = c("TP53","BAX","BCL2"),group.by = "orig.ident") + RotatedAxis()

gene_list <- read.csv("gene_list2.csv")
genes <- gene_list$Genes

VlnPlot(DA_TH_NR4A2,features = genes,ncol = 4) + RotatedAxis()

FeaturePlot(DA_TH_NR4A2,features = genes,ncol = 4,pt.size = 0.001) + RotatedAxis()

VlnPlot(DA_TH_NR4A2,features = c("CALB1"),split.by = "orig.ident" )

VlnPlot(DA_TH_NR4A2,features = c("CALB1"),group.by  = "orig.ident" )

VlnPlot(DA_TH_NR4A2,features = c("LMO3"),split.by = "orig.ident" )

VlnPlot(DA_TH_NR4A2,features = c("LMO3"),group.by = "orig.ident" )


VlnPlot(DA_TH_NR4A2,features = c("CALB1"),group.by = "orig.ident") + RotatedAxis()
VlnPlot(DA_TH_NR4A2,features = c("CALB1"),split.by = "orig.ident") + RotatedAxis()

VlnPlot(DA_TH_NR4A2,features = c("LMO3"),group.by = "orig.ident") + RotatedAxis()
VlnPlot(DA_TH_NR4A2,features = c("LMO3"),split.by = "orig.ident") + RotatedAxis()

saveRDS(DA_TH_NR4A2, file = "./DA_TH_NR4A2.rds")
