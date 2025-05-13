#Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

## anterior P7 snRNA-seq data loading ##
MC1.counts <- Read10X((data.dir = "data/3886/filtered_feature_bc_matrix/"))
MC1 <- CreateSeuratObject(counts = MC1.counts, project = "MC1")
MC1@meta.data$Sample <- "MAC1"
MC1@meta.data$Sex <- "Male"

MC2.counts <- Read10X((data.dir = "data/3887/filtered_feature_bc_matrix/"))
MC2 <- CreateSeuratObject(counts = MC2.counts, project = "MC2")
MC2@meta.data$Sample <- "MAC2"
MC2@meta.data$Sex <- "Male"

FC1.counts <- Read10X((data.dir = "data/4207/filtered_feature_bc_matrix/"))
FC1 <- CreateSeuratObject(counts = FC1.counts, project = "FC1")
FC1@meta.data$Sample <- "FAC1"
FC1@meta.data$Sex <- "Female"

# merge data
anterior <- merge(MC1, y = c(MC2, FC1), add.cell.ids = c("MAC1", "MAC2","FAC1"))

# Get mitochondrial reads in data
anterior[["percent.mt"]] <- PercentageFeatureSet(anterior, pattern = "^Mt-")

# Filtering based on nCount and nFeature
anterior.filter <- subset(anterior, subset = nFeature_RNA >= 400 & nFeature_RNA < quantile(nFeature_RNA, 0.98) & nCount_RNA < quantile(nCount_RNA, 0.98) & percent.mt < 15)

# normalization and clustering/visualization
anterior.filter <- NormalizeData(anterior.filter)
anterior.filter <- FindVariableFeatures(anterior.filter, selection.method = "vst", nfeatures = 2000)
anterior.filter <- ScaleData(anterior.filter)
anterior.filter <- RunPCA(anterior.filter)
anterior.filter <- FindNeighbors(anterior.filter, dims = 1:30, reduction = "pca")
anterior.filter <- FindClusters(anterior.filter, resolution = 0.3, cluster.name = "before_doublet")
anterior.filter <- RunUMAP(anterior.filter, dims = 1:30, reduction = "pca", reduction.name = "umap")
DimPlot(anterior.filter, reduction = "umap", label = T)
VlnPlot(anterior.filter, features = c("nFeature_RNA", "percent.mt"), pt.size = 0)

# remove cluster 0 and 3 due to low feauture count #######
anterior.filter <- subset(anterior.filter, idents = setdiff(levels(Idents(anterior.filter)), c("0","3")))

# repeat normalization and clustering analysis
anterior.filter <- NormalizeData(anterior.filter)
anterior.filter <- FindVariableFeatures(anterior.filter, selection.method = "vst", nfeatures = 2000)
anterior.filter <- ScaleData(anterior.filter)
anterior.filter <- RunPCA(anterior.filter)
anterior.filter <- FindNeighbors(anterior.filter, dims = 1:30, reduction = "pca")
anterior.filter <- FindClusters(anterior.filter, resolution = 0.4, cluster.name = "after_removal")
anterior.filter <- RunUMAP(anterior.filter, dims = 1:30, reduction = "pca", reduction.name = "umap")
DimPlot(anterior.filter, reduction = "umap", label = T, split.by = "Sample")
DimPlot(anterior.filter, reduction = "umap", label = T)
anterior <- JoinLayers(anterior.filter)

# cell type annotation
anterior.filter.label <- RenameIdents(anterior.filter, `4` ="D1MSN",`1` ="D2MSN",`17` ="CIN",
                                      `7` ="CART_N",`11` ="Sst_inNa",`18` ="Sst_inNa",`6` ="Sst_inNb",`19` ="inNc",`10` ="inNd",`13` ="NPC",`0` ="exLayer2/3", 
                                      `8` ="exLayer4",`5` ="exLayer5/6",`15` ="PN",`2` ="Astrocyte",`12` ="Astrocyte",`16` ="OL",`3` ="OPC",`9` ="Microglia",
                                      `20` ="PSC",`14` ="EC")

anterior.filter.label <- readRDS("out/manuscript/12325/anterior.adjacent.control.anno.new.rds")

## fig. 2a dimplot
DimPlot(anterior.filter.label, reduction = "umap", label = T)

# fig. 2c dotplot 
markers.to.plot <- c("Gad2","Satb2","Drd1","Drd2","Cartpt","Isl1","Tll1","Pde5a","Hap1","Sst","Crabp1","Tshz1","Top2a","Cux2","Nr4a3","Tle4","Ndst4","Slc1a2","Mbp","Cspg4","C1qb",
                     "Kdr","Col1a1")
levels <- c("D1MSN","D2MSN", "CART_N","CIN","Sst_inNa","Sst_inNb","inNc", "inNd", "NPC","exLayer2/3","exLayer4","exLayer5/6","PN","Astrocyte","OL","OPC","Microglia","EC","PSC"
)
anterior.filter.label@active.ident <- factor(anterior.filter.label@active.ident, levels = levels)
DotPlot(anterior.filter.label, features = markers.to.plot,cols = c("lightgrey", "red"))

###### Posterior #########
# posterior P7 snRNA-seq data loading
PMC1.counts <- Read10X((data.dir = "data/4203/filtered_feature_bc_matrix/"))
PMC1 <- CreateSeuratObject(counts = PMC1.counts, project = "PMC1")
PMC1@meta.data$Sample <- "MPC1"
PMC1@meta.data$Group <- "Control"

PMC2.counts <- Read10X((data.dir = "data/4204/filtered_feature_bc_matrix/"))
PMC2 <- CreateSeuratObject(counts = PMC2.counts, project = "PMC2")
PMC2@meta.data$Sample <- "MPC2"
PMC2@meta.data$Group <- "Control"

PFC.counts <- Read10X((data.dir = "data/4209/filtered_feature_bc_matrix/"))
PFC <- CreateSeuratObject(counts = PFC.counts, project = "PFC")
PFC@meta.data$Sample <- "PFC"
PFC@meta.data$Group <- "Control"

# merge data
posterior <- merge(PMC1, y = c(PMC2, PFC), add.cell.ids = c("PMC1", "PMC2","PFC"))

# Get mitochondrial reads in data
posterior[["percent.mt"]] <- PercentageFeatureSet(posterior, pattern = "^Mt-")

# Filtering based on nCount and nFeature
posterior.filter <- subset(posterior, subset = nFeature_RNA >= 400 & nFeature_RNA < quantile(nFeature_RNA, 0.98) & nCount_RNA < quantile(nCount_RNA, 0.98) & percent.mt < 15)

# normalization and clustering/visualization
posterior.filter <- NormalizeData(posterior.filter)
posterior.filter <- FindVariableFeatures(posterior.filter, selection.method = "vst", nfeatures = 2000)
posterior.filter <- ScaleData(posterior.filter)
posterior.filter <- RunPCA(posterior.filter)
posterior.filter <- FindNeighbors(posterior.filter, dims = 1:30, reduction = "pca")
posterior.filter <- FindClusters(posterior.filter, resolution = 0.3, cluster.name = "before_doublet")
posterior.filter <- RunUMAP(posterior.filter, dims = 1:30, reduction = "pca", reduction.name = "umap")
DimPlot(posterior.filter, reduction = "umap", label = T)
VlnPlot(posterior.filter, features = c("nFeature_RNA", "percent.mt"), pt.size = 0)

#### remove cluster 0 due to low feature count #######
posterior.filter <- subset(posterior.filter, idents = setdiff(levels(Idents(posterior.filter)), c("0","7")))

# repeat normalization and clustering analysis
posterior <- NormalizeData(posterior)
posterior <- FindVariableFeatures(posterior, selection.method = "vst", nfeatures = 2000)
posterior <- ScaleData(posterior)
posterior <- RunPCA(posterior)
posterior <- FindNeighbors(posterior, dims = 1:30, reduction = "pca")
posterior <- FindClusters(posterior, resolution = 0.4, cluster.name = "after_removal")
posterior <- RunUMAP(posterior, dims = 1:30, reduction = "pca", reduction.name = "umap")
DimPlot(posterior, reduction = "umap", label = T, split.by = "Sample")
DimPlot(posterior, reduction = "umap", label = T)

# annotation
posterior.label <- RenameIdents(posterior, `1` ="MSN",`19` ="MSN",`5` ="Sst_inNa",`2` ="Sst_inNb", 
                                `18` ="CIN",`23` ="inNc",
                                `10` ="inNd",`9` ="RN",`11` ="HPNa",`16` ="HPNb",
                                `17` ="NPC", `6` ="exLayer4",`4` ="exLayer2/3",`21` ="exLayer4",`7` ="exLayer5/6",`14` ="PN",
                                `22` ="Astrocyte",`0` ="Astrocyte",`15` ="OL",`3` ="OPC",
                                `13` ="Microglia",`12` ="EC",`20` ="PSC",`24` ="CRs",`8` ="immatureN")

saveRDS(posterior.label, file="out/manuscript/12325/posterior.adjacent.control.anno.rds")
posterior.label <- readRDS("out/manuscript/12325/posterior.adjacent.control.anno.rds")

## fig. 2d dotplot 
markers.to.plot <- c("Gad2","Satb2","Rarb","Isl1","Pde5a","Hap1","Sst","Crabp1","Tshz1","Lhx9","Sema3c","Neurod6","Top2a","Cux2","Nr4a3","Tle4","Ndst4","Slc1a2","Mbp","Cspg4","C1qb",
                     "Kdr","Col1a1","Tp73") #"Nxph1","Ptchd4","Cdh18",
levels <- c("MSN","CIN","18","Sst_inNa","Sst_inNb","inNc","inNd", "RN","HPNa","HPNb","NPC","6","exLayer2/3","exLayer4","exLayer5/6","PN","Astrocyte","OL","OPC","Microglia","EC","PSC","CRs","immatureN")
posterior.label@active.ident <- factor(posterior.label@active.ident, levels = levels)
DotPlot(posterior.label, features = markers.to.plot,cols = c("lightgrey", "red"))

## fig. 2b dimplot
DimPlot(posterior.label, reduction = "umap", label = T)


