sessionInfo()
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

## 3 months snRNA-seq data loading ##
MC.counts <- Read10X((data.dir = "data/GEX/MC/filtered_feature_bc_matrix/"))
MC <- CreateSeuratObject(counts = MC.counts, project = "3month")
MC@meta.data$Sample <- "MC"
MC@meta.data$Group <- "Control"
MC@meta.data$age <- "3month"

ME.counts <- Read10X((data.dir = "data/GEX/ME/filtered_feature_bc_matrix/"))
ME <- CreateSeuratObject(counts = ME.counts, project = "3month")
ME@meta.data$Sample <- "ME"
ME@meta.data$Group <- "Ecig"
ME@meta.data$age <- "3month"

#load the dataset#
FC.counts <- Read10X((data.dir = "data/GEX/FC/filtered_feature_bc_matrix/"))
FC <- CreateSeuratObject(counts = FC.counts, project = "3month")
FC@meta.data$Sample <- "FC"
FC@meta.data$Group <- "Control"
FC@meta.data$age <- "3month"

FE.counts <- Read10X((data.dir = "data/GEX/FE/filtered_feature_bc_matrix/"))
FE <- CreateSeuratObject(counts = FE.counts, project = "3month")
FE@meta.data$Sample <- "FE"
FE@meta.data$Group <- "Ecig"
FE@meta.data$age <- "3month"

# merge data
combined <- merge(MC, y = c(ME, FC, FE), add.cell.ids = c("MC", "ME","FC", "FE"))
head(combined@meta.data)

# Get mitochondrial reads in data
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^Mt-")

# Filtering based on nCount and nFeature
combined <- subset(combined, subset = nFeature_RNA >= 300 & nFeature_RNA < quantile(nFeature_RNA, 0.98) & nCount_RNA < quantile(nCount_RNA, 0.98) & percent.mt < 15)

# normalization and clustering/visualization
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca")
combined <- FindClusters(combined, resolution = 0.3, cluster.name = "before_integrate")
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", reduction.name = "umap")
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
DimPlot(combined, reduction = "umap", label = T)

#### remove cluster 0 and 2 due to low feature count #######
combined <- subset(combined, idents = setdiff(levels(Idents(combined)), c("0","1")))

# repeat normalization and clustering analysis
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca")
combined <- FindClusters(combined, resolution = 0.2, cluster.name = "before_integrate")
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", reduction.name = "umap")
DimPlot(combined, reduction = "umap", label = T)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

### annotation ###
combined.label <- RenameIdents(combined, `14` ="PSC",`6` ="EC",`0` ="OL",`10` ="Microglia",`16` ="Microglia",
                               `8` ="OPC",`2` ="Astrocyte",`15` ="Astrocyte",`18` ="Astrocyte",
                               `1` ="exLayer2/3",`13` ="exLayer4",`7` ="exLayer5/6",
                               `11` ="PNa",`4` ="PNb",`19` ="PNb",`5` ="immatureN",
                               `3` ="Sst_ina",`9` ="Sst_inb",
                               `17` ="inCGE",`12` ="inLGE")

saveRDS(combined.label, file="data/combine.annotated.12.2.24.rds")
combined.label <- readRDS("data/combine.annotated.12.2.24.rds")
summary(combined.label$nCount_RNA)
# fig.6d dimplot
DimPlot(combined.label, reduction = "umap", label = T)

# fig.6e dotplot
markers.to.plot <- c("Slc17a7","Gad2","Satb2","Col1a1","Kdr","Mog","C1qb","Cspg4",
                     "Aqp4","Rarb","Adarb2","Sst","Pde5a","Tle4","Fezf2","Cux2",
                     "Ntf3","Ndst4")
levels <- c("PSC","EC","OL","Microglia","OPC","Astrocyte","inLGE","inCGE",
            "Sst_inb","Sst_ina","exLayer5/6","exLayer4","exLayer2/3","5","PNb","PNa","immatureN")

combined.label@active.ident <- factor(combined.label@active.ident, levels = levels)
DotPlot(combined.label, features = markers.to.plot,cols = c("lightgrey", "red")) + RotatedAxis()

### DEG ######
Idents(combined.label) <- "Group"
combined.label$cluster.stim <- paste(Idents(combined.label), combined.label$cell.type, sep = "_")
Idents(combined.label) <- "cluster.stim"

########
## DEG between 3 month male e-cig and control######
Idents(combined.label) <- "Group"
library(dplyr)
combined.label$group.sample.stim <- paste(Idents(combined.label), combined.label$Sample, combined.label$cell.type, sep = "_")
Idents(combined.label) <- "group.sample.stim"

cell <- read.csv("data/celltype.csv", header = 1)
cell <- as.data.frame(cell)
cell <- split(cell, seq(nrow(cell)))
cell[[1]]

degs <- FindMarkers(combined.label, ident.1 = paste0("Ecig_ME_", cell[[1]]),
                    ident.2 = paste0("Control_MC_", cell[[1]]), test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$cluster <- as.character(cell[[1]])
degs$gene <- row.names(degs)
df <- degs

for (i in 2:16)
{
  degs <- FindMarkers(combined.label, ident.1 = paste0("Ecig_ME_", cell[[i]]),
                      ident.2 = paste0("Control_MC_", cell[[i]]), test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$cluster <- as.character(cell[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}
saveRDS(df, file="out/final/3mon.male.deg.df.rds")

# filter deg by p.adj
deg.sig <- df[df$p_val_adj < 0.05, ]
write.csv(deg.sig, file = "out/final/3mon.male.degs.sig.csv")
table(deg.sig$cluster)
deg.3mon.male <- read.csv("out/final/3mon.male.degs.sig.csv", header = T, row.names = 1)
########### DEG female ###############
## DEG between 3 month female e-cig and control######
Idents(combined.label) <- "group.sample.stim"

degs <- FindMarkers(combined.label, ident.1 = paste0("Ecig_FE_", cell[[1]]),
                    ident.2 = paste0("Control_FC_", cell[[1]]), test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$cluster <- as.character(cell[[1]])
degs$gene <- row.names(degs)
df <- degs

for (i in 2:16)
{
  degs <- FindMarkers(combined.label, ident.1 = paste0("Ecig_FE_", cell[[i]]),
                      ident.2 = paste0("Control_FC_", cell[[i]]), test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$cluster <- as.character(cell[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}
saveRDS(df, file="out/final/3mon.female.deg.df.rds")

# filter deg by p.adj
deg.sig <- df[df$p_val_adj < 0.05, ]
write.csv(deg.sig, file = "out/final/3mon.female.degs.sig.csv")
table(deg.sig$cluster)



