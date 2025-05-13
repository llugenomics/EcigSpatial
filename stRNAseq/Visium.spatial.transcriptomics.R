## Visium spatial transcriptomics data loading ##
#Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(glmGamPoi)
library(SCpubr)
library(stringr)

### create function for image rotation ##
rotimat=function(foo, rotation){
  if(!is.matrix(foo)){
    cat("Input is not a matrix")
    return(foo)
  }
  if(!(rotation %in% c("180","Hf", "Vf","R90", "L90"))){
    cat("Rotation should be either L90, R90, 180, Hf or Vf\n")
    return(foo)
  }
  if(rotation == "180"){
    foo <- foo %>%
      .[, dim(.)[2]:1] %>%
      .[dim(.)[1]:1,]
  }
  if(rotation == "Hf"){
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  if(rotation == "Vf"){
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "L90"){
    foo = t(foo)
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "R90"){
    foo = t(foo)
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  return(foo)
}

rotateSeuratImage = function(seuratVisumObject, slide = "MC67A", rotation="Vf"){
  if(!(rotation %in% c("180","Hf","Vf", "L90", "R90"))){
    cat("Rotation should be either 180, L90, R90, Hf or Vf\n")
    return(NULL)
  }else{
    seurat.visium = seuratVisumObject
    ori.array = (seurat.visium@images)[[slide]]@image
    img.dim = dim(ori.array)[1:2]/(seurat.visium@images)[[slide]]@scale.factors$lowres
    new.mx <- c()  
    # transform the image array
    for (rgb_idx in 1:3){
      each.mx <- ori.array[,,rgb_idx]
      each.mx.trans <- rotimat(each.mx, rotation)
      new.mx <- c(new.mx, list(each.mx.trans))
    }
    
    # construct new rgb image array
    new.X.dim <- dim(each.mx.trans)[1]
    new.Y.dim <- dim(each.mx.trans)[2]
    new.array <- array(c(new.mx[[1]],
                         new.mx[[2]],
                         new.mx[[3]]), 
                       dim = c(new.X.dim, new.Y.dim, 3))
    
    #swap old image with new image
    seurat.visium@images[[slide]]@image <- new.array
    
    ## step4: change the tissue pixel-spot index
    img.index <- (seurat.visium@images)[[slide]]@coordinates
    
    #swap index
    if(rotation == "Hf"){
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "Vf"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
    }
    
    if(rotation == "180"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "L90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[2]-img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.index$imagerow
    }
    
    if(rotation == "R90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[1]-img.index$imagerow
    }
    
    return(seurat.visium)
  }  
}

####load data#####
data.dir = "data/visium/"
samples <- c("P94A1","P95A1","P100A1", "P80A1","P81A1","P81A2","P83A1", "P112A1", "P113A1") # male anterior visium data
samples <- c("P93A1","P96A1","P101A1", "P74A1","P77A1","P82A2","P108A1", "P111A1") # female anterior visium data
samples <- c("P94P1","P95P1","P100P1", "P80P1","P81P1","P81P2","P83P1", "P112P1", "P113P1") # male posterior visium data
samples <- c("P93P1","P96P1","P101P1", "P74P1","P77P1","P82P2","P108P1", "P111P1") # female posterior visium data

#Select sample
curr.sample <- samples[1] # P94A1 as example here, done for all xx samples.
visium <- Load10X_Spatial(data.dir = paste0(data.dir, curr.sample), filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial", 
                             slice = curr.sample, 
                             filter.matrix = TRUE, 
                             to.upper = FALSE, 
                             image = NULL)
# Add metadata information
visium$treat <- "Control"
visium$sex <- "Male"
visium$ID <- "MC94A"

# Get mitochondrial reads in data
visium[["percent.mt"]] <- PercentageFeatureSet(visium, "^Mt-")

# Rotate the slide image as coronal right hemisphere
visium=rotateSeuratImage(visium,rotation = "L90", slide = curr.sample)
SpatialFeaturePlot(visium, features = "nCount_Spatial")

# Filtering based on nCount and nFeature
visium.filter <- subset(visium, subset = nFeature_Spatial > 500 & nCount_Spatial > 200 & percent.mt < 15)

# normalization
visium.sct <- SCTransform(visium.filter, 
                       assay = "Spatial", 
                       vst.flavor = "v2", verbose = FALSE)
MC1.sct <- visium.sct

# merge data
male.merge_sct <- merge(x = MC1.sct, y = list(MC2.sct, MC3_sct, ME1_sct, ME2_sct, 
                                              ME3_sct, ME4_sct, ME5_sct, ME6_sct), 
                        add.cell.ids = c("MC1", "MC2","MC3", "ME1","ME2","ME3", "ME4", "ME5","ME6"))

# Joint dimensional reduction on RNA expression data
VariableFeatures(male.merge_sct) <- c(VariableFeatures(MC1_sct), VariableFeatures(MC2_sct), VariableFeatures(MC3_sct), 
                                      VariableFeatures(ME1_sct), VariableFeatures(ME2_sct), VariableFeatures(ME3_sct),
                                      VariableFeatures(ME4_sct), VariableFeatures(ME5_sct), VariableFeatures(ME6_sct))
male.merge_sct <- RunPCA(male.merge_sct, verbose = FALSE)
male.merge_sct <- RunUMAP(male.merge_sct, reduction = "pca", dims = 1:30, verbose = FALSE)
male.merge_sct <- FindNeighbors(male.merge_sct, reduction = "pca", dims = 1:30, verbose = FALSE)
male.merge_sct <- FindClusters(male.merge_sct, resolution = 0.3, verbose = FALSE)

# Visualize
DimPlot(male.merge_sct, reduction = "umap", group.by = c("ident", "treat"))

# integration using harmony 
male.merge_sct_harmony <- IntegrateLayers(object = male.merge_sct, method = HarmonyIntegration,
                                          orig.reduction = "pca", normalization.method="SCT", new.reduction = "harmony",
                                          assay = "SCT", verbose = FALSE)
male.merge_sct_harmony <- FindNeighbors(male.merge_sct_harmony, reduction = "harmony", dims = 1:30)
male.merge_sct_harmony <- FindClusters(male.merge_sct_harmony, resolution = 0.4, cluster.name = "harmony_clusters")

# annotation regions 
male.merge_sct_harmony.label <- RenameIdents(male.merge_sct_harmony, `0` ="CPu",  `15` ="CPu",`13` ="CPu",
                                             `3` ="cc",`1` ="NAc",`9` ="Layer2/3",`2` ="Septal",`12` ="Pir",`5` ="Layer5/6",`11` ="LV",`7` ="Layer1",
                                             `6` ="ACC", `16` ="ACC",`4` ="Layer4",
                                             `8` ="AIC",`10` ="Sparse", `14` ="Layer5/6")
# add region and group information to the metadata
male.merge_sct_harmony.label$region <- Idents(male.merge_sct_harmony.label)
Idents(male.merge_sct_harmony.label) <- "treat"
male.merge_sct_harmony.label$cluster.stim <- paste(Idents(male.merge_sct_harmony.label), male.merge_sct_harmony.label$region, sep = "_")
saveRDS(male.merge_sct_harmony.label, file = "output/male.anterior.visium.rds")

########## fig.1e. male anterior dotplot #####
########## Suppl.Fig.2.MA.Dimplot / vlnplot #####
male.anterior.visium <- readRDS("output/final.final/visium.data/male.anterior.visium.rds") # filter datapoint
Idents(male.anterior.visium) <- "region"
markers.to.plot <- c("Slc22a6","Cux2","Tox","Cryab","Mbp","Drd1","Cartpt","Nr4a2",
                     "Egln3","Zic2","Tshz2","Dnah12") 
levels <- c("Layer1","Layer2/3", "Layer4","Layer5/6", 
            "cc","CPu","NAc","AIC","Pir","Septal", "ACC",
            "LV","Sparse")
male.anterior.visium@active.ident <- factor(male.anterior.visium@active.ident, levels = levels)
DotPlot(male.anterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))

pdf("output/final.final/figures/Fig.1e.MA.dotplot.pdf", width = 7, height = 4)
DotPlot(male.anterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))+RotatedAxis()
dev.off()
pdf("output/final.final/figures/Suppl.Fig.2.MA.Dimplot.pdf", width = 6, height = 4)
DimPlot(male.anterior.visium, reduction = "umap", group.by = "ident", label = T)
dev.off()
pdf("output/final.final/figures/Suppl.Fig.2.MA.region.gene.vlnplot.pdf", width = 7, height = 4)
VlnPlot(male.anterior.visium, features = "nFeature_Spatial", pt.size = 0)
dev.off()
pdf("output/final.final/figures/Suppl.Fig.2.MA.Dimplot.treat.pdf", width = 6, height = 4)
DimPlot(male.anterior.visium, reduction = "umap", group.by = "treat")
dev.off()

#### fig.1c.MA Dimplot ######
Idents(male.anterior.visium) <- "region"
levels <- c("CPu",  "Layer2/3","cc",
            "NAc","Layer4","Sparse","Layer5/6","Layer1",
            "ACC","AIC", "Septal",
            "LV", "Pir")
male.anterior.visium@active.ident <- factor(male.anterior.visium@active.ident, levels = levels)
pdf("output/final.final/figures/Fig1.MA.SpatialDimplot.anno.pdf", width = 15, height = 15)
SpatialDimPlot(male.anterior.visium, label = T, label.size = 3, pt.size.factor = 2, ncol = 3)
dev.off()

########## female anterior #####
########## Suppl.Fig.2.FA dotplot / Dimplot / vlnplot #####
female.anterior.visium <- readRDS("output/final.final/visium.data/female.anterior.visium.rds") # filter datapoint
Idents(female.anterior.visium) <- "region"
levels <- c("CPu",  "Layer2/3","cc",
            "NAc","Layer4","Sparse","Layer5/6","Layer1",
            "ACC","AIC", "Septal",
            "LV", "Pir")
female.anterior.visium@active.ident <- factor(female.anterior.visium@active.ident, levels = levels)

pdf("output/final.final/figures/Fig1.FA.SpatialDimplot.anno.pdf", width = 15, height = 15)
SpatialDimPlot(female.anterior.visium, label = T, label.size = 3, pt.size.factor = 2, ncol = 3)
dev.off()

# female anterior dotplot 
Idents(female.anterior.visium) <- "region"
markers.to.plot <- c("Slc22a6","Cux2","Tox","Cryab","Mbp","Drd1","Cartpt","Nr4a2",
                     "Egln3","Zic2","Tshz2","Dnah12") #"Sulf1",
levels <- c("Layer1","Layer2/3", "Layer4","Layer5/6", 
            "cc","CPu","NAc","AIC","Pir","Septal", "ACC",
            "LV","Sparse")
female.anterior.visium@active.ident <- factor(female.anterior.visium@active.ident, levels = levels)
DotPlot(female.anterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))

pdf("output/final.final/figures/Suppl.Fig.2.FA.dotplot.pdf", width = 7, height = 4)
DotPlot(female.anterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))+RotatedAxis()
dev.off()
pdf("output/final.final/figures/Suppl.Fig.2.FA.Dimplot.pdf", width = 6, height = 4)
DimPlot(female.anterior.visium, reduction = "umap", group.by = "ident", label = T)
dev.off()
pdf("output/final.final/figures/Suppl.Fig.2.FA.region.gene.vlnplot.pdf", width = 7, height = 4)
VlnPlot(female.anterior.visium, features = "nFeature_Spatial", pt.size = 0)
dev.off()

########## male posterior #####
########## Fig.1f.MP dotplot / Dimplot / supple.fig2.vlnplot #####
male.posterior.visium <- readRDS("output/final.final/visium.data/male.posterior.visium.rds") 
Idents(male.posterior.visium) <- "region"
levels <- c("HY", "Sparse", "cc", "CPu","Layer2/3","TH", "Layer5/6", "RSP","HP","Layer1","Amy", "LV","Layer4")
male.posterior.visium@active.ident <- factor(male.posterior.visium@active.ident, levels = levels)

pdf("output/final.final/figures/Fig1.MP.SpatialDimplot.anno.pdf", width = 15, height = 15)
SpatialDimPlot(male.posterior.visium, label = T, label.size = 3, pt.size.factor = 2, ncol = 3)
dev.off()
pdf("output/final.final/figures/Fig1.MP.SpatialDimplot.anno.no.label.pdf", width = 15, height = 15)
SpatialDimPlot(male.posterior.visium, label.size = 3, pt.size.factor = 2, ncol = 3)
dev.off()

# fig.1f.male posterior dotplot 
Idents(male.posterior.visium) <- "region"
markers.to.plot <- c("Slc22a6","Cux2","Etv1","Tle4","Mbp","Drd1","Tshz2","Neurod1","Tcf7l2","Ecel1",
                     "Nrp2",
                     "Dnah12") #"Sulf1","Tgfb2","Pbx3",
levels <- c("Layer1","Layer2/3", "Layer4","Layer5/6",
            "cc","CPu","RSP","HP","TH","HY","Amy",
            "LV","Sparse")
male.posterior.visium@active.ident <- factor(male.posterior.visium@active.ident, levels = levels)
DotPlot(male.posterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))

pdf("output/final.final/figures/Fig.1f.MP.dotplot.pdf", width = 7, height = 4)
DotPlot(male.posterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))+RotatedAxis()
dev.off()
pdf("output/final.final/figures/Suppl.Fig.2.MP.Dimplot.pdf", width = 6, height = 4)
DimPlot(male.posterior.visium, reduction = "umap", group.by = "ident", label = T)
dev.off()
pdf("output/final.final/figures/Suppl.Fig.2.MP.region.gene.vlnplot.pdf", width = 7, height = 4)
VlnPlot(male.posterior.visium, features = "nFeature_Spatial", pt.size = 0)
dev.off()

########## female posterior #####
########## Suppl.Fig.2.FP dotplot / Dimplot / vlnplot #####
female.posterior.visium <- readRDS("output/final.final/visium.data/female.posterior.visium.rds") 
levels <- c("HY", "Sparse", "cc", "CPu","Layer2/3","TH", "Layer5/6", "RSP","HP","Layer1","Amy", "LV","Layer4")
female.posterior.visium@active.ident <- factor(female.posterior.visium@active.ident, levels = levels)
Idents(female.posterior.visium) <- "region"
pdf("output/final.final/figures/Fig1.FP.SpatialDimplot.anno.pdf", width = 15, height = 15)
SpatialDimPlot(female.posterior.visium, label = T, label.size = 3, pt.size.factor = 2, ncol = 3)
dev.off()
pdf("output/final.final/figures/Fig1.FP.SpatialDimplot.anno.no.label.pdf.pdf", width = 15, height = 15)
SpatialDimPlot(female.posterior.visium, label.size = 3, pt.size.factor = 2, ncol = 3)
dev.off()

# female posterior dotplot 
markers.to.plot <- c("Slc22a6","Cux2","Etv1","Tle4","Mbp","Drd1","Tshz2","Neurod1","Tcf7l2","Ecel1",
                     "Nrp2",
                     "Dnah12") #"Sulf1","Tgfb2","Pbx3",
levels <- c("Layer1","Layer2/3", "Layer4","Layer5/6",
            "cc","CPu","RSP","HP","TH","HY","Amy",
            "LV","Sparse")
female.posterior.visium@active.ident <- factor(female.posterior.visium@active.ident, levels = levels)
DotPlot(female.posterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))

pdf("output/final.final/figures/supple.fig.2.FP.dotplot.pdf", width = 7, height = 4)
DotPlot(female.posterior.visium, features = markers.to.plot,cols = c("lightgrey", "red"))+RotatedAxis()
dev.off()

##### fig.1g upset plot for region specific marker genes on anterior control sections
##### region specific marker genes in male anterior####
Idents(male.anterior.visium) <- "treat"
MAC <- subset(x = male.anterior.visium,idents = "Control")
Idents(MAC) <- "region"
region <- read.csv("data/MA.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
MAC <- PrepSCTFindMarkers(MAC, assay = "SCT", verbose = TRUE)
# FinderMarkers for each subregion
degs <- FindMarkers(MAC, ident.1 = region[[1]],
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs
for (i in 2:13)
{
  degs <- FindMarkers(MAC, ident.1 = region[[i]],
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}

# filter deg by p.adj
deg.sig <- df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 1, ]
dim(deg.sig) #5482
write.csv(df, file = "output/final.final/table/fig1.df.MA.region.specific.csv")
write.csv(deg.sig, file = "output/final.final/table/fig1.MA.region.specific.degs.sig.csv")

##### region specific marker genes in female anterior####
female.anterior.visium <- readRDS("output/final.final/visium.data/female.anterior.visium.rds") # filter datapoint
Idents(female.anterior.visium) <- "treat"
FAC <- subset(x = female.anterior.visium,idents = "Control")
Idents(FAC) <- "region"
region <- read.csv("data/MA.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
region[[1]]
FAC <- PrepSCTFindMarkers(FAC, assay = "SCT", verbose = TRUE)
# FinderMarkers for each subregion
degs <- FindMarkers(FAC, ident.1 = region[[1]],
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs
for (i in 2:13)
{
  degs <- FindMarkers(FAC, ident.1 = region[[i]],
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}

# filter deg by p.adj
deg.sig <- df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 1, ]
dim(deg.sig) #5703
write.csv(df, file = "output/final.final/table/fig1.df.FA.region.specific.csv")
write.csv(deg.sig, file = "output/final.final/table/fig1.FA.region.specific.degs.sig.csv")

##################################################################
## common region specific marker genes in male and female anterior
library(dplyr)
library(purrr)
MA <- read.csv("output/final.final/table/fig1.df.MA.region.specific.csv", row.names = 1)
FA <- read.csv("output/final.final/table/fig1.df.FA.region.specific.csv", row.names = 1)
# filter region specific marker genes with 1.5 fold increased positive expression
MA <- MA[MA$avg_log2FC > 0 & MA$p_val_adj < 0.05 & abs(MA$avg_log2FC) > 0.6, ]
FA <- FA[FA$avg_log2FC > 0 & FA$p_val_adj < 0.05 & abs(FA$avg_log2FC) > 0.6, ]
# common region specific marker genes between male and female anterior
common_rows <- reduce(list(MA, FA), inner_join, by = c("gene", "region"))
table(common_rows$region)
write.csv(common_rows, file = "output/final.final/table/fig1g.MA.FA.region.specific.common.degs.csv")

## Upset plotting fig.1g
library(UpSetR)
common_A <- read.csv("output/final.final/table/fig1.MA.FA.region.specific.common.degs.csv", row.names = 1)
CPu <- common_A[common_A$region=="CPu",]
Septal <- common_A[common_A$region=="Septal",]
cc <- common_A[common_A$region=="cc",]
ACC <- common_A[common_A$region=="ACC",]
AIC <- common_A[common_A$region=="AIC",]
Layer2.3 <- common_A[common_A$region=="Layer2/3",]
Layer4 <- common_A[common_A$region=="Layer4",]
Layer5.6 <- common_A[common_A$region=="Layer5/6",]
NAc <- common_A[common_A$region=="NAc",]
Pir <- common_A[common_A$region=="Pir",]
LV <- common_A[common_A$cluster=="LV",]

deginput <- list(CPu =CPu$gene, Septal= Septal$gene,
                 ACC = ACC$gene,  Pir = Pir$gene,Layer2.3 = Layer2.3$gene, 
                 Layer4 = Layer4$gene, Layer5.6 = Layer5.6$gene, NAc = NAc$gene)

png("output/final.final/figures/fig1g.Anterior.region.specific.png", width = 5000, height = 3500, res = 600) 
upset(fromList(deginput),
      sets = c("CPu", "Septal",  "ACC","Pir",
               "Layer2.3", "Layer4", "Layer5.6", "NAc"),
      mainbar.y.label = "Region-specific genes (Anterior)",
      sets.x.label = "Regions",
      order.by = "freq",
      nintersects = 12,
      main.bar.color = "grey20", 
      sets.bar.color = "seagreen3", 
      matrix.color = "grey40",
      text.scale = c(2.2, 1.6, 1.5, 1.6, 1.6, 1.6))
dev.off()

############# posterior 
######## MP ###########
male.posterior.visium <- readRDS("output/final.final/visium.data/male.posterior.visium.rds") 
Idents(male.posterior.visium) <- "treat"
MPC <- subset(x = male.posterior.visium,idents = "Control")
Idents(MPC) <- "region"
region <- read.csv("data/MP.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
MPC <- PrepSCTFindMarkers(MPC, assay = "SCT", verbose = TRUE)
# FinderMarkers for each subregion
degs <- FindMarkers(MPC, ident.1 = region[[1]],
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs

for (i in 2:13)
{
  degs <- FindMarkers(MPC, ident.1 = region[[i]],
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}

# filter deg by p.adj
deg.sig <- df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 1, ]
dim(deg.sig) #6215
saveRDS(df, file = "output/final.final/table/fig1.df.MP.region.specific.rds")
write.csv(df, file = "output/final.final/table/fig1.df.MP.region.specific.csv")
write.csv(deg.sig, file = "output/final.final/table/fig1.MP.region.specific.degs.sig.csv")

##### FP ###
female.posterior.visium <- readRDS("output/final.final/visium.data/female.posterior.visium.rds") 
Idents(female.posterior.visium) <- "treat"
FPC <- subset(x = female.posterior.visium,idents = "Control")
Idents(FPC) <- "region"
region <- read.csv("data/MP.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
FPC <- PrepSCTFindMarkers(FPC, assay = "SCT", verbose = TRUE)
# FinderMarkers for each subregion
degs <- FindMarkers(FPC, ident.1 = region[[1]],
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs

for (i in 2:13)
{
  degs <- FindMarkers(FPC, ident.1 = region[[i]],
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}

# filter deg by p.adj
deg.sig <- df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 1, ]
dim(deg.sig) #6085
saveRDS(df, file = "output/final.final/table/fig1.df.FP.region.specific.rds")
write.csv(df, file = "output/final.final/table/fig1.df.FP.region.specific.csv")
write.csv(deg.sig, file = "output/final.final/table/fig1.FP.region.specific.degs.sig.csv")

##################################################################
## common region specific marker genes in male and female posterior
library(dplyr)
library(purrr)
MP <- read.csv("output/final.final/table/fig1.df.MP.region.specific.csv", row.names = 1)
FP <- read.csv("output/final.final/table/fig1.df.FP.region.specific.csv", row.names = 1)
# filter region specific marker genes with 1.5 fold increased positive expression
MP <- MP[MP$avg_log2FC > 0 & MP$p_val_adj < 0.05 & abs(MP$avg_log2FC) > 0.6, ]
FP <- FP[FP$avg_log2FC > 0 & FP$p_val_adj < 0.05 & abs(FP$avg_log2FC) > 0.6, ]
# common region specific marker genes between male and female anterior
common_P <- reduce(list(MP, FP), inner_join, by = c("gene", "region"))
write.csv(common_P, file = "output/final.final/table/fig1h.MP.FP.region.specific.common.degs.csv")

## Upset plotting fig.1h
######## Posterior region specific visium degs #########
table(common_P$region)

CPu <- common_P[common_P$region=="CPu",]
Amy <- common_P[common_P$region=="Amy",]
cc <- common_P[common_P$region=="cc",]
HP <- common_P[common_P$region=="HP",]
HY <- common_P[common_P$region=="HY",]
Layer2.3 <- common_P[common_P$region=="Layer2/3",]
Layer4 <- common_P[common_P$region=="Layer4",]
Layer5.6 <- common_P[common_P$region=="Layer5/6",]
RSP <- common_P[common_P$region=="RSP",]
TH <- common_P[common_P$region=="TH",]
LV <- common_P[common_P$region=="LV",]

deginput <- list(CPu =CPu$gene, 
                 cc = cc$gene, HP = HP$gene, HY = HY$gene, Layer2.3 = Layer2.3$gene, 
                 Layer4 = Layer4$gene, Layer5.6 = Layer5.6$gene, TH = TH$gene, LV = LV$gene)

png("output/final.final/figures/fig1h.Posterior.region.specific.png", width = 5000, height = 3500, res = 600) 
upset(fromList(deginput),
      sets = c("CPu", "HP", "cc",  "HY",
               "Layer2.3", "Layer4", "Layer5.6", "TH","LV"),
      mainbar.y.label = "Region-specific genes (Posterior)",
      sets.x.label = "Regions",
      order.by = "freq",
      nintersects = 12,
      main.bar.color = "grey20", 
      sets.bar.color = "seagreen3", 
      matrix.color = "grey40",
      text.scale = c(2.2, 1.6, 1.5, 1.6, 1.6, 1.6))
dev.off()

##### DEGs in male anterior visium data
male.anterior.visium <- readRDS("output/final.final/visium.data/male.anterior.visium.rds") # filter datapoint
Idents(male.anterior.visium) <- "cluster.stim"
region <- read.csv("data/MA.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
# Finder degs between ecig and control for each subregion
degs <- FindMarkers(male.anterior.visium, ident.1 = paste0("Ecig_", region[[1]]),
                    ident.2 = paste0("Control_", region[[1]]),
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs
for (i in 2:13)
{
  degs <- FindMarkers(male.anterior.visium, ident.1 = paste0("Ecig_", region[[i]]),
                      ident.2 = paste0("Control_", region[[i]]),
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}
degs.sig <- df[df$p_val_adj < 0.05, ]
write.csv(degs.sig, file = "output/final.final/table/fig3.MA.degs.sig.csv")

#########################################
##### DEGs in female anterior visium data
female.anterior.visium <- readRDS("output/final.final/visium.data/female.anterior.visium.rds") # filter datapoint
Idents(female.anterior.visium) <- "cluster.stim"
region <- read.csv("data/MA.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
# Finder degs between ecig and control for each subregion
degs <- FindMarkers(female.anterior.visium, ident.1 = paste0("Ecig_", region[[1]]),
                    ident.2 = paste0("Control_", region[[1]]),
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs
for (i in 2:13)
{
  degs <- FindMarkers(female.anterior.visium, ident.1 = paste0("Ecig_", region[[i]]),
                      ident.2 = paste0("Control_", region[[i]]),
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}
degs.sig <- df[df$p_val_adj < 0.05, ]
write.csv(degs.sig, file = "output/final.final/table/fig3.FA.degs.sig.csv")

#########################################
##### DEGs in male posterior visium data
male.posterior.visium <- readRDS("output/final.final/visium.data/male.posterior.visium.rds") 
Idents(male.posterior.visium) <- "treat"
male.posterior.visium$cluster.stim <- paste(Idents(male.posterior.visium), male.posterior.visium$region, sep = "_")
head(male.posterior.visium@meta.data)
Idents(male.posterior.visium) <- "cluster.stim"
region <- read.csv("data/MP.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
# Finder degs between ecig and control for each subregion
degs <- FindMarkers(male.posterior.visium, ident.1 = paste0("Ecig_", region[[1]]),
                    ident.2 = paste0("Control_", region[[1]]),
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs
for (i in 2:13)
{
  degs <- FindMarkers(male.posterior.visium, ident.1 = paste0("Ecig_", region[[i]]),
                      ident.2 = paste0("Control_", region[[i]]),
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}
degs.sig <- df[df$p_val_adj < 0.05, ]
write.csv(degs.sig, file = "output/final.final/table/fig3.MP.degs.sig.csv")

#########################################
##### DEGs in female posterior visium data
female.posterior.visium <- readRDS("output/final.final/visium.data/female.posterior.visium.rds") 
Idents(female.posterior.visium) <- "treat"
female.posterior.visium$cluster.stim <- paste(Idents(female.posterior.visium), female.posterior.visium$region, sep = "_")
head(female.posterior.visium@meta.data)
Idents(female.posterior.visium) <- "cluster.stim"
region <- read.csv("data/MP.region.csv", header = F)
colnames(region) <- "Region"
region <- split(region, seq(nrow(region)))
# Finder degs between ecig and control for each subregion
degs <- FindMarkers(female.posterior.visium, ident.1 = paste0("Ecig_", region[[1]]),
                    ident.2 = paste0("Control_", region[[1]]),
                    test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$region <- as.character(region[[1]])
degs$gene <- row.names(degs)
df <- degs
for (i in 2:13)
{
  degs <- FindMarkers(female.posterior.visium, ident.1 = paste0("Ecig_", region[[i]]),
                      ident.2 = paste0("Control_", region[[i]]),
                      test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$region <- as.character(region[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}
degs.sig <- df[df$p_val_adj < 0.05, ]
write.csv(degs.sig, file = "output/final.final/table/fig3.FP.degs.sig.csv")

### fig.3a bargraph of MA FA common degs
MA <- readRDS("output/final.final/visium.data/degs/male.anterior.degs.total.rds")
FA <- readRDS("output/final.final/visium.data/degs/female.anterior.degs.total.rds")
MP <- readRDS("output/final.final/visium.data/degs/male.posterior.degs.total.rds")
FP <- readRDS("output/final.final/visium.data/degs/female.posterior.degs.total.rds")
MA <- MA[MA$p_val_adj < 0.05, ]
FA <- FA[FA$p_val_adj < 0.05, ]
MP <- MP[MP$p_val_adj < 0.05, ]
FP <- FP[FP$p_val_adj < 0.05, ]
# common genes from male and female anterior visium data, respectively
common_rows <- reduce(list(MA, FA), inner_join, by = c("Gene", "cluster"))
MA <- common_rows[,c(2,5,6,7)]
colnames(MA) <- c("avg_log2FC", "p_val_adj", "cluster","Gene")
FA <- common_rows[,c(6,7,9,12)]
colnames(FA) <- c("cluster","Gene","avg_log2FC", "p_val_adj")
MA$Sex <- "Male"
FA$Sex <- "Female"
A <- rbind(MA, FA)
# common genes from male and female posterior visium data, respectively
common_rows <- reduce(list(MP, FP), inner_join, by = c("Gene", "cluster"))
MP <- common_rows[,c(2,5,6,7)]
colnames(MP) <- c("avg_log2FC", "p_val_adj", "cluster","Gene")
FP <- common_rows[,c(6,7,9,12)]
colnames(FP) <- c("cluster","Gene","avg_log2FC", "p_val_adj")
MP$Sex <- "Male"
FP$Sex <- "Female"
P <- rbind(MP, FP)

################################
library(stringr)
deg <- A # calculate for anterior
deg <- P # calculate for posterior

deg_up <- deg[deg$avg_log2FC > 0,]
deg_up$Sex_cluster <- paste0(deg_up$cluster, "_", deg_up$Sex)
deg.count <- as.data.frame(table(deg_up$Sex_cluster))
deg.count$dir <- "up"
deg_down <- deg[deg$avg_log2FC < 0,]
deg_down$Sex_cluster <- paste0(deg_down$cluster, "_", deg_down$Sex)

deg.count <- as.data.frame(table(deg_up$Sex_cluster))
deg.count$dir <- "up"
deg.count2 <- as.data.frame(table(deg_down$Sex_cluster))
deg.count2$dir <- "down"
deg.count2$Freq <- deg.count2$Freq * -1
deg.count <- rbind(deg.count, deg.count2)

deg.count[c('Region', 'Sex')] <- str_split_fixed(deg.count$Var1, '_', 2)
deg.count$col <- paste0(deg.count$dir, "_", deg.count$Sex)
saveRDS(deg.count, file = "output/final.final/visium.data/Fig3a.deg.bar.rds")
saveRDS(deg.count, file = "output/final.final/visium.data/Fig3b.deg.bar.rds")

# plot fig3a bargraph
deg.count$Var1 <- factor(deg.count$Var1, levels = c("Pir_Female", "Pir_Male","LV_Female", "LV_Male","AIC_Female", "AIC_Male",
                                                    "Layer4_Female", "Layer4_Male","Septal_Female", "Septal_Male",
                                                    "ACC_Female", "ACC_Male","NAc_Female", "NAc_Male","Layer2/3_Female", "Layer2/3_Male",
                                                    "Layer5/6_Female", "Layer5/6_Male","cc_Female", "cc_Male","CPu_Female", "CPu_Male"))
pdf("output/final.final/figures/fig3a.MA.FA.common.deg.bar.pdf", width = 4.5, height = 5)
ggplot(deg.count, aes(x=Var1, y=Freq, fill = col)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  scale_y_continuous(labels = abs) +
  theme_minimal() +
  labs(title = "Number of DEGs", x = "Region", y = "Count") +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  ) +
  scale_fill_manual(values = c("up_Male" = "red", "down_Male" = "blue","up_Female" = "salmon", "down_Female" = "steelblue"))
dev.off()

# plot fig3b bargraph
deg.count$Var1 <- factor(deg.count$Var1, levels = c("LV_Female", "LV_Male","RSP_Female", "RSP_Male","HP_Female", "HP_Male",
                                                    "Amy_Female", "Amy_Male","Layer4_Female", "Layer4_Male",
                                                    "cc_Female", "cc_Male","Layer5/6_Female", "Layer5/6_Male","CPu_Female", "CPu_Male",
                                                    "TH_Female", "TH_Male","Layer2/3_Female", "Layer2/3_Male","HY_Female", "HY_Male"))

pdf("output/final.final/figures/fig3b.MP.FP.common.deg.bar.pdf", width = 4.5, height = 5)
ggplot(deg.count, aes(x=Var1, y=Freq, fill = col)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  scale_y_continuous(labels = abs) +
  theme_minimal() +
  labs(title = "Number of DEGs", x = "Region", y = "Count") +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  ) +
  scale_fill_manual(values = c("up_Male" = "red", "down_Male" = "blue","up_Female" = "salmon", "down_Female" = "steelblue"))
dev.off()

##### fig.3c.upset plot
### anterior regions ####
MA <- readRDS("output/final.final/visium.data/degs/male.anterior.degs.total.rds")
FA <- readRDS("output/final.final/visium.data/degs/female.anterior.degs.total.rds")
MA <- MA[MA$p_val_adj < 0.05, ]
FA <- FA[FA$p_val_adj < 0.05, ]
common_rows <- reduce(list(MA, FA), inner_join, by = c("Gene", "cluster"))
CPu <- common_rows[common_rows$cluster=="CPu",]
cc <- common_rows[common_rows$cluster=="cc",]
Pir <- common_rows[common_rows$cluster=="Pir",]
AIC <- common_rows[common_rows$cluster=="AIC",]
Layer2.3 <- common_rows[common_rows$cluster=="Layer2/3",]
Layer4 <- common_rows[common_rows$cluster=="Layer4",]
Layer5.6 <- common_rows[common_rows$cluster=="Layer5/6",]
Septal <- common_rows[common_rows$cluster=="Septal",]
NAc <- common_rows[common_rows$cluster=="NAc",]
LV <- common_rows[common_rows$cluster=="LV",]
ACC <- common_rows[common_rows$cluster=="ACC",]
deginput <- list(CPu =CPu$Gene, 
                 NAc = NAc$Gene, Layer2.3 = Layer2.3$Gene, 
                 Layer4 = Layer4$Gene, Layer5.6 = Layer5.6$Gene, 
                 Septal = Septal$Gene,ACC = ACC$Gene, cc = cc$Gene)
# fig.3c
png("output/final.final/figures/fig3c.Anterior.upset.seven.regions.png", width = 4000, height = 3000, res = 600) 
upset(fromList(deginput),
      sets = c("CPu", "NAc",
               "Layer2.3", "Layer4", "Layer5.6", "Septal", "ACC","cc"),
      mainbar.y.label = "No. of DEGs (Anterior)",
      sets.x.label = "Regions",
      order.by = "freq",
      nintersects = 13,
      main.bar.color = "grey20", 
      sets.bar.color = "seagreen3", 
      matrix.color = "grey40",
      text.scale = c(2.1, 1.6, 1.6, 1.4, 1.4, 1.8)) # #CCCCFF #shade.color = "green"
dev.off()

##### fig.3d.upset plot
### posterior regions ####
MP <- readRDS("output/final.final/visium.data/degs/male.posterior.degs.total.rds")
FP <- readRDS("output/final.final/visium.data/degs/female.posterior.degs.total.rds")
MP <- MP[MP$p_val_adj < 0.05, ]
FP <- FP[FP$p_val_adj < 0.05, ]
common_rows <- reduce(list(MP, FP), inner_join, by = c("Gene", "cluster"))
CPu <- common_rows[common_rows$cluster=="CPu",]
Amy <- common_rows[common_rows$cluster=="Amy",]
cc <- common_rows[common_rows$cluster=="cc",]
HP <- common_rows[common_rows$cluster=="HP",]
HY <- common_rows[common_rows$cluster=="HY",]
Layer2.3 <- common_rows[common_rows$cluster=="Layer2/3",]
Layer4 <- common_rows[common_rows$cluster=="Layer4",]
Layer5.6 <- common_rows[common_rows$cluster=="Layer5/6",]
RSP <- common_rows[common_rows$cluster=="RSP",]
TH <- common_rows[common_rows$cluster=="TH",]
LV <- common_rows[common_rows$cluster=="LV",]
# fig.3d
deginput <- list(CPu =CPu$Gene, Amy= Amy$Gene,
                 cc = cc$Gene, HP = HP$Gene, HY = HY$Gene, Layer2.3 = Layer2.3$Gene, 
                 Layer4 = Layer4$Gene, Layer5.6 = Layer5.6$Gene,TH = TH$Gene)
png("output/final.final/figures/fig3d.posterior.upset.seven.regions.png", width = 4000, height = 3500, res = 600) 
upset(fromList(deginput),
      sets = c("CPu", "HP", "cc", "Amy", "HY",
               "Layer2.3", "Layer4", "Layer5.6", "TH"),
      mainbar.y.label = "No. of DEGs (Posterior)",
      sets.x.label = "Regions",
      order.by = "freq",
      nintersects = 13,
      main.bar.color = "grey20", 
      sets.bar.color = "seagreen3", 
      matrix.color = "grey40",
      text.scale = c(2.1, 1.6, 1.6, 1.4, 1.4, 1.8)) # #CCCCFF #shade.color = "green"
dev.off()











