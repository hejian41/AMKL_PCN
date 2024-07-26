# 1. load data
counts_pre <- Read10X("../rawdata/HSCT0314/", gene.column = 1, cell.column = 1, unique.features = T)
counts_post <- Read10X("../rawdata/HSCT0323/", gene.column = 1, cell.column = 1, unique.features = T)
dim(counts_pre);dim(counts_post)
colnames(counts_pre) <- paste("pre",colnames(counts_pre), sep = "_")
colnames(counts_post) <- paste("post",colnames(counts_post), sep = "_")

seu_pre <- CreateSeuratObject(counts = counts_pre, min.cells = 5, min.features = 200)
seu_post <- CreateSeuratObject(counts = counts_post, min.cells = 5, min.features = 200)
seu <- merge(seu_pre, seu_post)
dim(counts_pre);dim(counts_post);dim(seu)
rm(seu_pre, seu_post)

seu_list <- SplitObject(seu, split.by = "orig.ident")

seu_list <- DoubletFinder_multisamples(seu_list = seu_list)

table(seu_list$pre$doubFind_res)
table(seu_list$post$doubFind_res)

seu <- merge(seu_list$pre, seu_list$post)
rm(seu_list, counts_pre, counts_post)
table(seu$doubFind_res)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, features = c("nCount_RNA","nFeature_RNA","percent.mt"), group.by = "orig.ident", pt.size = 0)
ncol(seu)
seu <- subset(seu, percent.mt < 10)
seu <- subset(seu, doubFind_res == "Singlet")

saveRDS(seu, file = "./seu_all.Rds")


# 2. Figure 1
## 2.1 pre-data
seu <- readRDS("./seu_all.Rds")
seu_pre <- subset(seu, orig.ident == "pre")
seu_pre <- NormalizeData(seu_pre, scale.factor = 1e4)
seu_pre <- FindVariableFeatures(seu_pre, nfeatures = 2000)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu_pre <- CellCycleScoring(seu_pre, s.features = s.genes, g2m.features = g2m.genes)
seu_pre <- ScaleData(seu_pre, vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

seu_pre <- RunPCA(seu_pre, npcs = 50)
seu_pre <- RunUMAP(seu_pre, dims = 1:20)

seu_pre <- FindNeighbors(seu_pre, dims = 1:20)
seu_pre <- FindClusters(seu_pre, resolution = 0.25)


DimPlot(seu_pre, reduction = "umap", group.by = "seurat_clusters", label = T, cols = louvain_colors, 
        repel = T, label.size = 4, raster = T)

Idents(seu_pre) <- "seurat_clusters"
seu_pre <- RenameIdents(seu_pre, '0' = "MK", '1' = "T_CD8A", '2' = "T_IL7R", '3' = "T_CD247", 
                        '4' = "Ery", '5' = "Plasma", '6' = "Mono", '7' = "B", '8' = "Mac")
seu_pre$names <- Idents(seu_pre)
seu_pre$names <- factor(seu_pre$names, levels = c("T_IL7R", "T_CD8A", "T_CD247", "B", "Plasma", 
                                                  "Mono", "Mac", "MK","Ery"))
ph_colors$names <- setNames(c(louvain_colors[1:3], louvain_colors[c(5:7, 10, 8:9)]), levels(seu_pre$names))


Idents(seu_pre) <- "names"
Degs_pre <- FindAllMarkers(seu_pre, min.pct = 0.25, logfc.threshold = 0.5, only.pos = T)
Degs_pre_sig <- subset(Degs_pre, p_val_adj < 0.05)
table(Degs_pre_sig$cluster)
write.csv(Degs_pre_sig, file = "Degs_cluster.csv")

saveRDS(seu_pre, file = "./rds_files/seu_pre.Rds")


## 2.2 plot Fig1A,B
library(SCP)
fig1a <- ClassDimPlot(seu_pre, group.by = c("names"), reduction = "UMAP", 
                      theme_use = "theme_blank", label = T, pt.size = 1,
                      raster = T, label_repel = T, title = "4,325 cells")
ggsave(filename = "./Fig1a.pdf", plot = fig1a, width = 6, height = 5, dpi = 300)



marker_genes = c('IL7R', "TCF7", "TRAC","IL32", "CD2", "CD3G","CD8A", "GZMK", 
                 "CCR5","EOMES","GZMA", "CCL5", "NKG7","CD247","GZMB",
                 "CD22","MS4A1","CD19","BLK","IGKC","IGLC2","JCHAIN", "TNFRSF17", "SDC1",
                 "CSF3R","CD14","CD68","CSF1R", "FCGR3A", "ITGA2B","GFI1B", "GP1BB","GATA1",
                 "HBB","HBA1","ALAS2")

ht <- GroupHeatmap(
  srt = seu_pre,
  features = marker_genes,
  group.by = c("names"),
  show_row_names = TRUE, #show_column_names = T, 
  add_dot = TRUE, add_reticle = TRUE
)
fig1b <- ht$plot

ggsave(filename = "./Fig1b.pdf", plot = fig1b, width = 7, height = 10, dpi = 300)


## 2.3 reference projection
## 2.3.1 preprocessing 
library(SeuratDisk)
Convert("./HCA_immune/BM_adata.h5ad", dest = "h5seurat", overwrite = F)
seu_hBM <- LoadH5Seurat("./HCA_immune/BM_adata.h5seurat",meta.data = F, misc=F)

library(anndata)
anndata <- read_h5ad("./HCA_immune/BM_adata.h5ad")
anno_hBM <- anndata$obs
anno_hBM <- as.data.frame(anno_hBM)
colnames(anno_hBM)
anno_hBM <- anno_hBM[, c(4,19,31)]
tmp <- reshape2::colsplit(anno_hBM$cell_names, "-", c("name","batch"))
anno_hBM[,c("name2","batch")] <- tmp

seu_hBM$batch <- anno_hBM$batch
seu_hBM$barcode <- anno_hBM$cell_names

## take the sample
sample_cells <- c()
for (i in sort(unique(seu_hBM$batch))) {
  sample_tmp <- sample(rownames(subset(seu_hBM@meta.data, batch == i)), 
                       size = ceiling(60000*nrow(subset(seu_hBM@meta.data, 
                                                        batch == i))/nrow(seu_hBM@meta.data)))
  sample_cells <- c(sample_cells, sample_tmp)
}

seu_hBM2 <- seu_hBM[, sample_cells]
table(seu_hBM2$batch)
table(seu_hBM$batch)

saveRDS(sample_cells, file = "./HCA_immune/sample_cell_IDs.Rds")

seu_hBM2 <- PercentageFeatureSet(seu_hBM2, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(seu_hBM2, features = "percent.mt", group.by = "batch")
VlnPlot(seu_hBM2, features = "nFeature_RNA", group.by = "batch")
seu_hBM2 <- subset(seu_hBM2, subset = percent.mt < 10)
seu_hBM2 <- subset(seu_hBM2, subset = nFeature_RNA > 500)

seu_hBM2_mnn <- NormalizeData(seu_hBM2)
seu_hBM2_mnn <- FindVariableFeatures(seu_hBM2_mnn, nfeatures = 3000)
seu_hBM2_mnn <- RunFastMNN(object.list = SplitObject(seu_hBM2_mnn, split.by = "batch"))

seu_hBM2_mnn <- RunUMAP(seu_hBM2_mnn, reduction = "mnn", dims = 1:30, return.model = T)


DimPlot(seu_hBM2_mnn, group.by = c("batch"), ncol = 2, reduction = "umap", cols = mycolors, label = F)


seu_hBM2_mnn <- FindNeighbors(seu_hBM2_mnn, reduction = "mnn", dims = 1:30)
seu_hBM2_mnn <- FindClusters(seu_hBM2_mnn, resolution = 0.4)

Idents(seu_hBM2_mnn) <- "seurat_clusters"

seu_mnn_c7 <- subset(seu_hBM2_mnn, idents = 7)
seu_mnn_c7 <- FindNeighbors(seu_mnn_c7, dims = 1:30, k.param = 30, reduction = "mnn")
seu_mnn_c7 <- FindClusters(seu_mnn_c7, resolution = 0.2, algorithm = 1)
seu_hBM2_mnn$sub_cluster <- as.character(Idents(seu_hBM2_mnn))
seu_hBM2_mnn$sub_cluster[Cells(seu_mnn_c7)] <- paste("c7",Idents(seu_mnn_c7), sep = "_")

seu_mnn_c13 <- subset(seu_hBM2_mnn, idents = 13)
seu_mnn_c13 <- FindNeighbors(seu_mnn_c13, dims = 1:30, k.param = 30, reduction = "mnn")
seu_mnn_c13 <- FindClusters(seu_mnn_c13, resolution = 0.2, algorithm = 1)
seu_hBM2_mnn$sub_cluster[Cells(seu_mnn_c13)] <- paste("c13",Idents(seu_mnn_c13), sep = "_")


DimPlot(seu_hBM2_mnn, group.by = c("sub_cluster"), ncol = 1, reduction = "umap", 
        cols = hej_colors, label = T)

anno_hBM <- seu_hBM2_mnn@meta.data

library(future)
supportsMulticore()
plan(strategy = "multicore", workers = 10)
Idents(seu_hBM2_mnn) <- "sub_cluster"
Idents(seu_hBM2_mnn) <- "cluster"
allmarkers_hBM <- FindAllMarkers(seu_hBM2_mnn, logfc.threshold = 0.4, min.pct = 0.2, only.pos = T)
allmarkers_hBM_sig <- subset(allmarkers_hBM, p_val_adj < 0.05)
table(allmarkers_hBM_sig$cluster)
write.csv(allmarkers_hBM_sig, file = "./HCA_immune/Degs_hBM.csv")

### annotation acoording to the DEGs
# 0: Naive_T_Cell
# 1: cytotoxic_T_cell
# 2: CD14+Monocyte
# 3: B_Cell
# 4: Naive_B_Cell
# 5: Erythrocyte
# 6: cDC
# c7_0: HSPC
# c7_1: Myeloblast
# c7_2: MEP
# 8: pre_B_cell
# 9: CD16+Monocyte
# 10: pDC
# 11: Cyc_B_cell
# 12: plasma_cell
# c13_0: NK_cell
# c13_1: MK
# 14: Memory_plasma_cell
# 15: CD14+Monocyte


seu_hBM2_mnn <- RenameIdents(seu_hBM2_mnn, "0" = "Naive_T_Cell", "1" = "Cytotoxic_T_cell", "2" = "CD14+Monocyte", "3" = "B_Cell", 
                             "4" = "Naive_B_Cell", "5" = "Erythrocyte", "6" = "cDC", "c7_0" = "HSPC", "c7_1" = "Myeloblast", 
                             "c7_2" = "MEP", "8" = "pre_B_cell", "9" = "CD16+Monocyte", "10" = "pDC", "11" = "Cyc_B_cell", 
                             "12" = "Plasma_cell", "c13_0" = "NK_cell", "c13_1" = "MK", "14" = "Memory_plasma_cell", 
                             "15" = "CD14+Monocyte")
seu_hBM2_mnn$cluster <- Idents(seu_hBM2_mnn)
seu_hBM2_mnn$cluster <- factor(seu_hBM2_mnn$cluster, 
                               levels = c("HSPC","Myeloblast","cDC","CD14+Monocyte","CD16+Monocyte","pDC",
                                          "MEP","MK","Erythrocyte",
                                          "Naive_T_Cell", "Cytotoxic_T_cell","Cyc_B_cell", "pre_B_cell", "Naive_B_Cell",
                                          "B_Cell", "Plasma_cell","Memory_plasma_cell","NK_cell"))

DimPlot(seu_hBM2_mnn, group.by = c("cluster"), ncol = 1, reduction = "umap",cols = hej_colors, label = T)
fig1c_1 <- DimPlot(seu_hBM2_mnn, group.by = c("cluster"), ncol = 1, reduction = "umap",cols = hej_colors, 
                   label = T, raster = T)
ggsave(filename = "./Fig1c_1.pdf", plot = fig1c_1, width = 8, height = 5, dpi = 300)


### 2.3.2 projection
features <- Reduce(union, list(VariableFeatures(seu_pre), VariableFeatures(seu_hBM2_mnn)))
seu_pre_mew <- RunKNNMap(srt_query = seu_pre, srt_ref = seu_hBM2_mnn, 
                         ref_group = "cluster", features = features,
                         ref_umap = "umap", k = 30)

fig1c <- ProjectionPlot(
  srt_query = seu_pre_mew, srt_ref = seu_hBM2_mnn, pt.size = 0.5,
  query_group = "names", ref_group = "cluster", stroke.highlight = 0.5,
  #query_param = list(palette = "Paired", raster = T), 
  ref_param = list(palette = "Set2", raster = T)
)

ggsave(filename = "./Fig1c.pdf", plot = fig1c, width = 8, height = 6, dpi = 300)

### 2.4 proportion
library(RColorBrewer)
colors <- brewer.pal(n = 9, name = "Paired")
ph_colors$names <- setNames(c(colors[1:4], colors[7:8], colors[c(5:6,9)]), levels(seu_pre$names))
mypie_plot(pbmc = seu_pre, group.by = "names", col.group = ph_colors$names, 
           split.by = "orig.ident")

p1 <- mypie_plot(pbmc = seu_pre, group.by = "names", col.group = c(colors[1:4], colors[7:8], colors[c(5:6,9)]), 
                 split.by = "orig.ident")

colors <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(18)
seu_hBM2_mnn$exp <- "nBM"
length(unique(seu_hBM2_mnn$cluster))
p2 <- mypie_plot(pbmc = seu_hBM2_mnn, group.by = "cluster", col.group = colors, 
                 split.by = "exp")
p2
prop.table(table(seu_hBM2_mnn$cluster))
fig1c_2 <- p1 + p2
fig1c_2

ggsave(filename = "./Fig1c_2.pdf", plot = fig1c_2, width = 16, height = 7, dpi = 300)


## 2.5 megakaryoblastic data
### 2.5.1 import data
count_18m <- Read10X("../rawdata/2022Frot/M7/run_count_A330/outs/filtered_feature_bc_matrix/")
count_19m <- Read10X("../rawdata/2022Frot/M7/run_count_M002//outs/filtered_feature_bc_matrix/")
count_49m <- Read10X("../rawdata/2022Frot/M7/M704/outs/filtered_feature_bc_matrix/")

colnames(count_18m) <- paste("18m", colnames(count_18m), sep = "_")
colnames(count_19m) <- paste("19m", colnames(count_19m), sep = "_")
colnames(count_49m) <- paste("49m", colnames(count_49m), sep = "_")

meta_mkb <- data.frame(row.names = c(colnames(count_18m), 
                                     colnames(count_19m), 
                                     colnames(count_49m)),
                       sample = c(rep("18m", ncol(count_18m)),
                                  rep("19m", ncol(count_19m)),
                                  rep("49m", ncol(count_49m))))

seu_mkb <- CreateSeuratObject(counts = cbind(count_18m, count_19m,count_49m),
                              meta.data = meta_mkb, min.cells = 3)


seu_list <- SplitObject(seu_mkb, split.by = "orig.ident")

seu_list <- DoubletFinder_multisamples(seu_list = seu_list)

table(seu_list$`18m`$doubFind_res)
table(seu_list$`19m`$doubFind_res)
table(seu_list$`49m`$doubFind_res)

seu_mkb <- merge(seu_list$`18m`, list(seu_list$`19m`, seu_list$`49m`))
rm(seu_list, count_18m, count_19m, count_49m)
table(seu_mkb$doubFind_res)

seu_mkb[["percent.mt"]] <- PercentageFeatureSet(seu_mkb, pattern = "^MT-")
VlnPlot(seu_mkb, features = c("nCount_RNA","nFeature_RNA","percent.mt"), 
        group.by = "orig.ident", pt.size = 0)
ncol(seu_mkb)
seu_mkb <- subset(seu_mkb, percent.mt < 10)
seu_mkb <- subset(seu_mkb, doubFind_res == "Singlet")

### 2.5.2 preprocess
seu_mkb <- NormalizeData(seu_mkb)
seu_mkb <- FindVariableFeatures(seu_mkb, nfeatures = 3000)
seu_mkb <- RunFastMNN(object.list = SplitObject(seu_mkb, split.by = "sample"))

seu_mkb <- RunUMAP(seu_mkb, reduction = "mnn", dims = 1:30, return.model = T)

DimPlot(seu_mkb, group.by = c("sample"), ncol = 2, reduction = "umap", cols = mycolors, label = F)


seu_mkb <- FindNeighbors(seu_mkb, reduction = "mnn", dims = 1:30)
seu_mkb <- FindClusters(seu_mkb, resolution = 0.2)

Idents(seu_mkb) <- "seurat_clusters"


DimPlot(seu_mkb, group.by = c("seurat_clusters"), ncol = 1, reduction = "umap", 
        cols = hej_colors, label = T)

library(future)
supportsMulticore()
plan(strategy = "multicore", workers = 10)
Idents(seu_mkb) <- "seurat_clusters"
Idents(seu_mkb) <- "cluster"
allmarkers_mkb <- FindAllMarkers(seu_mkb, logfc.threshold = 0.5, min.pct = 0.2, only.pos = T)
allmarkers_mkb_sig <- subset(allmarkers_mkb, p_val_adj < 0.05)
table(allmarkers_mkb_sig$cluster)
write.csv(allmarkers_mkb_sig, file = "../rawdata/2022Frot/Degs_mkb.csv")

# 0: MK_MYL4
# 1: Naive_T_Cell
# 2: Myeloblast
# 3: cDC
# 4: CD14+Monocyte
# 5: B_Cell
# 6: Neutrophil 
# 7: MK_GP1BA
# 8: cytotoxic_T_cell
# 9: cytotoxic_T_cell
# 10: Erythrocyte
# 11: Erythrocyte
# 12: plasma_cell
# 13: pre_B_cell
# 14: Memory_plasma_cell

seu_mkb <- RenameIdents(seu_mkb, "0" = "MK", "1" = "Naive_T_Cell", "2" = "Myeloblast", "3" = "cDC", 
                        "4" = "Monocyte", "5" = "B_Cell", "6" = "Neutrophil", "7" = "MK",  
                        "8" = "Cytotoxic_T_cell", "9" = "Cytotoxic_T_cell", "10" = "Erythrocyte", 
                        "11" = "Erythrocyte", "12" = "Plasma_cell", "13" = "pre_B_cell", "14" = "Memory_plasma_cell")
seu_mkb$cluster <- Idents(seu_mkb)
seu_mkb$cluster <- factor(seu_mkb$cluster, 
                          levels = c("Myeloblast","cDC","Monocyte","Neutrophil","MK","Erythrocyte",
                                     "Naive_T_Cell", "Cytotoxic_T_cell","pre_B_cell", 
                                     "B_Cell", "Plasma_cell","Memory_plasma_cell"))

DimPlot(seu_mkb, group.by = c("cluster"), ncol = 1, reduction = "umap",cols = hej_colors, label = T)


### 2.5.3 projection
features <- Reduce(union, list(VariableFeatures(seu_pre), VariableFeatures(seu_mkb)))
seu_pre_mew <- RunKNNMap(srt_query = seu_pre, srt_ref = seu_mkb, 
                         ref_group = "cluster", features = features,
                         ref_umap = "umap", k = 30)

fig1c_4 <- ProjectionPlot(
  srt_query = seu_pre_mew, srt_ref = seu_mkb, pt.size = 0.5,
  query_group = "names", ref_group = "cluster", stroke.highlight = 0.5,
  query_param = list(palette = "Paired"), 
  ref_param = list(palette = "Set2", raster = T)
)
fig1c_4
ggsave(filename = "./Fig1c_4.pdf", plot = fig1c_4, width = 8, height = 6, dpi = 300)


