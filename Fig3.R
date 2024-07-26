# 1.GSVA for hallmark genesets
library(GSEABase)
hm.sets <- getGmt("./gmt_files/h.all.v2022.1.Hs.symbols.gmt")
length(hm.sets)
names(hm.sets)
annot_gsve <- seu_pre@meta.data
table(annot_gsve$names)
exprs_gsva <- GetAssayData(seu_pre, slot= "data")
exprs_gsva <- exprs_gsva[Matrix::rowSums(exprs_gsva > 0) > 5, ]
dim(exprs_gsva)

library(GSVA)
datasets.hm <- gsva(as.matrix(exprs_gsva), hm.sets, method="gsva", parallel.sz = 10, verbose = TRUE)

datasets.hm1 <- MinMax(datasets.hm, min = -0.5, max = 0.5)
library(pheatmap)
ph_colors$names <- setNames(c(louvain_colors[1:3], louvain_colors[c(5:7, 10, 8:9)]), levels(seu_pre$names))

### determine the significant geneset
datasets.hm2 <- datasets.hm[c(4,8,28,35,2), ]
datasets.hm2 <- as.data.frame(t(datasets.hm2))
datasets.hm2$cluster <- annot_gsve$names

datasets.hm2 <- reshape2::melt(datasets.hm2, id = "cluster")
colnames(datasets.hm2) <- c("cluster", "Geneset", "expression")

fig2a <- ggviolin(datasets.hm2, x="cluster", y = "expression", fill ="cluster", trim = T,
                  add = "boxplot", add.params = list(fill = "grey90", width = 0.1)) + 
  facet_wrap(facets = "Geneset",scales = "free_y", ncol = 1)+
  scale_fill_manual(values = ph_colors$names) +
  labs(y = "Relative signature score") 
fig2a
ggsave(filename = "./Fig2a_1.pdf", plot = fig2a, width = 4, height = 12, dpi = 300)


# 2 compare with RPMM
## 2.1 load RRMM data 
dir <- "./rawdata/2021_NC/files/"

files <- list.files(path = dir, pattern = "*.csv.gz", recursive = T)

store_csv = paste(dir,"combined.csv")
library(readr)
library(data.table)
df <- fread(file = file.path(dir,files[1]),encoding = 'UTF-8', data.table = T)
rs <- matrix(nrow = nrow(df), ncol = 0)
rownames(exprs) <- df$gene

for(i in 1:length(files)){
  df = fread(file = file.path(dir,files[i]),encoding = 'UTF-8')
  df = df[, -1]
  df <- as.matrix(df)              
  df <- as(df, 'dgCMatrix')
  exprs <- cbind(exprs, df)
}


meta <- read.csv("./rawdata/2021_NC/GSE161801_K43R_metadata_table.csv.gz", row.names = 1)

seu_nc <- CreateSeuratObject(counts = exprs, meta.data = meta, project = "2021_NC")

## 2.2 singleR prediction
library(SingleR)
library(BiocParallel)
library(future)
require(ggplotify)

seu_WBM <- subset(seu_nc, sorting == "WBM")
scRNA_pred <- SingleR(test = seu_pre@assays$RNA@data,
                      ref = seu_WBM@assays$RNA@data,
                      labels = seu_WBM$celltype_1,
                      de.method="wilcox",
                      BPPARAM=MulticoreParam(20))
ph_colors_SR <- list(Clusters = ph_colors$names, 
                     Labels = setNames(hej_colors[1:34], unique(seu_WBM$celltype_1))) 
p1 <- plotScoreHeatmap(scRNA_pred, clusters = seu_pre$names, show_colnames = F, cluster_cols = F, 
                       annotation_colors = ph_colors_SR)
p1 = as.ggplot(p1)
p1

tab <- table(Assigned=scRNA_pred$pruned.labels, Cluster=seu_pre$names)

p2 <- pheatmap(log2(tab+1), color=colorRampPalette(c("white", "red"))(101), cluster_cols = F)
fig2b = as.ggplot(p2)

# (p1 | p2) + plot_layout(widths = c(2,1))
ggsave(filename = "./Fig2b.pdf", plot = fig2b, width = 5, height = 8, dpi = 300)


## 2.3 Plasma integration
plasma_pre <- subset(seu_pre, names == "Plasma")

Myeloma <- subset(seu_nc, timepoint %in% c("pre"))
Myeloma <- subset(Myeloma, celltype_1 %in% c("Myeloma", "PC"))
Myeloma <- subset(Myeloma, sorting %in% c("CD138_pos"))

Myeloma$exp <- "Myeloma"

plasma_nBM$exp <- "normal"
plasma_nBM$celltype_1 <- "plasma"
plasma_nBM$PID_sample_new <- "normal"
plasma_pre$exp <- "pre"
plasma_pre$tumor_clone <- "pre"
plasma_pre$celltype_1 <- "plasmaAMKL"
plasma_pre$PID_sample_new <- "pre"

seu_merge_plsma <- merge(plasma_pre, list(Myeloma))

seu_merge_plsma <- NormalizeData(seu_merge_plsma, scale.factor = 1e4)

seu_merge_plsma <- FindVariableFeatures(seu_merge_plsma)
seu_merge_plsma <- RunFastMNN(object.list = SplitObject(seu_merge_plsma, split.by = "exp"))
seu_merge_plsma <- RunUMAP(seu_merge_plsma, reduction = "mnn", dims = 1:10, n.components = 2, 
                           n.neighbors = 30, return.model = T)


DimPlot(seu_merge_plsma, group.by = c("celltype_1"), ncol = 1, reduction = "umap", cols = mycolors, label = F)
fig2c
ggsave(filename = "./Fig2c.pdf", plot = fig2c, width = 6, height = 5, dpi = 300)


## 2.4 +1q evaluation
genes_1q <- readxl::read_xlsx("./41467_2021_26951_MOESM4_ESM.xlsx", sheet = 1)
seu_merge_plsma <- AddModuleScore(seu_merge_plsma, features = list(genes_1q$Gene), name = "score_1q")

fig2c_2 <- ggviolin(ggdata, x="celltype_1", y = "score_1q1", fill ="celltype_1", trim = T,
                    add = "boxplot", add.params = list(fill = "grey90", width = 0.1))+
  # scale_fill_manual(values = c("#DC8CC3","#E31A1C")) +
  labs(y = "+1q score") + 
  stat_compare_means(comparisons = list(c("plasmaAMKL", "Myeloma"), 
                                        c("PC", "plasmaAMKL"),
                                        c("Myeloma", "PC")), label = "p.format")
fig2c_2
ggsave(filename = "./Fig2c_3.pdf", plot = fig2c_2, width = 5, height = 4, dpi = 300)


## 2.5 CNV analysis
ibrary(infercnv)
library(AnnoProbe)
source("./inferCNV_heatmap.R")
seu_merge_plsma@meta.data$celltype_1 <- factor(seu_merge_plsma@meta.data$celltype_1, levels = c("PC", "Myeloma", "plasma"))
groupinfo <- seu_merge_plsma@meta.data[, "celltype_1", drop = F]
table(groupinfo$celltype_1)
DefaultAssay(seu_merge_plsma) <- "RNA"
count.matrix <- GetAssayData(seu_merge_plsma, slot = "count")
count.matrix <- count.matrix[rowSums(count.matrix) > 10, ]
geneInfor <- annoGene(rownames(count.matrix), "SYMBOL",'human')
colnames(geneInfor)
geneInfor <- geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
rownames(geneInfor) <- geneInfor$SYMBOL
geneInfor <- geneInfor[, -1]
head(geneInfor)## 这里可以去除性染色体# 也可以把染色体排序方式改变
table(geneInfor$chr)
geneInfor$chr <- factor(geneInfor$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
geneInfor <- geneInfor[order(geneInfor$chr), ]
geneInfor <- na.omit(geneInfor)


infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=count.matrix[rownames(geneInfor),],
                                     gene_order_file=geneInfor,
                                     annotations_file=groupinfo,
                                     ref_group_names=c("PC"),
                                     delim = "\t",
                                     max_cells_per_group = NULL,
                                     min_max_counts_per_cell = c(100, +Inf))
range(colSums(infercnv_obj@expr.data))

infercnv_obj2 <- infercnv::run(infercnv_obj, cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir= "./infercnv_output_new", # dir is auto-created for storing outputs
                               cluster_by_groups=F,   # cluster  
                               HMM = FALSE,
                               denoise=TRUE,
                               output_format = "pdf",
                               hclust_method="ward.D2", plot_steps=T)
saveRDS(infercnv_obj2, file = "./infercnv.Rds")





