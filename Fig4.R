# 1. CD38 expression
VlnPlot(seu_pre, group.by = "names", features = c("CD38"))
Fig3a <- ExpStatPlot(
  srt = seu_pre, group.by = "names", 
  features = c("CD38"))

Fig3a <- FeatureStatPlot(
  srt = seu, group.by = "clusters", 
  add_box = T,
  split.by = "orig.ident",
  palcolor = c("#F6B352","#9055A2"),
  stat.by = c("CD38"))

ggsave(filename = "./Fig3a.pdf", plot = Fig3a, width = 6, height = 3, dpi = 300)

# 2 post(anti-CD38) data integration
seu <- NormalizeData(seu, scale.factor = 1e4)
seu <- FindVariableFeatures(seu, nfeatures = 2000)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
seu <- ScaleData(seu, vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

seu <- RunPCA(seu, npcs = 50)

seu <- RunUMAP(seu, dims = 1:40)

seu <- FindNeighbors(seu, dims = 1:30, k.param = 20)
seu <- FindClusters(seu, resolution = 0.25)

DimPlot(seu, reduction = "umap", group.by = c( "orig.ident","seurat_clusters"), label = T, cols = louvain_colors)

Idents(seu) <- "seurat_clusters"
seu <- RenameIdents(seu, '0' = "T_CD8A", '1' = "MK", '2' = "Plasma", '3' = "T_IL7R", '4' = "Ery", '5' = "T_CD247", 
                    '6' = "Mono_CD14", '7' = "Ery", '8' = "B",'9' = "Mono_CD16", '10' = "Pre_B")
seu$clusters <- Idents(seu)
seu$clusters <- factor(seu$clusters, levels = c("T_IL7R", "T_CD8A", "T_CD247", "Pre_B","B", "Plasma", 
                                                "Mono_CD14", "Mono_CD16","MK","Ery"))
seu$orig.ident <- factor(seu$orig.ident, levels = c("pre","post"))
table(seu$orig.ident, seu$clusters)
saveRDS(seu, file = "./rds_files/seu_merge.Rds")



#3. MK comparison
## 3.1 signature gene
seu_MK <- subset(seu, clusters  == "MK")
Idents(seu_MK) <- "orig.ident"
degs_mk <- FindAllMarkers(seu_MK, logfc.threshold = 0.5, min.pct = 0.25, only.pos = T)
degs_mk_sig <- subset(degs_mk, p_val_adj < 0.01)
table(degs_mk_sig$cluster)

seu_MK <- AnnotateFeatures(seu_MK, species = "Homo_sapiens", db = c("TF", "SP"))

ht <- FeatureHeatmap(
  srt = seu_MK, group.by = "orig.ident", features = degs_mk_sig$gene, feature_split = degs_mk_sig$cluster,
  anno_keys = F, anno_features = F, anno_terms = F, nlabel = 30,
  feature_annotation = c("TF", "SP"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 6, width = 5
)

fig3f_1 <- as.ggplot(ht$plot)
fig3f_1

# fig3f_2 <- FeatureDimPlot(
#   srt = seu_MK, features = c("ANXA2","APOE",  "PKM", "MYLK"), ncol = 1,
#   reduction = "umap", theme_use = "theme_blank"
# )
# fig3f <- (fig3f_1 | fig3f_2) + plot_layout(widths = c(3, 1))
fig3f <- FeatureStatPlot(
  srt = seu_MK, group.by = "orig.ident",  
  comparisons = list(c("pre", "post")), add_box = TRUE,  palcolor = c("#F6B352","#9055A2"), 
  stat.by = c("ANXA2","APOE",  "PKM", "MYLK"), ncol = 2,sig_label = "p.format")
ggsave(filename = "./Fig3a.pdf", plot = Fig3a, width = 6, height = 3, dpi = 300)
(fig3f_1 | fig3f) + plot_layout(widths = c(6,1))


seu_MK <- AddModuleScore(seu_MK, features = list(MK_diff_genes,Blood_coagulation_genes, Cell_Prolif_genes), 
                         name = c("MK_diff",  "Blood_coagulation", "Cell_Prolif"))

MK_score <- FetchData(seu_MK, vars = c("orig.ident", "MK_diff1",  "Blood_coagulation2", "Cell_Prolif3"))
MK_score <- reshape2::melt(MK_score, id = "orig.ident")
colnames(MK_score) <- c("stage", "Geneset", "expression")
MK_score$stage <- factor(MK_score$stage, levels = c("pre", "post"))

fig3f_2 <- ggviolin(MK_score, x="stage", y = "expression", fill ="stage", trim = T,
                    add = "boxplot", add.params = list(fill = "grey90", width = 0.1)) + 
  facet_wrap(facets = "Geneset",scales = "free_y", ncol = 3)+
  scale_fill_manual(values = c("#F6B352","#9055A2")) +
  labs(y = "Relative signature score") + 
  stat_compare_means(comparisons = list(c("pre", "post")), label = "p.format")
fig3f_2
ggsave(filename = "./Fig3f_2.pdf", plot = last_plot(), width = 6, height = 3, dpi = 300)


## 3.2 GSVA
library(GSEABase)
hm.sets <- getGmt("~/common_data/gmt_files/h.all.v7.5.1.symbols.gmt")
length(hm.sets)
names(hm.sets)
annot_gsve <- seu_MK@meta.data
table(annot_gsve$orig.ident)
exprs_gsva <- GetAssayData(seu_MK, slot= "data")
exprs_gsva <- exprs_gsva[Matrix::rowSums(exprs_gsva > 0) > 5, ]
dim(exprs_gsva)

# Run GSVA 
library(GSVA)
datasets.hm <- gsva(as.matrix(exprs_gsva), hm.sets, method="gsva", parallel.sz = 20, verbose = TRUE)
datasets.hm <- datasets.hm[c(4,8,28,35,2), ]
datasets.hm2 <- as.data.frame(t(datasets.hm))

datasets.hm2$cluster <- annot_gsve$orig.ident

datasets.hm2 <- reshape2::melt(datasets.hm2, id = "cluster")
colnames(datasets.hm2) <- c("cluster", "Geneset", "expression")

fig3g_1<- ggviolin(datasets.hm2, x="cluster", y = "expression", fill ="cluster", trim = T,
                   add = "boxplot", add.params = list(fill = "grey90", width = 0.1)) +
  facet_wrap(facets = "Geneset",scales = "free_y", ncol = 2)+
  scale_fill_manual(values = ph_colors$orig.ident) +
  labs(y = "Relative signature score")  +
  stat_compare_means(comparisons = list(c("pre", "post")), label = "p.format")

fig3g_1 <- ggviolin(datasets.hm2, x="cluster", y = "expression", fill ="cluster", trim = T,
                    add = "boxplot", add.params = list(fill = "grey90", width = 0.1)) + 
  facet_wrap(facets = "Geneset",scales = "free_y", ncol = 5)+
  scale_fill_manual(values = c("#F6B352","#9055A2")) +
  labs(y = "Relative signature score") + 
  stat_compare_means(comparisons = list(c("pre", "post")), label = "p.format")

fig3g_1
ggsave(filename = "./Fig3g_1.pdf", plot = fig3g_1, width = 10, height = 3, dpi = 300)

## 3.3 GSEA pathway
seu_MK <- RunDEtest(srt = seu_MK, group_by = "orig.ident", fc.threshold = 1, only.pos = FALSE,
                    BPPARAM=BiocParallel::MulticoreParam(20))
VolcanoPlot(srt = seu_MK, group_by = "orig.ident")
seu_MK <- RunGSEA(
  srt = seu_MK, group_by = "orig.ident", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.01", BPPARAM=BiocParallel::MulticoreParam(20)
)

fig3g_2 <- GSEAPlot(srt = seu_MK, group_by = "orig.ident", group_use = "post", topTerm = 4)

fig3g_3 <- GSEAPlot(srt = seu_MK, group_by = "orig.ident", group_use = "post", geneSetID = "GO:0006323")
fig3g_3


# 4. Myeloma subset
## 4.1 signautr genes
seu_plasma <- subset(seu, clusters  == "Plasma")

Idents(seu_plasma) <- "orig.ident"
degs_mye <- FindAllMarkers(seu_plasma, min.pct = 0.25, logfc.threshold = 0.5, only.pos = T)

degs_mye_sig <- subset(degs_mye, p_val_adj < 0.01)
table(degs_mye_sig$cluster)

seu_plasma <- AnnotateFeatures(seu_plasma, species = "Homo_sapiens", db = c("TF", "SP"))

ht <- FeatureHeatmap(
  srt = seu_plasma, group.by = "orig.ident", features = degs_mye_sig$gene, feature_split = degs_mye_sig$cluster,
  anno_keys = F, anno_features = F, anno_terms = F, nlabel = 30,
  feature_annotation = c("TF", "SP"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 6, width = 5
)

fig3g <- as.ggplot(ht$plot)
fig3g


## 4.2 GSVA
library(GSEABase)
hm.sets <- getGmt("~/common_data/gmt_files/h.all.v7.5.1.symbols.gmt")
length(hm.sets)
names(hm.sets)
annot_gsve <- seu_plasma@meta.data
table(annot_gsve$orig.ident)
exprs_gsva <- GetAssayData(seu_plasma, slot= "data")
exprs_gsva <- exprs_gsva[Matrix::rowSums(exprs_gsva > 0) > 5, ]
dim(exprs_gsva)

# Run GSVA 
library(GSVA)
datasets.hm <- gsva(as.matrix(exprs_gsva), hm.sets, method="gsva", parallel.sz = 20, verbose = TRUE)
dim(datasets.hm)
range(datasets.hm)
write.csv(datasets.hm, file = "Hallmark_GSVA_scores_myeloma.csv", quote = F)
datasets.hm <- read.csv("Hallmark_GSVA_scores_myeloma.csv", row.names = 1)

rownames(datasets.hm)
datasets.hm <- datasets.hm[c(4,8,28,35,2), ]
datasets.hm2 <- as.data.frame(t(datasets.hm))

datasets.hm2$cluster <- annot_gsve$orig.ident

datasets.hm2 <- reshape2::melt(datasets.hm2, id = "cluster")
colnames(datasets.hm2) <- c("cluster", "Geneset", "expression")

fig3h_1<- ggviolin(datasets.hm2, x="cluster", y = "expression", fill ="cluster", trim = T,
                   add = "boxplot", add.params = list(fill = "grey90", width = 0.1)) +
  facet_wrap(facets = "Geneset",scales = "free_y", ncol = 5)+
  scale_fill_manual(values = ph_colors$orig.ident) +
  labs(y = "Relative signature score")  +
  stat_compare_means(comparisons = list(c("pre", "post")), label = "p.format")


## 4.3 GSEA pathway
seu_plasma <- RunDEtest(srt = seu_plasma, group_by = "orig.ident", fc.threshold = 1, only.pos = FALSE,
                        BPPARAM=BiocParallel::MulticoreParam(20))
VolcanoPlot(srt = seu_plasma, group_by = "orig.ident")
seu_plasma <- RunGSEA(
  srt = seu_plasma, group_by = "orig.ident", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.01", BPPARAM=BiocParallel::MulticoreParam(20), GO_simplify = T,
)

fig3i_1 <- GSEAPlot(srt = seu_plasma, group_by = "orig.ident", group_use = "post", topTerm = 4)

fig3i_2 <- GSEAPlot(srt = seu_plasma, group_by = "orig.ident", group_use = "post", geneSetID = "GO:0008152")

ggsave(filename = "./Fig3i.pdf", plot = fig3i, width = 14, height = 4, dpi = 300)




## 5. cell interaction
library(celltalker)
## Split dataset
seu_cci <- subset(seu, clusters != "Pre_B")
seu_cci$cell_type <- as.character(seu_cci$clusters)
seu_cci$cell_type <- gsub("T_", "", seu_cci$cell_type)
ser_split <- SplitObject(seu_cci, split="orig.ident")

## 5.1 pre
pre_interactions <- celltalk(input_object=ser_split[["pre"]],
                             metadata_grouping="cell_type",
                             ligand_receptor_pairs=ramilowski_pairs,
                             number_cells_required=20,
                             min_expression=50,
                             max_expression=20000,
                             scramble_times=10)

## Check out the interactions - pre
pre_interactions$p_val_adj <- p.adjust(pre_interactions$p_val,method="fdr")
pre_interactions_filter <- subset(pre_interactions, p_val_adj<0.01)

## 5.2 post
post_interactions <- celltalk(input_object=ser_split[["post"]],
                              metadata_grouping="cell_type",
                              ligand_receptor_pairs=ramilowski_pairs,
                              number_cells_required=20,
                              min_expression=50,
                              max_expression=20000,
                              scramble_times=10)
## Check out the interactions - post
post_interactions$p_val_adj <- p.adjust(post_interactions$p_val,method="fdr")
post_interactions_filter <- subset(post_interactions, p_val_adj<0.01)

write.csv(pre_interactions_filter, file = "./pre_dataset_CCI.csv")
write.csv(post_interactions_filter, file = "./post_dataset_CCI.csv")


# 5.3 Identify the top 10 interactions for cluster MK or Plasma
top_stats_pre <- pre_interactions_filter %>%
  filter(cell_type1 %in% c("MK", "Plasma") | cell_type2 %in% c("MK", "Plasma") ) %>%
  group_by(cell_type1) %>%
  top_n(10,interact_ratio) %>%
  ungroup()
# Assign colors to cell types
all_cell_types <- unique(seu[["clusters"]][,1])

# Suppress messages to silence the circlize functions
circos_plot(ligand_receptor_frame=top_stats_pre,
            cell_group_colors=ph_colors$clusters,
            ligand_color="lightblue",
            receptor_color="orange",
            cex_outer=0.5,
            cex_inner=0.35)

top_stats_post <- post_interactions_filter %>%
  filter(cell_type1 %in% c("MK", "Plasma") | cell_type2 %in% c("MK", "Plasma") ) %>%
  group_by(cell_type1) %>%
  top_n(10,interact_ratio) %>%
  ungroup()
# Suppress messages to silence the circlize functions
suppressMessages(
  circos_plot(ligand_receptor_frame=top_stats_post,
              cell_group_colors=ph_colors$clusters,
              ligand_color="lightblue",
              receptor_color="orange",
              cex_outer=0.5,
              cex_inner=0.35)
)





## 5.4 compare

seu_cci <- subset(seu, clusters %in% c("MK", "Plasma", "Mono_CD14", "Mono_CD16"))
seu_cci$cell_type <- as.character(seu_cci$clusters)
seu_cci$group <- as.character(seu_cci$orig.ident)
seu_cci$cell_type <- gsub("Mono_", "", seu_cci$cell_type)
table(seu_cci$cell_type, seu_cci$group)

# prepare l_r data
human_lr <- read.table("~/common_data/Human-2022-Dimitrov-LR-pairs.txt", sep = ",", header = T)
human_lr <- human_lr[, 3:4]
tmp <- colsplit(human_lr$target_genesymbol, "_", c("r1", "r2"))
human_lr$target_genesymbol <- tmp$r1
human_lr$pair <- paste0(human_lr$source_genesymbol, "_", human_lr$target_genesymbol)
colnames(human_lr) <- c("ligand","receptor","pair")
human_lr <- human_lr[!duplicated(human_lr$pair), ]

## calculate CCI
interactions_list <- list()
for (i in c("pre", "post")) {
  seu_tmp <- subset(seu_cci, group == i)
  interactions <- celltalk(input_object=seu_tmp,
                           metadata_grouping="cell_type",
                           ligand_receptor_pairs= human_lr,
                           number_cells_required=10,
                           min_expression=10,
                           max_expression=20000,
                           scramble_times=10)
  #interactions$p_val_adj <- p.adjust(interactions$p_val,method="fdr")
  #interactions <- subset(interactions, p_val_adj < 0.1)
  interactions_list[[i]] <- interactions
}

interactions_merge <- data.frame()
for (i in c("pre", "post")) {
  interactions_df <- interactions_list[[i]]
  interactions_df$category <- i
  interactions_merge <- rbind(interactions_merge, interactions_df)
}

View(interactions_merge)
interactions_merge <- subset(interactions_merge, p_val < 0.05)
interactions_merge <- interactions_merge[interactions_merge$interaction_pairs %in% c("CD14_MK", "MK_CD14",
                                                                                     "CD16_MK", "MK_CD16",
                                                                                     "CD14_Plasma", "Plasma_CD14",
                                                                                     "CD16_Plasma", "Plasma_CD16"), ]
ggdata <- interactions_merge[, c("interaction_pairs", "category")]
ggdata$interaction_pairs <- gsub("CD14_MK|MK_CD14", "CD14-MK",ggdata$interaction_pairs)
ggdata$interaction_pairs <- gsub("CD16_MK|MK_CD16", "CD16-MK",ggdata$interaction_pairs)
ggdata$interaction_pairs <- gsub("CD14_Plasma|Plasma_CD14", "CD14-Plasma",ggdata$interaction_pairs)
ggdata$interaction_pairs <- gsub("CD16_Plasma|Plasma_CD16", "CD16-Plasma",ggdata$interaction_pairs)
ggdata <- as.data.frame(table(ggdata$category, ggdata$interaction_pairs))
colnames(ggdata) <- c("group", "pair", "Freq")
ggdata$group <- factor(ggdata$group, levels = c("pre","post"))


### L-R counts 
ggplot(ggdata,aes(x = pair,y = Freq))+
  geom_bar(stat = 'identity', aes(fill = group),position = position_dodge(0.9)) +
  theme_bw() + scale_fill_manual(values = c("#F6B352","#9055A2")) +
  labs(y = "No. of L-R pairs between clusters", x = "Cell Pair")
ggsave(filename = "./Fig3i_1.pdf", plot = last_plot(), width = 4, height = 3, dpi = 300)

write.csv(interactions_merge, file = "CCI_Mono_MK_Plasma.csv")

# Create circos plots
## Identify the top 10 interactions for Ery or Mac
top_list <- list()
for (i in c("pre", "post")) {
  interactions <- subset(interactions_merge, category == i)
  top_stats <- interactions %>%
    top_n(100,interact_ratio) %>%
    ungroup()
  top_list[[i]] <- top_stats
}

View(top_list$pre)


# Suppress messages to silence the circlize functions
circos_plot(ligand_receptor_frame=top_list$post[1:30,],
            cell_group_colors=ph_colors$cluster,
            ligand_color="lightblue",
            receptor_color="orange",
            cex_outer=0.8,
            cex_inner=0.3)


top_stats_merge <- data.frame()
for (i in c("pre", "post")) {
  top_stats_df <- top_list[[i]]
  top_stats_df$category <- i
  top_stats_merge <- rbind(top_stats_merge, top_stats_df)
}


### Visulization for significant L-R pairs
selected <- c( "TNF_FFAR2","S100A8_TLR4", "GZMB_CHRM3","TNC_CNTN1","BDNF_DRD4","PAM_DPP4", 
               "ZG16B_TLR4", "ZG16B_TLR2",
               "PRG4_TLR2", "TNC_ITGA9",  "SCT_RAMP2",  "SCT_TSHR","TNC_SDC1",
               "CXCL8_SDC1","CXCL2_DPP4", "CXCL2_CXCR2", "BMP6_AMHR2", "APOE_LSR","WNT5A_PTPRK","EFNB2_EPHB1")
ggdata <- top_stats_merge
ggdata <- top_stats_merge[top_stats_merge$interaction %in% selected, ]


library(grafify)
ggdata$category <- factor(ggdata$category, levels = c("pre","post"))
ggplot(ggdata, aes(x = interaction_pairs, y = interaction)) + 
  geom_point(aes(size = interact_ratio, color = p_val)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,colour = "black"))+
  facet_wrap(~category) + scale_y_discrete(limits = rev(selected))+
  scale_color_gradientn(colours = rev(colorRampPalette(c("grey60","orange","red","brown"))(11)))
ggsave(filename = "./Fig3i_2.pdf", plot = last_plot(), width = 6, height = 5, dpi = 300)



