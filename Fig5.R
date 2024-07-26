# 1. FAC proportion
prop_orig <- readxl::read_xlsx("./巨核白血病流式比例.xlsx", sheet = 2)
colnames(prop_orig) <- c("type", "pre","post", "AZA_Vin")
prop_orig[, 2:4] <- apply(prop_orig[, 2:4], 2, function(x) x / sum(x))
prop_melt <- reshape2::melt(prop_orig, ide = "type")
prop_melt$type <- factor(prop_melt$type, levels = c("Lymphocyte", "pre-B", "Plasma", "pMyeloid", 
                                                    "Neutrophil","Monocyte", "MK", "Erythrocyte"))

colors <- RColorBrewer::brewer.pal(n = 9, name = "Paired")
ph_colors_names <- setNames(c(colors[1:4], colors[7:8], colors[c(5:6,9)]), 
                            c("T_IL7R", "T_CD8A", "T_CD247", "B", "Plasma", 
                              "Mono", "Mac", "MK","Ery"))

prop_melt <- prop_melt %>% dplyr::group_by(variable) %>% dplyr::arrange(desc(type)) %>%
  dplyr::mutate(ypos = cumsum(value) - 0.5*value)

p <- ggplot(prop_melt, aes(x = variable, group = type)) +
  geom_bar(aes(y=value, fill=type), stat = "identity", width = 0.8)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(mapping = aes(x = variable, y = ypos, label = paste0(sprintf("%1.1f", 100*value),"%")))+
  scale_fill_manual(values = as.character(ph_colors_names[c(1,4,5,7,3,6,8,9)])) +  
  cowplot::theme_cowplot()+
  scale_x_discrete(limit = rev(c("pre", "post", "AZA_Vin")))+
  labs(x = "group", y = "Flowcytometry proportion")+
  ggtitle(label = "Fig. 3d")+
  coord_flip()

p 

p <- ggplot(data=prop_melt, mapping=aes(x=1, y = value, fill=type))+
  geom_bar(stat="identity", width=1,position='stack', size=1)+
  coord_polar("y", start=0)+
  facet_wrap(~variable)+
  scale_fill_manual(values=c(colors[1:4], colors[7:8], colors[c(5:6,9)]))+
  geom_text_repel(stat="identity",
                  aes(y=value, 
                      label = paste0(round(value*100, 2), "%")), size=3, 
                  position=position_stack(vjust = 0.5)) +
  theme_minimal()+ 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 12))


ggsave(filename = "./Fig3d_2.pdf", plot = p, width = 15, height = 3.5, dpi = 300)

# 2. Dynamics
df_ratio <- data.frame(group = c("pre","post","post_A"), 
                       rat_MK = c(20.83,6.15,12.3),
                       rat_plasma = c(1.47, 15.3, 0))
df_ratio <- reshape2::melt(df_ratio, id = "group")
df_ratio$group <- factor(df_ratio$group, levels = c("pre","post","post_A"))

ggplot(df_ratio, aes(x = group, y = value, col = variable)) + 
  geom_point(size = 5)+geom_line(aes(group = variable), size = 1.5)
ggsave(filename = "./Fig4.pdf", plot = last_plot(), width = 5, height = 3.4, dpi = 300)
