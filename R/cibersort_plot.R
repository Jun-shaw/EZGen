


cibersort_plot <- function(TME.res.list=TME.res.list,
                           clust.res.list=clust.res.list,
                           palette.col=NA,
                           heatmap.col=c("#5770A6", 'white',"#CE5C69"))
{


library(dplyr)
library(tidyr)
library(forcats)
library(ggpubr)
library(stringr)
library(ComplexHeatmap)
library(tidyHeatmap)

  if (!dir.exists('TME_plot')) {
    dir.create('TME_plot')
  }
cibersort_long <- TME.res.list$CIBERSORT %>%
  select(`P-value.CIBERSORT`,Correlation.CIBERSORT, RMSE.CIBERSORT,ID,everything()) %>%
  pivot_longer(- c(1:4),names_to = "cell_type",values_to = "fraction") %>%
  dplyr::mutate(cell_type = gsub("_CIBERSORT","",cell_type),
                cell_type = gsub("_"," ",cell_type))
if (is.na(palette.col)) {
  palette.col <- c("#5770A6", "#CE5C69", "#A281B1", "#678A74", "#F7C5A8", "#FFEADB", "#FFBABA", "#7A4579", "#D56073", "#EC9E69", "#FADCAA", "#D79ABC", "#BAABDA", "#9FDFCD", "#DCFFCC", "#63B7AF", "#347474", "#35495E", "#EE8572", "#6C5B7C", "#C06C84", "#F67280", "#F8B595", "#F5E8C7", "#818FB4", "#435585", "#363062", "#967E76", "#D7C0AE", "#EEE3CB", "#9BABB8")
}

#条状图
pdf(file=paste0('TME_plot/',"Bar_chart.pdf"),width = 12,height = 8)
p1 <- cibersort_long %>%
  ggplot(aes(ID,fraction))+
  geom_bar(stat = "identity",position = "stack",aes(fill=cell_type))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette.col,name=NULL)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")
print(p1)
dev.off()
#总箱式图
pdf(file=paste0('TME_plot/',"cell_box_plot.pdf"),width = 12,height = 8)
p2 <- ggplot(cibersort_long,aes(fct_reorder(cell_type, fraction),fraction,fill = cell_type)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2,aes(color=cell_type))+
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = palette.col)
print(p2)
dev.off()
#分组箱式图

# 确保 clust.res.list$iClusterBayes$clust.res 的 samID 列是行名
rownames(clust.res.list$iClusterBayes$clust.res) <- clust.res.list$iClusterBayes$clust.res$samID

# 使用 match 函数将 clust 列添加到 cibersort_long
cibersort_long$sub_type <- clust.res.list$iClusterBayes$clust.res[cibersort_long$ID, "clust"]
cibersort_long$sub_type <- paste0("CS", cibersort_long$sub_type)
cibersort_long$cell_type <- gsub("\\.CIBERSORT", "", cibersort_long$cell_type)

pdf(file=paste0('TME_plot/',"subtype_box_plot.pdf"),width = 12,height = 8)
p3 <- ggplot(cibersort_long, aes(fct_reorder(cell_type, fraction), fraction)) +
  geom_boxplot(aes(fill = sub_type), outlier.shape = 21, color = "black") +
  scale_fill_manual(values = palette.col[c(1, 2)]) +
  theme_bw() +
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black", size = 12)) +
  stat_compare_means(aes(group = sub_type, label = ..p.signif..),
                     method = "kruskal.test", label.y = 0.4)
print(p3)
dev.off()
#heatmap
pdf(file=paste0('TME_plot/',"heatmap_plot.pdf"),width = 12,height = 8)
p4 <- heatmap(.data = cibersort_long,
              .row = cell_type,
              .column = ID,
              .value = fraction,
              scale = "column",
              palette_value = circlize::colorRamp2(seq(-2, 2, length.out = 3), heatmap.col),
              show_column_names=F,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 7),
              column_title_gp = gpar(fontsize = 7),
              row_title_gp = gpar(fontsize = 7)
) %>% add_tile(sub_type)
print(p4)
dev.off()
return(list(p1=p1,p2=p2,p3=p3,p4=p4))
}
