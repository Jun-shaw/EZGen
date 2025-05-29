#' @export
#' @import clusterProfiler
#' @import msigdbr
#' @import GSVA
#' @import ggplot2
#' @import GSEABase
#' @import tidyverse
#' @import ggthemes
#' @import ggprism
#' @import limma
#' @import stringr
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @references
#' Yu G, Wang L, Han Y, He Q (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5):284-287.
#' @examples # There is no example and please refer to vignette.

subtype_GSVA <- function(DEGs.res       = DEGs.res,
                         data.base      = 'GO',#c('GO','KEGG')
                         subtype.list   =c('CS1','CS2')#"exp/ctrl"
                         )

{
  library(clusterProfiler)
  library(msigdbr)
  library(GSVA)
  library(ggplot2)
  library(GSEABase)
  library(tidyverse)
  library(ggthemes)
  library(ggprism)
  library(limma)
  library(stringr)

  if (!dir.exists('EA_plot')) {
    dir.create('EA_plot')
  }
##KEGG
KEGG_df_all <-  msigdbr(species = "Homo sapiens", # Homo sapiens or Mus musculus
                        category = "C2",
                        subcategory = "CP:KEGG")
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) ##按照gs_name给gene_symbol分组

##GO
GO_df_all <- msigdbr(species = "Homo sapiens",
                     category = "C5")
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
go_list <- split(GO_df$gene_symbol, GO_df$gs_name) ##按照gs_name给gene_symbol分组

####  GSVA  ####
#GSVA算法需要处理logCPM, logRPKM,logTPM数据或counts数据的矩阵####
#dat <- as.matrix(counts)
#dat <- as.matrix(log2(edgeR::cpm(counts))+1)
dat <- as.matrix(log2(DEGs.res$expr+1))

if (data.base=='GO') {geneset <- go_list}else{geneset <- kegg_list}

gsva_mat <- gsva(expr=dat,
                 gset.idx.list=geneset,
                 kcdf="Gaussian",
                 verbose=T,
                 parallel.sz = parallel::detectCores())

write.csv(gsva_mat,paste0('EA_plot/',"gsva_",data.base,"_matrix.csv"))

group_list <- DEGs.res[["clust.res"]]
design <- model.matrix(~0+factor(group_list$clust))
colnames(design) <- levels(factor(group_list$clust))
rownames(design) <- colnames(gsva_mat)
contrast.matrix <- makeContrasts(contrasts=paste0(subtype.list[1],'-',subtype.list[2]),  #"exp/ctrl"
                                 levels = design)

fit1 <- lmFit(gsva_mat,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正

summary(decideTests(efit,lfc=0, p.value=0.5)) #统计查看差异结果
tempOutput <- topTable(efit, coef=paste0(subtype.list[1],'-',subtype.list[2]), n=Inf)
degs <- na.omit(tempOutput)
write.csv(degs,paste0('EA_plot/',"gsva_",data.base,"_result.csv"))

p_cutoff=0.5
Diff <- rbind(subset(degs,logFC>0)[1:10,], subset(degs,logFC<0)[1:10,]) #选择上下调前20通路
dat_plot <- data.frame(id  = row.names(Diff),
                       p   = Diff$P.Value,
                       lgfc= Diff$logFC)
dat_plot$group <- ifelse(dat_plot$lgfc>0 ,1,-1)    # 将上调设为组1，下调设为组-1
dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$group # 将上调-log10p设置为正，下调-log10p设置为负
dat_plot$id <- str_replace(dat_plot$id, "KEGG_","");dat_plot$id[1:10]
dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= p_cutoff,
                                    ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
                             levels=c('Up','Down','Not'))

dat_plot <- dat_plot %>% arrange(lg_p)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

## 设置不同标签数量
low1 <- dat_plot %>% filter(lg_p < log10(p_cutoff)) %>% nrow()
low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
high0 <- dat_plot %>% filter(lg_p < -log10(p_cutoff)) %>% nrow()
high1 <- nrow(dat_plot)

p1 <- ggplot(data = dat_plot,aes(x = id, y = lg_p,
                                 fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#CE5C69','Not'='#cccccc','Down'='#5770A6')) +
  geom_hline(yintercept = c(-log10(p_cutoff),log10(p_cutoff)),color = 'white',size = 0.5,lty='dashed') +
  xlab('') +
  ylab('-log10(P value) of GSVA score') +
  guides(fill="none")+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black',size = 2) + #黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey',size = 2) + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey',size = 2) + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black',size = 2) # 黑色标签
ggsave("EA_plot/GSVA_barplot_pvalue.pdf", p1, width = 10, height = 15)
p1
}
