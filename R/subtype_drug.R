
subtype_drug <- function(expr  =DEGs.res$expr,
                         group =DEGs.res$group,
                         calc  =T,
                         col   =c("#CE5C69","#5770A6"),
                         subtype.list = c('CS1','CS3'),
                         gene  ='SLC7A11',
                         drug  ='oxaliplatin'
    )
{
  library(limma)
  library(oncoPredict)
  library(parallel)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ggpubr)
  library(Hmisc)
  library(ggstatsplot)

CTRP2_Expr <- readRDS(file="inst/DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
CTRP2_Res <- readRDS(file = "inst/DataFiles/Training Data/CTRP2_Res.rds")

expr <- as.matrix(expr)


if (calc==T) {
#设置参数
batchCorrect<-"eb"
powerTransformPhenotype<-TRUE
removeLowVaryingGenes<-0.2
removeLowVaringGenesFrom<-"homogenizeData"
minNumSamples=10
selection<- 1
printOutput=TRUE
pcr=FALSE
report_pc=FALSEcc=FALSE
rsq=FALSE
percent=80

#计算并自动输出
calcPhenotype(trainingExprData=CTRP2_Expr,
              trainingPtype=CTRP2_Res,
              testExprData=expr,
              batchCorrect=batchCorrect,
              powerTransformPhenotype=powerTransformPhenotype,
              removeLowVaryingGenes=removeLowVaryingGenes,
              minNumSamples=minNumSamples,
              selection=selection,
              printOutput=printOutput,
              pcr=pcr,
              removeLowVaringGenesFrom=removeLowVaringGenesFrom,
              report_pc=report_pc,
              percent=percent,
              rsq=rsq)
}
#p1---分组柱状图
if (!dir.exists('drug_plot')) {
  dir.create('drug_plot')
}
res <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
tumor_matrix <- data.frame(t(expr))
tumor_matrix <- tumor_matrix[order(tumor_matrix$SLC7A11), ]# 排序
tumor_matrix <- data.frame(t(tumor_matrix))
Data <- data.frame(row.names = colnames(tumor_matrix), group = group)# 创建分组数据框
res <- as.data.frame(res[,'oxaliplatin', drop = FALSE])
rownames(res_oxaliplatin) <- rownames(res) # 保留行名
cam <- cbind(res,Data)

colnames(cam)[2]<- gene
colnames(cam)[1]<- paste0(drug,'(Ic50)')

cam[[colnames(cam)[1]]] <- as.numeric(as.character(cam[[colnames(cam)[1]]]))


p1 <- ggplot(cam, aes(x = !!sym(gene), y = !!sym(colnames(cam)[1]), color = !!sym(gene))) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1)) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 18)) +
  scale_x_discrete(labels = subtype.list) +
  theme_minimal(base_size = 12, base_family = "sans") +
  scale_color_manual(values = col) +
  geom_signif(comparisons = list(subtype.list), test = "t.test", map_signif_level = TRUE)

print(p1)
ggsave(filename = paste0('drug_plot/',gene, '_',drug,'_Bar_chart.pdf'), plot = p1, width = 6, height = 8)

#p2---散点图
res <- read.csv('calcPhenotype_Output/DrugPredictions.csv')
names(res)[1] <- "sample"

fpkm_gene_drug <- tumor_matrix %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(gene) %>%
  rownames_to_column("sample") %>%
  inner_join(res)
fpkm_gene_drug[1:4,1:4]

p2 <- ggscatterstats(
  data = fpkm_gene_drug,
  x = !!sym(gene),
  y = !!sym(drug),
  xlab = gene,
  ylab = drug,
  smooth.line.args = list(linewidth = 1.5, color = 'black', method = "lm", formula = y ~x),
  point.args = list(size = 3, alpha = 0.4, stroke = 0,color ='#A281B1'),
  marginal = TRUE,
  marginal.type = "densigram",
  margins = "both",
  xsidehistogram.args = list(fill = col[[1]], color = "white", na.rm = TRUE), # 分别设置颜色
  ysidehistogram.args = list(fill = col[[2]], color = "white", na.rm = TRUE), # 分别设置颜色
  title = paste0("Relationship between",gene," and ",drug),
  messages = FALSE
)
  print(p2)
  ggsave(filename = paste0('drug_plot/',gene, '_',drug,'_Scatter_plot.pdf'), plot = p2, width = 6, height = 8)

}


