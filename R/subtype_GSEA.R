#' @export
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import stringr
#' @import ggplot2
#' @import enrichplot
#' @import pathview
#' @import DOSE
#' @import topGO
#' @references
#' Yu G, Wang L, Han Y, He Q (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5):284-287.
#' @examples # There is no example and please refer to vignette.

subtype_GSEA <- function(DEGs.res       = DEGs.res,
                         data.base      = 'GO',#c('GO','KEGG')
                         geneset.KEGG   = NA,
                         geneset.GO     = NA,
                         Pvalue         = 0.5,
                         col            = '#CE5C69')

{
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  library(ggplot2)
  library(enrichplot)
  library(pathview)
  library(DOSE)
  library(topGO)
  if (!dir.exists('EA_plot')) {
    dir.create('EA_plot')
  }

  DEGs.dt <- DEGs.res$DEGs
  colnames(DEGs.dt)[c(1,3)] <- c("SYMBOL","logFC") # Modify column names to prepare for the next merge
  Gene <- bitr(DEGs.dt$SYMBOL,
               fromType = "SYMBOL", # Type of input data
               toType = c("ENTREZID"), # Type of data to convert to
               OrgDb = org.Hs.eg.db) # Species
  DEGs <- merge(DEGs.dt,Gene,by = "SYMBOL")

  DEGs1<-DEGs[,c(1,2,9)]
  dega<-DEGs1[,c(3,2)]
  geneList = dega[,2]
  names(geneList) = as.character(dega[,1])
  geneList = sort(geneList, decreasing = TRUE)

  if (data.base=='KEGG') {

  KEGG <- enrichKEGG(gene = Gene$ENTREZID,
                     organism = "hsa",
                     keyType = "kegg", #KEGG数据库
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)

  kk2 <- gseKEGG(geneList = geneList,
                 organism = 'hsa',
                 keyType = "kegg",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = Pvalue,
                 pAdjustMethod = "none",
                 by = "fgsea")
  kegg_result <- as.data.frame(kk2)
  rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
  af=as.data.frame(kk2@result)
  if (purrr::is_empty(rownames(af))) {
    stop('\n=> The GSEA results from the KEGG database are empty; please adjust the p-value.')
  }
  View(af)
  write.csv(af,file="EA_plot/KEGG_GSEA.xls",sep="\t",quote=F,col.names=T)

  # 结果可视化--分别取GSEA结果的前5个后5个展示 -----------------------------------------------
  pdf(file="GSEA_KEGG_top5.pdf",width = 8,height = 8)
  View(as.data.frame(kk2@result))
  num=5
  p <- gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
  print(p)
  dev.off()
  p

  # GSEA 结果表达分布的脊线图 ---------------------------------------------------------
  #跑分和预排序是可视化 GSEA 结果的传统方法,可视化基因集的分布和富集分数
  pdf(file="EA_plot/GSEA_KEGG_ridge.pdf",width = 8,height = 8)
  p1 <-ridgeplot(kk2,
            showCategory = 20,
            fill = "p.adjust",
            core_enrichment = TRUE,
            label_format = 32)
  print(p1)
  dev.off()

  # 单独展示某一个条目 ---------------------------------------------------------------
  if (!is.na(geneset.KEGG)) {
    pdf(file="EA_plot/GSEA_KEGG_single.pdf",width = 8,height = 8)
    p2 <- gseaplot2(kk2,
                    title = "GSEA KEGG",  #设置标题
                    geneset.KEGG, #绘制hsa04658通路的结果，通路名称与编号对应
                    color=col, #线条颜色
                    base_size = 20, #基础字体的大小
                    subplots = 1:3,
                    pvalue_table = T) # 显示p值
    print(p2)
    dev.off()
    }
  }

  if (data.base=='GO') {
  gsea_go <- gseGO(geneList     = geneList,#根据LogFC排序后的基因列表
                   OrgDb        = org.Hs.eg.db,
                   ont          = "ALL",#GO分析的模块
                   minGSSize    = 10,#最小基因集的基因数
                   maxGSSize    = 500,#最大基因集的基因数
                   pvalueCutoff = Pvalue,#p值的阈值
                   verbose      = FALSE)#是否输出提示信息
  bf <- as.data.frame(gsea_go)
  if (purrr::is_empty(rownames(bf))) {
    stop('\n=> The GSEA results from the GO database are empty; please adjust the p-value.')
  }
  View(bf)
  write.table(bf,file="EA_plot/GO_GSEA.xls",sep="\t",quote=F,col.names=T)
  #绘图
  if (!is.na(geneset.GO)){
    pdf(file="EA_plot/GSEA_GO_single.pdf",width = 8,height = 8)
      p <- gseaplot2(gsea_go,
                     geneSetID=geneset.GO,
                     title = "GSEA GO",
                     color = col,
                     base_size = 11,
                     rel_heights = c(1.5, 0.5, 1),
                     subplots = 1:3,
                     pvalue_table = TRUE,
                     ES_geom = "line")
      print(p)
    dev.off()
    }
  }
}
