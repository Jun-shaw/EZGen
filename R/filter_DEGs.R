#' Calculating the area under the curve after developing the category predictive model
#'
#' @param res.by.ML.Dev.Pred.Category.Sig  Output of function ML.Dev.Pred.Category.Sig
#' @param cohort.for.cal A data frame with the 'ID' and 'Var' as the first two columns. Starting in the fourth column are the variables that contain variables of the model you want to build. The second column 'Var' only contains 'Y' or 'N'.
#' @import DESeq2
#' @import limma
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @return A data frame containing the AUC of each predictive model.
#' @export
#'
#' @examples no examples
#'
filter_DEGs<-function(expr.data ='',
                      data.type ='mRNA',
                      deg.method='DESeq2',
                      Pvalue    =0.05,
                      log2FC    =2,
                      TCGA      =T,
                      tumor.num =17,
                      normal.num=13,
                      color1    ='#CE5C69',
                      color2    ='#5770A6',
                      title     ='TCGA',
                      key.genes =NA)
{
  library(DESeq2)#read count数据即可
  library(limma)
  library(dplyr)
  # 读取整理表达矩阵 ------------------------------------------------------------------
  #读取文件------------------------------------------------------------------
  df <- read.csv(expr.data)
  expr <- df %>%
    group_by(across(1)) %>%  # 根据第一列进行分组
    summarize(across(everything(), max, .names = "{col}"))  # 计算每列的最大值
  expr <- as.data.frame(expr)

  row.names(expr) <- expr[,1]
  expr <-expr[,-1]

  expr <- expr %>%
    mutate(across(everything(), as.numeric))
  expr <- ceiling(expr)#取整数
  expr <- expr[rowMeans(expr) > 1, ]  # 保留行均值大于1的行（排除不重要的极端值基因）
  expr[expr < 0] <- 0
  # 创建分组信息 -------------------------------------------------------------------------
  if (TCGA==T) {


  Tumor <- grep('01A',colnames(expr))
  Normal <- grep('11A',colnames(expr))
  Tumor_matrix <- expr[,Tumor]#提取肿瘤样本组矩阵
  Normal_matrix <- expr[,Normal]#提取正常样本组矩阵
  expr <- cbind(Tumor_matrix,Normal_matrix)#根据行合并矩阵
  group <- factor(c(rep("Tumor",times=length(Tumor)),rep("Normal",times=length(Normal))))#创建分组（因子变量）
  Data <- data.frame(row.names = colnames(expr), group = group)#创建分组数据框
  write.csv(expr,paste0(data.type, "_matrix.csv"),row.names = TRUE)#重新写入整理后的矩阵文件
  message("\n=> Complete the basic configuration.")
  }else {
    group <- factor(c(rep("Tumor",times=tumor.num),rep("Normal",normal.num)))
    Data <- data.frame(row.names = colnames(expr), group = group)
    write.csv(expr,paste0(data.type, "_matrix.csv"),row.names = TRUE)#重新写入整理后的矩阵文件
    message("\n=> Complete the basic configuration.")
  }
  if (deg.method=='DESeq2'){dds <- DESeqDataSetFromMatrix(countData = expr,
                                                       colData = Data,
                                                       design = ~ group)
  #第二步：开始差异分析
  dds2 <- DESeq(dds)
  res <- results(dds2, contrast=c("group", "Tumor", "Normal"))#后者为对照组
  res <- res[order(res$pvalue),]#按P值从小到大排序

  my_result <- as.data.frame(res)#转成容易查看的数据框
  my_result <- na.omit(my_result)#删除倍数为0的值

  #第三步：保存差异分析的结果
  my_result$Gene_symbol<-rownames(my_result)
  my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
  rownames(my_result) <- NULL
  my_result$padj[my_result$padj == 0] <- 9.9e-300
  write.csv(my_result,file=paste0(data.type, "_deseq2.csv"))#写入
  }
  if (deg.method=='limma') {
    design <- model.matrix(~0+factor(group))
    colnames(design) <- levels(factor(group))
    rownames(design) <- colnames(expr)

    # #构建比较矩阵——contrast -------------------------------------------------------
    contrast.matrix <- makeContrasts(Tumor-Normal,levels = design)

    # #线性拟合模型构建 ---------------------------------------------------------------
    fit <- lmFit(expr,design) #非线性最小二乘法
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
    my_result <- topTable(fit2, coef = 1,n = Inf)
    my_result$Gene_symbol <- row.names(my_result)
    my_result <- setNames(my_result, c("log2FoldChange", 'AveExpr','t',"pvalue", "padj",'B','Gene_symbol'))
    my_result$padj[my_result$padj == 0] <- 9.9e-300
    write.csv(my_result,file=paste0(data.type, "_limma.csv"))#写入
  }
  # DEGs的筛选 ------------------------------------------------------------------
  my_result$regulate <- ifelse(my_result$padj > Pvalue, "unchanged",
                               ifelse(my_result$log2FoldChange > log2FC, "up-regulated",
                                      ifelse(my_result$log2FoldChange < -log2FC, "down-regulated", "unchanged")))
  #可以把上调基因和下调基因取出放在一块
  DEGs <-subset(my_result, padj < Pvalue & abs(log2FoldChange) > log2FC)
  upgene <- DEGs[DEGs$regulate=='up-regulated',]
  downgene <- DEGs[DEGs$regulate=='down-regulated',]

  write.csv(DEGs,file= paste0(data.type, "_DEGs.csv"))#写入
  message("\n=> Visualization of the data.")
  pdf(file=paste0(data.type, "_Volcano Plot.pdf"),width = 10,height = 8)
  if (is.character(key.genes) && length(key.genes) > 0) {
    library(ggrepel)
    my_result$log10padj <- -log10(my_result$padj)#生成新的一列v
    non_zero_values <- my_result$v[my_result$v != Inf]# 提取padj列中所有非零值
    my_result$v[my_result$v == 0] <- sample(non_zero_values, sum(my_result$v == 0), replace = TRUE)# 替换padj列中的0值
    pdf(file=paste0(data.type, "_Volcano Plot.pdf"),width = 10,height = 8)
    vp <- ggscatter(my_result,
                    x = "log2FoldChange",
                    y = "log10padj",
                    ylab = "-log10(adjust p-value)",
                    size = 2,
                    color = "regulate",
                    palette = c(color2,'#DFE0DF',color1)) +#P值分界线
      geom_vline(xintercept = c(-log2FC,log2FC),lty=4,col ="gray",lwd=0.8)+ #FC分界线
      geom_hline(yintercept=-log10(Pvalue),lty=2,col = "gray",lwd=0.6)+#P值分界线 #可以去掉
      geom_text_repel(data = subset(my_result, Gene_symbol %in% key.genes),
                      aes(label = Gene_symbol),
                      color = "black",
                      box.padding = 0.5,
                      point.padding = 1,
                      segment.color = "black",
                      show.legend = FALSE, max.overlaps = 3000 )
    print(vp)
    dev.off()
    print(vp)}else{
      library(ggplot2)
      vp <-ggplot(data=my_result, aes(x=log2FoldChange, y=-log10(padj),color=regulate)) +
        geom_point(shape = 16, size=2) +
        theme_set(theme_set(theme_bw(base_size=20))) +
        xlab("log2 fold change") + #X轴标题
        ylab("-log10 p-value") + #Y轴标题
        theme(plot.title = element_text(size=15,hjust = 2.5)) +
        theme_classic()+
        scale_colour_manual(values = c(color2,'#DFE0DF',color1))+#颜色自定义
        geom_vline(xintercept = c(-log2FC,log2FC),lty=4,col ="gray",lwd=0.8)+ #FC分界线
        geom_hline(yintercept=-log10(Pvalue),lty=2,col = "gray",lwd=0.6)+#P值分界线
        labs(title=title)+#标题
        annotate("text",x=upgene$log2FoldChange[1:3],y=(-log10(upgene$padj[1:3])),label=upgene$Gene_symbol[1:3], size=5.0)+
        annotate("text",x=downgene$log2FoldChange[1:3],y=(-log10(downgene$padj[1:3])),label=downgene$Gene_symbol[1:3], size=5.0)
      print(vp)
      dev.off()
      print(vp)
  }
  }

