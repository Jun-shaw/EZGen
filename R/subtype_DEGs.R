#' Calculating the area under the curve after developing the category predictive model
#'
#' @param res.by.ML.Dev.Pred.Category.Sig  Output of function ML.Dev.Pred.Category.Sig
#' @param cohort.for.cal A data frame with the 'ID' and 'Var' as the first two columns. Starting in the fourth column are the variables that contain variables of the model you want to build. The second column 'Var' only contains 'Y' or 'N'.
#' @import DESeq2
#' @import limma
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import ggpubr
#' @return A data frame containing the AUC of each predictive model.
#' @export
#'
#' @examples
#'
subtype_DEGs<-function(clust.res.list='',
                      method         ='iClusterBayes', #("IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering")
                      expr.data           ='TCGA_matrix.csv',
                      subtype.list   =c('CS1','CS3'),
                      subtype.ctrl   ='CS2',
                      deg.method     ='DESeq2',
                      Pvalue         =0.05,
                      log2FC         =0.3,
                      col1           ='#CE5C69',
                      col2           ='#5770A6',
                      title          ='CS1 vs CS3',
                      key.genes      =NA)
  {
  library(DESeq2)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  df <- read.csv(expr.data)
  expr <- df %>%
    group_by(across(1)) %>%  # 根据第一列进行分组
    summarize(across(everything(), max, .names = "{col}"))  # 计算每列的最大值
  expr <- as.data.frame(expr)

  row.names(expr) <- expr[,1]
  expr <-expr[,-1]
  expr <- expr %>%
    mutate(across(everything(), as.numeric))
  expr <- ceiling(expr)
  expr <- expr[rowMeans(expr) > 1, ]
  expr[expr < 0] <- 0
  # 创建分组信息 -------------------------------------------------------------------------
  if (method=='all') {
    Data <-  clust.res.list$clust.res
    message(paste0("\n=> Select the comprehensive clustering machine learning algorithm."))}
  if (method=='IntNMF') {
    Data <-  clust.res.list$IntNMF$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='CIMLR') {
    Data <-  clust.res.list$CIMLR$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='PINSPlus') {
    Data <-  clust.res.list$PINSPlus$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='NEMO') {
    Data <-  clust.res.list$NEMO$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='COCA') {
    Data <-  clust.res.list$COCA$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='MoCluster') {
    Data <-  clust.res.list$MoCluster$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='LRAcluster') {
    Data <-  clust.res.list$LRAcluster$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='iClusterBayes') {
    Data <-  clust.res.list$iClusterBayes$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='SNF') {
    Data <-  clust.res.list$SNF$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}
  if (method=='ConsensusClustering'){
    Data <-  clust.res.list$ConsensusClustering$clust.res
    message(paste0("\n=> Select the ",method , " clustering machine learning algorithm."))}

  Data <- subset(Data, clust %in% as.numeric(gsub("\\D", "", subtype.list)))
  group <- factor(paste0("CS", Data$clust))
  Data <- data.frame(row.names = row.names(Data), group = group)
  expr <- expr[, colnames(expr) %in% rownames(Data), drop = FALSE]
  data_row_order <- rownames(Data)
  expr <- expr[, data_row_order, drop = FALSE]
  message("\n=> Complete the basic configuration.")

  if (deg.method=='DESeq2'){
    dds <- DESeqDataSetFromMatrix(countData = expr,
                                 colData = Data,
                                 design = ~ group)
  #第二步：开始差异分析
  dds2 <- DESeq(dds)
  res <- results(dds2, contrast=c("group", subtype.list[subtype.list != subtype.ctrl], subtype.ctrl))#后者为对照组
  res <- res[order(res$pvalue),]#按P值从小到大排序

  my_result <- as.data.frame(res)#转成容易查看的数据框
  my_result <- na.omit(my_result)#删除倍数为0的值

  #第三步：保存差异分析的结果
  my_result$Gene_symbol<-rownames(my_result)
  my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
  rownames(my_result) <- NULL
  my_result$padj[my_result$padj == 0] <- 9.9e-300
  write.csv(my_result,file=paste0(subtype.list[subtype.list != subtype.ctrl],'_vs_',subtype.ctrl,"_mRNA_deseq2.csv"))#写入
  }
  if (deg.method=='limma') {
    design <- model.matrix(~0+factor(group$group))
    colnames(design) <- levels(factor(group$group))
    rownames(design) <- colnames(expr)

    # #构建比较矩阵——contrast -------------------------------------------------------
    contrast.matrix <- makeContrasts(Tumor-Normal,levels = design)

    # #线性拟合模型构建 ---------------------------------------------------------------
    fit <- lmFit(expr,design) #非线性最小二乘法
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
    my_result <- topTable(fit2, coef = 1,n = Inf)
    my_result$Gene_symbol <- row.names(my_result)
    setnames(my_result, c("logFC", "P.Value", "adj.P.Val"), c("log2FoldChange", "pvalue", "padj"))
    my_result$padj[my_result$padj == 0] <- 9.9e-300
    write.csv(my_result,file=paste0(subtype.list[subtype.list != subtype.ctrl],'_vs_',subtype.ctrl,"_mRNA_limma.csv"))#写入
  }
  # DEGs的筛选 ------------------------------------------------------------------
  my_result$regulate <- ifelse(my_result$padj > Pvalue, "unchanged",
                               ifelse(my_result$log2FoldChange > log2FC, "up-regulated",
                                      ifelse(my_result$log2FoldChange < -log2FC, "down-regulated", "unchanged")))
  #可以把上调基因和下调基因取出放在一块
  DEGs <-subset(my_result, padj < Pvalue & abs(log2FoldChange) > log2FC)
  upgene <- DEGs[DEGs$regulate=='up-regulated',]
  downgene <- DEGs[DEGs$regulate=='down-regulated',]
  write.csv(DEGs,file= paste0(subtype.list[subtype.list != subtype.ctrl],'_vs_',subtype.ctrl,"_mRNA_DEGs.csv"))#写入
  message("\n=> Visualization of the data.")
  pdf(file=paste0(subtype.list[subtype.list != subtype.ctrl],'_vs_',subtype.ctrl,"_Volcano_Plot.pdf"),width = 8,height = 8)
  if (is.character(key.genes) && length(key.genes) > 0) {
    my_result$log10padj <- -log10(my_result$padj)#生成新的一列v
    non_zero_values <- my_result$v[my_result$v != Inf]# 提取padj列中所有非零值
    my_result$v[my_result$v == 0] <- sample(non_zero_values, sum(my_result$v == 0), replace = TRUE)# 替换padj列中的0值
    vp <- ggscatter(my_result,
                    x = "log2FoldChange",
                    y = "log10padj",
                    ylab = "-log10(adjust p-value)",
                    size = 2,
                    color = "regulate",
                    palette = c(col2,'#DFE0DF',col1)) +#P值分界线
      geom_vline(xintercept = c(-log2FC,log2FC),lty=4,col ="gray",lwd=0.8)+ #FC分界线
      geom_hline(yintercept=-log10(Pvalue),lty=2,col = "gray",lwd=0.6)+#P值分界线 #可以去掉
      geom_text_repel(data = subset(my_result, Gene_symbol %in% key.genes),
                      aes(label = Gene_symbol),
                      color = "black",
                      box.padding = 0.5,
                      point.padding = 1,
                      segment.color = "black",
                      show.legend = FALSE, max.overlaps = 3000 )
    plot(vp)
    dev.off()
  }else{
    vp <-ggplot(data=my_result, aes(x=log2FoldChange, y=-log10(padj),color=regulate)) +
      geom_point(shape = 16, size=2) +
      theme_set(theme_set(theme_bw(base_size=20))) +
      xlab("log2 fold change") + #X轴标题
      ylab("-log10 p-value") + #Y轴标题
      theme(plot.title = element_text(size=15,hjust = 2.5)) +
      theme_classic()+
      scale_colour_manual(values = c(col2,'#DFE0DF',col1))+#颜色自定义
      geom_vline(xintercept = c(-log2FC,log2FC),lty=4,col ="gray",lwd=0.8)+ #FC分界线
      geom_hline(yintercept=-log10(Pvalue),lty=2,col = "gray",lwd=0.6)+#P值分界线
      labs(title='CS3 vs. CS1')+#标题
      annotate("text",x=upgene$log2FoldChange[1:3],y=(-log10(upgene$padj[1:3])),label=upgene$Gene_symbol[1:3], size=5.0)+
      annotate("text",x=downgene$log2FoldChange[1:3],y=(-log10(downgene$padj[1:3])),label=downgene$Gene_symbol[1:3], size=5.0)
    plot(vp)
    dev.off()
  }
  clust.res <- data.frame(
    samID = row.names(Data),  # 使用 Data 的行名作为 samID 列
    clust = group,            # 将 group 向量作为 clust 列
    row.names = row.names(Data)  # 设置行名
  )
  return(list(vp = vp, expr=expr, clust.res= clust.res, res = my_result, DEGs = DEGs))
}

