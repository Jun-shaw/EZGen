#' @name subtype_heatmap
#' @param method A string to select each methon.
#' @param genelist.title A string vector to assign titles for each geneset.
#' @param legend.name A string vector to assign legend title for each subdata.
#' @param clust.dend A dendrogram object returned returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param show.col.dend A logical vector to indicate if showing the dendrogram for column at the top of heatmap.
#' @param show.colnames A logical vector to indicate if showing the names for column at the bottom of heatmap.
#' @param show.row.dend A logical vector to indicate if showing the dendrogram for row of each subdata.
#' @param show.rownames A logical vector to indicate if showing the names for row of each subdata.
#' @param clust.dist.row A string vector to assign distance method for clustering each subdata at feature dimension.
#' @param clust.method.row A string vector to assign clustering method for clustering each subdata at feature dimension.
#' @param clust.col A string vector storing colors for annotating each subtype at the top of heatmap.
#' @param color A list of string vectors storing colors for each subheatmap of subdata.
#' @param annCol A data.frame storing annotation information for samples with exact the same sample order with data parameter.
#' @param annColors A list of string vectors for colors matched with annCol.
#' @param width An integer value to indicate the width for each subheatmap with unit of cm.
#' @param height An integer value to indicate the height for each subheatmap with unit of cm.
#' @param file.save.path A string value to indicate the output path for storing the comprehensive heatmap.
#' @param file.name A string value to indicate the name of the comprehensive heatmap.
#' @return A pdf of multi-omics comprehensive heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap rowAnnotation anno_mark draw ht_opt %v%
#' @importFrom ClassDiscovery distanceMatrix
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom circlize colorRamp2
#' @importFrom dplyr %>%
#' @references Lu, X., Meng, J., Zhou, Y., Jiang, L., and Yan, F. (2020). MOVICS: an R package for multi-omics integration and visualization in cancer subtyping. Bioinformatics, btaa1018. [doi.org/10.1093/bioinformatics/btaa1018]
#' @references Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#' @export
#' @examples # There is no example and please refer to vignette.

subtype_heatmap <- function(data=data.list,
                            method='iClusterBayes',
                            genelist1.name='Disulfidptosis',
                            genelist2.name='Glucose_starvation',
                            feat.num=10,
                            col1 = c("#00FF00", "#008000", "#000000", "#800000", "#FF0000"),
                            col2 = c("#5770A6", "white"  , "#CE5C69"),
                            clust.col = c("#5770A6", "#A281B1", "#CE5C69", "#BDD5EA", "#FFA5AB", "#011627","#023E8A","#9D4EDD"),
                            genelist.title=c("Disulfidptosis","Glucose starvation"),
                            legend.name   = c("mRNA","mRNA "),
                            clust.dend    = NULL, # 不显示树状图
                            show.rownames = c(F,F), # 基因名显示设置
                            show.colnames = FALSE, # 样本名显示设置
                            annCol        = NULL, # 不对样本进行注释
                            annColors     = NULL, # 不设置注释颜色
                            width         = 10, # 每个子热图的宽度
                            height        = 5, # 每个子热图的高度
                            file.save.path=getwd(),
                            file.name      = "subtype heatmap")
{
library(ComplexHeatmap)
clust.dist.row   = c("pearson","pearson")
clust.method.row = c("ward.D","ward.D")
show.col.dend    = TRUE
show.row.dend    = c(TRUE, TRUE)
n_dat <- length(data.list)
halfwidth  = c(1,1);centerFlag = c(T,T);scaleFlag  = c(T,T)
# 数据标准化处理
# check data
if(is.null(names(data))){
  names(data) <- sprintf("dat%s", 1:length(data))
}
stdata <- list()
for (i in 1:n_dat) {
  tmp <- t(scale(t(data[[i]]), center = centerFlag[i], scale = scaleFlag[i]))
  if (!is.na(halfwidth[i])) {
    tmp[tmp > halfwidth[i]] = halfwidth[i]
    tmp[tmp < (-halfwidth[i])] = -halfwidth[i]
  }
  stdata[[names(data)[i]]] <- tmp
}


# 根据选定的方法获取特征和聚类结果
if (method == 'iClusterBayes') {
  feat <- clust.res.list$iClusterBayes$feat.res
  clust.res <- clust.res.list$iClusterBayes$clust.res
}
else if (method == 'MoCluster') {
  feat <- clust.res.list$MoCluster$feat.res
  clust.res <- clust.res.list$MoCluster$clust.res
}
else if (method == 'CIMLR') {
  feat <- clust.res.list$CIMLR$feat.res
  clust.res <- clust.res.list$CIMLR$clust.res
}

feat1  <- feat[which(feat$dataset == genelist1.name),][1:feat.num,"feature"]
feat2  <- feat[which(feat$dataset == genelist2.name),][1:feat.num,"feature"]
annRow <- list(feat1, feat2)
col.list   <- list(col1, col2)

ht_opt$message = FALSE
defaultW <- getOption("warn")
options(warn = -1)

colvec <- clust.col[1:length(unique(clust.res$clust))]
names(colvec) <- paste0("CS",unique(clust.res$clust))

if(!is.null(annCol) & !is.null(annColors)) {

  annCol <- annCol[colnames(stdata[[1]]), , drop = FALSE]
  annCol$Subtype <- paste0("CS",clust.res[colnames(stdata[[1]]),"clust"])
  annColors[["Subtype"]] <- colvec

  if(is.null(clust.dend)) {
    clust.res <- clust.res[order(clust.res$clust),]
    annCol <- annCol[clust.res$samID, , drop = FALSE]
  }

  ha <- ComplexHeatmap::HeatmapAnnotation(df     = annCol,
                                          col    = annColors,
                                         border = FALSE)
} else {
  annCol <- data.frame("Subtype" = paste0("CS",clust.res[colnames(stdata[[1]]),"clust"]),
                       row.names = colnames(stdata[[1]]),
                       stringsAsFactors = FALSE)
  annColors <- list("Subtype" = colvec)

  if(is.null(clust.dend)) {
    clust.res <- clust.res[order(clust.res$clust),]
    annCol <- annCol[clust.res$samID,,drop = FALSE]
  }

  ha <- ComplexHeatmap::HeatmapAnnotation(df     = annCol,
                                          col    = annColors,
                                          border = FALSE)
}

# 生成热图
ht <- list()
for (i in 1:n_dat) {

  hcg <- hclust(ClassDiscovery::distanceMatrix(as.matrix(t(stdata[[i]])), clust.dist.row[i]), clust.method.row[i])

  if(is.null(annRow[[i]][1])) {
    rowlab <- ""
    rowlab.index <- 0
  } else if (is.na(annRow[[i]][1])) {
    rowlab <- ""
    rowlab.index <- 0
  } else {
    rowlab <- intersect(rownames(stdata[[i]]),annRow[[i]])
    rowlab.index <- match(rowlab, rownames(stdata[[i]]))
  }


if(is.null(clust.dend)) {
  dt <- lapply(stdata, function(x) x[,clust.res$samID])
ht[[i]] <-  ComplexHeatmap::Heatmap(matrix               = as.matrix(dt[[i]]),
                                    row_title            = genelist.title[i],
                                    name                 = legend.name[i],
                                    cluster_columns      = F,
                                    cluster_rows         = hcg,
                                    show_column_dend     = show.col.dend,
                                    show_column_names    = show.colnames,
                                    show_row_dend        = show.row.dend[i],
                                    show_row_names       = show.rownames[i],
                                    col                  = grDevices::colorRampPalette(col.list[[i]])(64),
                                    top_annotation       = switch((i == 1) + 1, NULL, ha),
                                    width                = grid::unit(width, "cm"),
                                    height               = grid::unit(height, "cm"),
                                    heatmap_legend_param = list(at     = pretty(range(dt[[i]])),
                                                                labels = pretty(range(dt[[i]]))),
                                    right_annotation     = ComplexHeatmap::rowAnnotation(link =
                                                                                           anno_mark(at         = rowlab.index,
                                                                                                     labels     = rowlab,
                                                                                                     which      = "row",
                                                                                                     lines_gp   = grid::gpar(fontsize = 5),
                                                                                                     link_width = grid::unit(3, "mm"),
                                                                                                     padding    = grid::unit(0.8, "mm"),
                                                                                                     labels_gp  = grid::gpar(fontsize = 7))))


}else {
  ht[[i]] <-  ComplexHeatmap::Heatmap(matrix               = as.matrix(dt[[i]]),
                                      row_title            = genelist.title[i],
                                      name                 = legend.name[i],
                                      cluster_columns      = F,
                                      cluster_rows         = hcg,
                                      show_column_dend     = show.col.dend,
                                      show_column_names    = show.colnames,
                                      show_row_dend        = show.row.dend[i],
                                      show_row_names       = show.rownames[i],
                                      col                  = grDevices::colorRampPalette(col.list[[i]])(64),
                                      top_annotation       = switch((i == 1) + 1, NULL, ha),
                                      width                = grid::unit(width, "cm"),
                                      height               = grid::unit(height, "cm"),
                                      heatmap_legend_param = list(at     = pretty(range(dt[[i]])),
                                                                  labels = pretty(range(dt[[i]]))),
                                      right_annotation     = ComplexHeatmap::rowAnnotation(link =
                                                                                             anno_mark(at         = rowlab.index,
                                                                                                       labels     = rowlab,
                                                                                                       which      = "row",
                                                                                                       lines_gp   = grid::gpar(fontsize = 5),
                                                                                                       link_width = grid::unit(3, "mm"),
                                                                                                       padding    = grid::unit(0.8, "mm"),
                                                                                                       labels_gp  = grid::gpar(fontsize = 7))))
  }
}
# 将子热图合并
if(n_dat == 1){
  ht_list <- ht[[1]]
}
if(n_dat == 2){
  ht_list <- ht[[1]] %v% ht[[2]]
}
# 输出到 PDF 文件
outFile <- file.path(file.save.path,paste0(file.name,".pdf"))
if(is.null(annCol)) {
  pdf(outFile, width = width, height = height * n_dat/2)
} else {
  pdf(outFile, width = width, height = height * n_dat/1.5)
}
draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right") # output to pdf
invisible(dev.off())

draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right") # output to screen

options(warn = defaultW)


}
