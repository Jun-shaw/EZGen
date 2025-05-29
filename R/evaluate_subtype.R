#' 评估亚型
#'
#' 通过开发类别预测模型后计算曲线下面积
#'
#' @param expr 包含表达数据的文件路径
#' @param genelist1 包含基因列表1的文件路径
#' @param genelist2 包含基因列表2的文件路径
#' @param clust.num 聚类数目
#'
#' @return 包含每个预测模型的AUC的数据框
#' @export
#'
#' @examples
#' \dontrun{
#' evaluate_subtype(expr = "matrix.csv",
#'                genelist1 = "genelist1.csv",
#'                genelist2 = "genelist2.csv",
#'                clust.num = 2:8,
#'                center = TRUE,
#'                scale = TRUE,
#'                filepath = getwd(),
#'                filename ='optimal_number_cluster')
#' }
#'
evaluate_subtype <- function(expr=NULL,
                           genelist1.file='',
                           genelist2.file='',
                           genelist1.name='',
                           genelist2.name='',
                           clust.num=2:8,
                           center = TRUE,
                           scale = TRUE,
                           filepath = getwd(),
                           filename ='optimal_number_cluster')
{
  df <- expr
  genelist1.file <- read.csv(genelist1.file)
  genelist1.file <- genelist1.file[,1]
  df1 <- df[rownames(df) %in% genelist1.file,]
  data1 <- df1

  genelist2.file <- read.csv(genelist2.file)
  genelist2.file <- genelist2.file[,1]
  df2 <- df[rownames(df) %in% genelist2.file,]
  data2 <- df2

  data.list <- list(data1)
  names(data.list) <- c(genelist1.name)
  assign("data.list", data.list, envir = .GlobalEnv)
  saveRDS(data.list, file = "data.list.rds")

  dat <- lapply(data.list, function(dd) {
    if (!all(dd >= 0))
      dd <- pmax(dd + abs(min(dd)), 0) + .Machine$double.eps
    dd <- dd/max(dd)
    return(dd %>% as.matrix)
  })
  dat <- lapply(dat, function(x) t(x) + .Machine$double.eps)

  message("calculating Cluster Prediction Index...")
  optk1 <- IntNMF::nmf.opt.k(dat = dat, n.runs = 5, n.fold = 5,
                             k.range = clust.num, st.count = 10, maxiter = 100,
                             make.plot = FALSE)
  optk1 <- as.data.frame(optk1)
  message("calculating Gap-statistics...")
  moas <- data.list %>% mogsa::mbpca(ncomp = 2, k = "all",
                                       method = "globalScore", center = center, scale = scale,
                                       moa = TRUE, svd.solver = "fast", maxiter = 1000, verbose = FALSE)
  gap <- mogsa::moGap(moas, K.max = max(clust.num), cluster = "hclust",
                      plot = FALSE)
  optk2 <- as.data.frame(gap$Tab)[-1, ]
  N.clust <- as.numeric(which.max(apply(optk1, 1, mean) + optk2$gap)) + 1
  if (length(N.clust) == 0) {
    message("--fail to define the optimal cluster number!")
    N.clust <- "null"
  }
  outFig <- paste0(filename, ".pdf")
  par(bty = "o", mgp = c(1.9, 0.33, 0), mar = c(3.1, 3.1, 2.1,
                                                3.1) + 0.1, las = 1, tcl = -0.25)
  plot(NULL, NULL, xlim = c(min(clust.num), max(clust.num)),
       ylim = c(0, 1), xlab = "Number of Clusters",
       ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "#EAE9E9", border = FALSE)
  grid(col = "white", lty = 1, lwd = 1.5)
  points(clust.num, apply(optk1, 1, mean), pch = 19, col = ggplot2::alpha("#5770A6"),
         cex = 1.5)
  lines(clust.num, apply(optk1, 1, mean), col = "#5770A6",
        lwd = 2, lty = 4)
  mtext("Cluster Prediction Index", side = 2, line = 2, cex = 1.5,
        col = "#5770A6", las = 3)
  par(new = TRUE, xpd = FALSE)
  plot(NULL, NULL, xlim = c(min(clust.num), max(clust.num)),
       ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(clust.num, optk2$gap, pch = 19, col = ggplot2::alpha("#CE5C69",
                                                                0.8), cex = 1.5)
  lines(clust.num, optk2$gap, col = "#CE5C69", lwd = 2, lty = 4)
  axis(side = 4, at = seq(0, 1, 0.2), labels = c("0.0", "0.2",
                                                 "0.4", "0.6", "0.8", "1.0"))
  mtext("Gap-statistics", side = 4, line = 2, las = 3, cex = 1.5,
        col = "#CE5C69")
  invisible(dev.copy2pdf(file = file.path(filepath, outFig),
                         width = 5.5, height = 5))
  message("visualization done...")
  if (N.clust > 1) {
    message(paste0("--the imputed optimal cluster number is ", N.clust))
  }
  return(list(N.clust = N.clust, CPI = optk1, Gapk = optk2))
  opt.clust.num <- N.clust

}
