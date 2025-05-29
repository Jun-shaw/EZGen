#' @name
#' @title title
#' @describeIn Calculating the area under the curve after developing the category predictive model
#'
#' @import dplyr
#' @import corrplot
#'
#' @param expr  Output of function ML.Dev.Pred.Category.Sig
#' @param group A data frame with the 'ID' and 'Var' as the first two columns. Starting in the fourth column are the variables that contain variables of the model you want to build. The second column 'Var' only contains 'Y' or 'N'.
#' @param gene
#' @param subtype.list
#' @param clust.col
#'
#' @return A data frame containing the AUC of each predictive model.
#' @export
#'
#' @examples
#'
cor_plot <- function(expr=DEGs.res$expr,
                     group=DEGs.res$clust.res,
                     gene.list=c('AR','KLK3','TMPRSS2','SPDEF','CHGA','SYP','ENO2','PEG10','SPIC','NEUROD1','AMIGO2','SLC7A11'),
                     subtype.list  = c('CS1','CS3'),
                     col=c("#5770A6",'white',"#CE5C69"),
                     file.name='')
{
  library(corrplot)
  library(dplyr)
  # 读取需要处理的数据 ---------------------------------------------------------------
  if (length(expr)==1){
    expr.dt <- read.csv(expr)
    expr.dt <- aggregate(. ~ expr.dt[[1]], data = expr.dt, FUN = max)
    row.names(expr.dt) <- expr.dt[,1]
    expr.dt <-expr.dt[,-c(1,2)]

  }else{expr.dt <- expr}

  expr.dt <- expr.dt %>%
    mutate(across(everything(), as.numeric))
  expr.use <- t(log2(expr.dt+1))

  if(is.na(group)){
    missing.genes <- gene.list[!gene.list %in% colnames(expr.use)]
    if (length(missing.genes) > 0) {

      stop(paste0("\n=> fail to match gene ",missing.genes))
    }
    expr.use <- data.frame(expr.use[,gene.list])
  }else{
  if (!grepl("^CS",group$clust[1])) {
    group$clust <- factor(paste0("CS", group$clust))
  }

    sample <- group[group$clust%in%subtype.list,]
    sample.id <- sample$samID
    missing.samples <- sample.id[!sample.id %in% row.names(expr.use)]
    if (length(missing.samples) > 0) {
      stop("\n=> fail to match sample id between expr and group.")
    }
    expr.use <- data.frame(expr.use[sample.id,gene.list])
  }

  expr.use <- na.omit(expr.use)
  cor.res <- cor (expr.use, method="pearson")
  test.res = cor.mtest(cor.res, method="pearson",conf.level = 0.95)
  # 开始绘图 --------------------------------------------------------------------
  if (!dir.exists('expr_plot')) {
    dir.create('expr_plot')
  }
  addcol <- colorRampPalette(col)
  pdf(paste0('expr_plot/',file.name,'_cor_plot.pdf'), width = 8, height = 8)
  corrplot(cor.res, method = "pie" ,type = "upper",col = addcol(100),
           #method = c("circle", "square", "ellipse", "number", "shade", "color", "pie"),
           tl.col = "black", tl.cex = 1.2, tl.srt = 45,tl.pos = "lt",
           p.mat = test.res$p, diag = T,
           sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
           insig = 'label_sig', pch.col = 'grey20')

  corrplot(cor.res, method = "number", type = "lower",col = addcol(100),
           tl.col = "n", tl.cex = 0.8, tl.pos = "n",add = T)
  dev.off()
}
