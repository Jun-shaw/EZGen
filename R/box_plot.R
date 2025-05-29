#' @name
#' @title title
#' @describeIn Calculating the area under the curve after developing the category predictive model
#'
#' @import ggplot2
#' @import ggpubr
#' @import Hmisc
#' @import patchwork
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

box_plot <- function(expr='ML_data/WCDT-MCRPC_mrna_expr_fpkm.csv',
                     group=clust.res.list$iClusterBayes$clust.res,
                     gene.list=c('AR','KLK3','TMPRSS2','SPDEF','CHGA','SYP','ENO2','PEG10','SPIC','NEUROD1','CMTM2','AMIGO2'),
                     subtype.list  = c('CS1','CS3'),
                     clust.col=c("#5770A6","#CE5C69"))
{
# 读取需要处理的数据 ---------------------------------------------------------------
  library(ggplot2)
  library(ggpubr)
  library(Hmisc)
  library(patchwork)

  if (length(expr)==1){
    expr.dt <- read.csv(expr)
    expr.dt <- aggregate(. ~ expr.dt[[1]], data = expr.dt, FUN = max)
    row.names(expr.dt) <- expr.dt[,1]
    expr.dt <-expr.dt[,-c(1,2)]

  }else{expr.dt <- expr}

  expr.dt <- expr.dt %>%
    mutate(across(everything(), as.numeric))
  expr.use <- t(log2(expr.dt+1))

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
plot.list <- list()
for (i in 1:length(gene.list)) {
gene <- gene.list[i]
gene.expr <- data.frame(expr.use[,gene])
gene.expr <- cbind(sample$clust,gene.expr)
colnames(gene.expr)[1] <- "Subtype"
colnames(gene.expr)[2] <- gene
gene.expr <- na.omit(gene.expr)

# 开始绘图 --------------------------------------------------------------------
if (!dir.exists('expr_plot')) {
  dir.create('expr_plot')
}

p <- ggplot(gene.expr, aes(x = Subtype, y = !!sym(gene), color = Subtype)) +
  geom_boxplot(size = 1.5) +  # Adjust the thickness of the boxplot
  geom_point(position = position_jitter(width = 0.1)) +
  theme(text = element_text(size = 12),  # Adjust text size and axis tithttp://127.0.0.1:20903/graphics/plot_zoom_png?width=1111&height=1078le size
        axis.title = element_text(size = 18)) +
  scale_x_discrete(labels = subtype.list) +  # Adjust x-axis labels
  labs(y = paste0(gene,' expression log2(fpkm+1)')) +  # Set title for y-axis
  theme_minimal(base_size = 12, base_family = "sans") +
  scale_color_manual(values = clust.col) +
  geom_signif(comparisons = list(subtype.list), test = "wilcox.test", map_signif_level = TRUE)
plot.list[[gene]] <- p
ggsave(filename = paste0('expr_plot/',gene, '_expression.pdf'), plot = p, width = 6, height = 8)
  }

wrap_plots(plot.list, byrow = T, nrow = 3)
}
