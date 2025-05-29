#' Calculating the area under the curve after developing the category predictive model
#'
#' @param res.by.ML.Dev.Pred.Category.Sig  Output of function ML.Dev.Pred.Category.Sig
#' @param cohort.for.cal A data frame with the 'ID' and 'Var' as the first two columns. Starting in the fourth column are the variables that contain variables of the model you want to build. The second column 'Var' only contains 'Y' or 'N'.
#'
#' @return A data frame containing the AUC of each predictive model.
#' @export
#'
#' @examples
#'

filter_merge <- function(DEGs1=data1,
                         DEGs2=data2,
                         genelist=genes,
                         color=c('#CE5C69','#5770A6','#E0A980'),
                         labels=c("DEGs of GSE6752", "DEGs of GSE48403 ",'Disulfidptosis-related genes'),
                         save.path='G:/importance/undergraduated/12_R包开发/ctDNA/VENE.TIF')
  {
  library (VennDiagram)
  library(openxlsx)
  set1=t(DEGs1)
  set2=t(DEGs2)
  set3=t(genelist)
venn.diagram(x=list(set1,set2,set3),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=color, #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=color, # 填充色 配色https://www.58pic.com/
             category.names = labels , #标签名
             cat.dist = 0.05, # 标签距离圆圈的远近
             cat.pos = c(-0, -0, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 1, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "text",  # 标签位置, outer内;text 外
             output=TRUE,
             filename=save.path,# 文件保存
             imagetype="tiff",  # 类型（tiff png svg）
             resolution = 600,  # 分辨率
             compression = "lzw"# 压缩算法

)
intersection_1_2 <- intersect(set1, set2)
intersection_1_3 <- intersect(set1, set3)
intersection_2_3 <- intersect(set2, set3)
intersection_1_2_3 <- intersect(intersection_1_2, set3)
write.csv(intersection_1_2,"candidate_genelist.csv")
write.csv(intersection_1_2_3,"popular_genelist.csv")
}
