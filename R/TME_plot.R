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
#多轨热图绘制#


TME_plot <-function(TME.res.list=TME.res.list,
                    group =DEGs.res$group,
                        color2=c( "#663366", "#ffffff","#ff9933"),
                        color3=c( "#cccc00", "#ffffff", "#66cc99"),
                        color4=NA,
                        threshold1=0.05,
                        threshold2=0.01,
                        category=T,
                        category_colors=c('#CC88B0', '#DBE0ED', '#87B5B2', '#F4CEB4', '#F1DFA4', '#998DB7')
)
{
  #假设有个热图的矩阵数据（这里仅为一组重复两次以作示范）
  library(ComplexHeatmap)
  library(data.table)
  library(tidyverse)
  library(circlize)
  library(dplyr)
  library(RColorBrewer)
  print('处理数据......')
  data <- TME.res.list
  data <- lapply(data, function(mat) mat[, -1])

  data <- lapply(names(data), function(name) {
    mat <- data[[name]]  # 获取当前矩阵
    new_row <- matrix(name, nrow = 1, ncol = ncol(mat))  # 创建新行
    colnames(new_row) =colnames(mat)  # 设置新行的列名
    mat <- rbind(mat, new_row)  # 将新行添加到矩阵
    return(mat)  # 返回修改后的矩阵
  })


  df=data %>% as.data.frame()#数据转化为数据框

  df <- t(df)
  colnames(df)[ncol(df)] <- "category"

  #如果矩阵数据分组，可用split参数来指定分类变量


  df <- as.data.frame(df)

  group$ID <- row.names(group)

  subtype1 <- c(group$ID[group$group == unique(group$group)[1]],"category")
  subtype2 <- c(group$ID[group$group == unique(group$group)[2]],"category")

  category.list <- unique(ann_row)

  cir1 <- df %>% select(all_of(subtype1))
  cir2 <- df %>% select(all_of(subtype2))

  cir.res.list1 <- list()
  cir.res.list2 <- list()

  for (i in 1:nrow(category.list)) {
    current_category <- category.list[i,]
    cir1_i <- cir1[cir1$category ==current_category,]
    cir1_i <- cir1_i[,-ncol(cir1_i)]
    cir1_i[] <- lapply(cir1_i, function(x) as.numeric(as.character(x)))
    cir.res.list1[[i]] <- cir1_i
    names(cir.res.list1)[[i]] <- as.character(current_category)
  }

  for (i in 1:nrow(category.list)) {
    current_category <- category.list[i, ]
    cir2_i <- cir2[cir1$category ==current_category,]
    cir2_i <- cir2_i[,-ncol(cir2_i)]
    cir2_i[] <- lapply(cir2_i, function(x) as.numeric(as.character(x)))
    cir.res.list2[[i]] <- cir2_i
    names(cir.res.list2)[[i]] <- as.character(current_category)
  }

  if(category==T){
    ann_row =data.frame(pathway=df$category)#对行进行注释，用于后续的热图分裂

    row.names(ann_row) = rownames(cir1)

    ann_row <- as.matrix(ann_row)}#在circlize函数中，需要为matrix
  cir_group <- as.matrix(df[, ncol(df)])

  get_stars <- function(value) {
    if (value < threshold2) {
      stars <- "**"
    } else if (value < threshold1) {
      stars <- "*"
    } else {
      stars <- ""
    }
    return(stars)
  }

  cir1 <- cir1[,-27]
  cir2 <- cir2[,-15]

  cir1[] <- lapply(cir1, function(x) as.numeric(as.character(x)))
  cir2[] <- lapply(cir2, function(x) as.numeric(as.character(x)))

  cir1$category <- ann_row
  cir2$category <- ann_row

  print('开始生成热图.......')


  #定义热图颜色梯度：
  color1=c("#5770A6","#A281B1","#CE5C69");range1=c(0, 0.5, 1)
  color2=c("#25849d","#d4e9f1","#9871e0");range2=c(0, 5, 10)
  color3=c("#ffd984","#85a543","#368838");range3=c(0, 5, 10)
  color4=c("#5770A6","#A281B1","#CE5C69");range4=c(0, 0.5, 1)
  color5=c("#b39943","#ffc3ad","#fc8c79");range5=c(0, 0.5, 2)
  color6=c("#ff6f7c","#ff4ca2","#cb49c7");range6=c(-4000, 0, 4000)
  color7=c("#00aca0","#00cd96","#5bec7e");range7=c(0, 0.5, 1)

  mycol1=colorRamp2(range1,color1)#设置legend颜色，
  mycol2=colorRamp2(range2,color2)
  mycol3=colorRamp2(range3,color3)
  mycol4=colorRamp2(range4,color4)
  mycol5=colorRamp2(range5,color5)
  mycol6=colorRamp2(range6,color6)
  mycol7=colorRamp2(range7,color7)

  # 确保没有 NULL 值
  color1 <- ifelse(is.null(color1), "black", color1)

  mycol <- c(mycol1, mycol2,mycol3,mycol4,mycol5,mycol6,mycol7)

  # 创建命名向量将这两者配对
  row_colors <- setNames(mycol, category.list)

  circos.clear()
  circos.par(gap.after=c(2,2,2,2,2,2,2))#circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息

  circos.heatmap(cir1,

                 col= mycol1, #用行注释分裂热图=

                 split=ann_row,

                 rownames.col="black",

                 show.sector.labels = F,

                 track.height = 0.1, #轨道的高度，数值越大圆环越粗

                 rownames.side="outside",

                 rownames.cex=0.5,#字体大小

                 rownames.font=1.2,#字体粗细

                 bg.border="black", #背景边缘颜色

                 dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈

                 cluster=F,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类

                 dend.track.height=0.1,#调整行聚类树的高度

                 dend.callback=function(dend,m,si) { #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称

                   color_branches(dend,k=10,col=1:10) #color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色

                 }

  )





    circos.heatmap(cir2,

                   col = mycol1,

                   split=ann_row,

                   track.height = 0.1,

                   bg.border="black", #背景边缘颜色

                   rownames.cex=0.3)#加入第二个热图

  if(category==T){
    # 创建行名分类信息与颜色的配对
    group_names <-  unique(df$category)
    group_colors <- category_colors

    # 创建命名向量将这两者配对
    row_colors <- setNames(group_colors, group_names)
    circos.heatmap(cir_group,

                   col = row_colors,

                   split=ann_row,

                   track.height = 0.05,

                   bg.border="black", #背景边缘颜色

                   rownames.cex=0.3,

                   show.sector.labels = F,)
  }

  #添加列名#

  #第一个环形列名

  circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){

    if(CELL_META$sector.numeric.index==6){# the last sector

      cn=colnames(cir1)

      n=length(cn)

      circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(1,"mm"),#x坐标

                  (1:n)*1.2+5,#调整y坐标,行距+距离中心距(1:n)*1.2+5,

                  cn,cex=1,adj=c(0,0),facing="inside")

    }

  },bg.border=NA)

  print('添加图例......')
  library(circlize)
  library(gridBase)

  for (i in 1:outcome.num){
    lg_Exp=Legend(title=outcome[i],col_fun=mycol[[i]],direction = c("vertical"))
    assign(paste("lg_Exp", i, sep = ""), lg_Exp)
  }

  # 创建图例

  circle_size= unit(0.6,"snpc")

  h= dev.size()

  lgd_list= packLegend(lg_Exp1,lg_Exp2,lg_Exp3, max_height = unit(2*h,"inch"))
  if (outcome.num>1){lgd_list= packLegend(lg_Exp1,lg_Exp2, max_height = unit(2*h,"inch"))}
  if (outcome.num>2){lgd_list= packLegend(lg_Exp1,lg_Exp2,lg_Exp3, max_height = unit(2*h,"inch"))}
  if (outcome.num>3){lgd_list= packLegend(lg_Exp1,lg_Exp2,lg_Exp3,lg_Exp4, max_height = unit(2*h,"inch"))}
  draw(lgd_list, x = circle_size, just ="midle")

}

