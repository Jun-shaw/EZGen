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

create_nomogram <- function( cli.data=cli.data,
                             rs.data=NA,
                             cont.var=NULL,
                             cate.var=NULL,
                             Pvalue=0.10,
                             rate.range=c(1,3,5),
                             col=c("#5770A6", "#A281B1", "#CE5C69"))
{
  library(timeROC)
  library(rms) # 拟合模型
  library(regplot) # 添加个案的列线图
  library(DynNom) # 动态列线图
  library(VRPM) # 彩色条带式静态列线图
  library(survival) # 生存函数
  library(tidyr)
  library(forestplot)
  library(grid)
  library(checkmate)
  library(abind)
  library(dplyr)
  library(pROC)
  library(ggplot2)
  library(ggDCA)
  library(tidymodels)
  library(tidyverse)

  if (!dir.exists('nomogram_plot')) {
    dir.create('nomogram_plot')
  }
  message("\n=> Cox regression analysis.")
  #时间变量的处理
  units(mydata$time) <- "month"#将时间单位设置为月
  df <- mydata
  df<-na.omit(df)
  colnames(df)
  ####1. 设置相关参数
  ##告诉程序那些是连续变量，那些为分类变量
  cont_covariates <- cont.var
  cate_covariates <- cate.var
  ##分类变量转换成因子
  ###注意，对于有些变量，需要使用levels告诉程序这个变量值的顺序，如：Stage I、Stage II、 Stage III
  for (i in cate_covariates) {
    df[,i] <- factor(df[,i])
  }

  y<- Surv(time=df$time,event=df$status==1)#1为事件发生

  ############################start from here#########
  ############################start from here#########
  ############################start from here#########
  ############################start from here#########
  ####2. 单因素COX分析
  Uni_cox_model<- function(x){
    surv <- as.formula(paste0 ("y~",x))
    cox <- coxph(surv,data=df)
    sum <- summary(cox)
    HR <- round(sum$coefficients[,2],3)
    PValue <- round(sum$coefficients[,5],3) %>% sprintf("%1.3f",.)
    lower <- round(sum$conf.int[,3],3)
    upper <- round(sum$conf.int[,4],3)
    subchar <- rownames(sum$coefficients)
    HRa <- paste0(sprintf("%1.3f",HR) ," [",sprintf("%1.3f",lower),", ",sprintf("%1.3f",upper),"]")
    Uni_cox_model<- data.frame('Characteristics' = paste(x,"_",subchar,sep=""),
                               'HRa'=HRa,
                               'PValue' = PValue,
                               'HR' = HR,
                               'lower' = lower,
                               'upper' = upper
    )
    return(Uni_cox_model)
  }
  # Uni_cox_model("age")
  #转换成数据框，并转置
  univ_results  <- do.call(rbind,lapply(c(cont_covariates,cate_covariates) %>% na.omit ,Uni_cox_model))
  #最后，将P值=0的变为p<0.0001
  univ_results$PValue <- as.numeric(as.character(univ_results$PValue))
  univ_results$PValue[univ_results$PValue < 0.001] <-"<0.001"
  names(univ_results) <- c("Variants","Hazard Ratio (95%CI) ","P-value ","HR ","Lower ","Upper ")

  ########3. 多因素分析
  #1.提取单因素p<0.05变量
  univ2mul <- univ_results$Variants[univ_results$`P-value `< Pvalue]

  #多因素模型建立

  # mul_Variants <- do.call(c,lapply(univ2mul,function(x) strsplit(x,split="_")[[1]][1])) %>%
  #  unique()
  mul_Variants <- stringr::str_extract(univ2mul,paste0(c(cont_covariates,cate_covariates) %>% na.omit,
                                                       collapse = "|"))
  mul_cox_model<- as.formula(paste0 ("y~",
                                     paste0(mul_Variants,
                                            collapse = "+")))
  mul_cox <- coxph(mul_cox_model,data=df)
  cox_sum <- summary(mul_cox)

  #提取多因素回归的信息
  HR<- round(cox_sum$coefficients[,2],3)
  PValue<- round(cox_sum$coefficients[,5],3)
  lower<-round(cox_sum$conf.int[,3],3)
  upper<-round(cox_sum$conf.int[,4],3)

  #多因素结果优化并成表：mul_cox1
  HRa = paste0(sprintf("%1.3f",HR) ," [",sprintf("%1.3f",lower),", ",sprintf("%1.3f",upper),"]")

  mul_results<- data.frame("HRa"=HRa,"PValue"=PValue,"HR"=HR,"lower"=lower,"upper"=upper) %>%
    tibble::rownames_to_column(var = "Variants")
  mul_results$PValue[mul_results$PValue < 0.001] <-"<0.001"
  colnames(mul_results)=c("Variants","Hazard Ratio (95%CI)","P-value","HR","Lower","Upper")


  #### 4. 合并、整理单、多因素COX结果
  # univ_results$Variants <- do.call(c,lapply(univ_results$Variants,function(x) strsplit(x,split="_")[[1]][2]))
  univ_results$Variants <- stringr::str_remove(univ_results$Variants,paste0(paste0(c(cont_covariates,cate_covariates),"_") %>% na.omit,
                                                                            collapse = "|"))

  final_results<-merge(univ_results,mul_results,by="Variants",all=TRUE)

  #### 5. 格式化、美化
  #分类变量各分类样本数计算
  data_sum <- do.call(rbind,lapply(cate_covariates,function(x){
    ddd <- data.frame(Var1 = x,Freq = NA,Variants="")
    ddd <- table(df[,x]) %>% as.data.frame() %>% dplyr::mutate(Variants = paste0(x,.data$Var1)) %>%
      rbind(ddd,.)

  }))
  #添加进去连续变量，如果有
  data_sum <- rbind(data.frame(Var1 = cont_covariates,
                               Freq = apply(df[cont_covariates], 2, function(x) length(na.omit(x))),
                               Variants=cont_covariates),
                    data_sum)

  #如需调整输出顺序，可以在此处调整之后再进行下一步，也可以导出excel之后整理
  final_results2 <- left_join(data_sum,final_results,by="Variants",keep=T)
  final_results2 <- final_results2[-c(3:4)]
  final_results2$`Hazard Ratio (95%CI) `[which(final_results2$Var1 %in% c(cate_covariates))+1] = 1
  final_results2$`Hazard Ratio (95%CI)`[which(final_results2$Var1 %in%
                                                intersect(cate_covariates,mul_Variants))  +1] =1

  index <- which(final_results2$Var1 %in% c(cate_covariates))
  colnames(final_results2)[1:2] <- c("characteristic","N")
  message("\n=> Save the cox regression analysis result.")
  write.csv(final_results2,"cox.res.csv",row.names = F)

  #### 5. Forest plot for univariate cox results
  message("\n=> Forest plot for univariate cox results.")
  univ_results2 <- rbind(c(colnames(final_results2)[1:4],NA,NA,NA),final_results2[1:7])
  is.sum <- c(T,rep(F,nrow(univ_results2)-1))
  is.sum[which(univ_results2$characteristic %in% c(cate_covariates))] <- TRUE
  pdf('nomogram_plot/uni_cox.pdf', width = 8, height = 8)
  p1 <- forestplot(univ_results2[,colnames(univ_results2)[c(1,3,4)]],
                   mean=as.numeric(univ_results2[,"HR "]),
                   lower=as.numeric(univ_results2[,"Lower "]),
                   upper=as.numeric(univ_results2[,"Upper "]),
                   is.summary= is.sum,
                   hrzl_lines=list("1" = gpar(lty=1,lwd=1),
                                   "2" = gpar(lty=2)),
                   txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                  ticks=gpar(cex=1.1),
                                  xlab=gpar(cex = 1.2),
                                  title=gpar(cex = 1.2)),
                   ##fpColors函数设置颜色
                   col=fpColors(box=col[1], lines=col[1], zero = "gray50"),
                   #箱线图中基准线的位置
                   zero=1,
                   cex=0.9, lineheight = "auto",
                   colgap=unit(8,"mm"),xlog=T,
                   #箱子大小，线的宽度
                   lwd.ci=2,boxsize=0.5,
                   #箱线图两端添加小竖线，高度
                   ci.vertices=TRUE, ci.vertices.height = 0.2,
                   graph.pos = 4)
  print(p1)
  dev.off()
  #### 6. Forest plot for multivariate cox results
  message("\n=> Forest plot for multivariate cox results")
  mul_results2 <- rbind(c(colnames(final_results2)[c(1,8:9)],NA,NA,NA),final_results2[c(1,8:12)])
  is.sum <- c(T,rep(F,nrow(univ_results2)-1))
  is.sum[which(mul_results2$characteristic %in% c(cate_covariates))] <- TRUE
  pdf('nomogram_plot/mul_cox.pdf', width = 8, height = 8)
  p2 <- forestplot(mul_results2[,colnames(mul_results2)[c(1:3)]],
                   mean=as.numeric( mul_results2[,"HR"]),
                   lower=as.numeric(mul_results2[,"Lower"]),
                   upper=as.numeric(mul_results2[,"Upper"]),
                   is.summary= is.sum,
                   hrzl_lines=list("1" = gpar(lty=1,lwd=1),
                                   "2" = gpar(lty=2)),
                   txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                  ticks=gpar(cex=1.1),
                                  xlab=gpar(cex = 1.2),
                                  title=gpar(cex = 1.2)),
                   ##fpColors函数设置颜色
                   col=fpColors(box=col[2], lines=col[2], zero = "gray50"),
                   #箱线图中基准线的位置
                   zero=1,
                   cex=0.9, lineheight = "auto",
                   colgap=unit(8,"mm"),xlog=T,
                   #箱子大小，线的宽度
                   lwd.ci=2,boxsize=0.5,
                   #箱线图两端添加小竖线，高度
                   ci.vertices=TRUE, ci.vertices.height = 0.2,
                   graph.pos = 4)
  print(p2)
  dev.off()
  # 开始构建nomogram模型 ------------------------------------------------------------
  message("\n=> Nomogram model plot")
  #数据打包

  Variants <- unique(mul_Variants)

  #拟合模型
  surv_object <- with(mydata, Surv(time, status==1))#构建回归分析因变量

  formula <- as.formula(paste("surv_object ~", paste(Variants, collapse = " + ")))

  model <-  cph(formula, x = TRUE ,y = TRUE, surv = TRUE, data = mydata)

  # 设置不同节点的生存函数 -------------------------------------------------------------
  surv <- Survival(model)#拟合生存函数

  surv1 <- function(x)surv(12*rate.range[1],lp=x)#一年生存函数
  surv2 <- function(x)surv(12*rate.range[2],lp=x)#三年生存函数
  surv3 <- function(x)surv(12*rate.range[3],lp=x)#五年生存函数


  #初始版本的列线图
  Nomogram_1 <- nomogram(model,fun = list(surv1,surv2,surv3),lp=F,#模型，#要放入的生存函数，#风险预测轴，T或F
                         funlabel = c(paste0(rate.range[1],' year survival rate'),paste0(rate.range[2],' years survival rate'),paste0(rate.range[3],' years survival rate')),#风险预测轴的名称
                         maxscale = 100,fun.at = c(0.1,seq(0.1,0.9,by=0.1),0.90))#分数轴总分，#风险预测轴概率取值
  #调整一些参数
  pdf('nomogram_plot/Nomogram_1.pdf', width = 8, height = 8)
  p3 <- plot(Nomogram_1,
             xfrac = .35,#变量与图形的占比（调整变量与坐标抽距离）
             cex.var = 1,#变量字体加粗
             cex.axis = 1,#数轴：字体的大小
             tcl = -0.5,#数轴：刻度的长度
             lmgp = 0.3,#数轴：文字与刻度的距离
             label.every = 2,#数轴：刻度下的文字，1=连续显示，2=隔一个显示一个
             col.grid = gray(c(0.8,0.95)))
  print(p3)
  dev.off()
  # 添加个案的列线图 ----------------------------------------------------------------
  formula_2 <- as.formula(paste("Surv(time, status == 1) ~", paste(Variants, collapse = " + ")))

  Nomogram_2 <- cph(formula_2, x = TRUE ,y = TRUE, surv = TRUE, data = mydata)

  # 使用 regplot 函数创建一个逻辑回归图，其中包括密度图和箱线图
  #pdf('nomogram_plot/Nomogram_2.pdf', width = 8, height = 8)
  #p4 <- regplot(Nomogram_2, # Nomogram_2参数指定使用刚才创建的cox回归模型
  #              plots = c("density","boxes") ,# plots 参数指定要绘制的图形类型
  #              observation = mydata[1,],# 指定在图形中展示的样本，选择第几行
  #              center = T, # center 参数指定是否将图形的中心点设置为零，对齐变量
  #             points = TRUE,# points 参数指定是否绘制每个观测值的点
  #              droplines = F,# droplines 参数指定是否绘制垂直于 x 轴的线
  #              title = "individual nomogram",# title 参数指定图形的标题
  #              odds = F, # odds 参数指定是否显示比值或者概率
  #              interval="confidence",# interval 参数指定使用置信区间还是预测区间
  #              rank = NULL,# rank 参数指定是否对观测值进行排名
  #              clickable = F,# 是否启用交互式模式，这个大家一定要换成T，玩一下
  #              dencol = "#5770A6",# dencol 参数指定密度图的颜色
  #              boxcol = "#5770A6"# boxcol 参数指定箱线图的颜色
  #)
  #print(p4)
  #dev.off()
  # 模型评估与验证 -----------------------------------------------------------------
  # C-INDEX曲线 -------------------------------------------------------------------
  message("\n=> ROC curve")
  pred <- predict(model,mydata,type="lp")# 使用模型预测数据
  colnames(mydata)
  ROC_table <- data.frame(time = mydata[,"time"],status = mydata[,"status"],score = pred)
  # 计算不同时间点的时间C-INDEX曲线
  time_roc_res <- timeROC(T = ROC_table$time,
                          delta = ROC_table$status,
                          marker = ROC_table$score,
                          cause = 1,
                          weighting="marginal",
                          times = c(12*rate.range[1], 12*rate.range[2],12*rate.range[3]),
                          ROC = TRUE,
                          iid = TRUE
  )
  time_ROC_df <- data.frame(TP_1year = time_roc_res$TP[, 1],
                            FP_1year = time_roc_res$FP[, 1],
                            TP_3year = time_roc_res$TP[, 2],
                            FP_3year = time_roc_res$FP[, 2],
                            TP_5year = time_roc_res$TP[, 3],
                            FP_5year = time_roc_res$FP[, 3]
  )

  year1 <- paste0(rate.range[1],' year')
  year2 <- paste0(rate.range[2],' years')
  year3 <- paste0(rate.range[3],' years')

  df <- data.frame(
    Time = rep(c(year1,year2,year3), each=nrow(time_ROC_df)),
    FP = c(time_ROC_df$FP_1year, time_ROC_df$FP_3year, time_ROC_df$FP_5year),
    TP = c(time_ROC_df$TP_1year, time_ROC_df$TP_3year, time_ROC_df$TP_5year)
  )
  pdf('nomogram_plot/time_ROC.pdf', width = 8, height = 8)
  p5 <- ggplot(df, aes(x = FP, y = TP, color = Time)) +
    geom_smooth(se = FALSE, size = 1.2, method = 'loess',na.rm = TRUE) +  # 平滑曲线
    geom_abline(slope = 1, intercept = 0, color = "grey10", linetype = 2) +
    scale_color_manual(values = c("#5770A6", "#A281B1", "#CE5C69"),
                       name = NULL,
                       labels = c(paste0("AUC at ", year1, ": ", round(time_roc_res$AUC[[1]], 2)),
                                  paste0("AUC at ", year2, ": ", round(time_roc_res$AUC[[2]], 2)),
                                  paste0("AUC at ", year3, ": ", round(time_roc_res$AUC[[3]], 2)))) +
    coord_fixed(ratio = 1) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    scale_y_continuous(limits = c(0, 1)) +  # 设置 y 轴范围
    theme_minimal(base_size = 14, base_family = "sans") +
    theme(legend.position = c(0.7, 0.15),
          panel.border = element_rect(fill = NA),
          axis.text = element_text(color = "black"))
  print(p5)
  dev.off()

  # 校准曲线 --------------------------------------------------------------------
  message("\n=> Calibration curve")
  plot_error <- function(x, y, sd, len = 1, col = "black") {
    len <- len * 0.05
    arrows(x0 = x, y0 = y, x1 = x, y1 = y - sd*y, col = col, angle = 90, length = len)
    arrows(x0 = x, y0 = y, x1 = x, y1 = y + sd*y, col = col, angle = 90, length = len)
  }
  pdf('nomogram_plot/cal_curve.pdf', width = 8, height = 8)
  par(mar=c(5,5,2,1))
  plot(x=1,type = "n",
       xlim = c(0,1),
       ylim = c(0,1),
       xaxs = "i",
       yaxs = "i",
       xlab = "Predicted Probability",
       ylab="Observed  Probability",
       legend =FALSE,
       subtitles = FALSE,
       cex=1.5,
       cex.axis=1.5,
       cex.lab=1.5)

  cal.list <- list()
  for (i in seq_along(rate.range)) {
    model.name <- paste0("model", i)
    cal.name <- paste0("cal", i)

    assign(model.name, cph(formula_2, x = TRUE, y = TRUE, surv = TRUE, data = mydata, time.inc = 12 * rate.range[i]))

    assign(cal.name, calibrate(get(model.name),  # 使用动态生成的模型名称
                               cmethod = 'KM',
                               method = 'boot',
                               u = 12 * rate.range[i],
                               m = nrow(cli.data) / 3,
                               B = 200))  # 执行200次自助法计算置信区间
    x<-get(cal.name)[,c("mean.predicted")] #x1表示模型中的预测值
    y<-get(cal.name)[,c("KM")]  #y1为km
    sd<-get(cal.name)[,c("std.err")] #模型中的标准误传递给sd1

    points(x,y,
           type = "o",
           pch = i,
           col = col[i],
           lty = 1,
           lwd = 2)
    plot_error(x,y,
               sd=sd,
               col=col[i])

    cal.list[[cal.name]] <- get(cal.name)
  }
  abline(0,1,lty=3,lwd=2,col="black")
  legend(0.01,0.95,
         legend=c(paste0(year1,' OS'),paste0(year2,' OS'),paste0(year3,' OS')),
         lty = c(1,1),
         lwd = c(2,2),
         pch = c(1,2,3),
         col = col,
         horiz = TRUE,
         bty="n",
         cex=1.2)
  dev.off()

}

