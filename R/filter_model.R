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

filter_model <- function(model1=best.model.res1,
                         model2=best.model.res2,
                         model1.name='CoxBoost',
                         model2.name='RSF',
                         rs.data=best.model.rs,
                         col=c('#5770A6','#CE5C69'))
{
  library(survival)
  library(survminer)
  library(ggplot2)
  library(patchwork)
  library(ggrisk)
  library(pROC)
  if (!dir.exists('ML_plot')) {
    dir.create('ML_plot')
  }
if (model1.name=='RSF') {
  pdf(paste0('ML_plot/',model1.name,"_plot.pdf"),width = 8, height = 8)
  plot(model1)
  dev.off()
}

  if (model2.name=='RSF') {
    pdf(paste0('ML_plot/',model2.name,"_plot.pdf"),width = 8, height = 8)
    plot(model2)
    dev.off()
  }
  ##################### ggrisk ####################

  for (i in 1:length(rs.data)) {
    message(paste0('ggrisk analysis dataset ', i))

    rs.dt <- rs.data[[i]]

    if (length(unique(rs.dt$status)) < 2) {
      warning(paste("Dataset", names(rs.data)[[i]], "does not have two unique status values. Skipping."))
      next  # 跳过当前迭代
    }

    formula <- as.formula(paste("Surv(time, status) ~", 'RS'))
    cox.model <- coxph(formula, data = rs.dt)
    pdf(paste0('ML_plot/', names(rs.data)[[i]], "_RS_plot.pdf"), width = 8, height = 8)
    two_scatter(cox.model,
                cutoff.value = 'roc',
                code.0 = 'Still Alive',
                code.1 = 'Already Dead',
                code.highrisk = 'High Risk',
                code.lowrisk = 'Low Risk',
                title.A.ylab = 'Risk Score',
                title.B.ylab = 'OS (months)',
                title.A.legend = 'Risk Group',
                title.B.legend = 'Status',
                size.ylab.title = 14,
                size.Atext = 11,
                size.Btext = 11,
                size.points = 2,
                size.dashline = 1,
                size.cutoff = 5,
                size.legendtitle = 13,
                size.legendtext = 12,
                color.A = c(low = col[1], high = col[2]),
                color.B = c(code.0 = col[1], code.1 = col[2]),
                vjust.A.ylab = 1,
                vjust.B.ylab = 2,
                family = 'sans',
                expand.x = 3
    )
    dev.off()
  }

    ##################### ROC #######################
  for (i in 1:length(rs.data)) {
    message(paste0('ROC analysis dataset ', i))
    rs.dt <- rs.data[[i]]
    if (length(unique(rs.dt$status)) < 2) {
      warning(paste("Dataset", names(rs.data)[[i]], "does not have two unique status values. Skipping."))
      next
    }
    pdf(paste0('ML_plot/',names(rs.data)[[i]],"_roc_plot.pdf"), width = 8, height = 8)
    roc.res<- roc(rs.dt$status, rs.dt$RS)
    plot(smooth(roc.res), col="#CE5C69", legacy.axes=T, print.auc=T, print.thres=T)##显示AUC面积

    dev.off()
}
    ##################### KM #######################
  for (i in 1:length(rs.data)) {
      message(paste0('KM analysis dataset ', i))
    rs.dt <- rs.data[[i]]
    if (length(unique(rs.dt$status)) < 2) {
      warning(paste("Dataset", names(rs.data)[[i]], "does not have two unique status values. Skipping."))
      next
    }
    res.cut <- surv_cutpoint(rs.dt, time = "time", event = "status",variables = 'RS')
    highRS.data <- rs.dt[rs.dt$RS >= res.cut[["cutpoint"]][["cutpoint"]],]
    lowRS.data<- rs.dt[rs.dt$RS<res.cut[["cutpoint"]][["cutpoint"]],]
    highRS.data$group <- '1'
    lowRS.data$group <- '0'
    rs.dt <- rbind(highRS.data,lowRS.data)
    rs.dt$group <- as.factor(rs.dt$group)
    fit<-survfit(Surv(time,status)~group,data = rs.dt)

   km.plot <-  ggsurvplot(fit,
               data = rs.dt,
               surv.median.line = "hv", # Add medians survival
               size = 1, # change line size
               cex.lab=2,
               break.time.by = floor(max(rs.dt$time)/10), # break X axis in time intervals by 500.
               xlim = c(0,max(rs.dt$time)),
               axis.title.x =element_text(size=5),
               axis.title.y = element_text(size=5),
               palette = col,# custom color palettes
               conf.int = T, # Add confidence interval
               pval = TRUE, # Add p-value
               risk.table = T, # Add risk table
               xlab = "Follow-up months", # customize X axis label.
               ylab="Survival probability ",
               risk.table.col = "strata",# Risk table color by groups
               legend.labs =  c("Low risk","High risk"), # Change legend labels
               risk.table.height = 0.3, # Useful to change when you have multiple groups
               ggtheme = theme_bw())# Add title

    pdf(paste0('ML_plot/',names(rs.data)[[i]],"_KM_plot.pdf"), width = 8, height = 8)
    print(km.plot)
    dev.off()
  }
}
