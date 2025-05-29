#' @name subtype_KM
#' @title Comparison of prognosis by Kaplan-Meier survival curve
#' @description This function calculates Kaplan-meier estimator and generate survival curves with log-rank test to detect prognostic difference among identified subtypes. If more than 2 subtypes are identified, pair-wise comparisons will be performed with an additional table printed on the survival curve.
#' @param clust.res.list An object returned by `subtype_harmony()`
#' @param cli.data A data.frame with rownames of observations and with at least two columns of `time` for survival time and `status` for survival status (0: censoring; 1: event)
#' @param clust.col A string vector storing colors for each subtype.
#' @param p.adjust.method A string value for indicating method for adjusting p values (see \link[stats]{p.adjust}). Allowed values include one of c(`holm`, `hochberg`, `hommel`, `bonferroni`, `BH`, `BY`, `fdr`, `none`); "BH" by default.
#' @param file.name A string value to indicate the output path for storing the kaplan-meier curve.
#' @param file.save.path A string value to indicate the name of the kaplan-meier curve.
#' @param surv.median.line A string value for drawing a horizontal/vertical line at median survival. Allowed values include one of c(`none`, `hv`, `h`, `v`). v: vertical, h:horizontal; "none" by default.
#' @return A figure of multi-omics Kaplan-Meier curve (.pdf) and a list with the following components:
#'
#'         \code{fitd}       an object returned by \link[survival]{survdiff}.
#'
#'         \code{fid}        an object returned by \link[survival]{survfit}.
#'
#'         \code{overall.p}  a nominal p.value calculated by Kaplain-Meier estimator with log-rank test.
#'
#'         \code{pairwise.p} an object of class "pairwise.htest" which is a list containing the p values (see \link[survminer]{pairwise_survdiff}); (\emph{only returned when more than 2 subtypes are identified}).
#' @import survival
#' @import survminer
#' @import ggplot2
#' @import ggpmisc
#' @importFrom grDevices pdf dev.off pdf.options
#' @importFrom ggpp geom_table
#' @importFrom tibble tibble
#' @export
#' @examples # There is no example and please refer to vignette.
subtype_KM <- function(clust.res.list   = harmony.res,
                       method           = 'all', #("IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering")
                       cli.data         = cli.data,
                       subtype.list     = c('CS1','CS3'),
                       cli.limit        = NULL,
                       break.time.by    = 12,
                       clust.col        = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB"),
                       p.adjust.method  = "BH",
                       surv.median.line = "none",
                       file.save.path   = getwd(),
                       file.name        = NULL)
  {
  if(!all(is.element(c("time","status"), colnames(cli.data)))) {
    stop("\n=> fail to find variables of time and status.")
  }
if (method=='all') {
  # get common samples
  if (is.na(clust.res.list$clust.res$clust[1])) {message('please input harmoney result')}else{
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                              clust.res.list$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))}
}
if (method=='IntNMF') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$IntNMF$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$IntNMF$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
}
if (method=='CIMLR') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$CIMLR$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$CIMLR$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='PINSPlus') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$PINSPlus$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$PINSPlus$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='NEMO') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$NEMO$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$NEMO$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='COCA') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$COCA$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$COCA$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='MoCluster') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$MoCluster$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$MoCluster$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='LRAcluster') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$LRAcluster$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$LRAcluster$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='iClusterBayes') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$iClusterBayes$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$iClusterBayes$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='SNF') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$SNF$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$SNF$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }
if (method=='ConsensusClustering') {
    # get common samples
    comsam <- intersect(rownames(cli.data),rownames(clust.res.list$ConsensusClustering$clust.res))
    cli.res <- cbind.data.frame(cli.data[comsam,c("time","status")],
                                clust.res.list$ConsensusClustering$clust.res[comsam, "clust", drop = FALSE])
    message(paste0("\n=> a total of ",length(comsam), " samples are identified."))
  }

  subtype.list.num <- as.numeric(gsub("CS", "", subtype.list))
  clust.num <- length(unique(cli.res$clust))
  if (!is.numeric(subtype.list.num) || length(subtype.list.num) != 2) {
    stop("\n=> The 'subtype.list' must consist of two subtypes, such as CS1 and CS2.")
  }
  if (max(subtype.list.num) >max(clust.num)) {
    stop("\n=> The 'subtype.list' must be less than the maximum number of subtypes.")
  }
    cli.res <- subset(cli.res, clust %in% subtype.list.num)
    clust.col <- clust.col[subtype.list.num]

  # remove missing data if possible
  if(sum(c(is.na(cli.res$time), is.na(cli.res$status))) > 0) {
    message("\n=> removed missing values.")
    cli.res <- as.data.frame(na.omit(cli.res))
    message(paste0("\n=> leaving ",nrow(cli.res), " observations."))
  }

  if(is.null(cli.limit)) {
    xlim = c(0, max(cli.res$time))
  } else {
    message(paste0("\n=> cut survival curve up to ",cli.limit," Months"))
    xlim = c(0, cli.limit)
  }

  # basic survival analysis
  fit <- survfit(Surv(time, status)~ clust,
                 data      = cli.res,
                 type      = "kaplan-meier",
                 error     = "greenwood",
                 conf.type = "plain",
                 na.action = na.exclude)

  # hack strata for better survival curve
  names(fit$strata) <- gsub("Subtype=", "", names(fit$strata))
  # kaplan-meier curve
  p <- suppressWarnings(ggsurvplot(fit,
             data = cli.res,
             surv.median.line = surv.median.line, # Add medians survival
             size = 1, # change line size
             cex.lab=2,
             break.time.by = 36, # break X axis in time intervals by 500.
             xlim = xlim,
             axis.title.x =element_text(size=5),
             axis.title.y = element_text(size=5),
             palette = clust.col[1:clust.num],# custom color palettes
             conf.int = T, # Add confidence interval
             pval = TRUE, # Add p-value
             risk.table = T, # Add risk table
             xlab = "Follow-up months", # customize X axis label.
             ylab="Survival probability",
             risk.table.col = "strata",# Risk table color by groups
             risk.table.height = 0.3, # Useful to change when you have multiple groups
             ggtheme = theme_bw()))# Add title

  # output curve to pdf
  if(!is.null(file.name)) {
    outFile <- file.path(file.save.path,paste0(file.name,".pdf"))
  } else {
    outFile <- file.path(file.save.path,paste0("km_curve_",harmony.res$method,".pdf"))
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(outFile, width = 8, height = 8)
  # ggsave(outFile, width = 6, height = 7)
  print(p)
  dev.off()

  # output curve to screen
  print(p)

  return(list(fit = fit,  overall.p = p.val))
}
