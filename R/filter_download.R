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
filter_Download<-function(project='',
                          type='mrna')
{
  library(TCGAbiolinks)
  library(easyTCGA)
  if(type=='mrna'){getmrnaexpr(project)}
  if(type=='mirna'){getmirnaexpr(project)}
  if(type=='met'){query.met<- GDCquery(project=project,
                                     data.category ="DNA Methylation",
                                     data.type='Methylation Beta Value',
                                     platform ="Illumina Human Methylation 450")
  GDCprepare(
    query = query.met,save = TRUE,save.filename =paste0(project,'_DNAmet.rda'),
    summarizedExperiment=TRUE)}

}
