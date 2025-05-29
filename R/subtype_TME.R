#' @name subtype_TME
#' @import xCell
#' @import MCPcounter
#' @import EPIC
#' @import CIBERSORT
#' @import estimate
#' @import TIMER
#' @import quantiseqr
#' #' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `xCell`, `MCPcounter`, `EPIC`,`CIBERSORT`, `estimate`, `TIMER`, `ips`, `quantiseq`, `svr`,`lsei`.
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#'#' @references
#' 1. Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457.
#' 2. Vegesna R, Kim H, Torres-Garcia W, …, Verhaak R. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4, 2612.
#' 3. Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., …, Sopper, S. (2019). Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34.
#' 4. Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174.
#' 5. P. Charoentong et al., Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Reports 18, 248-262 (2017).
#' 6. Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218.
#' 7. Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220.
#' 8. Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476.
#'
#' @author Junxiao Shen
#' @export
subtype_TME <- function(clust.res.list = clust.res.list,
                        clust.method ='iClusterBayes',#("IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering")
                        expr = expr, # 基因表达矩阵，可以是一个数值型矩阵或数据框，行名为基因，列名为样本名
                        subtype.list = c('CS1','CS2'),
                        method = "CIBERSORT", # c('xCell','MCPcounter', 'EPIC',  'CIBERSORT', 'ips', 'estimate', 'svr', 'lsei', 'TIMER', 'quantiseq')
                        perm = 100, # 统计分析的置换次数（建议≥100次）影响'CIBERSORT'方法。
                        reference=NA,
                        indications=NULL)
{
  message(paste0("\n=> Select the ",subtype.list[1], ' and ' ,subtype.list[2]))
  subtype.list.num <- as.numeric(gsub("CS", "", subtype.list))

  if (clust.method=='all') {
    Data <-  subset(clust.res.list$clust.res,clust.res.list$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the comprehensive clustering machine learning algorithm."))}
  if (clust.method=='IntNMF') {
    Data <-  subset(clust.res.list$IntNMF$clust.res,clust.res.list$IntNMF$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='CIMLR') {
    Data <-  subset(clust.res.list$CIMLR$clust.res,clust.res.list$CIMLR$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='PINSPlus') {
    Data <-  subset(clust.res.list$PINSPlus$clust.res,clust.res.list$PINSPlus$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='NEMO') {
    Data <-  subset(clust.res.list$NEMO$clust.res,clust.res.list$NEMO$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='COCA') {
    Data <-  subset(clust.res.list$COCA$clust.res,clust.res.list$COCA$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='MoCluster') {
    Data <-  subset(clust.res.list$MoCluster$clust.res,clust.res.list$MoCluster$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='LRAcluster') {
    Data <-  subset(clust.res.list$LRAcluster$clust.res,clust.res.list$LRAcluster$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='iClusterBayes') {
    Data <-  subset(clust.res.list$iClusterBayes$clust.res,clust.res.list$iClusterBayes$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='SNF') {
    Data <-  subset(clust.res.list$SNF$clust.res,clust.res.list$SNF$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}
  if (clust.method=='ConsensusClustering'){
    Data <-  subset(clust.res.list$ConsensusClustering$clust.res,clust.res.list$ConsensusClustering$clust.res$clust %in% subtype.list.num)
    message(paste0("\n=> Select the ",clust.method , " clustering machine learning algorithm."))}

  expr <- t(expr)
  expr <- subset(expr,row.names(expr) %in% Data$samID)
  expr <- t(expr)

  TME.harmony <- function(res){
    res<-as.data.frame(t(res))
    colnames(res)<-gsub(colnames(res),pattern = "\\ ",replacement = "\\.")
    colnames(res)<-paste0(colnames(res),".",method)
    res$ID<-rownames(res)
    res <- res[, c("ID", setdiff(names(res), "ID"))]
    return(res)
  }

  #' Deconvolute using CIBERSORT
  if (method=='CIBERSORT'){
    library(CIBERSORT)
    library(limma)
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
    data('CIBERSORT_LM22')
    CIBERSORT_LM22 <- as.matrix(CIBERSORT_LM22)
    expr=normalizeBetweenArrays(expr)
    res<-cibersort(sig_matrix = CIBERSORT_LM22,
                   mixture_file = expr,
                   perm = perm,
                   QN = T
    )
    res<-as.data.frame(t(res))
    res <- TME.harmony(res = res)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(res)

  }
  #' Deconvolute using svr
  if (method=="svr"){
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))

    expr<-as.data.frame(expr)
    res<-CIBERSORT(sig_matrix = reference,
                   mixture_file = expr,
                   perm = perm,
                   QN = T,
                   absolute = F,)
    res<-as.data.frame(t(res))
    res <- TME.harmony(res = res)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(res)

  }
  #' Deconvolute using lsei
  if (method=="lsei"){
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
    expr <- as.matrix(expr)
    reference <- as.matrix(reference)
    # scale_reference
    scale_reference <- T
    if (scale_reference){
      reference <- (reference - mean(reference))/sd(as.vector(reference))
    }
    Ymedian <- max(median(expr),1)
    # common expr
    common <- intersect(rownames(expr), rownames(reference))
    expr <- expr[match(common, rownames(expr)), ]
    reference <- reference[match(common, rownames(reference)), ]
    # deconvolution
    output <- matrix()
    Numofx <- ncol(reference)
    AA <- reference
    EE <- rep(1, Numofx)
    FF <- 1
    GG <- diag(nrow=Numofx)
    HH <- rep(0, Numofx)
    out.all <- c()
    itor <- 1
    samples <- ncol(expr)
    while (itor <= samples){
      BB <- expr[, itor]
      BB <- (BB - mean(BB))/sd(BB)
      out <- lsei(AA, BB, EE, FF, GG, HH)
      out.all <- rbind(out.all, out$X)
      itor <- itor + 1
    }
    rownames(out.all) <- colnames(expr)
    res <- out.all
    res<-as.data.frame(t(res))
    res <- TME.harmony(res = res)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(res)

  }
  #' Deconvolute using xCell
  if (method=='xCell') {
    library(xCell)
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
    expr <- log2(expr+1)
    res<- xCellAnalysis(expr,rnaseq = T)
    res <- TME.harmony(res = res)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(res)

  }
  #' Deconvolute using MCPcounter
  if (method=='MCPcounter'){
    library(MCPcounter)
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
    expr <- log2(expr+1)
    data('MCPcounter_probesets');data("MCPcounter_genes")
    res<-MCPcounter.estimate(expr,
                             featuresType = "HUGO_symbols",
                             probesets= MCPcounter_probesets,
                             genes= MCPcounter_genes)
    res <- TME.harmony(res = res)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(res)

  }
  #' Deconvolute using EPIC
  if (method=='EPIC'){
    library(EPIC)
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
    expr <- log2(expr+1)
    data('TRef')
    out <- EPIC(bulk = expr, reference = TRef)
    res<-as.data.frame(t(out$cellFractions))
    res <- TME.harmony(res = res)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(res)

  }
  #' Deconvolute using TIMER
  if (method=="TIMER"){
    if (length(indications) == 0) {
      warning("indications is empty. Skipping this TIMER method.")
      return(NULL)  # 或者使用 next 来跳过
    }
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
    indications <- tolower(indications)
    assert("indications fit to mixture matrix", length(indications) == ncol(expr))
    args <- new.env()
    args$outdir <- tempdir()
    args$batch <- tempfile()
    lapply(unique(indications), function(ind) {
      tmp_file <- tempfile()
      tmp_mat <- expr[, indications == ind, drop = FALSE] %>% as_tibble(rownames = "gene_symbol")
      write_tsv(tmp_mat, tmp_file)
      cat(paste0(tmp_file, ",", ind, "\n"), file = args$batch, append = TRUE)
    })
    # reorder results to be consistent with input matrix
    results <- deconvolute_timer.default(args)[, make.names(colnames(expr))]

    colnames(results) <- colnames(expr)
    results <- TME.harmony(results)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(results)

  }
  #' Deconvolute using ips
  if (method=='ips'){
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
      sample_names<-colnames(expr)
      ####################################################
      ##
      ##   This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram from "EXPR.txt" and "IPS_genes.txt"
      ##   (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
      ##   Version 1.0 08.07.2016
      ##   Needs packages ggplot2,grid,gridExtra
      ##
      ####################################################
      ################################################
      ipsmap<- function (x) {
        # if(is.na(x)) x<-0
        if (x<=0) {
          ips<-0
        } else {
          if (x>=3) {
            ips<-10
          } else {
            ips<-round(x*10/3, digits=0)
          }
        }
        return(ips)
      }
      mapcolors<-function (x) {
        za<-NULL
        if (x>=3) {
          za=1000
        } else {
          if (x<=-3) {
            za=1
          } else {
            za=round(166.5*x+500.5,digits=0)
          }
        }
        return(my_palette[za])
      }
      mapbw<-function (x) {
        za2<-NULL
        if (x>=2) {
          za2=1000
        } else {
          if (x<=-2) {
            za2=1
          } else {
            za2=round(249.75*x+500.5,digits=0)
          }
        }
        return(my_palette2[za2])
      }
      ## Assign colors
      my_palette <- colorRampPalette(c("#5770A6", "white", "#CE5C69"))(n = 1000)
      my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)

      data('ips_gene_set')
      IPSG<-ips_gene_set

      IPSG<-IPSG[IPSG$GENE%in%rownames(expr),]
      unique_ips_genes<-as.vector(unique(IPSG$NAME))

      IPS<-NULL
      MHC<-NULL
      CP<-NULL
      EC<-NULL
      SC<-NULL
      AZ<-NULL

      # Gene names in expression file
      GVEC<-row.names(expr)
      # Genes names in IPS genes file
      VEC<-as.vector(IPSG$GENE)
      # Match IPS genes with genes in expression file
      ind<-which(is.na(match(VEC,GVEC)))
      # List genes missing or differently named
      MISSING_GENES<-VEC[ind]
      dat<-IPSG[ind,]
      if (length(MISSING_GENES)>0) {
        cat("differently named or missing genes: ",MISSING_GENES,"\n")
      }
      # for (x in 1:length(ind)) {
      #   print(IPSG[ind,])
      # }
      for (i in 1:ncol(expr)) {
        GE<-expr[,i]
        mGE<-mean(GE,na.rm=TRUE)
        sGE<-sd(GE,na.rm=TRUE)
        Z1<-(expr[as.vector(IPSG$GENE),i]-mGE)/sGE
        W1<-IPSG$WEIGHT
        WEIGHT<-NULL
        MIG<-NULL
        k<-1
        for (gen in unique_ips_genes) {
          MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
          WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
          k<-k+1
        }
        WG<-MIG*WEIGHT
        MHC[i]<-mean(WG[1:10],na.rm=TRUE)
        CP[i]<-mean(WG[11:20],na.rm=TRUE)
        EC[i]<-mean(WG[21:24],na.rm=TRUE)
        SC[i]<-mean(WG[25:26],na.rm=TRUE)
        AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i],na.rm = TRUE)
        # print(paste0(">>> Processing sample ", i))
        # if(is.na(AZ[i])) {
        #   print(paste0(">>> ", i," Sample with error"))
        #   AZ[i]<-0
        # }
        IPS[i]<-ipsmap(AZ[i])
        res<-data.frame(ID=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
        res <- res[, c("ID", setdiff(names(res), "ID"))]
        message(paste0("\n=> ",method , " done."))
        message('file saved.')
        return(res)
      }

  }
  #' Deconvolute using estimate
  if (method=='estimate'){

    if (!requireNamespace("estimate", quietly = TRUE)) {install.packages("estimate", repos = "http://R-Forge.R-project.org")}
    library(estimate)
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))
    expr<-as.data.frame(expr)
    expr<-tibble::rownames_to_column(expr,var = "symbol")
    sampleData<-paste0(NA,"-expr.txt")
    write.table(expr,sampleData,sep = "\t",row.names = F,quote = F)
    filterCommonGenes(input.f= sampleData,
                      output.f= paste0(NA,"_Tumor_purity.gct"),
                      id="GeneSymbol")
    file.remove(paste0(NA,"-expr.txt"))
    estimateScore(input.ds = paste0(NA,"_Tumor_purity.gct"),
                  output.ds= paste0(NA,"_Tumor_estimate_score.gct"),
                  platform= "affymetrix")
    file.remove(paste0(NA,"_Tumor_purity.gct"))
    scores=read.table(paste0(NA,"_Tumor_estimate_score.gct"),skip = 2,header = T)
    file.remove(paste0(NA,"_Tumor_estimate_score.gct"))
    rownames(scores)=scores[,1]
    scores=t(scores[,3:ncol(scores)])
    colnames(scores)<-paste0(colnames(scores),"_estimate")
    scores<-tibble::rownames_to_column(as.data.frame(scores),var = "ID")
    scores$ID<-gsub(scores$ID,pattern = "\\.",replacement = "-")
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(scores)

  }
  #' Deconvolute using quantiseq
  if (method=="quantiseq"){
    if (!requireNamespace("quantiseqr", quietly = TRUE)) {BiocManager::install("quantiseqr")}
    library(quantiseqr)
    message(paste0("\n=> Select the ",method , " deconvolution machine learning algorithm."))

    res <- run_quantiseq(
      expression_data = expr,
      is_arraydata = F,
      is_tumordata = T,
      scale_mRNA = F)
    res <- res[,-1]
    res<-as.data.frame(t(res))
    res  <- TME.harmony(res)
    message(paste0("\n=> ",method , " done."))
    message('file saved.')
    return(res)

  }

  if (method=='all') {
    result_list <- list()

    result_list[['CIBERSORT']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr, subtype.list = subtype.list,method='CIBERSORT', perm=perm)
    result_list[['xCell']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr, subtype.list = subtype.list, method='xCell')
    result_list[['MCPcounter']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr, subtype.list = subtype.list,method='MCPcounter')
    result_list[['EPIC']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr, subtype.list = subtype.list,method='EPIC')
    result_list[['TIMER']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr, subtype.list = subtype.list,method='TIMER')
    result_list[['ips']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr,subtype.list = subtype.list, method='ips')
    result_list[['estimate']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr, subtype.list = subtype.list ,method='estimate')
    result_list[['quantiseq']] <- subtype_TME(clust.res.list=clust.res.list,clust.method=clust.method,expr, subtype.list =subtype.list ,method='quantiseq')

    return(result_list)
  }
}


