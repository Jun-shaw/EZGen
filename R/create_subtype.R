#' @name create_subtype
#' @title Utilize 10 clustering machine learning algorithms for cluster analysis
#' @param data Matrix file containing two sets of genes. .
#' @param clust.num Number of clusters.
#' @examples # There is no example and please refer to vignette.
#' @export
#' @return A list of results returned by each specified algorithms.
#' @import IntNMF
#' @import CIMLR
#' @import PINSPlus
#' @import coca
#' @import mogsa
#' @import iClusterPlus
#' @import SNFtool
#' @import ConsensusClusterPlus
#' @importFrom dplyr %>%
#' @importFrom PINSPlus PerturbationClustering
#' @importFrom vegan vegdist
#' @details
#'
#' @references
#' Pierre-Jean M, Deleuze J F, Le Floch E, et al. Clustering and variable selection evaluation of 13 unsupervised methods for multi-omics data integration[J]. Briefings in Bioinformatics, 2019.
#'
#' intNMF:
#' Chalise P, Fridley BL. Integrative clustering of multi-level omic data based on non-negative matrix factorization algorithm. PLoS One. 2017;12(5):e0176278.
#'
#' CIMLR:
#' Ramazzotti D, Lal A, Wang B, Batzoglou S, Sidow A. Multi-omic tumor data reveal diversity of molecular mechanisms that correlate with survival. Nat Commun. 2018;9(1):4453.
#'
#' PINSPlus:
#' Nguyen H, Shrestha S, Draghici S, Nguyen T. PINSPlus: a tool for tumor subtype discovery in integrated genomic data. Bioinformatics. 2019;35(16):2843-2846.
#'
#' NEMO:
#' Rappoport N, Shamir R. NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics. 2019;35(18):3348-3356.
#'
#' COCA:
#' Hoadley KA, Yau C, Wolf DM, et al. Multiplatform analysis of 12 cancer types reveals molecular classification within and across tissues of origin. Cell. 2014;158(4):929-944.
#'
#' Mocluster:
#' Meng C, Helm D, Frejno M, Kuster B. moCluster: Identifying Joint Patterns Across Multiple Omics Data Sets. J Proteome Res. 2016;15(3):755-765.
#'
#' LRAcluster:
#' Wu D, Wang D, Zhang MQ, Gu J. Fast dimension reduction and integrative clustering of multi-omics data using low-rank approximation: application to cancer molecular classification. BMC Genomics. 2015;16:1022.
#'
#' iClusterBayes:
#' Mo Q, Shen R, Guo C, Vannucci M, Chan KS, Hilsenbeck SG. A fully Bayesian latent variable model for integrative clustering analysis of multi-type omics data. Biostatistics. 2018;19(1):71-86.
#'
#' SNF:
#' Wang B, Mezlini AM, Demir F, et al. Similarity network fusion for aggregating data types on a genomic scale. Nat Methods. 2014;11(3):333-337.
#'
#' ConsensusClustering:
#' Monti S, Tamayo P, Mesirov J, et al. Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. Machine Learning. 2003;52:91-118.
#'
#'
#'
create_subtype <- function(data= data.list,
                           clust.num=opt.res$N.clust,
                           methodslist = list("IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering"))
{
  type        = rep("gaussian", length(data))
  clust.res.list <- list()
  for (method in unlist(methodslist)){

    dat <- lapply(data, function (dd){
      if (!all(dd >= 0)) dd <- pmax(dd + abs(min(dd)), 0) + .Machine$double.eps
      dd <- dd/max(dd)
      return(dd %>% as.matrix)
    })

    if (method=="IntNMF"){
      library(dplyr)
      message(paste0(method," start..."))
      dat <- lapply(data, function (dd){
        if (!all(dd >= 0)) dd <- pmax(dd + abs(min(dd)), 0) + .Machine$double.eps
        dd <- dd/max(dd)
        return(dd %>% as.matrix)
      })

      dat <- lapply(dat, function(x) t(x) + .Machine$double.eps)

      result.intNMF <- dat %>% IntNMF::nmf.mnnals(k = clust.num)
      clust.intNMF <- result.intNMF$clusters

      clustres <- data.frame(samID = colnames(data[[1]]),
                             clust = as.numeric(clust.intNMF),
                             row.names = colnames(data[[1]]),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]
      clust.res.list[[method]] <- list(fit = result.intNMF, clust.res = clustres, method = "IntNMF")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))
    }

    if (method=="CIMLR"){
      library(CIMLR)
      message(paste0(method," start..."))
      fit <- CIMLR(data,
                   c= clust.num,
                   cores.ratio = 0)
      input_dat <- do.call(rbind,lapply(seq(along = data), function(x){
        ddd <- data[[x]]
        rownames(ddd) <- paste(rownames(ddd), names(data)[x], sep = "+")
        ddd
      }))
      ranks <- CIMLR_Feature_Ranking(A = fit$S, X = input_dat)
      ranks$names <- rownames(input_dat)[ranks$aggR]
      fit$selectfeatures <- ranks
      clustres <- data.frame(samID = colnames(data[[1]]),
                             clust = fit$y$cluster,
                             row.names = colnames(data[[1]]),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]

      message("clustering done...")

      f <- sapply(strsplit(ranks$name, "+",fixed = TRUE), "[",1)
      d <- sapply(strsplit(ranks$name, "+",fixed = TRUE), "[",2)

      featres <- data.frame(feature = f,
                            dataset = d,
                            pvalue = ranks$pval,
                            stringsAsFactors = FALSE)
      feat.res <- NULL
      for (d in unique(featres$dataset)) {
        tmp <- featres[which(featres$dataset == d),]
        tmp <- tmp[order(tmp$pvalue, decreasing = FALSE),]
        tmp$rank <- 1:nrow(tmp)
        feat.res <- rbind.data.frame(feat.res,tmp)
      }
      message("feature selection done...")
      clust.res.list[[method]] <-  list(fit = fit, clust.res = clustres, feat.res = feat.res, method = "CIMLR")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))
    }

    if (method=="PINSPlus"){
      library(PINSPlus)
      message(paste0(method," start..."))
      d <- do.call(rbind, data)
      df <- t(d)
      fit <- PerturbationClustering(data             = df,
                                    kMin             = clust.num,
                                    kMax             = clust.num,
                                    clusteringMethod = 'kmeans',#c("kmeans", "hclust", "pam")
                                    iterMin          = 50,
                                    iterMax          = 500,
                                    verbose          = TRUE)

      clustres <- data.frame(samID = rownames(df),
                             clust = fit$cluster,
                             row.names = rownames(df),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]
      clust.res.list[[method]] <- list(fit = fit, clust.res = clustres, method = "PINSPlus")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))
    }

    if (method=="NEMO"){
      affinityMatrix <- function(Diff,K=20,sigma=0.5) {
        N = nrow(Diff)
        Diff = (Diff + t(Diff)) / 2
        diag(Diff) = 0;
        sortedColumns = as.matrix(t(apply(Diff,2,sort)))
        finiteMean <- function(x) { mean(x[is.finite(x)]) }
        means = apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;

        avg <- function(x,y) ((x+y)/2)
        Sig = outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
        Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
        densities = dnorm(Diff,0,sigma*Sig,log = FALSE)

        W = (densities + t(densities)) / 2
        return(W)
      }

      dist2 <- function(X,C) {
        ndata = nrow(X)
        ncentres = nrow(C)
        sumsqX = rowSums(X^2)
        sumsqC = rowSums(C^2)
        XC = 2 * (X %*% t(C))

        res = matrix(rep(sumsqX,times=ncentres),ndata,ncentres) + t(matrix(rep(sumsqC,times=ndata),ncentres,ndata)) - XC
        res[res < 0] = 0
        return(res)
      }

      spectralClustering = SNFtool::spectralClustering
      nemo.num.clusters <- function(W, NUMC=2:15) {
        if (min(NUMC) == 1) {
          warning("Note that we always assume there are more than one cluster.")
          NUMC = NUMC[NUMC > 1]
        }
        W = (W + t(W))/2
        diag(W) = 0
        if (length(NUMC) > 0) {
          degs = rowSums(W)
          degs[degs == 0] = .Machine$double.eps
          D = diag(degs)
          L = D - W
          Di = diag(1/sqrt(degs))
          L = Di %*% L %*% Di
          print(dim(L))
          eigs = eigen(L)
          eigs_order = sort(eigs$values, index.return = TRUE)$ix
          eigs$values = eigs$values[eigs_order]
          eigs$vectors = eigs$vectors[, eigs_order]
          eigengap = abs(diff(eigs$values))
          eigengap = (1:length(eigengap)) * eigengap

          t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = TRUE)$ix
          return(NUMC[t1[1]])
        }
      }

      nemo.affinity.graph <- function(raw.data, k = NA, NUM.NEIGHBORS.RATIO = 6) {
        if (is.na(k)) {
          k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]]) / NUM.NEIGHBORS.RATIO)))
        } else if (length(k) == 1) {
          k = rep(k, length(raw.data))
        }
        sim.data = lapply(1:length(raw.data), function(i) {affinityMatrix(dist2(as.matrix(t(raw.data[[i]])),
                                                                                as.matrix(t(raw.data[[i]]))), k[i], 0.5)})
        affinity.per.omic = lapply(1:length(raw.data), function(i) {
          sim.datum = sim.data[[i]]
          non.sym.knn = apply(sim.datum, 1, function(sim.row) {
            returned.row = sim.row
            threshold = sort(sim.row, decreasing = TRUE)[k[i]]
            returned.row[sim.row < threshold] = 0
            row.sum = sum(returned.row)
            returned.row[sim.row >= threshold] = returned.row[sim.row >= threshold] / row.sum
            return(returned.row)
          })
          sym.knn = non.sym.knn + t(non.sym.knn)
          return(sym.knn)
        })
        patient.names = Reduce(union, lapply(raw.data, colnames))
        num.patients = length(patient.names)
        returned.affinity.matrix = matrix(0, ncol = num.patients, nrow=num.patients)
        rownames(returned.affinity.matrix) = patient.names
        colnames(returned.affinity.matrix) = patient.names

        shared.omic.count = matrix(0, ncol = num.patients, nrow=num.patients)
        rownames(shared.omic.count) = patient.names
        colnames(shared.omic.count) = patient.names

        for (j in 1:length(raw.data)) {
          curr.omic.patients = colnames(raw.data[[j]])
          returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, curr.omic.patients]
          shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
        }

        final.ret = returned.affinity.matrix / shared.omic.count
        lower.tri.ret = final.ret[lower.tri(final.ret)]
        final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])

        return(final.ret)
      }

      nemo.clustering <- function(omics.list, num.clusters = NULL, num.neighbors = NA) {
        if (is.null(num.clusters)) {
          num.clusters = NA
        }

        graph = nemo.affinity.graph(omics.list, k = num.neighbors)
        if (is.na(num.clusters)) {
          num.clusters = nemo.num.clusters(graph)
        }
        clustering = spectralClustering(graph, num.clusters)
        names(clustering) = colnames(graph)
        return(clustering)
      }
      message(paste0(method," start..."))
      fit <- nemo.clustering(omics.list    = data,
                             num.clusters  = clust.num,
                             num.neighbors = NA)

      clustres <- data.frame(samID = colnames(data[[1]]),
                             clust = as.numeric(fit),
                             row.names = colnames(data[[1]]),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]
      clust.res.list[[method]] <- list(fit = fit, clust.res = clustres, method = "NEMO")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))
    }

    if (method=="COCA"){
      library(vegan)
      message(paste0(method," start..."))
      dt <- lapply(data, t)
      outputBuildMOC <- coca::buildMOC(dt,
                                       M         = length(dt),
                                       K         = clust.num,
                                       methods   = 'hclust',
                                       distances = 'euclidean')

      moc <- outputBuildMOC$moc
      datasetIndicator <- outputBuildMOC$datasetIndicator

      hcs <- hclust(vegdist(as.matrix(moc), method = "jaccard"), "ward.D")
      coca <- cutree(hcs,clust.num)

      clustres <- data.frame(samID = rownames(dt[[1]]),
                             clust = as.numeric(coca),
                             row.names = rownames(dt[[1]]),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]
      clust.res.list[[method]] <- list(fit = outputBuildMOC, clust.res = clustres, clust.dend = hcs, method = "COCA")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))
    }

    if (method=="MoCluster"){
      library(mogsa)
      message(paste0(method," start..."))
      moas <- data %>% mogsa::mbpca(ncomp      = clust.num,
                                    k          = 10,
                                    method     = 'globalScore',
                                    option     = "lambda1",
                                    center     = T,
                                    scale      = T,
                                    moa        = TRUE,
                                    svd.solver = "fast",
                                    maxiter    = 1000,
                                    verbose    = FALSE)

      scrs <- moas %>% moaScore
      dist <- scrs %>% dist
      clust.dend <- hclust(dist, method = "ward.D")

      clustres <- data.frame(samID = colnames(data[[1]]),
                             clust = cutree(clust.dend,k = clust.num),
                             row.names = colnames(data[[1]]),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]
      message("clustering done...")

      featres <- moas@loading[which(moas@loading[,1] != 0),]
      f <- sub('_[^_]*$', '', rownames(featres))
      d <- sub('.*_', '', rownames(featres))
      featres <- data.frame(feature = f,
                            dataset = d,
                            load = featres[,1],
                            stringsAsFactors = FALSE)
      feat.res <- NULL
      for (d in unique(featres$dataset)) {
        tmp <- featres[which(featres$dataset == d),]
        feat.res <- rbind.data.frame(feat.res,tmp)
      }
      message("feature selection done...")
      clust.res.list[[method]] <- list(fit = moas, clust.res = clustres, feat.res = feat.res, clust.dend = clust.dend, method = "MoCluster")
      message(paste0(method," done..."))
    }

    if (method=="LRAcluster"){
      message(paste0(method," start..."))
      LRAcluster <- function(data, types, dimension = 2, names = as.character(1:length(data)))
      {
        #--------#
        # binary #
        #--------#
        epsilon.binary<-2.0
        check.binary.row<-function(arr)
        {
          if (sum(!is.na(arr))==0)
          {
            return (F)
          }
          else
          {
            idx<-!is.na(arr)
            if (sum(arr[idx])==0 || sum(arr[idx])==sum(idx))
            {
              return (F)
            }
            else
            {
              return (T)
            }
          }
        }
        check.binary<-function(mat,name)
        {
          index<-apply(mat,1,check.binary.row)
          n<-sum(!index)
          if (n>0)
          {
            w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
            warning(w)
          }
          mat_c<-mat[index,]
          rownames(mat_c)<-rownames(mat)[index]
          colnames(mat_c)<-colnames(mat)
          return (mat_c)
        }

        base.binary.row<-function(arr)
        {
          idx<-!is.na(arr)
          n<-sum(idx)
          m<-sum(arr[idx])
          return (log(m/(n-m)))
        }

        base.binary<-function(mat)
        {
          mat_b<-matrix(0,nrow(mat),ncol(mat))
          ar_b<-apply(mat,1,base.binary.row)
          mat_b[1:nrow(mat_b),]<-ar_b
          return (mat_b)
        }

        update.binary<-function(mat,mat_b,mat_now,eps)
        {
          mat_p<-mat_b+mat_now
          mat_u<-matrix(0,nrow(mat),ncol(mat))
          idx1<-!is.na(mat) & mat==1
          idx0<-!is.na(mat) & mat==0
          index<-is.na(mat)
          arr<-exp(mat_p)
          mat_u[index]<-mat_now[index]
          mat_u[idx1]<-mat_now[idx1]+eps*epsilon.binary/(1.0+arr[idx1])
          mat_u[idx0]<-mat_now[idx0]-eps*epsilon.binary*arr[idx0]/(1.0+arr[idx0])
          return (mat_u)
        }

        stop.binary<-function(mat,mat_b,mat_now,mat_u)
        {
          index<-!is.na(mat)
          mn<-mat_b+mat_now
          mu<-mat_b+mat_u
          arn<-exp(mn)
          aru<-exp(mu)
          idx1<-!is.na(mat) & mat==1
          idx0<-!is.na(mat) & mat==0
          lgn<-sum(log(arn[idx1]/(1+arn[idx1])))+sum(log(1/(1+arn[idx0])))
          lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
          return (lgu-lgn)
        }

        LL.binary<-function(mat,mat_b,mat_u)
        {
          index<-!is.na(mat)
          mu<-mat_b+mat_u
          aru<-exp(mu)
          idx1<-!is.na(mat) & mat==1
          idx0<-!is.na(mat) & mat==0
          lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
          return (lgu)
        }

        LLmax.binary<-function(mat)
        {
          return (0)
        }

        LLmin.binary<-function(mat,mat_b)
        {
          index<-!is.na(mat)
          aru<-exp(mat_b)
          idx1<-!is.na(mat) & mat==1
          idx0<-!is.na(mat) & mat==0
          lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
          return (lgu)
        }

        binary_type_base <- function( data,dimension=2 ,name="test")
        {
          data<-check.binary(data,name)
          data_b<-base.binary(data)
          data_now<-matrix(0,nrow(data),ncol(data))
          data_u<-update.binary(data,data_b,data_now)
          data_u<-nuclear_approximation(data_u,dimension)
          while (T)
          {
            thr<-stop.binary(data,data_b,data_now,data_u)
            message(thr)
            if (thr<0.2)
            {
              break
            }
            data_now<-data_u
            data_u<-update.binary(data,data_b,data_now)
            data_u<-nuclear_approximation(data_u,dimension)
          }
          return (data_now)
        }

        #----------#
        # gaussian #
        #----------#

        epsilon.gaussian=0.5

        check.gaussian.row<-function(arr)
        {
          if (sum(!is.na(arr))==0)
          {
            return (F)
          }
          else
          {
            return (T)
          }
        }
        check.gaussian<-function(mat,name)
        {
          index<-array(T,nrow(mat))
          for(i in 1:nrow(mat))
          {
            if (sum(is.na(mat[i,])==ncol(mat)))
            {
              war<-paste("Warning: ",name,"'s ",as.character(i)," line is all NA. Delete this line",sep="")
              warning(war)
              index[i]<-F
            }
          }
          mat_c<-mat[index,]
          rownames(mat_c)<-rownames(mat)[index]
          colnames(mat_c)<-colnames(mat)
          return (mat_c)
        }

        base.gaussian.row<-function(arr)
        {
          idx<-!is.na(arr)
          return (mean(arr[idx]))
        }

        base.gaussian<-function(mat)
        {
          mat_b<-matrix(0,nrow(mat),ncol(mat))
          ar_b<-apply(mat,1,base.gaussian.row)
          mat_b[1:nrow(mat_b),]<-ar_b
          return (mat_b)
        }

        update.gaussian<-function(mat,mat_b,mat_now,eps)
        {
          mat_p<-mat_b+mat_now
          mat_u<-matrix(0,nrow(mat),ncol(mat))
          index<-!is.na(mat)
          mat_u[index]<-mat_now[index]+eps*epsilon.gaussian*(mat[index]-mat_p[index])
          index<-is.na(mat)
          mat_u[index]<-mat_now[index]
          return (mat_u)
        }

        stop.gaussian<-function(mat,mat_b,mat_now,mat_u)
        {
          index<-!is.na(mat)
          mn<-mat_b+mat_now
          mu<-mat_b+mat_u
          ren<-mat[index]-mn[index]
          reu<-mat[index]-mu[index]
          lgn<- -0.5*sum(ren*ren)
          lgu<- -0.5*sum(reu*reu)
          return (lgu-lgn)
        }

        LL.gaussian<-function(mat,mat_b,mat_u)
        {
          index<-!is.na(mat)
          mu<-mat_b+mat_u
          reu<-mat[index]-mu[index]
          lgu<- -0.5*sum(reu*reu)
          return (lgu)
        }

        LLmax.gaussian<-function(mat)
        {
          return (0.0)
        }

        LLmin.gaussian<-function(mat,mat_b)
        {
          index<-!is.na(mat)
          reu<-mat[index]-mat_b[index]
          lgu<- -0.5*sum(reu*reu)
          return (lgu)
        }

        gaussian_base<-function(data,dimension=2,name="test")
        {
          data<-check.gaussian(data,name)
          data_b<-base.gaussian(data)
          data_now<-matrix(0,nrow(data),ncol(data))
          data_u<-update.gaussian(data,data_b,data_now)
          data_u<-nuclear_approximation(data_u,dimension)
          while(T)
          {
            thr<-stop.gaussian(data,data_b,data_now,data_u)
            message(thr)
            if (thr<0.2)
            {
              break
            }
            data_now<-data_u
            data_u<-update.gaussian(data,data_b,data_now)
            data_u<-nuclear_approximation(data_u,dimension)
          }
          return (data_now)
        }

        #---------#
        # poisson #
        #---------#

        epsilon.poisson<-0.5

        check.poisson.row<-function(arr)
        {
          if (sum(!is.na(arr))==0)
          {
            return (F)
          }
          else
          {
            idx<-!is.na(arr)
            if (sum(arr[idx]<0)>0)
            {
              return (F)
            }
            else
            {
              return (T)
            }
          }
        }

        check.poisson<-function(mat,name)
        {
          w<-paste(name," is poisson type. Add 1 to all counts",sep="")
          message(w)
          index<-apply(mat,1,check.poisson.row)
          n<-sum(!index)
          if (n>0)
          {
            w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
            warning(w)
          }
          mat_c<-mat[index,]+1
          rownames(mat_c)<-rownames(mat)[index]
          colnames(mat_c)<-colnames(mat)
          return (mat_c)
        }

        base.poisson.row<-function(arr)
        {
          idx<-!is.na(arr)
          m<-sum(log(arr[idx]))
          n<-sum(idx)
          return(m/n)
        }

        base.poisson<-function(mat)
        {
          mat_b<-matrix(0,nrow(mat),ncol(mat))
          ar_b<-apply(mat,1,base.poisson.row)
          mat_b[1:nrow(mat_b),]<-ar_b
          return (mat_b)
        }

        update.poisson<-function(mat,mat_b,mat_now,eps)
        {
          mat_p<-mat_b+mat_now
          mat_u<-matrix(0,nrow(mat),ncol(mat))
          index<-!is.na(mat)
          mat_u[index]<-mat_now[index]+eps*epsilon.poisson*(log(mat[index])-mat_p[index])
          index<-is.na(mat)
          mat_u[index]<-mat_now[index]
          return (mat_u)
        }

        stop.poisson<-function(mat,mat_b,mat_now,mat_u)
        {
          index<-!is.na(mat)
          mn<-mat_b+mat_now
          mu<-mat_b+mat_u
          lgn<-sum(mat[index]*mn[index]-exp(mn[index]))
          lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
          return (lgu-lgn)
        }

        LL.poisson<-function(mat,mat_b,mat_u)
        {
          index<-!is.na(mat)
          mu<-mat_b+mat_u
          lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
          return (lgu)
        }

        LLmax.poisson<-function(mat)
        {
          index<-!is.na(mat)
          lgu<-sum(mat[index]*log(mat[index])-mat[index])
          return (lgu)
        }

        LLmin.poisson<-function(mat,mat_b)
        {
          index<-!is.na(mat)
          lgu<-sum(mat[index]*mat_b[index]-exp(mat_b[index]))
          return (lgu)
        }

        poisson_type_base<-function(data,dimension=2,name="test")
        {
          data<-check.poisson(data,name)
          data_b<-base.poisson(data)
          data_now<-matrix(0,nrow(data),ncol(data))
          data_u<-update.poisson(data,data_b,data_now)
          data_u<-nuclear_approximation(data_u,dimension)
          while(T)
          {
            thr<-stop.poisson(data,data_b,data_now,data_u)
            message(thr)
            if (thr<0.2)
            {
              break
            }
            data_now<-data_u
            data_u<-update.poisson(data,data_b,data_now)
            data_u<-nuclear_approximation(data_u,dimension)
          }
          return (data_now)
        }

        #----#
        # na #
        #----#

        nuclear_approximation<-function(mat,dimension)
        {
          svd<-svd(mat,nu=0,nv=0)
          if (dimension<length(svd$d))
          {
            lambda<-svd$d[dimension+1]
            svd<-svd(mat,nu=dimension,nv=dimension)
            indexh<-svd$d>lambda
            indexm<-svd$d<lambda
            dia<-array(svd$d,length(svd$d))
            dia[indexh]<-dia[indexh]-lambda
            dia[indexm]<-0
            mat_low<-svd$u%*%diag(c(dia[1:dimension],0))[1:dimension,1:dimension]%*%t(svd$v)
          }
          else
          {
            mat_low<-mat
          }
          return (mat_low)
        }

        #------------#
        # LRAcluster #
        #------------#
        check.matrix.element<-function(x)
        {
          if (!is.matrix(x))
          {
            return (T)
          }
          else
          {
            return (F)
          }
        }

        ncol.element<-function(x)
        {
          return (ncol(x))
        }

        nrow.element<-function(x)
        {
          return (nrow(x))
        }

        check<-function(mat,type,name)
        {
          if (type=="binary")
          {
            return (check.binary(mat,name))
          }
          else if (type=="gaussian")
          {
            return (check.gaussian(mat,name))
          }
          else if (type=="poisson")
          {
            return (check.poisson(mat,name))
          }
          else
          {
            e<-paste("unknown type ",type,sep="")
            stop(e)
          }
        }

        eps<-0.0
        if (!is.list(data))
        {
          stop("the input data must be a list!")
        }
        c<-sapply(data,check.matrix.element)
        if (sum(c)>0)
        {
          stop("each element of input list must be a matrix!")
        }
        c<-sapply(data,ncol.element)
        if (length(levels(factor(c)))>1)
        {
          stop("each element of input list must have the same column number!")
        }
        if (length(data)!=length(types))
        {
          stop("data and types must be the same length!")
        }
        nSample<-c[1]
        loglmin<-0
        loglmax<-0
        loglu<-0.0
        nData<-length(data)
        for (i in 1:nData)
        {
          data[[i]]<-check(data[[i]],types[[i]],names[[i]])
        }
        nGeneArr<-sapply(data,nrow.element)
        nGene<-sum(nGeneArr)
        indexData<-list()
        k=1
        for(i in 1:nData)
        {
          indexData[[i]]<- (k):(k+nGeneArr[i]-1)
          k<-k+nGeneArr[i]
        }
        base<-matrix(0,nGene,nSample)
        now<-matrix(0,nGene,nSample)
        update<-matrix(0,nGene,nSample)
        thr<-array(0,nData)
        for (i in 1:nData)
        {
          if (types[[i]]=="binary")
          {
            base[indexData[[i]],]<-base.binary(data[[i]])
            loglmin<-loglmin+LLmin.binary(data[[i]],base[indexData[[i]],])
            loglmax<-loglmax+LLmax.binary(data[[i]])
          }
          else if (types[[i]]=="gaussian")
          {
            base[indexData[[i]],]<-base.gaussian(data[[i]])
            loglmin<-loglmin+LLmin.gaussian(data[[i]],base[indexData[[i]],])
            loglmax<-loglmax+LLmax.gaussian(data[[i]])
          }
          else if (types[[i]]=="poisson")
          {
            base[indexData[[i]],]<-base.poisson(data[[i]])
            loglmin<-loglmin+LLmin.poisson(data[[i]],base[indexData[[i]],])
            loglmax<-loglmax+LLmax.poisson(data[[i]])
          }
        }
        for (i in 1:nData)
        {
          if (types[[i]]=="binary")
          {
            update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
          }
          else if (types[[i]]=="gaussian")
          {
            update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
          }
          else if (types[[i]]=="poisson")
          {
            update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
          }
        }
        update<-nuclear_approximation(update,dimension)
        nIter<-0
        thres<-array(Inf,3)
        epsN<-array(Inf,2)
        while(T)
        {
          for (i in 1:nData)
          {
            if (types[[i]]=="binary")
            {
              thr[i]<-stop.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
            }
            else if (types[[i]]=="gaussian")
            {
              thr[i]<-stop.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
            }
            else if (types[[i]]=="poisson")
            {
              thr[i]<-stop.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
            }
          }
          nIter<-nIter+1
          thres[1]<-thres[2]
          thres[2]<-thres[3]
          thres[3]<-sum(thr)
          epsN[1]<-epsN[2]
          epsN[2]<-eps
          if (nIter>5)
          {
            if (runif(1)<thres[1]*thres[3]/(thres[2]*thres[2]+thres[1]*thres[3]))
            {
              eps<-epsN[1]+0.05*runif(1)-0.025
            }
            else
            {
              eps<-epsN[2]+0.05*runif(1)-0.025
            }
            if (eps< -0.7)
            {
              eps<- 0
              epsN<-c(0,0)
            }
            if (eps > 1.4)
            {
              eps<-0
              epsN<-c(0,0)
            }
          }
          if (sum(thr)<nData*0.2)
          {
            break
          }
          now<-update
          for (i in 1:nData)
          {
            if (types[[i]]=="binary")
            {
              update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
            }
            else if (types[[i]]=="gaussian")
            {
              update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
            }
            else if (types[[i]]=="poisson")
            {
              update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
            }
          }
          update<-nuclear_approximation(update,dimension)
        }
        for (i in 1:nData)
        {
          if (types[[i]]=="binary")
          {
            loglu<-loglu+LL.binary(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
          }
          else if (types[[i]]=="gaussian")
          {
            loglu<-loglu+LL.gaussian(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
          }
          else if (types[[i]]=="poisson")
          {
            loglu<-loglu+LL.poisson(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
          }
        }
        sv<-svd(update,nu=0,nv=dimension)
        coordinate<-diag(c(sv$d[1:dimension],0))[1:dimension,1:dimension]%*%t(sv$v)
        colnames(coordinate)<-colnames(data[[1]])
        rownames(coordinate)<-paste("PC ",as.character(1:dimension),sep="")
        ratio<-(loglu-loglmin)/(loglmax-loglmin)
        return (list("coordinate"=coordinate,"potential"=ratio))
      }
      data <- lapply(data, as.matrix)
      fit <- LRAcluster(data, dimension = clust.num, types = as.list(type))
      dist <- fit$coordinate %>% t %>% dist
      clust.dend <- hclust(dist, method = "ward.D")

      clustres <- data.frame(samID = colnames(data[[1]]),
                             clust = cutree(clust.dend,k = clust.num),
                             row.names = colnames(data[[1]]),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]
      clust.res.list[[method]] <- list(fit = fit, clust.res = clustres, clust.dend = clust.dend, method = "LRAcluster")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))
    }

    if (method=="iClusterBayes"){
      library(iClusterPlus)
      message(paste0(method," start..."))
      data.backup <- data
      dt <- lapply(data, t)
      res <- iClusterBayes(dt1         = dt[[1]],
                           dt2         = dt[[2]],
                           type        = type,
                           K           = clust.num - 1,
                           n.burnin    = 18000,
                           n.draw      = 12000,
                           prior.gamma = rep(0.5,length(data)),
                           sdev        = 0.05,
                           thin        = 3)
      message("clustering done...")
      clustres <- data.frame(samID = rownames(dt[[1]]),
                             clust = res$clusters,
                             row.names = rownames(dt[[1]]),
                             stringsAsFactors = FALSE)
      clustres <- clustres[order(clustres$clust),]


      featres <- data.frame(feature = as.character(unlist(lapply(data.backup, rownames))),
                            dataset = rep(names(data),as.numeric(sapply(data.backup, function(x) length(rownames(x))))),
                            post.prob = unlist(res$beta.pp),
                            stringsAsFactors = FALSE)
      feat.res <- NULL
      for (d in unique(featres$dataset)) {
        tmp <- featres[which(featres$dataset == d),]
        tmp <- tmp[order(tmp$post.prob, decreasing = TRUE),]
        tmp$rank <- 1:nrow(tmp)
        feat.res <- rbind.data.frame(feat.res,tmp)

        clust.res.list[[method]] <- list(fit = res, clust.res = clustres, feat.res = feat.res, method = "iClusterBayes")

      }
      message("feature selection done...")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))

    }

    if (method=="SNF") {
    message(paste0(method," start..."))
    dt <- lapply(data, t)

    dat <- lapply(dt, function (dd){
      dd <- dd %>% as.matrix
      W <- dd %>% SNFtool::dist2(dd) %>% SNFtool::affinityMatrix(K = 30, sigma = 0.5)
    })
    W <-  SNFtool::SNF(Wall = dat,
                       K    = 30,
                       t    = 20)
    clust.SNF = W %>% SNFtool::spectralClustering(clust.num)

    clustres <- data.frame(samID = rownames(dt[[1]]),
                           clust = clust.SNF,
                           row.names = rownames(dt[[1]]),
                           stringsAsFactors = FALSE)
    clustres <- clustres[order(clustres$clust),]
    clust.res.list[[method]] <- list(fit = W, clust.res = clustres, method = "SNF")
    print(clust.res.list[[method]])
    message(paste0(method," done..."))
    }

    if (method=="ConsensusClustering"){
      library(ConsensusClusterPlus)
      message(paste0(method," start..."))
      d <- do.call(rbind, data)
      fit <-  ConsensusClusterPlus(d            = as.matrix(d),
                                   maxK         = ifelse(clust.num == 2, 3, clust.num), # cannot set as 2
                                   reps         = 500,
                                   pItem        = 0.8,
                                   pFeature     = 0.8,
                                   clusterAlg   = "hc",
                                   innerLinkage = "ward.D",
                                   finalLinkage = "ward.D",
                                   distance     = "pearson",
                                   seed         = 20000709,
                                   verbose      = F,
                                   plot         = NULL,
                                   writeTable   = F,
                                   title        = 'consensuscluster')
      res <- fit[[clust.num]]

      clustres <- data.frame(samID = colnames(data[[1]]),
                             clust = as.numeric(res$consensusClass),
                             row.names = colnames(data[[1]]),
                             stringsAsFactors = F)
      clustres <- clustres[order(clustres$clust),]
      clust.res.list[[method]] <- list(fit = fit, clust.res = clustres, clust.dend = res$consensusTree, method = "ConsensusClustering")
      print(clust.res.list[[method]])
      message(paste0(method," done..."))
    }

    assign("clust.res.list", clust.res.list, envir = .GlobalEnv)
    message("file saved...")
    saveRDS(clust.res.list,file = "clust.res.list.rds")
  }

}
