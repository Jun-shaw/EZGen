
create_model <- function(prog.data.list=prog.data.list,
                         train.data.pos=1,
                         gene.list=gene.list,
                         unicox_pcutoff=0.05,
                         top.num = 100,
                         method='all',
                         hm.col=c("#5770A6", "#FFFFFF", "#CE5C69"),
                         cohort.col=c('#CE5C69','#5770A6','#A281B1'))
{
  if(method=='out_RSF'){
    res <- create_model_outRSF(prog.data.list=prog.data.list,
                         train.data.pos=train.data.pos,
                         gene.list=gene.list,
                         unicox_pcutoff=unicox_pcutoff,
                         top.num = top.num,
                         hm.col=hm.col,
                         cohort.col=cohort.col)
    return(res)
  }
  if(method=='all'){
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(data.table)
  library(ggplot2)
  library(ggsci)
  library(survival)
  library(randomForestSRC)
  library(glmnet)
  library(plsRcox)
  library(superpc)
  library(gbm)
  library(CoxBoost)
  library(survivalsvm)
  library(dplyr)
  library(tibble)
  library(BART)
  library(ggbreak)
  library(tidyr)
  library(ggbreak)
  library(edgeR)
  library(limma)
  library(survival)
  library(survminer)
  library(stringi)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(beepr)
  library(pheatmap)
  library(data.table)
  library(ggsignif)
  library(RColorBrewer)
  library(future.apply)
  library(gplots)
  library(DESeq2)
  library(ggrepel)
  library(Rcpp)
  library(survivalsvm)
  library(dplyr)
  library(rms)
  library(pec)
  library(ggDCA)
  library(glmnet)
  library(foreign)
  library(regplot)
  library(randomForestSRC)
  library(timeROC)
  library(tidyr)
  library(tibble)
  library(caret)
  library(regplot)
  library(gbm)
  library(tidyverse)
  library(gbm)
  library(obliqueRSF)
  library(remotes)
  library(aorsf)
  library(party)
  library(partykit)
  library(xgboost)

  if (!dir.exists('ML_plot')) {
    dir.create('ML_plot')
  }
  harmony <- function(dat=rs)
  {rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]}

  options(stringsAsFactors = FALSE)

  for (i in  1:length(prog.data.list)) {
    colnames(prog.data.list[[i]])[1:3]=c('ID','time','status')
    prog.data.list[[i]]$time=as.numeric(prog.data.list[[i]]$time)
    prog.data.list[[i]]$status=as.numeric(prog.data.list[[i]]$status)
  }


  result <- data.frame()
  rf_nodesize <- 5
  seed <- 20000709

  cox.res.list <- fliter_unicox (prog.data.list=prog.data.list,
                                      train.data.pos=train.data.pos,
                                      gene.list=gene.list,
                                      Pvalue=unicox_pcutoff)
  train <- cox.res.list$train.vali.list[[train.data.pos]]
  ml.res <- list()
  rs.res <- list()
#### 1-1 RSF ###################
  message(paste0("\n=> 1-1 RSF", " machine learning algorithm."))
  set.seed(seed)
  fit <- rfsrc(Surv(time,status)~.,data = train,
               ntree = 1000,nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  best <- which.min(fit$err.rate)
  set.seed(seed)
  fit <- rfsrc(Surv(time,status)~.,data = train,
               ntree = best,nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(cox.res.list$train.vali.list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- 'RSF'
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs

#### 1-2 RSF+Enet ###############
  message(paste0("\n=> 1-2 RSF+Enet", " machine learning algorithm."))

  #基因重要性
  vi <- data.frame(imp=vimp.rfsrc(fit)$importance)
  vi$imp <- (vi$imp-min(vi$imp))/(max(vi$imp)-min(vi$imp))
  vi$ID <- rownames(vi)

  #基因重要性可视化
  pdf("ML_plot/rsf_highgene(rsf).pdf")
  ggplot(vi,aes(imp,reorder(ID,imp)))+
    geom_bar(stat = 'identity',fill='#CE5C69',color='black',width=0.7)+
    geom_vline(xintercept = 0.01,color='grey50',linetype=2)+
    labs(x='Relative importance by Random Forest',y=NULL)+
    theme_bw(base_rect_size = 1.5)+
    theme(axis.text.x = element_text(size = 11,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.title = element_text(size=13,color='black'),
          legend.text = element_text(size=12,color='black'),
          legend.title = element_text(size=13,color='black'))+
    scale_y_discrete(expand = c(0.03,0.03))+
    scale_x_continuous(expand = c(0.01,0.01))
  dev.off()
  #提取重要性大于0.01的基因
  rid <- rownames(vi)[vi$imp>0.01]
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply(cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  x1 <- as.matrix(train2[,rid])
  x2 <- as.matrix(Surv(train2$time,train2$status))

  #利用循环探索Enet最佳模型
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
    rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
    harmony(rs)
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('RSF + Enet','[α=',alpha,']')
    result <- rbind(result,cc)
    ml.res [[unique(cc$Model)]]<- fit
    rs.res [[unique(cc$Model)]]<- rs
  }
  #利用交叉验证探索Enet最佳模型
  set.seed(seed)
  modelexp=as.matrix(trainlist2[[train.data.pos]][,c(3:ncol(trainlist2[[train.data.pos]]))])
  modelstat=Surv(trainlist2[[train.data.pos]]$time,trainlist2[[train.data.pos]]$status)
  Enetmodel <- glmnet(modelexp,modelstat,family = 'cox',nfolds=10)
  Enetmodel_cv<-cv.glmnet(modelexp,modelstat,family = 'cox',nfolds=10)
  #建立最优模型
  model_fit<-glmnet(modelexp,modelstat,family = 'cox',nfolds=10,keep=T,lambda = Enetmodel_cv$lambda.min)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model_fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=model_fit$lambda.min)))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + Enet','[lambda=',round(Enetmodel_cv$lambda.min,3),']')
  result <- rbind(result,cc)
  ml.res [[paste0('RSF + Enet','[lambda=',round(Enetmodel_cv$lambda.min,3),']')]]<- model_fit
  rs.res [[paste0('RSF + Enet','[lambda=',round(Enetmodel_cv$lambda.min,3),']')]]<- rs
#### 1-3 RSF+stepcox #######
message(paste0("\n=> 1-3 RSF+stepcox", " machine learning algorithm."))
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time,status)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + StepCox','[',direction,']')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}
#### 1-4 RSF+CoxBoost ##########
message(paste0("\n=> 1-4 RSF+CoxBoost", " machine learning algorithm."))
set.seed(seed)
#计算最佳penalty
pen <- optimCoxBoostPenalty(train2[,'time'],train2[,'status'],as.matrix(train2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
#计算最佳stepno
cv.res <- cv.CoxBoost(train2[,'time'],train2[,'status'],as.matrix(train2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
#构建最佳模型
fit <- CoxBoost(train2[,'time'],train2[,'status'],as.matrix(train2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 1-5 RSF+plsRcox ###########
message(paste0("\n=> 1-5 RSF+plsRcox", " machine learning algorithm."))
set.seed(seed)
#建立专用矩阵
model_exp=data.frame(train2[,-c(1:2)])
model_time=train2$time
model_stat=train2$status
#建立模型
model<-plsRcox(model_exp,time = model_time,event = model_stat,nt=10)
#进行交叉验证
cv.model<-cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = T)
#构建最优模型
model<-plsRcox(model_exp,
               time = model_time,
               event = model_stat,
               nt=cv.model$lambda.min5,
               alpha.pvals.expli = 0.05,
               sparse = T,
               pvals.expli = T)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model,type="lp",newdata=x[,-c(1,2)])))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + plsRcox')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs
#### 1-6 RSF+superpc ############
message(paste0("\n=> 1-6 RSF+superpc", " machine learning algorithm."))
data <- list(x=t(train2[,-c(1,2)]),y=train2$time,censoring.status=train2$status,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5)
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$time,censoring.status=w$status,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + SuperPC')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs
#### 1-7 RSF+gbm ################
message(paste0("\n=> 1-7 RSF+gbm", " machine learning algorithm."))
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,
           data = train2,
           distribution = 'coxph',
           n.minobsinnode = 10,
           n.cores = 1,
           n.trees = 1000,
           shrinkage = 0.005,
           interaction.depth = 2,
           cv.folds = 5)
#构建最优模型
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,
           data = train2,
           distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,
           n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + GBM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs
#### 1-8 RSF+survivalsvm ########
message(paste0("\n=> 1-8 RSF+survivalsvm", " machine learning algorithm."))
set.seed(seed)
fit = survivalsvm(Surv(time,status)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + survival-SVM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs
#### 1-9 RSF+Ridge ##############
#message(paste0("\n=> 1-9 RSF+Ridge", " machine learning algorithm."))
#set.seed(seed)
#modelexp=as.matrix(train2[,c(3:ncol(train2))])
#利用循环探索Ridge最佳模型
#for (alpha in seq(0,1,0.1)) {
#  alpha <- 0.1
#  sum(is.na(modelexp))
#  sum(is.na(train2$status))
#  str(modelexp)
#  str(train2$status)
#  train2$status <- as.factor(train2$status)
 # modelexp <-

 # model <- glmnet(modelexp,train2$status,family = 'binomial',alpha = alpha,nfolds=10)
 # model_cv<-cv.glmnet(modelexp,train2$status,family = 'binomial',alpha =alpha,nfolds=10)
 # fit<-glmnet(modelexp,train2$status,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
 # rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
 # harmony(rs)
#  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
#    rownames_to_column('ID')
#  cc$Model <- paste0('RSF + Ridge','[α=',alpha,']')
#  result <- rbind(result,cc)
 # ml.res [[unique(cc$Model)]]<- fit
 # rs.res [[unique(cc$Model)]]<- rs
#}

#### 1-10 RSF+obliqueRSF ########
message(paste0("\n=> 1-10 RSF+obliqueRSF", " machine learning algorithm."))
set.seed(seed)
model<-orsf(data = train2,n_tree = 100,formula = Surv(time,status)~.)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + obliqueRSF')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs
#### 1-11 RSF+xgboost ###########
message(paste0("\n=> 1-11 RSF+xgboost", " machine learning algorithm."))
set.seed(seed)
#建立专用矩阵
model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$time)
#构建参数
object<-list(bojective="surivival:cox",
             booster="gbtree",
             eval_metric="cox-nloglik",
             eta=0.01,
             max_depth=3,
             subsample=1,
             colsample_bytree=1,
             gamma=0.5)
#构建模型
model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + xgboost')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 1-12 RSF+CForest###########
message(paste0("\n=> 1-12 RSF+CForest", " machine learning algorithm."))
set.seed(seed)
model<-party::cforest(Surv(time,status)~.,data=train2,controls = cforest_unbiased(ntree=50))
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CForest')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 1-13 RSF+CTree##############
message(paste0("\n=> 1-13 RSF+CTree", " machine learning algorithm."))
set.seed(seed)
#建立模型
model<-ctree(Surv(time,status)~.,data=train2)
#计算变量重要性
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CTree')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-1 Enet ##########
message(paste0("\n=> 2-1 Enet", " machine learning algorithm."))
modelexp=as.matrix(train[,c(3:ncol(train))])
modelstat=Surv(train$time,train$status)
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,modelstat,family = 'cox',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,modelstat,family = 'cox',alpha = alpha,nfolds=10)
  fit <-glmnet(modelexp,modelstat,family = 'cox',alpha =alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(cox.res.list$train.vali.list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}


#### 2-2 Lasso+RSF ###########
message(paste0("\n=> 2-2 Lasso+RSF", " machine learning algorithm."))
set.seed(seed)
fit = cv.glmnet(modelexp, modelstat,family = "cox")
coef.min = coef(fit, s = "lambda.min")
rid <- coef.min@Dimnames[[1]]

train2 <- train[,c('time','status',rid)]
trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

set.seed(seed)
fit <- rfsrc(Surv(time,status)~.,data = train2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(time,status)~.,data = train2,
             ntree = best,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + RSF')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

##### 2-3 Lasso+StepCox ########
message(paste0("\n=> 2-3 Lasso+StepCox", " machine learning algorithm."))

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time,status)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + StepCox','[',direction,']')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}


#### 2-4 Lasso+CoxBoost #########
message(paste0("\n=> 2-4 Lasso+CoxBoost", " machine learning algorithm."))

set.seed(seed)
#计算最佳penalty
modelpen<-optimCoxBoostPenalty(time = train2$time,
                               status = train2$status,
                               as.matrix(train2[,-c(1:2)]),
                               trace = T,
                               parallel = T)
#计算最佳stepno
cvmodel<-cv.CoxBoost(time = train2$time,
                     status = train2$status,
                     as.matrix(train2[,-c(1:2)]),
                     maxstepno = 100,
                     K = 3,
                     type = "verweij",
                     penalty=modelpen$penalty)
#构建CoxBoost模型
fit<-CoxBoost(time = train2$time,
              status = train2$status,
              as.matrix(train2[,-c(1:2)]),
              stepno = cvmodel$optimal.step,
              penalty = modelpen$penalty)

rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-5 Lasso+plsRcox #########
message(paste0("\n=> 2-5 Lasso+plsRcox", " machine learning algorithm."))

set.seed(seed)
#建立专用矩阵
model_exp=data.frame(train2[,-c(1:2)])
model_time=train2$time
model_stat=train2$status
#建立模型
model<-plsRcox(model_exp,time = model_time,event = model_stat,nt=10)
#进行交叉验证
cv.model<-cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = F)
#构建最优模型
cv.plsRcox.res=cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = F)
fit <- plsRcox(model_exp,
               time = model_time,
               event = model_stat,
               nt=cv.model$lambda.min5,
               alpha.pvals.expli = 0.05,
               sparse = T,
               pvals.expli = T)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + plsRcox')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-6 Lasso+superpc #########
message(paste0("\n=> 2-6 Lasso+superpc", " machine learning algorithm."))

data <- list(x=t(train2[,-c(1,2)]),y=train2$time,censoring.status=train2$status,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5)
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$time,censoring.status=w$status,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + SuperPC')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-7 Lasso+gbm ##############
message(paste0("\n=> 2-7 Lasso+gbm", " machine learning algorithm."))

set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,data = train2,distribution = 'coxph',
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,data = train2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + GBM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-8 Lasso+survivalsvm #########
message(paste0("\n=> 2-8 Lasso+survivalsvm", " machine learning algorithm."))

set.seed(seed)
fit = survivalsvm(Surv(time,status)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + survival-SVM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-9 Lasso+Ridge #########
message(paste0("\n=> 2-9 Lasso+Ridge", " machine learning algorithm."))

set.seed(seed)
modelexp=as.matrix(train2[,c(3:ncol(train2))])
#利用循环探索Ridge最佳模型
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,train2$status,family = 'binomial',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,train2$status,family = 'binomial',alpha =alpha,nfolds=10)
  fit<-glmnet(modelexp,train2$status,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + Ridge','[α=',alpha,']')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 2-10 Lasso+obliqueRSF ########
message(paste0("\n=> 2-10 Lasso+obliqueRSF", " machine learning algorithm."))

set.seed(seed)
model<-orsf(data = train2,n_tree = 100,formula = Surv(time,status)~.)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + obliqueRSF')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-11 Lasso+xgboost #########
message(paste0("\n=> 2-11 Lasso+xgboost", " machine learning algorithm."))

set.seed(seed)
#建立专用矩阵
model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$time)
#构建参数
object<-list(bojective="surivival:cox",
             booster="gbtree",
             eval_metric="cox-nloglik",
             eta=0.01,
             max_depth=3,
             subsample=1,
             colsample_bytree=1,
             gamma=0.5)
#构建模型
model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + xgboost')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-12 Lasso+CForest#########
message(paste0("\n=> 2-12 Lasso+CForest", " machine learning algorithm."))

set.seed(seed)
model<-party::cforest(Surv(time,status)~.,data=train2,controls = cforest_unbiased(ntree=50))
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CForest')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 2-13 Lasso+CTree############
message(paste0("\n=> 2-13 Lasso+CTree", " machine learning algorithm."))

set.seed(seed)
model<-ctree(Surv(time,status)~.,data=train2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CTree')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 3-1 StepCox ################
message(paste0("\n=> 3-1 StepCox", " machine learning algorithm."))

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rs <- lapply(cox.res.list$train.vali.list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-2 StepCox+RSF ############
message(paste0("\n=> 3-2 StepCox+RSF", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))

  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply(cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  set.seed(seed)
  fit <- rfsrc(Surv(time,status)~.,data = train2,
               ntree = 1000,nodesize = rf_nodesize,##该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  best <- which.min(fit$err.rate)
  set.seed(seed)
  fit <- rfsrc(Surv(time,status)~.,data = train2,
               ntree = best,nodesize = rf_nodesize,##该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + RSF')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-3 StepCox+Enet ##########
message(paste0("\n=> 3-3 StepCox+Enet", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply(cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})
  x1 <- as.matrix(train2[,rid])
  x2 <- as.matrix(Surv(train2$time,train2$status))

  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
    rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
    harmony(rs)
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox','[',direction,']',' + Enet','[α=',alpha,']')
    result <- rbind(result,cc)
    ml.res [[unique(cc$Model)]]<- fit
    rs.res [[unique(cc$Model)]]<- rs
  }
}

#### 3-4 StepCox+CoxBoost ########
message(paste0("\n=> 3-4 StepCox+CoxBoost", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  set.seed(seed)
  pen <- optimCoxBoostPenalty(train2[,'time'],train2[,'status'],as.matrix(train2[,-c(1,2)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(train2[,'time'],train2[,'status'],as.matrix(train2[,-c(1,2)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(train2[,'time'],train2[,'status'],as.matrix(train2[,-c(1,2)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + CoxBoost')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-5 StepCox+plsRcox ########
message(paste0("\n=> 3-5 StepCox+plsRcox", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  set.seed(seed)
  cv.plsRcox.res=cv.plsRcox(list(x=train2[,rid],time=train2$time,status=train2$status),nt=10,nfold = 10,verbose = F)
  fit <- plsRcox(train2[,rid],time=train2$time,event=train2$status,nt=as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + plsRcox')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-6 StepCox+superpc ########
message(paste0("\n=> 3-6 StepCox+superpc", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  data <- list(x=t(train2[,-c(1,2)]),y=train2$time,censoring.status=train2$status,featurenames=colnames(train2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5)
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                       n.fold = 5,
                       n.components=3,
                       min.features=1,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  rs <- lapply(trainlist2,function(w){
    test <- list(x=t(w[,-c(1,2)]),y=w$time,censoring.status=w$status,featurenames=colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2],RS=rr)
    return(rr2)
  })
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + SuperPC')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-7 StepCox+gbm ########
message(paste0("\n=> 3-7 StepCox+gbm", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  set.seed(seed)
  fit <- gbm(formula = Surv(time,status)~.,data = train2,distribution = 'coxph',
             n.trees = 1000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 5,n.cores = 1)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(time,status)~.,data = train2,distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 5,n.cores = 1)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + GBM')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-8 StepCox+survival-SVM #########
message(paste0("\n=> 3-8 StepCox+survival-SVM", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  #direction='both'
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  fit = survivalsvm(Surv(time,status)~., data= train2, gamma.mu = 1)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + survival-SVM')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-9 StepCox+Ridge ########
message(paste0("\n=> 3-9 StepCox+Ridge", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})
  set.seed(seed)
  modelexp=as.matrix(train2[,c(3:ncol(train2))])
  #利用循环探索Ridge最佳模型
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    model <- glmnet(modelexp,train2$status,family = 'binomial',alpha = alpha,nfolds=10)
    model_cv<-cv.glmnet(modelexp,train2$status,family = 'binomial',alpha =alpha,nfolds=10)
    fit<-glmnet(modelexp,train2$status,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
    rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
    harmony(rs)
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox','[',direction,']','+ Ridge','[α=',alpha,']')
    result <- rbind(result,cc)
    ml.res [[unique(cc$Model)]]<- fit
    rs.res [[unique(cc$Model)]]<- rs
  }
}

#### 3-10 StepCox+obliqueRSF ##########
message(paste0("\n=> 3-10 StepCox+obliqueRSF", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})
  model<-orsf(data = train2,n_tree = 100,formula = Surv(time,status)~.)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']', '+ obliqueRSF')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-11 StepCox+xgboost #########
message(paste0("\n=> 3-11 StepCox+xgboost", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  #direction='both'
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

  #建立专用矩阵
  model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$time)
  #构建参数
  object<-list(bojective="surivival:cox",
               booster="gbtree",
               eval_metric="cox-nloglik",
               eta=0.01,
               max_depth=3,
               subsample=1,
               colsample_bytree=1,
               gamma=0.5)
  #构建模型
  model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']', '+ xgboost')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-12 StepCox+CForest #########
message(paste0("\n=> 3-12 StepCox+CForest", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})
  model<-party::cforest(Surv(time,status)~.,data=train2,controls = cforest_unbiased(ntree=50))
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + CForest')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 3-13 StepCox+CTree #########
message(paste0("\n=> 3-13 StepCox+CTree", " machine learning algorithm."))

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(time,status)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('time','status',rid)]
  trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})
  model<-ctree(Surv(time,status)~.,data=train2)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + CTree')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 4-1 CoxBoost #########
message(paste0("\n=> 4-1 CoxBoost", " machine learning algorithm."))

set.seed(seed)
pen <- optimCoxBoostPenalty(train[,'time'],train[,'status'],as.matrix(train[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train[,'time'],train[,'status'],as.matrix(train[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(train[,'time'],train[,'status'],as.matrix(train[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply( cox.res.list$train.vali.list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 4-2 CoxBoost+Enet ########
message(paste0("\n=> 4-2 CoxBoost+Enet", " machine learning algorithm."))

rid <- names(coef(fit)[which(coef(fit)!=0)])
train2 <- train[,c('time','status',rid)]
trainlist2 <- lapply( cox.res.list$train.vali.list,function(x){x[,c('time','status',rid)]})

x1 <- as.matrix(train2[,rid])
x2 <- as.matrix(Surv(train2$time,train2$status))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 4-3 CoxBoost+stepcox #######
message(paste0("\n=> 4-3 CoxBoost+stepcox", " machine learning algorithm."))

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time,status)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  harmony(rs)
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + StepCox','[',direction,']')
  result <- rbind(result,cc)
  ml.res [[unique(cc$Model)]]<- fit
  rs.res [[unique(cc$Model)]]<- rs
}

#### 4-4 CoxBoost+RSF ########
message(paste0("\n=> 4-4 CoxBoost+RSF", " machine learning algorithm."))

set.seed(seed)
fit <- rfsrc(Surv(time,status)~.,data = train2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(time,status)~.,data = train2,
             ntree = best,nodesize = rf_nodesize,##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'CoxBoost + RSF'
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 4-5 RSF+plsRcox ########
message(paste0("\n=> 4-5 RSF+plsRcox", " machine learning algorithm."))

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train2[,rid],time=train2$time,status=train2$status),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(train2[,rid],time=train2$time,event=train2$status,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + plsRcox')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 4-6 CoxBoost+superpc ########
message(paste0("\n=> 4-6 CoxBoost+superpc", " machine learning algorithm."))

data <- list(x=t(train2[,-c(1,2)]),y=train2$time,censoring.status=train2$status,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 5,
                     n.components=3,
                     min.features=1,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$time,censoring.status=w$status,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + SuperPC')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 4-7 CoxBoost+gbm #######
message(paste0("\n=> 4-7 CoxBoost+gbm", " machine learning algorithm."))

set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,data = train2,distribution = 'coxph',
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 3,n.cores = 1)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,data = train2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 3,n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + GBM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 4-8.CoxBoost+survivalsvm ########
message(paste0("\n=> 4-8.CoxBoost+survivalsvm", " machine learning algorithm."))

fit = survivalsvm(Surv(time,status)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + survival-SVM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 5 plsRcox########
message(paste0("\n=> 5 plsRcox", " machine learning algorithm."))

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train[,-c(1,2)],time=train$time,status=train$status),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(train[,-c(1,2)],time=train$time,event=train$status,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply( cox.res.list$train.vali.list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 6 superpc########
message(paste0("\n=> 6 superpc", " machine learning algorithm."))
data <- list(x=t(train[,-c(1,2)]),y=train$time,censoring.status=train$status,featurenames=colnames(train)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                     n.fold = 10,
                   n.components=3,
                    min.features=1,
                    max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                    compute.preval=TRUE)
rs <- lapply(cox.res.list$train.vali.list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$time,censoring.status=w$status,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 7 GBM ########
message(paste0("\n=> 7 GBM", " machine learning algorithm."))

set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,data = train,distribution = 'coxph',
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~.,data = train,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
rs <- lapply( cox.res.list$train.vali.list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### 8 survivalsvm ########
message(paste0("\n=> 8 survivalsvm", " machine learning algorithm."))

fit = survivalsvm(Surv(time,status)~., data= train, gamma.mu = 2)
rs <- lapply( cox.res.list$train.vali.list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
harmony(rs)
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(time,status)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')
result <- rbind(result,cc)
ml.res [[unique(cc$Model)]]<- fit
rs.res [[unique(cc$Model)]]<- rs

#### plot ########

result$Cindex=round(result$Cindex,3)
result2 <- result
result2=setDT(result2)  # 将数据框转换为data.table
result2.train=result2[result2$ID==names(cox.res.list$train.vali.list)[[train.data.pos]],]
result2.train=as.data.frame(result2.train)
result2.test=result2[result2$ID!=names(cox.res.list$train.vali.list)[[train.data.pos]],]
result2.test=as.data.frame(result2.test)
test.name <- unique(result2.test$ID)
Cindexnums=data.frame(Modle=result2.train$Model,
                      Train=result2.train$Cindex)
names(Cindexnums)[2] <- names(cox.res.list$train.vali.list)[[train.data.pos]]
for (i in 1:length(test.name)) {
  result2.test.i <- result2.test[result2.test$ID == test.name[[i]],]
  Cindextest <- data.frame(Test=result2.test.i$Cindex)
  names(Cindextest) <- test.name[[i]]
  Cindexnums <- cbind(Cindexnums,Cindextest)
}

Cindexnums[,-1] <- apply(Cindexnums[,-1], 2, as.numeric)
Cindexnums$All <- apply(Cindexnums[,2:(length(test.name)+1)], 1, mean)
# 根据C指数排序
Cindexnums <- Cindexnums[order(Cindexnums$All, decreasing = T),]
Cindexnums <- Cindexnums[1:top.num,]
#输出C指数结果
write.csv(Cindexnums,"ML_plot/out_Cindex.csv")
nums <- Cindexnums[, 2:(length(test.name)+2)]%>%as.matrix()
rownames(nums)=Cindexnums$Modle
##热图绘制
Cindex_mat=nums
# 计算每种算法在所有队列中平均C-index
avg_Cindex <- apply(Cindex_mat, 1, mean)
# 对各算法C-index由高到低排序
avg_Cindex <- sort(avg_Cindex, decreasing = T)
# 对C-index矩阵排序
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]
# 保留三位小数
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3))
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)
names(cohort.col) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = cohort.col),
                          show_annotation_name = F)

cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha,
              top_annotation = col_ha,
              col = hm.col,
              rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
              cluster_columns = FALSE, cluster_rows = FALSE, # 不进行聚类，无意义
              show_column_names = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path("ML_plot/c-index.pdf"), width = cellwidth * ncol(Cindex_mat) + 5, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())
return(list(unicox.res=cox.res.list[["candidate.genes"]],cindex.res = Cindexnums,ml.res=ml.res,rs.res=rs.res))
}
}
