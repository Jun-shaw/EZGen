fliter_unicox <- function(prog.data.list=prog.data.list,
                          train.data.pos=train.data.pos,
                          gene.list=gene.list,
                          Pvalue=unicox_pcutoff) #单因素回归的筛选阈值
{
# Replace '-' in column names with '.'

prog.data.list <- lapply(prog.data.list,function(x){
  colnames(x) = gsub('-','.',colnames(x))
  return(x)})
prog.data.list <- lapply(prog.data.list,function(x){
  colnames(x) = gsub('_','.',colnames(x))
  return(x)})
gene.list <-gsub('-','.',gene.list)
gene.list <- gsub("_",".",gene.list)
# Data preprocessing

# Matching candidate genes to genes in each cohort
common.genes = c('ID', 'time', 'status',gene.list)

for (i in 1:length(prog.data.list)) {
  common.genes = intersect(common.genes, colnames(prog.data.list[[i]]))
}

message(paste0('---the number of the all genes is ', length(gene.list),' ---'))
message(paste0('---the number of the common genes in all cohorts is ', length(common.genes)-3,' ---'))

returnIDtoRS = function(rs.table.list, rawtableID){

  for (i in names(rs.table.list)) {
    rs.table.list[[i]] $ID = rawtableID[[i]]$ID
    rs.table.list[[i]] = rs.table.list[[i]] %>% dplyr::select('ID', everything())
  }

  return(rs.table.list)
}

train.data <- prog.data.list[[train.data.pos]]

if(identical(c('ID','time','status'),colnames(train.data)[1:3])&
   length(gene.list)>0 &
   identical(c('ID','time','status'),common.genes[1:3]) &
   length(common.genes)>3 ){

  message('--- Data preprocessing ---')

  prog.data.list <- lapply(prog.data.list, function(x){x[, common.genes]})

  prog.data.list <- lapply(prog.data.list,function(x){
    x[,-c(1:3)] <- apply(x[,-c(1:3)],2,as.numeric)
    return(x)})

  prog.data.list <- lapply(prog.data.list,function(x){
    x[,c(1:2)] <- apply(x[,c(1:2)],2,as.factor)
    return(x)})

  prog.data.list <- lapply(prog.data.list,function(x){
    x[,c(2:3)] <- apply(x[,c(2:3)],2,as.numeric)
    return(x)})

  prog.data.list <- lapply(prog.data.list,function(x){
    x <- x[!is.na(x$time)&!is.na(x$status), ]
    return(x)})

  prog.data.list <- lapply(prog.data.list,function(x){
    x <- x[x$time > 0, ]
    return(x)})
  # use the mean replace the NA
  prog.data.list <- lapply(prog.data.list,function(x){
    x[,-c(1:3)] <-apply(x[,-c(1:3)],2,function(x){
      x[is.na(x)]=mean(x,na.rm = T)
      return(x)
    })

    return(x)})

  if(is.null(Pvalue)){
    Pvalue =0.05
  } else {
    Pvalue = Pvalue
  }

      display.progress = function (index, totalN, breakN=20) {
        if ( index %% ceiling(totalN/breakN)  ==0  ) {
          cat(paste(round(index*100/totalN), "% ", sep=""))
        }
      }

      ###############unicox#######
      print("Stating the univariable cox regression")

      train.data <- prog.data.list[[train.data.pos]]
      unicox <- data.frame()

      for(i in 1:ncol(train.data[,4:ncol(train.data)])){

        display.progress(index = i, totalN = ncol(train.data[,4:ncol(train.data)]))
        gene <- colnames(train.data[,4:ncol(train.data)])[i]
        tmp <- data.frame(expr = as.numeric(train.data[,4:ncol(train.data)][,i]),
                          futime = train.data$time,
                          fustat = train.data$status,
                          stringsAsFactors = F)
        cox <- coxph(Surv(futime, fustat) ~ expr, data = tmp)
        coxSummary <- summary(cox)
        unicox <- rbind.data.frame(unicox,
                                   data.frame(gene = gene,
                                              HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                              z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                              pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                              lower = as.numeric(coxSummary$conf.int[,3][1]),
                                              upper = as.numeric(coxSummary$conf.int[,4][1]),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
      }

      print("Finished the univariable cox regression")

      candidate.genes <- unicox[which(unicox$pvalue < Pvalue), "gene"]

}
    message(paste0('---the number of the final unicox filtered candidate genes is ', length(candidate.genes),' ---'))

    common.feature =c('ID', 'time', 'status',candidate.genes)
    colnames(train.data)[1:20]
    train.data =  as.data.frame(train.data)[,common.feature]
    # 假设 prog.data.list 是一个包含数据框的列表，common.feature 是一个包含特征名称的向量
    train.vali.list <- lapply(prog.data.list, function(x) { new_data <- x[, common.feature[-1]];rownames(new_data) <- x[,1];return(new_data)})
    return(list(train.data=train.data,train.vali.list=train.vali.list,candidate.genes=candidate.genes,unicox.res=unicox))

}
