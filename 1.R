library(devtools)
library(EZGen)
library(MOVICS)
install.packages('tibble')
usethis::use_r("heatmap_plot")
devtools::load_all()
options(error = NULL)
#filter_Download  <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
filter_Download(project = 'TCGA-CaAD',type='mrna')
#filter_DEGs  <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
expr <- read.csv("G:/importance/undergraduated/数据集整理/bulk data（clinical）/NEPC/NE.csv",row.names = 1)
cli <- read.csv("G:/importance/undergraduated/2_CRPC/第二版/NE亚型/bulk_group.csv")
cli <- subset(cli,cli$subtype%in%c('NEC-NE','EMT-NE'))
cli$subtype <- factor(cli$subtype,levels = c('NEC-NE','EMT-NE'))
ordered_sampleIDs <- cli$sampleID[order(match(cli$subtype, levels(cli$subtype)))]

expr <- expr[,colnames(expr)%in%ordered_sampleIDs]
expr <- expr[, match(ordered_sampleIDs, colnames(expr))]
table(cli$subtype)

filter_DEGs(expr.data ='',
            data.type ='mRNA',
            deg.method='DESeq2',
            Pvalue    =0.05,
            log2FC    =2,
            TCGA      =T,
            tumor.num =50,
            normal.num=49,
            color1    ='#CE5C69',
            color2    ='#5770A6',
            title     ='TCGA',
            key.genes =NA)

#filter_merge <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
filter_merge (DEGs1    = data1,
              DEGs2    = data2,
              genelist = genes,
              color    = c('#CE5C69','#5770A6','#E0A980'),
              labels   = c("DEGs of GSE6752", "DEGs of GSE48403 ",'Disulfidptosis-related genes'),
              save.path= 'G:/importance/undergraduated/12_R包开发/ctDNA/VENE.TIF')
#evaluate_subtype <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
library(readr)
expr <- read.csv('NE.csv',row.names = 1)
expr <- read_csv('PCaProfilter_matrix（未去批次效应）.csv')
colnames(expr)[1] <- 'gene'

group <- read.csv('cli_all.csv')
group <- subset(group,group$Molecular.Type=='NEPC')
NE.ID <- group$ID

expr <- expr[,colnames(expr)%in% c('gene',NE.ID)]
expr <- as.data.frame(expr)

row.names(expr) <- expr[,1]
expr <- expr[,-1]

expr <- read.csv("G:/importance/undergraduated/肾癌数据集/bulk/TCGA/TCGA.csv",row.names = 1)

opt.res <- evaluate_subtype(expr= expr,
                            genelist1.file= 'Palmitoylation.csv',
                            #genelist2.file= 'EMT.csv',
                            genelist1.name= 'Palmitoylation',
                            #genelist2.name= 'EMT_R',
                            clust.num= 2:5,
                            filename = 'optimal_number_cluster')
summary(data.list[[1]])
sum(is.na(data.list[[1]]))

dim(data.list[[1]])
colnames(data.list[[1]])
str(data.list[[1]])
#create_subtype <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
create_subtype(data        = data.list,
               clust.num   = 3,
               methodslist = list("IntNMF", "CIMLR", "PINSPlus", "NEMO",  "MoCluster", "LRAcluster", "ConsensusClustering"))


create_subtype(data        = data.list,
               clust.num   = 3,
               methodslist = list('iClusterBayes'))

#subtype_harmony <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
harmony.res <- subtype_harmony(clust.res.list = clust.res.list,
                               map.col        = c("#364888","#5770A6", "#A281B1",'#E0A980','#F0A8B9',"#CE5C69",'#b30c2a'),
                               clust.col      = c("#5770A6", "#A281B1", "#CE5C69", "#BDD5EA",'#E0A980'),
                               linkage        = "average",
                               showID         = FALSE,
                               file.name      = "homoney_heatmap",
                               width          = 5,
                               height         = 5.5)
saveRDS(harmony.res,'harmony.res.rds')

#subtype_heatmap <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
subtype_heatmap (data          = data.list,
                 method        = 'CIMLR',#('all',"IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering")
                 genelist1.name= 'Prol_E',
                 #genelist2.name= 'EMT_R',
                 feat.num      = 5,
                 col1          = c("#5770A6", 'white', "#CE5C69"),
                 #col2          = c("#5770A6", "white"  , "#CE5C69"),
                 clust.col     = c("#5770A6", "#A281B1", "#CE5C69", "#BDD5EA",'#E0A980'),
                 genelist.title= c("Palmitoylation relative genes"),
                 legend.name   = c("mRNA"),
                 clust.dend    = NULL, # Do not display the dendrogram
                 show.rownames = c(F, FALSE), # Gene name display settings
                 show.colnames = FALSE, # Sample name display settings
                 annCol        = NULL, # Do not annotate samples
                 annColors     = NULL, # Do not set annotation colors
                 width         = 10, # Width of each sub-heatmap
                 height        = 5, # Height of each sub-heatmap
                 file.save.path= getwd(),
                 file.name     = "subtype_heatmap")

# extract PAM50, pathologic stage and age for sample annotation
#annCol    <- surv.info[,c("PAM50", "pstage", "age"), drop = FALSE]

# generate corresponding colors for sample annotation
#annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age),
#                                                           median(annCol$age),
#                                                           max(annCol$age)),
#                                                colors = c("#0000AA", "#555555", "#AAAA00")),
#                  PAM50  = c("Basal" = "blue",
#                             "Her2"   = "red",
#                             "LumA"   = "yellow",
#                             "LumB"   = "green",
#                             "Normal" = "black"),
#                  pstage = c("T1"    = "green",
#                             "T2"    = "blue",
#                             "T3"    = "red",
#                             "T4"    = "yellow",
#                             "TX"    = "black"))

#subtype_KM <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
cli.data <- read.csv("G:/importance/undergraduated/肾癌数据集/bulk/TCGA/cli.csv",row.names = 1)
cli.data$time <- cli.data$time /30
subtype_KM (clust.res.list   =clust.res.list,
            method           = 'CIMLR', #('all',"IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering")
            cli.data         = cli.data,
            subtype.list     = c('CS1','CS3'),
            cli.limit        = NULL,
            break.time.by    = 12,
            clust.col        = c("#5770A6", "#A281B1", "#CE5C69", "#BDD5EA",'#E0A980'),
            p.adjust.method  = "BH",
            surv.median.line = "none",
            file.save.path   = getwd(),
            file.name        = 'subtype_KM')

#subtype_DEGs <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
DEGs.res  <- subtype_DEGs(clust.res.list= clust.res.list,
                          method        = 'CIMLR', #('all',"IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering")
                          expr.data     = 'G:/importance/undergraduated/肾癌数据集/bulk/TCGA/TCGA.csv',
                          subtype.list  = c('CS1','CS3'),
                          subtype.ctrl  = 'CS1',
                          deg.method    = 'DESeq2',
                          Pvalue        = 0.05,
                          log2FC        = 2,
                          col1          = '#CE5C69',
                          col2          = '#5770A6',
                          title         = 'CS3 vs CS1',
                          key.genes     = NULL)
saveRDS (DEGs.res,'DEGs.res.rds')
plot(DEGs.res$vp)
#subtype_EA <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
subtype_EA (DEGs.res = DEGs.res,
             data.base= 'KEGG',#c('GO','KEGG')
             method   = 'Taiji',#c('basic','circos','dendrogram','Taiji','bubble')
             lfc.col  = c('#CE5C69','#5770A6'),
             zsc.col  = c('#CE5C69', '#EDF5F7', '#5770A6'),
             cat.col  = c('#CE5C69','#A281B1','#5770A6'))

#subtype_GSVA <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
subtype_GSVA (DEGs.res       = DEGs.res,
              data.base      = 'GO',#c('GO','KEGG')
              subtype.list   = c('CS1','CS3'))#"exp/ctrl"

#subtype_GSEA <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
subtype_GSEA (DEGs.res       = DEGs.res,
              data.base      = 'KEGG',#c('GO','KEGG')
              geneset.KEGG   = 'hsa04514',
              geneset.GO     = 'GO:0005515',
              Pvalue         = 0.5,
              col            = '#CE5C69')

#subtype_TME <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
TME.res.list <- subtype_TME(clust.res.list = clust.res.list,
                            clust.method   = 'iClusterBayes',#('all',"IntNMF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "iClusterBayes", "SNF", "ConsensusClustering")
                            expr           = expr, # 基因表达矩阵，可以是一个数值型矩阵或数据框，行名为基因，列名为样本名
                            subtype.list   = c('CS1','CS3'),
                            method         = "all", # c('all','xCell','MCPcounter', 'EPIC',  'CIBERSORT', 'ips', 'estimate', 'svr', 'lsei', 'TIMER', 'quantiseq')
                            perm           = 100, # 统计分析的置换次数（建议≥100次）影响'CIBERSORT'方法。
                            reference      = NA,
                            indications    = NULL)
saveRDS(TME.res.list,'TME.res.list.rds')

#CIBERSORT_plot <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
cibersort.plot.list <- cibersort_plot(TME.res.list=TME.res.list,
                                       palette.col=NA,
                                       heatmap.col=c("#5770A6", 'white',"#CE5C69"))
cibersort.plot.list[4]

#box_plot <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
box_plot(expr         ='WCDT-MCRPC_mrna_expr_count.csv',
         group        =DEGs.res$clust.res,
         gene.list =c('AR','KLK3','KLK2','SPDEF','CHGA','SYP','ENO2','NCAM1','SPIC','NEUROD1','CMTM2','SLC7A11'),
         subtype.list = c('CS1','CS3'),
         clust.col    =c ("#5770A6","#CE5C69"))

#cor_plot <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
cor_plot (expr='ML_data/mRNA_matrix/GSE6752_mRNA_matrix.csv',
          group=NA,
          gene.list=c('AR','KLK3','STEAP4','SPDEF','CHGA','ENO2','SYP','NCAM1','AMIGO2','POLRMT','DUS3L'),
          subtype.list  = c('CS1','CS3'),
          col=c("#5770A6",'white',"#CE5C69"),
          file.name='GSE6752-CRPC')
dev.off()
#subtype_drug <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
subtype_drug(expr  =DEGs.res$expr,
             group =DEGs.res$group,
             col   =c("#CE5C69","#5770A6"),
             subtype.list = c('CS1','CS3'),
             gene  ='SLC7A11',
             drug  ='oxaliplatin')

#create_model <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
train_data <- read.csv('ML_data/WCDT_mCRPC.csv')
valid_data_1 <- read.csv('ML_data/MCTP_mCRPC.csv')
valid_data_2 <- read.csv('ML_data/SU2C_mCRPC.csv')
prog.data.list <- list(train_data,valid_data_1,valid_data_2)
names(prog.data.list) <- c("WCDT_MCRPC", "MCTP_2012",'SU2C_2019')
saveRDS (prog.data.list,'prog.data.list.rds')

gene.list <- DEGs.res$DEGs$Gene_symbol
gene.list <- union(gene.list, 'SLC7A11')
gene.list <- read.csv('NEPC/candidate_genelist.csv')
colnames(gene.list)[1] <- 'gene'

gene.list <- read.csv('G:/importance/undergraduated/2_CRPC/第二版/DEG、WGNA/candidate.csv')[,1]
gene.list <- read.csv('intersection_up_sample.csv')[,1]

prog.data <- prog.data.list[[1]]
prog.data <- as.data.frame(lapply(prog.data, function(x) pmax(x, 0)))


saveRDS(prog.data.list,'G:/importance/undergraduated/数据集整理/bulk data（clinical）/prog.data.list.rds')

prog.data.list <- readRDS("G:/importance/undergraduated/数据集整理/bulk data（clinical）/prog.data.list.rds")

prog.data.list[[10]] <- NULL
model.res.list <- create_model(prog.data.list=prog.data.list,
                               train.data.pos=1,
                               gene.list=NEP100,
                               unicox_pcutoff=0.1,
                               top.num = 100,
                               method='out_RSF',#c('all','out_RSF')
                               hm.col=c("#5770A6", "#FFFFFF", "#CE5C69"),
                               cohort.col=c('#b30c2a','#ce5c69','#e0a980','#f4c889','#bdd5a3','#519981','#8ba1c6','#5770a6','#a281b1'))
saveRDS (model.res.list,'model.res.list.rds')

#filter_model <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
#查看最佳模型
best.model <- model.res.list[["cindex.res"]][2,]$Modle
#读取最佳模型结果
best.model.res1 <- model.res.list[["ml.res"]][['CoxBoost']]
best.model.res2 <- model.res.list[["ml.res"]][[best.model]]
#读取riskscore
best.model.rs <- model.res.list[["rs.res"]][[best.model]]
#出图
filter_model (model1=best.model.res1,
              model2=best.model.res2,
              model1.name='CoxBoost',
              model2.name='RSF',
              rs.data=best.model.rs,
              col=c('#5770A6','#CE5C69'))

#create_nomogram <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
cli.data <- read.csv('TCGA_cli2.csv',row.names = 1)
cli.data$time <- cli.data$time /30
rs.data=best.model.rs$WCDT_MCRPC
cli.data <- tibble::rownames_to_column(cli.data, var = "ID")
rs.data <- tibble::rownames_to_column(rs.data, var = "ID")
rs.data <- rs.data[,c(1,4)]
# 使用 dplyr 进行合并
mydata <- full_join(cli.data, rs.data, by = "ID")
row.names(mydata) <- mydata[,1]
mydata <- mydata[,-1]
dd <- datadist(mydata)
options(datadist = "dd")

create_nomogram (cli.data=cli.data,
                 cont.var=c('RS'),
                 cate.var=c('age','meta','race'),
                 Pvalue=0.10,
                 rate.range=c(1,3,5),
                 col=c("#5770A6", "#A281B1", "#CE5C69"))

#DCA <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-
library(ggDCA)
library(rms)
m1 <- cph(Surv(time,status==1)~age+meta+RS  ,x=T,y=T,data=cli.data, surv=T)
m2 <- cph(Surv(time,status==1)~meta  ,x=T,y=T,data=cli.data, surv=T)
m3 <- cph(Surv(time,status==1)~RS,x=T,y=T,data=cli.data, surv=T)
m4 <- cph(Surv(time,status==1)~age,x=T,y=T,data=cli.data, surv=T)

d_train <- dca(m1,m2,m3,m4,
               model.names =c('Nomograms Model','Meta','RS','Age'))#修改图例名称
ggplot(d_train)
