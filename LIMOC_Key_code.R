
# set working path and creat directory
tumor.path <- "./analysis"; 
setwd(tumor.path) #create dir
res.path    <- file.path(tumor.path, "Results")
fig.path    <- file.path(tumor.path, "Figures")
data.path    <- file.path(tumor.path, "InputData")

if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }
if (!file.exists(data.path)) { dir.create(data.path) }

# set colors
clust.col <- c("#2EC4B6", "#E71D36", "#FF9F1C")

# load R package
library(MOVICS)
library(wateRmelon)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(consensusMIBC)
library(preprocessCore)
library(ggpubr)
library(survival)
library(survminer)
library(survRM2)
library(ggpmisc)
library(clusterProfiler)
library(enrichplot)
library(estimate)
library(gplots)
library(ClassDiscovery)
library(ComplexHeatmap)
library(GSVA)
library(ggalluvial)

# customized fucntion
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
}

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


# load data ---------------------------------------------------------------

load('LUAD.tcga.Rdata')

# extract mo.data
mo.data <- LUAD.tcga[1:4]

#-------------------------------------------------------#
# identify optimal clustering number (may take a while) #
optk.LUAD <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), # note: the last data is somatic mutation which is a binary matrix
                         try.N.clust = 2:6, # try cluster number from 2 to 6
                         fig.path    = fig.path,
                         fig.name    = "CLUSTER NUMBER OF TCGA-LUAD")

#-------------------------------------------------------------------#
# perform multi-omics integrative clustering with the 7 algorithms #
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "ConsensusClustering",  "CIMLR", "MoCluster"),
                         N.clust     = 3,
                         type        = c("gaussian", "gaussian", "gaussian", "gaussian", "binomial"))

save(moic.res.list, file = file.path(res.path,"moic.res.list.rda"))

load(file=file.path(res.path,"moic.res.list.rda"))

SNF <- moic.res.list$SNF$clust.res
PINSPlus <- moic.res.list$PINSPlus$clust.res
NEMO <- moic.res.list$NEMO$clust.res
COCA <- moic.res.list$COCA$clust.res
ConsensusClustering <- moic.res.list$ConsensusClustering$clust.res
CIMLR <- moic.res.list$CIMLR$clust.res
MoCluster <- moic.res.list$MoCluster$clust.res

moic.res <- cbind.data.frame(
  SNF = paste0("LIMOC",SNF$clust),
  PINSPlus = paste0("LIMOC",PINSPlus$clust),
  NEMO = paste0("LIMOC",NEMO$clust),
  COCA = paste0("LIMOC",COCA$clust),
  ConsensusClustering = paste0("LIMOC",ConsensusClustering$clust),
  CIMLR = paste0("LIMOC",CIMLR$clust),
  MoCluster = paste0("LIMOC",MoCluster$clust))



rownames(moic.res) <- rownames(SNF)

#---------------------------------------------------#
# get consensus clustering from different algorithm #
cmoic.LUAD <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               mapcolor = c("#000004FF", "#56106EFF", "#BB3754FF", "#F98C0AFF", "#FCFFA4FF"),
                               clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627",
                                             "#023E8A", "#9D4EDD"),
                               fig.path      = fig.path,
                               distance      = "euclidean",
                               linkage       = "ward.D2",
                               width         = 7)

cmoic.LUAD$clust.res

table(cmoic.LUAD$clust.res$clust)

# complexheatmap ----------------------------------------------------------

identical(rownames(moic.res),cmoic.LUAD$clust.res$samID)

moic.res$real <- paste0("LIMOC",cmoic.LUAD$clust.res$clust)

moic.res <- moic.res[cmoic.LUAD$clust.dend$order,]

moic.res2 <- moic.res

moic.res2$sample <- rownames(moic.res2)

## 聚类之后，样本顺序
cmoic.LUAD$clust.dend$order

#定义热图注释的颜色
#rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c(brewer.pal(9, "BuGn")[4:6],
         brewer.pal(9, "OrRd")[4:6],
         brewer.pal(9, "Blues")[4:6],
         "#FF857B","#EE6C90","#C25E9F",
         "#D8AC93", "#C98B69","#B15928" ,
         
         "#B49ECC", "#9778B8","#6A3D9A" , 
         "#F18C8D", "#EB5F61","#E31A1C" , 
         '#2EC4B6','#E71D36','#FF9F1C')

colorList=list()

j=0

for(cluster in colnames(moic.res[,1:(ncol(moic.res))])){
  cliLength=length(levels(factor(moic.res[,cluster])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(moic.res[,cluster]))
  cliCol["unknow"]="grey75"
  colorList[[cluster]]=cliCol
}



#绘制热图
ha=HeatmapAnnotation(df=moic.res, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(moic.res))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

draw(Hm, merge_legend=TRUE, heatmap_legend_side="bottom", annotation_legend_side="bottom")


#输出热图
pdf(file="new_heatmap.pdf", width=7, height=5)
draw(Hm, merge_legend=TRUE, heatmap_legend_side="bottom", annotation_legend_side="bottom")
dev.off()


# convert beta value to M value for stronger signal
indata <- mo.data
#indata$meth <- log2(indata$meth / (1 - indata$meth))


# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation



mRNA.col   <-  c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
meth.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")

CNA.col <- c("#6699CC", "white"  , "#FF3C38")

mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col,  meth.col,CNA.col, mut.col)


subtype_ID <- moic.res[,"real",drop=F]

subtype_ID$samID <- rownames(subtype_ID)

subtype_ID <- subtype_ID[,c(2,1)]

names(subtype_ID)[2] <- 'clust'


annCol    <- LUAD.tcga$surv.info[,c("age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage",'OS'), drop = FALSE]

names(annCol) <- c('age','gender','stage','OS')

annCol$stage <- gsub("(A|B|C)$","",annCol$stage)

annCol$stage <- ifelse(annCol$stage=='[Discrepancy]',NA,annCol$stage)
annCol$stage <- ifelse(annCol$stage=='',NA,annCol$stage)

annCol$OS <- ifelse(annCol$OS==1,'dead','alive')

annCol$Barcode <- substring(rownames(annCol),1,12)



# age ---------------------------------------------------------------------

median(annCol5$age,na.rm = T)

annCol5_1 <- annCol5[annCol5$age>66,]


fitd <- survdiff(Surv(futime, fustat) ~ cluster,
                 data      = annCol5_1,
                 na.action = na.exclude)

p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(futime, fustat) ~ cluster,
               data      = annCol5_1,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)


# kaplan-meier curve
p_movics_1 <- ggsurvplot(fit               = fit,
                         conf.int          = FALSE,
                         risk.table        = F,
                         risk.table.col    = "strata",
                         palette           =c("#2EC4B6", "#E71D36", "#FF9F1C"),
                         data              = annCol5_1,
                         size              = 1,
                         #xlim              = c(0,120),
                         #break.time.by     = 20,
                         legend.title      = "",
                         pval              = T, # 不计算p值，改为手动添加
                         surv.median.line  = "hv",
                         xlab              = "Time (month)",
                         ylab              = "Survival probability",
                         risk.table.y.text = FALSE)



p_movics_1

ggsave(file="TCGA_p_movics_1_age_big_66.pdf",width = 5,height = 5)



median(annCol5$age,na.rm = T)

annCol5_2 <- annCol5[annCol5$age<66,]


fitd <- survdiff(Surv(futime, fustat) ~ cluster,
                 data      = annCol5_2,
                 na.action = na.exclude)

p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(futime, fustat) ~ cluster,
               data      = annCol5_2,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

# kaplan-meier curve
p_movics_2 <- ggsurvplot(fit               = fit,
                         conf.int          = FALSE,
                         risk.table        = FALSE,
                         risk.table.col    = "strata",
                         palette           = c("#2EC4B6", "#E71D36", "#FF9F1C"),
                         data              = annCol5_2,
                         size              = 1,
                         #xlim              = c(0,120),
                         #break.time.by     = 20,
                         legend.title      = "",
                         pval              = TRUE, # 不计算p值，改为手动添加
                         surv.median.line  = "hv",
                         xlab              = "Time (month)",
                         ylab              = "Survival probability",
                         risk.table.y.text = FALSE)

p_movics_2

ggsave(file="TCGA_p_movics_2_age_small_66.pdf",width = 5,height = 5)




# Gender ------------------------------------------------------------------


annCol5_3 <- annCol5[annCol5$gender=='MALE',]


fitd <- survdiff(Surv(futime, fustat) ~ cluster,
                 data      = annCol5_3,
                 na.action = na.exclude)

p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(futime, fustat) ~ cluster,
               data      = annCol5_3,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)



# kaplan-meier curve
p_movics_3 <- ggsurvplot(fit               = fit,
                         conf.int          = FALSE,
                         risk.table        = F,
                         risk.table.col    = "strata",
                         palette           = c("#2EC4B6", "#E71D36", "#FF9F1C"),
                         data              = annCol5_3,
                         size              = 1,
                         #xlim              = c(0,120),
                         #break.time.by     = 20,
                         legend.title      = "",
                         pval              = T, # 不计算p值，改为手动添加
                         surv.median.line  = "hv",
                         xlab              = "Time (month)",
                         ylab              = "Survival probability",
                         risk.table.y.text = FALSE)


p_movics_3

ggsave(file="TCGA_p_movics_3_male.pdf",width = 5,height = 5)


annCol5_4 <- annCol5[annCol5$gender=='FEMALE',]


fitd <- survdiff(Surv(futime, fustat) ~ cluster,
                 data      = annCol5_4,
                 na.action = na.exclude)

p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(futime, fustat) ~ cluster,
               data      = annCol5_4,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

# kaplan-meier curve
p_movics_4 <- ggsurvplot(fit               = fit,
                         conf.int          = FALSE,
                         risk.table        = FALSE,
                         risk.table.col    = "strata",
                         palette           = c("#2EC4B6", "#E71D36", "#FF9F1C"),
                         data              = annCol5_4,
                         size              = 1,
                         #xlim              = c(0,120),
                         #break.time.by     = 20,
                         legend.title      = "",
                         pval              = T, # 不计算p值，改为手动添加
                         surv.median.line  = "hv",
                         xlab              = "Time (month)",
                         ylab              = "Survival probability",
                         risk.table.y.text = FALSE)

p_movics_4

ggsave(file="TCGA_p_movics_4_female.pdf",width = 5,height = 5)





#------------------#
# compare survival #

surv.LUAD <- compSurv(moic.res         = cmoic.LUAD,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      p.adjust.method  = "none",
                      clust.col        = clust.col,
                      fig.path         = fig.path,
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")


surv.info2 <- LUAD.tcga$surv.info[,c("PFI","PFI.time")]

colnames(surv.info2) <- c("fustat","futime")

surv.LUAD <- compSurv(moic.res         = cmoic.LUAD,
                      surv.info        = surv.info2,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(3,5), # estimate 5 and 10-year survival
                      p.adjust.method  = "none",
                      clust.col        = clust.col,
                      fig.path         = fig.path,
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC PFS")

#---------------------------------#
# mutational frequency comparison #
mut.LUAD2 <- compMut(moic.res     = cmoic.LUAD,
                     mut.matrix   = LUAD.tcga$mut, # binary somatic mutation matrix
                     doWord       = TRUE, # generate table in .docx format
                     doPlot       = TRUE, # draw OncoPrint
                     freq.cutoff  = 0.03, # keep those genes that mutated in at least 3% of samples
                     p.cutoff     = 0.05, # keep those genes with nominal p value < 0.25 to draw OncoPrint
                     p.adj.cutoff = 0.25, # keep those genes with adjusted p value < 1 to draw OncoPrint
                     innerclust   = TRUE, # perform clustering within each subtype
                     annCol       = annCol, # same annotation for heatmap
                     annColors    = annColors, # same annotation color for heatmap
                     width        = 8,
                     height       = 5,
                     fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                     tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION",
                     res.path     = res.path,
                     fig.path     = fig.path)

write.csv(mut.LUAD2,file = 'mut.LUAD2.csv')


# Barplot -----------------------------------------------------------------

names(mut.LUAD2)[1] <- 'Gene'

mut.LUAD2$pvalue <- as.numeric(mut.LUAD2$pvalue)

mut.LUAD3 <- mut.LUAD2[order(mut.LUAD2$pvalue),][1:20,]

mut.LUAD3$Gene <- factor(mut.LUAD3$Gene,
                         levels = mut.LUAD3$Gene[20:1])
# 绘制条形图
ggplot(data = mut.LUAD3, aes(x = -log(pvalue), y = Gene)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(x = "p-value", y = "Mutation Gene") +
  theme_bw()+
  theme(axis.text = element_text(color = "black")) 


ggsave('../Revised/TP53_Rank_TCGA.pdf',width = 5,height = 5)



#-------------------------------#
# compare total mutation burden #
tmb.LUAD <- compTMB(moic.res     = cmoic.LUAD,
                    maf          = LUAD.tcga$maf,
                    nonSyn       = NULL,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.path     = fig.path,
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

#---------------------------------#
# compare fraction genome altered #
# change column names of segment data
segment2 <- segment[,setdiff(colnames(segment),"num.mark")]
colnames(segment2) <- c("sample","chrom","start","end","value")

# compare FGA, FGG, and FGL
fga.LUAD <- compFGA(moic.res     = cmoic.LUAD,
                    segment      = segment2,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.3, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    clust.col    = clust.col,
                    width = 8,height = 3,
                    fig.path     = fig.path,
                    fig.name     = "BARPLOT OF FGA")


# customize the factor level for pstage
annCol$stage <- factor(annCol$stage, levels = c("Stage I","Stage II","Stage III","Stage IV"))
annCol
# agreement comparison (support up to 6 classifications include current subtype)
agree.brca <- compAgree(moic.res  = cmoic.LUAD,
                        subt2comp = annCol[,c("OS","stage")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")

