rm(list = ls())
options(stringsAsFactors = F)
library(org.Hs.eg.db)
library(GEOquery)
# stromal cells from tumor-draining lymph nodes
gset <- getGEO('GSE113249', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)       ## 平台文件  
a=gset[[1]] #
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度
head(dat)
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData
colnames(dat)=pd$title
dat1=dat

gset <- getGEO('GSE113250', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)       ## 平台文件  
a=gset[[1]] #
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度
head(dat)
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData
colnames(dat)=pd$title
dat2=dat

dat=cbind(dat1,dat2)

dat[1:4,1:4] 
boxplot(dat)
# GPL21163    Agilent-074809 SurePrint G3 Mouse GE v2 8x60K Microarray [Probe Name version]
library(GEOquery)
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL21163', destdir=".")
colnames(Table(gpl))  
head(Table(gpl)[,c(1,6)]) ## you need to check this , which column do you need
probe2gene=Table(gpl)[,c(1,6)] 
save(probe2gene,file='probe2gene.Rdata')

load(file='probe2gene.Rdata')
ids=probe2gene 
head(ids)
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat),]

dat[1:4,1:4]   
dat=dat[ids$probe_id,] 

ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
dim(dat)

save(dat,file = 'step1-output.Rdata')
load(file = 'step1-output.Rdata')



rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(GEOquery)
load(file = 'step1-output.Rdata')
colnames(dat)
dat[1:4,1:4]

library(pheatmap)
ccg=trimws(strsplit('Cdc25c, Bub1, Ttk,  Cdk1',',')[[1]])
mrg=trimws(strsplit('Vcam1, Arhgap5, Cxcr3, Ccr2',',')[[1]])
mat=dat[c(ccg,mrg),3:4]
mat
pheatmap(mat)
pheatmap(dat[rownames(dat)[grepl('^Ig',rownames(dat))],])



rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(GEOquery)
load(file = 'step1-output.Rdata')
colnames(dat)
dat[1:4,1:4]

rownames(dat)=toupper(rownames(dat))

m=dat[,2]+dat[,1]/2
logFC1=log2(dat[,2]/dat[,1]);fivenum(logFC1)
logFC2=dat[,2]-dat[,1];fivenum(logFC2)
plot(logFC1,m)
plot(logFC2,m)
plot(logFC1,logFC2)

# 这里的 logFC1 和 logFC2 差异很小
geneList=scale( logFC2 )
# geneList=scale( dat[,2]/dat[,1] )
names(geneList)=rownames(dat)
geneList=sort(geneList,decreasing = T)
head(geneList)

#browseVignettes("gskb")  
library(gskb) 
data(mm_GO)
(g1=mm_GO$GO_BP_MM_HUMORAL_IMMUNE_RESPONSE)
(g2=mm_GO$GO_BP_MM_LEUKOCYTE_CHEMOTAXIS)
g1=g1[g1 %in% rownames(dat)]
g2=g2[g2 %in% rownames(dat)]
library(pheatmap)
pheatmap(dat[g1 ,1:2])
pheatmap(dat[g2 ,1:2])

library(ggplot2)
library(clusterProfiler)
library(GSEABase) # BiocManager::install('GSEABase')
# geneset <- read.gmt('h.all.v6.2.symbols.gmt')
# Then we know geneset is a simple data.frame
# first column is name of geneset, second column is gene symbol
geneset=data.frame(ont=c(rep('HUMORAL_IMMUNE_RESPONSE',length(g1)),
                         rep('LEUKOCYTE_CHEMOTAXIS',length(g2))),
                   gene=c(g1,g2))
egmt <- GSEA(geneList, 
             TERM2GENE=geneset, 
             verbose=T)
head(egmt)
gseaplot(egmt,'HUMORAL_IMMUNE_RESPONSE') 
gseaplot(egmt,'LEUKOCYTE_CHEMOTAXIS') 


##################R package for Y叔################
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(GEOquery)
load(file = 'step1-output.Rdata')
colnames(dat)
dat[1:4,1:4]

rownames(dat)=toupper(rownames(dat))

m=dat[,2]+dat[,1]/2
logFC1=log2(dat[,2]/dat[,1]);fivenum(logFC1)
logFC2=dat[,2]-dat[,1];fivenum(logFC2)
plot(logFC1,m)
plot(logFC2,m)
plot(logFC1,logFC2)

geneList=scale( logFC2 )
names(geneList)=rownames(dat)
geneList=sort(geneList,decreasing = T)
head(geneList)

s2g=select(org.Hs.eg.db,names(geneList),
           'ENTREZID','SYMBOL')

head(s2g)
s2g=s2g[!is.na(s2g$ENTREZID),]
s2g=s2g[s2g$SYMBOL %in%  names(geneList),]
geneList=geneList[s2g$SYMBOL]  
names(geneList)=s2g$ENTREZID  
head(geneList)
tail(geneList)

library(ggplot2)
library(clusterProfiler)
library(GSEABase) 
## 这里比较耗时：
go_bp_gsea <- gseGO(geneList     = geneList, 
                    OrgDb='org.Hs.eg.db',
                    ont = 'BP',
                    nPerm        = 1000,
                    minGSSize    = 10,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)

gseaplot(go_bp_gsea, geneSetID = rownames(go_bp_gsea[1,]))
gseaplot(go_bp_gsea, geneSetID = "GO:0006959")
gseaplot(go_bp_gsea, geneSetID = "GO:0030595")
