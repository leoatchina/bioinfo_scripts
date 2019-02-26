rm(list = ls())
options(stringsAsFactors = F)
library(org.Hs.eg.db)
library(GEOquery)

# link https://www.ncbi.nlm.nih.gov/pubmed/30643287 
# stromal cells from tumor-draining lymph nodes
gset <-  getGEO('GSE113249', destdir="~/raw_data",
               AnnotGPL = F,     # 注释文件
               getGPL = F)       # 平台文件  
a=gset[[1]] #
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度
head(dat)
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData
colnames(dat)=pd$title
dat1=dat

gset <- getGEO('GSE113250', destdir="~/raw_data",
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
dim(dat)
boxplot(dat)
# GPL21163 Agilent-074809 SurePrint G3 Mouse GE v2 8x60K Microarray [Probe Name version], from article
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL21163', destdir="~/raw_data")
colnames(Table(gpl))   # Table is not somatic function
head(Table(gpl)[,c(1,6)]) ## you need to check this , which column do you need
probe2gene=Table(gpl)[,c(1,6)] 
save(probe2gene,file='~/data/probe2gene.Rdata')
save(dat, file = '~/data/dat.Rdata')



##### load Rdata, go on analysis
load(file='~/data/probe2gene.Rdata')
load(file='~/data/dat.Rdata')
ids=probe2gene 
head(ids)
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat),]
dat=dat[ids$probe_id,] 

ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果,所要要先排序
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
dim(dat)
save(dat,file = '~/data/step1-output.Rdata')



################ pheatmap
rm(list = ls())
load(file = '~/data/step1-output.Rdata')
colnames(dat)
dat[1:4,1:4]

library(pheatmap)
# 下面这两行，有点懒了, 另外这些基因从哪里来的
ccg=trimws(strsplit('Cdc25c, Bub1, Ttk,  Cdk1',',')[[1]])
mrg=trimws(strsplit('Vcam1, Arhgap5, Cxcr3, Ccr2',',')[[1]])
mat=dat[c(ccg,mrg),3:4]
mat
pheatmap(mat)
pheatmap(dat[rownames(dat)[grepl('^Ig',rownames(dat))],])



############### gsea
rm(list = ls())
load(file = '~/data/step1-output.Rdata')
rownames(dat)=toupper(rownames(dat))
m=dat[,2]+dat[,1]/2   # 为什么除以2
m2 = dat[,2]+dat[,1]
logFC1=log2(dat[,2]/dat[,1])
fivenum(logFC1)
logFC2=dat[,2]-dat[,1]
fivenum(logFC2)
par(mfrow=c(2,2))
# 形状基本一样
plot(logFC1,m)
plot(logFC2,m)
plot(logFC1,m2)
plot(logFC2,m2)
par(mfrow=c(1,1))
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
pheatmap(dat[g1 ,1:2])
pheatmap(dat[g2 ,1:2])





# GSEA原理
# 给定一个排序的基因表L和一个预先定义的基因集S (比如编码某个代谢通路的产物的基因, 基因组上物理位置相近的基因，或同一GO注释下的基因)，GSEA的目的是判断S里面的成员s在L里面是随机分布还是主要聚集在L的顶部或底部。这些基因排序的依据是其在不同表型状态下的表达差异，若研究的基因集S的成员显著聚集在L的顶部或底部，则说明此基因集成员对表型的差异有贡献，也是我们关注的基因集。
# GSEA计算中几个关键概念：
# 计算富集得分 (ES, enrichment score). ES反应基因集成员s在排序列表L的两端富集的程度。计算方式是，从基因集L的第一个基因开始，计算一个累计统计值。当遇到一个落在s里面的基因，则增加统计值。遇到一个不在s里面的基因，则降低统计值。每一步统计值增加或减少的幅度与基因的表达变化程度（更严格的是与基因和表型的关联度）是相关的。富集得分ES最后定义为最大的峰值。正值ES表示基因集在列表的顶部富集，负值ES表示基因集在列表的底部富集。
# 评估富集得分(ES)的显著性。通过基于表型而不改变基因之间关系的排列检验 (permutation test)计算观察到的富集得分(ES)出现的可能性。若样品量少，也可基于基因集做排列检验 (permutation test)，计算p-value。
# 多重假设检验矫正。首先对每个基因子集s计算得到的ES根据基因集的大小进行标准化得到Normalized Enrichment Score (NES)。随后针对NES计算假阳性率。（计算NES也有另外一种方法，是计算出的ES除以排列检验得到的所有ES的平均值）
# Leading-edge subset，对富集得分贡献最大的基因成员。
# GSEA分析

# 然后进行GSEA分析：
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



#
library(GO.db)  
ls("package:GO.db")  
library(org.Hs.eg.db)
eg2go=toTable(org.Hs.egGO2ALLEGS)
head(eg2go)
data=eg2go[eg2go$go_id=='GO:0030595',]
columns(org.Hs.eg.db)
g11=select(org.Hs.eg.db,data$gene_id,'SYMBOL', 'ENTREZID')[,2]
g11=g11[g11 %in% rownames(dat)]

data=eg2go[eg2go$go_id=='GO:0006959',]
columns(org.Hs.eg.db)
g22=select(org.Hs.eg.db,data$gene_id,'SYMBOL', 'ENTREZID')[,2]
g22=g22[g22 %in% rownames(dat)]
pheatmap(dat[g11 ,1:2])
pheatmap(dat[g22 ,1:2])
geneset=data.frame(ont=c(rep('HUMORAL_IMMUNE_RESPONSE',length(g11)),
                         rep('LEUKOCYTE_CHEMOTAXIS',length(g22))),
                   gene=c(g11,g22))
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             verbose=T)
head(egmt)
gseaplot(egmt,'HUMORAL_IMMUNE_RESPONSE') 
gseaplot(egmt,'LEUKOCYTE_CHEMOTAXIS') 


##################R package for Y叔################
rm(list = ls())  
options(stringsAsFactors = F)
load(file = '~/data/step1-output.Rdata')
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


















####################################
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(GEOquery)
library(ggplot2)
library(clusterProfiler)
library(GSEABase) 
# 这里载入前面下载到的GEO数据集，看前面推文即可
load(file = '~/data/step1-output.Rdata')
colnames(dat)
dat[1:4,1:4] 

m=dat[,2]+dat[,1]/2
logFC1=log2(dat[,2]/dat[,1]);fivenum(logFC1)
logFC2=dat[,2]-dat[,1];fivenum(logFC2) 
## 这里我们也挑选基因
table(logFC2>0.7)
logFC2=logFC2[logFC2>0.7] 

geneList=scale(logFC2)
names(geneList)=rownames(geneList)
geneList=sort(geneList,decreasing = T)
head(geneList)
# 这里需要注意的是 要 以 entrez ID来进行后续GSEA分析
# 如果是symbol，需要使用其它包，比如 GSEABase
library(org.Mm.eg.db)
s2g=select(org.Mm.eg.db,names(geneList),
           'ENTREZID','SYMBOL')

head(s2g)
s2g=s2g[!is.na(s2g$ENTREZID),]
s2g=s2g[s2g$SYMBOL %in%  names(geneList),]
geneList=geneList[s2g$SYMBOL]  
names(geneList)=s2g$ENTREZID  
head(geneList)
tail(geneList)

## 这里比较耗时：
  go_bp_gsea <- gseGO(geneList     = geneList, 
                      OrgDb='org.Mm.eg.db',
                      ont = 'BP',
                      nPerm        = 1000,
                      minGSSize    = 10,
                      pvalueCutoff = 0.9,
                      verbose      = FALSE)
  
  gseaplot(go_bp_gsea, geneSetID = rownames(go_bp_gsea[1,]))
  gseaplot(go_bp_gsea, geneSetID = "GO:0006959")
  gseaplot(go_bp_gsea, geneSetID = "GO:0030595")
  data=go_bp_gsea@result
  table(data$pvalue<0.01)
  