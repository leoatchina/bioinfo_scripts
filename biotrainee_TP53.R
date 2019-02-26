# Step0 Before starting your project --------------------------------------

## Remove everything in the working environment, not including loaded libraries.
rm(list = objects( all = TRUE ))
options( stringsAsFactors = FALSE )
if (capabilities( "libcurl" ) == T ) {
  options( download.file.method = "libcurl" )
}

if (!is.null( dev.list() )) dev.off() 

clearhistory = function() {
  write( "", file = ".blank" )
  loadhistory( ".blank" )
  unlink( ".blank" )
}
clearhistory()
local({
  options( repos  = "https://mirrors.ustc.edu.cn/CRAN/" )
  options( BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/" )
})
bioPackages = 
  c( 
    "dplyr", "stringi", "purrr", ## 
    "R.utils", "data.table", ## unzip and read table
    "GEOquery", ## download
    "FactoMineR", "factoextra", "ggfortify", ## PCA
    "pheatmap", ## heatmap
    "ggplot2", ## Volcano plot
    "limma", "DESeq2", "edgeR", ## DEG
    "clusterProfiler", "org.Hs.eg.db", ## annotation
    "pathview" ## kegg
  )
lapply( bioPackages, 
  function(bioPackage) {
    if ( !require( bioPackage, character.only = T ) ) {
      CRANpackages = available.packages()
      
      ## install packages by CRAN
      if ( bioPackage %in% rownames( CRANpackages ) ) {
        install.packages( bioPackage )
        
      }else{
        ## install packages by bioconductor
        ## R version >= 3.5 ===> BiocManager
        if ( as.character( sessionInfo()$R.version$minor ) >= 3.5 ) {
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install(bioPackage, update = TRUE, ask = FALSE)
          
        }else{
          ## R version < 3.5 ===> BiocInstaller
          if (!requireNamespace("BiocInstaller", quietly = TRUE))
            source( "https://bioconductor.org/biocLite.R" )
          BiocInstaller::biocLite( bioPackage, ask = FALSE)
        }
      }
    }
  }
)

# Step1 download TCGA dateset ---------------------------------------------
# 这个技巧很好
if (!file.exists( './data/TCGA-BRCA.htseq_counts.Rdata' )) {
  gzfile = "./raw_data/TCGA-BRCA.htseq_counts.tsv.gz"
  download.file("https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.htseq_counts.tsv.gz", 
                destfile = gzfile)
  library(R.utils)
  gunzip(gzfile, remove = F)
  library(data.table)
  raw_data = fread( "./raw_data/TCGA-BRCA.htseq_counts.tsv",
                     sep = '	', header = T) # tidyverse
  raw_data = as.data.frame( raw_data )
  raw_data[1:5, 1:6] 
  rownames( raw_data ) = raw_data[, 1]
  raw_data = raw_data[, -1]
  raw_data[1:5, 1:6]
  raw_data = 2^raw_data - 1
  raw_data = ceiling( raw_data )
  raw_data[1:5, 1:6]
  pick_row = apply( raw_data, 1, function(x){
    sum(x == 0) < 10
  })
  raw_data = raw_data[pick_row, ]
  dim(raw_data  )
  save( raw_data, file = './data/TCGA-BRCA.htseq_counts.Rdata' )
}else{
  load('./data/TCGA-BRCA.htseq_counts.Rdata')
}


# Step2 Grouping by special clinical information --------------------------

if (!file.exists( './raw_data/TCGA-BRCA.GDC_phenotype.tsv.gz' )) {
  gzfile = "./raw_data/TCGA-BRCA.GDC_phenotype.tsv.gz"
  download.file("https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.GDC_phenotype.tsv.gz", 
                destfile = gzfile)
  phenoData = read.table( gzfile,
                           header = T,
                           sep = '	',
                           quote = '' )
  save( phenoData, file = './data/TCGA-BRCA.GDC_phenotype.Rdata' )
}else{
  load('./data/TCGA-BRCA.GDC_phenotype.Rdata')
}

pheno_num = c()
invisible(
  lapply(1:ncol(phenoData), 
         function(col_num){
           ## Assume that the classification project is between 2 and 4
           if (1 < dim(table(phenoData[,col_num])) & 
               dim(table(phenoData[,col_num])) < 5) {
             pheno_num <= append(pheno_num, col_num, after = length(pheno_num))
           }
         }
  )
)
names(phenoData[, pheno_num])


## Category 3: TP53
if (!file.exists( './raw_data/TCGA-BRCA.mutect2_snv.tsv.gz' )) {
  gzfile = "./raw_data/TCGA-BRCA.mutect2_snv.tsv.gz"
  download.file("https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.mutect2_snv.tsv.gz", 
                destfile = gzfile)
  mutype_file = read.table( gzfile,
                             header = T,
                             sep = '	',
                             quote = '' )
  save( mutype_file, file = './data/TCGA-BRCA.mutect2_snv.Rdata' )
}else{
  load('./data/TCGA-BRCA.mutect2_snv.Rdata')
}

## Pick columns that contains 'tp53'
TP53 = mutype_file[mutype_file$gene == 'tp53' | mutype_file$gene == 'TP53',]
TP53_sample = unique( sort( TP53$Sample_ID ) )
tumor_sample = colnames(raw_data)[substr( colnames(raw_data),14,15) < 10]
TP53_sample = intersect(tumor_sample, TP53_sample)
noTP53_sample = setdiff(tumor_sample, TP53_sample)
save(TP53_sample, noTP53_sample, file = './data/sample_by_TP53.Rdata')


# Step3 Filt sample ------------------------------------------------

load('./data/TCGA-BRCA.htseq_counts.Rdata')
tp53_sample = c(TP53_sample, noTP53_sample)
AssayData = raw_data[, tp53_sample]
dim(AssayData)
group_list = c(rep('TP53', length(TP53_sample)),
                rep('NO_TP53', length(noTP53_sample)))
save(AssayData, group_list, file = './data/tnbc_tumor_TP53_AssayData.Rdata')



#　heatmap and volcano map
draw_heatmap = function(nrDEG, type){
  library( "pheatmap" )
  nrDEG_Z = nrDEG[ order( nrDEG$logFC ), ]
  nrDEG_F = nrDEG[ order( -nrDEG$logFC ), ]
  choose_gene = c( rownames( nrDEG_Z )[1:50], rownames( nrDEG_F )[1:50] )
  choose_matrix = AssayData[ choose_gene, ]
  choose_matrix = t( scale( t( choose_matrix ) ) )
  
  choose_matrix[choose_matrix > 2] = 2
  choose_matrix[choose_matrix < -2] = -2
  
  annotation_col = data.frame( CellType = factor( group_list ) )
  rownames( annotation_col ) = colnames( AssayData )
  filename = paste('./fig/', type, '_heatmap_top100_logFC.png',
                    sep = "", collapse = NULL)
  pheatmap( fontsize = 6, choose_matrix, annotation_col = annotation_col, 
            show_rownames = T, show_colnames = F,
            annotation_legend = T, cluster_cols = T, 
            filename = filename)
}

draw_volcano = function(nrDEG, type){
  library( "ggplot2" )
  logFC_cutoff = with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
  
  nrDEG$change = as.factor( ifelse( 
    nrDEG$P.Value < 0.01 & abs(nrDEG$logFC) > logFC_cutoff,
    ifelse( nrDEG$logFC > logFC_cutoff, 'UP', 'DOWN' ), 'NOT' ) )
  nrDEGfile = paste('./data/', type, '_nrDEG_by_logFC.Rdata',
                     sep = "", collapse = NULL)
  save( nrDEG, file = nrDEGfile )
  
  this_tile = paste0( 
    'Cutoff for logFC is ', round( logFC_cutoff, 3 ),
    'The number of up gene is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),
    'The number of down gene is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ) )
    volcano = ggplot(data = nrDEG, 
                   aes( x = logFC, y = -log10(P.Value), color = change)) +
  geom_point( alpha = 0.4, size = 1.75 ) +
  theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
  xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
  ggtitle( this_tile ) + 
  theme( plot.title = element_text( size = 15, hjust = 0.5 )) +
  scale_colour_manual( values = c('blue','black','red') )
  print( volcano )
  filename = paste('./fig/', type, '_volcano_logFC.png',
                  sep = "", collapse = NULL)
  ggsave( volcano, filename = filename )
  dev.off()
}
# Step2 Then run edgeR ----------------------------------------------------

library(edgeR)
## A list-based S4 class for storing read counts and associated information 
## from digital gene expression or sequencing technologies.
DGElist <- DGEList( counts = AssayData, group = factor(group_list) )
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
## Calculate Normalization Factors to Align Columns of a Count Matrix
DGElist <- calcNormFactors( DGElist )
DGElist$samples

design <- model.matrix( ~0 + factor(group_list) )
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(group_list))

## Estimate Common Dispersion for Negative Binomial GLMs
DGElist <- estimateGLMCommonDisp(DGElist, design)
## Estimate Trended Dispersion for Negative Binomial GLMs
DGElist <- estimateGLMTrendedDisp(DGElist, design)
## Empirical Bayes Tagwise Dispersions for Negative Binomial GLMs
DGElist <- estimateGLMTagwiseDisp(DGElist, design)

## glmFit fits genewise negative binomial glms, all with the same design matrix 
## but possibly different dispersions, offsets and weights
fit <- glmFit(DGElist, design)
## https://www.biostars.org/p/110861/
## glmLRT conducts likelihood ratio tests for one or more coefficients in the 
## linear model.
results <- glmLRT(fit, contrast = c(-1, 1)) 
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
head(nrDEG_edgeR)
colnames(nrDEG_edgeR)[4] <- c("P.Value") 
draw_heatmap(nrDEG = nrDEG_edgeR, type = 'edgeR')
draw_volcano(nrDEG = nrDEG_edgeR, type = 'edgeR')


# Step3 Lastly run voom from limma ----------------------------------------

library(limma)
## A list-based S4 class for storing read counts and associated information 
## from digital gene expression or sequencing technologies.
DGElist <- DGEList( counts = AssayData, group = factor(group_list) )
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
## Calculate Normalization Factors to Align Columns of a Count Matrix
DGElist <- calcNormFactors( DGElist )
DGElist$samples

design <- model.matrix( ~0 + factor(group_list) )
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(group_list))

## Transform RNA-Seq Data Ready for Linear Modelling
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
## Fit linear model for each gene given a series of arrays
fit <- lmFit(v, design)

## Construct the contrast matrix corresponding to specified contrasts of a set 
## of parameters.
cont.matrix <- makeContrasts(contrasts = c('TP53-NO_TP53'), levels = design)
## Given a linear model fit to microarray data, compute estimated coefficients 
## and standard errors for a given set of contrasts.
fit2 <- contrasts.fit(fit, cont.matrix)
## Empirical Bayes Statistics for Differential Expression
fit2 <- eBayes(fit2)

nrDEG_limma_voom = topTable(fit2, coef = 'TP53-NO_TP53', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
draw_heatmap(nrDEG = nrDEG_limma_voom, type = 'limma_voom')
draw_volcano(nrDEG = nrDEG_limma_voom, type = 'limma_voom')


# Step4 Compare three methods ---------------------------------------------

cor(na.omit(lf))
mi <- unique(c(rownames(nrDEG_edgeR),
               rownames(nrDEG_limma_voom)))
lf <- data.frame(edgeR = nrDEG_edgeR[mi, 1],
                 limma_voom = nrDEG_limma_voom[mi, 1])
cor(na.omit(lf))
##                edgeR     limma_voom
##    edgeR      1.0000000  0.9118619
##    limma_voom 0.9118619  1.0000000

## function of KEGG pathway
kegg_plot = function(type) {
  kk.up = enrichKEGG(   gene          =  gene_up    ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.8        ,
                         qvalueCutoff  =  0.8        )
  kk.down = enrichKEGG( gene          =  gene_down  ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.8        ,
                         qvalueCutoff  =  0.8        )
  library( "ggplot2" )
  kegg_down_dt = as.data.frame( kk.down )
  kegg_up_dt = as.data.frame( kk.up )
  down_kegg = kegg_down_dt[ kegg_down_dt$pvalue < 0.05, ]
  down_kegg$group = -1
  up_kegg = kegg_up_dt[ kegg_up_dt$pvalue < 0.05, ]
  up_kegg$group = 1
  dat = rbind( up_kegg, down_kegg )
  dat$pvalue = -log10( dat$pvalue )
  dat$pvalue = dat$pvalue * dat$group
  dat = dat[ order( dat$pvalue, decreasing = F ), ]
  g_kegg = ggplot( dat, 
                    aes(x = reorder( Description, order( pvalue, decreasing = F )), 
                        y = pvalue, fill = group)) + 
    geom_bar( stat = "identity" ) + 
    scale_fill_gradient( low = "blue", high = "red", guide = FALSE ) + 
    scale_x_discrete( name = "Pathway names" ) +
    scale_y_continuous( name = "log10P-value" ) +
    coord_flip() + theme_bw() + 
    theme( plot.title = element_text( hjust = 0.5 ), 
           axis.text.x = element_text(size = 10),
           axis.text.y = element_text(size = 7)) +
    ggtitle( "Pathway Enrichment" ) 
  print( g_kegg )
  filename = paste('./fig/kegg_up_down_', type, '.png', sep = "", collapse = NULL)
  ggsave( g_kegg, filename = filename )
}

## function of GO pathway
go_plot = function(type) {
  go_enrich_results = lapply( g_list, function(gene) {
    lapply( c( 'BP', 'MF', 'CC' ) , function(ont) {
      cat( paste( 'Now process', ont ) )
      ego = enrichGO( gene          =  gene,
                       universe      =  gene_all,
                       OrgDb         =  org.Hs.eg.db,
                       ont           =  ont ,
                       pAdjustMethod =  "BH",
                       pvalueCutoff  =  0.99,
                       qvalueCutoff  =  0.99,
                       readable      =  TRUE)
      print( head( ego ) )
      return( ego )
    })
  })
  gofilename = paste('./data/go_enrich_result', type, '.Rdata', 
                      sep = "", collapse = NULL)
  save( go_enrich_results, file = gofilename )
  
  n1 = c( 'gene_up', 'gene_down', 'gene_diff' )
  n2 = c( 'BP', 'MF', 'CC' ) 
  for (i in 1:3) {
    for (j in 1:3) {
      fn = paste0( './fig/dotplot_', n1[i], '_', n2[j], '_', type, '.png' )
      cat( paste0( fn, '' ) )
      png( fn, res = 150, width = 1080 )
      print( dotplot( go_enrich_results[[i]][[j]] ) )
      dev.off()
    }
  }
}

# Step1 annotation --------------------------------------------------------

library( "clusterProfiler" )
library( "org.Hs.eg.db" )
keytypes(org.Hs.eg.db)
library("stringr")

load( "./data/edgeR_nrDEG_by_logFC.Rdata" )
load( "./data/limma_voom_nrDEG_by_logFC.Rdata" )

## tans1: ENSEMBL2ENTREZID
table( nrDEG$change )
rownames( nrDEG ) = str_sub(rownames( nrDEG ), start = 1, end = 15)
nrDEG$ENSEMBL = rownames( nrDEG )
df = bitr( rownames( nrDEG ), fromType = "ENSEMBL", toType = c( "ENTREZID" ), 
            OrgDb = org.Hs.eg.db )
head( df )
nrDEG = merge( nrDEG, df, by = 'ENSEMBL' )
head( nrDEG )

gene_up = nrDEG[ nrDEG$change == 'UP', 'ENTREZID' ] 
gene_down = nrDEG[ nrDEG$change == 'DOWN', 'ENTREZID' ]
gene_diff = c( gene_up, gene_down )
gene_all = as.character(nrDEG[ ,'ENTREZID'] )
g_list = list( gene_up = gene_up, gene_down = gene_down, gene_diff = gene_diff)


# Step2 pathway analysis ------------------------------------------

kegg_plot("edgeR")
go_plot("edgeR")

kegg_plot("limma_voom")
go_plot("limma_voom")

library(pathview)
geneList = nrDEG$logFC
names( geneList ) = nrDEG$ENTREZID
geneList = sort( geneList, decreasing = T )
pathview( gene.data = geneList, 
          pathway.id = dat$ID, 
          species = "hsa", 
          limit = list(gene = 5, cpd = 1))


