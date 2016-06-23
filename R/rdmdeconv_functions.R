# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'







.setupBiocLite <- function(){
  if( !require(BiocInstaller) ){
    setRepositories()
    source("https://bioconductor.org/biocLite.R")
    biocLite("BiocInstaller")
    library(BiocInstaller)
  }
}


.setupDESeq <-function(){
  .setupBiocLite()
  if( !require(dplyr) ){
    biocLite('dplyr')
    library(dplyr)
  }
  if( !require(DESeq) ){
    biocLite('DESeq')
    library(DESeq)
  }
}


.setupCellMix<- function(){
  .setupBiocLite()
  if( !require(CellMix) ){
    biocLite('CellMix', siteRepos = 'http://web.cbio.uct.ac.za/~renaud/CRAN')
    library(CellMix)
  }
}

make_signatures <- function(proportions=proportions57,all_samples=epigenomes57.N){

  signatures <- sapply(rownames(proportions),function(rowname) {
    ss = proportions[rowname,]
    whichcols = colSums(ss)>0 #Important columns
    ss = ss[,whichcols]
    head(ss)
    sum = sum(ss)
    return (rowSums(as.matrix(all_samples[,whichcols]) %*% diag(as.vector(ss)) / sum ))
  })

  all_samples = all_samples[,which(colSums(all_samples,na.rm = TRUE) != 0 )]

  rownames(signatures) = rownames(all_samples)
  (filter_empty_data(signatures))
}

normalize_data_by_columns <- function(inp){
  out <- sapply(colnames(inp),function(cname){(inp[,cname]-min(inp[,cname]))/(max(inp[,cname])-min(inp[,cname]))})
  colnames(out) <- colnames(inp)
  rownames(out) <- rownames(inp)
  (out)
}

extract_markers <- function(signatures=make_signatures()){
  .setupDESeq()
  cds <- newCountDataSet(floor(signatures), conditions = colnames(signatures))
  cds <- estimateSizeFactors( cds )
  cds <- estimateDispersions(cds,method='blind',sharingMode='fit-only')
  fit1 <- fitNbinomGLMs( cds, count ~ condition )
  fit0 <- fitNbinomGLMs( cds, count ~ 1 )
  pvals <- nbinomGLMTest( fit1, fit0 )
  padj <- p.adjust( pvals, "BH" )
  
  markers <- signatures / rowSums(signatures)
  markers <- normalize_data_by_columns(filter_empty_data(markers))
  max_row <- apply(markers,1,max)
  which_max_row <-apply(markers,1,which.max)
  max_signatures_row <- apply(signatures,1,max)
  markers_df <- data.frame(rownames(signatures),
                           padj,
                           pvals,
                           max_row,
                           colnames(signatures)[which_max_row],
                           max_signatures_row)
  
  markers_df = markers_df[!is.na(markers_df$pvals),]
  
  colnames(markers_df) <- c("id","pvals","padj",'max','sig','sigval')
  (markers_df)
}


map_genes <- function(input){
  (sapply(strsplit(input,"\\|"),function(gn){as.character(gene_info[which(gene_info$V7==gn[1]),]$V1)}))
}

filter_empty_data <-function(input){
  (input[which(rowSums(input,na.rm = TRUE) != 0 ), which(colSums(input,na.rm = TRUE) != 0 )])
}

filter_data <- function(input,signatures){
  input <- filter_empty_data(input)
  (input[intersect(rownames(input),rownames(signatures)),])
}


deconv <- function(mix, signatures=make_signatures(), markers= rownames(signatures), method='lsfit'){
  if (method=='deconRNASeq'){
    (deconv_deconrnaseq(mix,signatures,markers))
  } else{
    .setupCellMix()
    mix <- filter_data(mix,signatures)
    signatures <- filter_data(signatures,mix)
    markers <- intersect(markers,rownames(signatures))
    gse <- ExpressionMix(data.matrix(mix))
    res <- ged(gse, signatures,subset=markers,method)
    (coef(res))    
  }
}

deconv_deconrnaseq <- function(mix, signatures=make_signatures(), markers= rownames(signatures)) {
  if( !require(DeconRNASeq) ){
    .setupBiocLite()
    biocLite("DeconRNASeq")
    library(DeconRNASeq)
  }
  mmix <- mix[markers,]
  msignatures <- signatures[markers,]
  res <- DeconRNASeq(mmix,
                     as.data.frame(msignatures),
                     proportions = NULL,
                     checksig=TRUE,
                     known.prop = FALSE,
                     use.scale = TRUE,
                     fig = TRUE)
  resall <- res$out.all
  rownames(resall) <- colnames(mix)
  (t(resall))
}

