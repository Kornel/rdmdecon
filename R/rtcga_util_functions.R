


plotRes <- function(res,cutoff=0.01){
  res_trimmed = res[which(rowMeans(res,na.rm = TRUE) >= cutoff ),]
  boxplot(t(res_trimmed),xlab="Tissue", ylab="Proportions",main="Estimated proportions",las=2,cex.axis=0.50)
}
map_genes_brca <- function(dataset,signatures){
  mapped_genes <- map_genes(colnames(dataset))
  without_duplicates = setdiff(1:length(mapped_genes),which(duplicated(mapped_genes))) # remove duplicates from mapped games
  mix = data.frame(t(dataset[,-1])[without_duplicates,]) #remove duplicates from mix
  colnames(mix) <- dataset[,1] #set column names based on first row from 
  rownames(mix) <- mapped_genes[without_duplicates] # set proper rownames without duplicates
  #Filter data (ex. available genes) performed by deconv anyway
  (filter_data(mix,signatures))
}

test_gene_mapping <- function(dataset,mix){
  test <- data.frame(colnames(dataset)[without_duplicates],rownames(mix))
  colnames(test) <- c('BRCA','Roadmap')
  
  #Data mapped
  head(test[order(test[2]),],15)
  
  #Gene info
  selgi <- gene_info[gene_info$V1 %in% rownames(mix),]
  head(selgi[order(selgi$V1),],15)
  
}
