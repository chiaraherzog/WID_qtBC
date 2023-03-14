# mQTL BC (WID-qtBC) index trained on entire discovery set
# Author: Chiara Herzog
# Date: 13 Aug 2021

index_mQTL <- function(beta){
  
  # Load coefficients 
  index <- readRDS('1-analysis-pipeline/2-output/index_coef.Rds')
  scale <- readRDS('1-analysis-pipeline/2-output/scale.Rds')
  
  # Compute index
  ind <- na.omit(match(names(index), rownames(beta)))
  b <- beta[ind,]
  
  ind <- na.omit(match(rownames(b), names(index)))
  w <- index[ind] 
  
  if(!identical(names(w), rownames(b))){
    stop('***names mismatch***')
  }
  
  index_mQTL <- colSums(b*w, na.rm = TRUE)
  index_mQTL <- -(index_mQTL - mean(scale)) / sd(scale)

  return(index_mQTL)
}
