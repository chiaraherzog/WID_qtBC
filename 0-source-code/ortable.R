# Odds ratio table for different strata
# Author: chiara herzog, chiara.herzog@uibk.ac.at
# Date: 22 Nov 2021
# Description: Types and strata should be ordered, with the reference type being level 1 and reference strata being 1 also.

ortable <- function(type, strata, typereference = 1, reference = 1) {
  
  require(epitools)
  
  # Create a DF for the outpout
  dat <- as.data.frame(matrix(nrow = length(levels(strata)), ncol = 6))
  colnames(dat) <- c("group", "Cases", "Controls", "OR", "95% CI", "p")
  
  types <- levels(type)
  # Run through strata and get fisher's exact p value. OR estimated by midp
  
  for (i in 1:length(levels(strata))){
    ref <- levels(strata)[reference]
    alt <- levels(strata)[i]
    
    tab <- as.table(matrix(c(sum(strata==alt & type==types[2]),
                             sum(strata==ref & type==types[2]),
                             sum(strata==alt & type== types[1]),
                             sum(strata==ref & type== types[1])),
                           nrow = 2))
    
    or <- oddsratio(tab) 
    
    dat$group[i] <-  alt
    dat$Cases[i] <- sum(type==types[2] & strata==alt)
    dat$Controls[i] <- sum(type==types[1] & strata==alt)
    dat$OR[i] <- ifelse(i == 1, "1 (Reference)",
                        signif(or$measure[2,1],2))
    dat$`95% CI`[i] <- ifelse(i == 1, "-",
                              paste0(signif(or$measure[2,2],2),"-", signif(or$measure[2,3],2)))
    dat$p[i] <- ifelse(i == 1, "-",
                       signif(or$p.value[2,2], 3))
  }
  
  return(dat)
}
