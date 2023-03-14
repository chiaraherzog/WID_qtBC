ortable_percentiles <- function(score, type, percentiles, scoretype){
  require(epitools)
  results <- data.frame(matrix(nrow = length(percentiles),
                               ncol = 10))
  
  colnames(results) <- c("percentile",
                         "values",
                         "n_case",
                         "n_control",
                         "or",
                         "ci",
                         "ci_lo",
                         "ci_hi",
                         "p", 
                         "scoretype")
  
  results$scoretype <- scoretype
  
  for (i in 1:length(percentiles)){
    
    cutoff_lo <- ifelse(i == 1,
                        min(score),
                        quantile(score, percentiles[i-1]/100))
    cutoff_hi <- ifelse(i == max(length(percentiles)),
                        max(score),
                        quantile(score, percentiles[i]/100))
    
    alt <- score >= cutoff_lo & score < cutoff_hi
    ref <- score < cutoff_lo | score >= cutoff_hi
    
    tab <- as.table(matrix(c(sum(type[alt]!="Control"),
                             sum(type[ref]!="Control"),
                             sum(type[alt]=="Control"),
                             sum(type[ref]=="Control")),
                           nrow = 2))
    
    or <- oddsratio(tab) 
    
    results$percentile[i] <- ifelse(i == 1,
                                    paste0("0-", percentiles[i], "%"),
                                    paste0(percentiles[i-1], "-", percentiles[i], "%"))
    results$values[i] <- paste0(signif(cutoff_lo, 2), "-", signif(cutoff_hi,2))
    results$n_case[i] <- sum(type[alt]!="Control")
    results$n_control[i] <- sum(type[alt]=="Control")
    results$or[i] <- signif(or$measure[2,1],2)
    results$ci[i] <- paste0(signif(or$measure[2,2],2), "-",
                            signif(or$measure[2,3],2))
    results$ci_lo[i] <- signif(or$measure[2,2],2)
    results$ci_hi[i] <- signif(or$measure[2,3],2)
    results$p[i] <- signif(or$p.value[2,2], 3)
    
  }
  
  return(results)
}
