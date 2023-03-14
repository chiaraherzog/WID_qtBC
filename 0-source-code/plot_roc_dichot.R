plot_roc_dichot <- function(type, index, dichotom_var, title1 = "all", title2 = "below 30", title3 = "above 30", style = "default",
                            direction = "<"){
  
  require(pROC)
  
  #--------------------------------#  
  # make annotation labels
  auc <- round(as.numeric(roc(type,index, quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type,index, quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type,index, quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno1 <- paste('AUC (', title1, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  lvls <- levels(dichotom_var)
  
  auc <- round(as.numeric(roc(type[dichotom_var==lvls[1]],index[dichotom_var==lvls[1]], quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type[dichotom_var==lvls[1]],index[dichotom_var==lvls[1]], quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type[dichotom_var==lvls[1]],index[dichotom_var==lvls[1]], quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno2 <- paste('AUC (', title2, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  auc <- round(as.numeric(roc(type[dichotom_var==lvls[2]],index[dichotom_var==lvls[2]], quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type[dichotom_var==lvls[2]],index[dichotom_var==lvls[2]], quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type[dichotom_var==lvls[2]],index[dichotom_var==lvls[2]], quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno3 <- paste('AUC (', title3, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  if(style == "lancet"){
    anno1 <- gsub("[.]", "·", anno1)
    anno2 <- gsub("[.]", "·", anno2)
    anno2 <- gsub("[.]", "·", anno2)
  }
  
  #--------------------------------#  
  
  roc1 <- roc(type, index, direction = direction)
  roc2 <- roc(type[dichotom_var==lvls[1]], index[dichotom_var==lvls[1]], direction = direction)
  roc3 <- roc(type[dichotom_var==lvls[2]], index[dichotom_var==lvls[2]], direction = direction)
  title1 <- as.character(title1)
  title2 <- as.character(title2)
  title3 <- as.character(title3)
  
  ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              size = 0.7) +
    geom_path(aes(x=1-roc2$specificities,
                  y=(roc2$sensitivities)),
              size = 0.7,
              linetype = "dotted") +
    geom_path(aes(x=1-roc3$specificities,
                  y=(roc3$sensitivities)),
              size = 0.7,
              linetype = "dashed") +
    annotate("segment", x = 0, y = 0,
             xend = 1, yend = 1,
             colour = "gray60") +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_minimal() +
    theme(plot.title = element_text(size=10),
          panel.grid = element_blank()) +
    annotate("segment",
               x = 0.4, xend = 0.5,
               y = 0.55, yend = 0.55,
               linetype = "solid", linewidth = 0.9) +
    annotate("segment",
               x = 0.4, xend = 0.5,
               y = 0.35, yend = 0.35,
               linetype = "dotted", linewidth = 0.9) +
    annotate("segment",
               x = 0.4, xend = 0.5,
               y = 0.15, yend = 0.15,
               linetype = "dashed", linewidth = 0.9) +
    annotate(geom='text',
             x=0.55,
             y=0.55,
             label=anno1,
             hjust = 0,
             colour='black',
             size=3.1)  +
    annotate(geom='text',
             x=0.55,
             y=0.35,
             hjust = 0,
             label=anno2,
             colour='black',
             size=3.1) +
    annotate(geom='text',
             x=0.55,
             y=0.15,
             hjust = 0,
             label=anno3,
             colour='black',
             size=3.1) 
  
}
