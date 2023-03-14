plot_two_rocs_ER <- function(type, index, er, group2 = "Group2",
                          direction = "<"){
  
  require(pROC)
  
  #--------------------------------#  
  # make annotation labels
  ind1 <- (er=="Neg" & !is.na(er)) | type == "Control"
  ind2 <- (er=="Pos" & !is.na(er)) | type == "Control"
  
  auc <- round(as.numeric(roc(type, index, quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type, index, quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type, index, quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno1 <- paste('AUC = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  auc <- round(as.numeric(roc(type[ind1], index[ind1], quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type[ind1], index[ind1], quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type[ind1], index[ind1], quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno2 <- paste('AUC (ER-) = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  auc <- round(as.numeric(roc(type[ind2], index[ind2], quiet=T, ci=T, direction = direction)$auc),digits=2)
  cil <- round(as.numeric(roc(type[ind2], index[ind2], quiet=T, ci=T, direction = direction)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type[ind2], index[ind2], quiet=T, ci=T, direction = direction)$ci[3]),digits=2)
  anno3 <- paste('AUC (ER+) = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  #--------------------------------#  
  
  roc1 <- roc(type, index, direction = direction)
  roc2 <- roc(type[ind1], index[ind1], direction = direction)
  roc3 <- roc(type[ind2], index[ind2], direction = direction)
  
  ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              linetype = "solid",
              size = 0.7) +
    geom_path(aes(x=1-roc2$specificities,
                  y=(roc2$sensitivities)),
              linetype = "dotted",
              size = 0.7) +
    geom_path(aes(x=1-roc3$specificities,
                  y=(roc3$sensitivities)),
              linetype = "dashed",
              size = 0.7) +
    annotate("segment", x = 0, y = 0,
             xend = 1, yend = 1,
             colour = "gray60") +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme(plot.title = element_text(size=10),
          panel.background = element_blank()) +
    annotate(geom = "segment",
             x = 0.45, xend = 0.55,
             y = 0.4, yend = 0.4,
             linewidth = 0.7,
             linetype = "solid") +
    annotate(geom = "segment",
             x = 0.45, xend = 0.55,
             y = 0.25, yend = 0.25,
             linewidth = 0.7,
             linetype = "dotted") +
    annotate(geom = "segment",
             x = 0.45, xend = 0.55,
             y = 0.1, yend = 0.1,
             linewidth = 0.7,
             linetype = "dashed") +
    annotate(geom='text',
             x=0.6,
             y=0.4,
             label=anno1,
             hjust = 0,
             size=2.8)  +
    annotate(geom='text',
             x=0.6,
             y=0.25,
             label=anno2,
             hjust = 0,
             size=2.8) +
    annotate(geom='text',
             x=0.6,
             y=0.1,
             hjust = 0,
             label=anno3,
             size=2.8)
  
}
