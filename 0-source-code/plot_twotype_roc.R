plot_twotype_roc <- function(type, index, title = NULL,
                             cols = c("gray20", "gray40"),
                             style = "default",
                             textcol = "white",
                             direction = "<",
                             size = 5){
  
  require(pROC)
  
  types <- as.character(levels(type))
  types <- types[types != "Control"]
  
  rocs <- vector(mode = "list", length = length(types))
  names(rocs) <- types
  ci <- vector(mode = "list", length = length(types))
  anno <- vector(mode = "list", length = length(types))
  
  for (t in types){
    i <- which(types == t)
    
    tmptype <- droplevels(type[type %in% c("Control", t)])
    tmpindex <- index[type %in% c("Control", t)]
    rocs[[i]] <- roc(tmptype, tmpindex, direction = direction)
    ci[[i]] <- ci.auc(rocs[[i]], conf.level = 0.95)
    
    anno[[i]] <- paste('AUC (', t, ') = ', round(rocs[[i]]$auc[[1]], 2),'\n(95% CI: ', round(ci[[i]][1], 2),'-',round(ci[[i]][3],2),')',sep='')
    
    
    if(style == "lancet"){
      anno[[i]] <- gsub("[.]", "·", anno[[i]])
    }
    
  }
  
  for (i in 1:length(types)){
    tmp <- data.frame(sensitivities = rocs[[i]]$sensitivities,
                      specificities = rocs[[i]]$specificities,
                      type = rep(types[i], length(rocs[[i]]$sensitivities)),
                      anno = rep(anno[[i]], length(rocs[[i]]$specificities)))
    
    if (i == 1){
      rocdf <- tmp
    } else {
      rocdf <- rbind(rocdf, tmp)
    }
  }
  
  anno.y <- vector("numeric", length = 4)
  anno.y <- case_when(length(types) == 1 ~ c(0.2, 0, 0, 0),
                      length(types) == 2 ~ c(0.35, 0.2, 0, 0),
                      length(types) == 3 ~ c(0.35, 0.25, 0.15, 0),
                      length(types) == 4 ~ c(0.45, 0.35, 0.25, 0.15))
  anno.y <- anno.y[1:length(types)]
  cols <- cols[1:length(types)]
  
  p <- rocdf |> 
    mutate(type = factor(type, levels = types)) |> 
    ggplot(aes(x = 1-specificities,
               y = sensitivities,
               colour = type)) +
    geom_abline(slope = 1, intercept = 0, size = 0.5, colour = "gray40") +
    geom_path(linewidth = 0.9) +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    annotate("label", 
             x = 0.70,
             y = anno.y,
             label = anno,
             fill = cols,
             colour = textcol,
             size = size) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_colour_manual(values = cols)
  
  
  return(p)
}
