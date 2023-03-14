# Permutation testing for mQTL CpGs -----------------------
# entire dataset is loaded; CpGs are randomly sampled to check what AUC would be at random. EGAD00010002079 and EGAD00010002081

# Discovery dataset
load("~/Dropbox/data/3c/beta_merged.Rdata")
load("1-analysis-pipeline/0-data/data.Rdata")
disco <- data |> 
  dplyr::filter(grepl("discovery", dataset))

beta_tr <- beta_merged[,match(disco$basename, colnames(beta_merged))]
identical(disco$basename, colnames(beta_tr))
rm(beta_merged);gc()

# Validation dataset
load("~/Dropbox/data/3c-ext-validation/beta_merged.Rdata")
val <- data |>  # load validation dataset to report on validation in dataset
  dplyr::filter(grepl("^validation", dataset))
# Test
beta_val <- beta_merged[,match(val$basename, colnames(beta_merged))]
identical(val$basename, colnames(beta_val))
rm(beta_merged);gc()

# Find overlapping CpGs
intersect <- intersect(rownames(beta_tr), rownames(beta_val))
beta_tr <- beta_tr[match(intersect, rownames(beta_tr)),]
beta_val <- beta_val[match(intersect, rownames(beta_val)),]
identical(rownames(beta_tr), rownames(beta_val))

# get type vals
type_tr <- as.factor(disco$type)
type_val <- as.factor(val$type)
rm(data, intersect, disco, val); gc()

# Run 10k permutations -----------------------
df <- data.frame(matrix(nrow = 10000, ncol = 8))
colnames(df) <- c("auc", "confint_l", "confint_hi", "cpgs", "n_index", "index_cpgs", "mqtl_cpgs", "mqtl_cpgs_index")
mqtl <- readxl::read_xlsx("1-analysis-pipeline/0-data/10549_2021_6185_MOESM1_ESM.xlsx", sheet = 6, skip = 3) |>  # get mQTL CpGs
  janitor::clean_names() |>
  dplyr::pull(cp_g_site)
intersect <- intersect(mqtl, rownames(beta_tr))
mQTL <- mqtl[match(intersect, mqtl)];rm(mqtl, intersect)

source("0-source-code/el-classifier.R")

cg <- nrow(beta_tr)
n <- 10000

set.seed(314153467)
pB <- txtProgressBar(min=1,max=n, width =50L, style = 3)
for (i in 1:n){
  # sample
  ind <- sample(1:cg, 704, replace = FALSE)
  tmp_tr <- beta_tr[ind,]
  tmp_val <- beta_val[ind,]
  
  # save subsets
  df[i,]$cpgs <- paste(rownames(tmp_tr), collapse = ",")
  df[i,]$mqtl_cpgs <- sum(rownames(tmp_tr) %in% mQTL)
  
  # train
  res <- el_classifier(tmp_tr,
                       type_tr,
                       tmp_val,
                       type_val,
                       alpha = 1)
  
  #hi-lo tr set
  means <- mean(res$tr_predictor[type_tr == "Control",]) - mean(res$tr_predictor[type_tr!="Control",])
  dir_tr <- ifelse(means > 0, "control bigger", "case bigger")
  
  means_val <- mean(res$val_predictor[type_val=="Control",]) - mean(res$val_predictor[type_val!="Control",])
  dir_val <- ifelse(means_val > 0, "control bigger", "case bigger")
  
  # save aucs
  if(dir_tr == dir_val){
    df[i,]$confint_l <- suppressMessages(round(as.numeric(ci.auc(type_val, as.numeric(res$val_predictor)))[1],2))
    df[i,]$confint_hi <- suppressMessages(round(as.numeric(ci.auc(type_val, as.numeric(res$val_predictor)))[3],2))
    df[i,]$auc <- suppressMessages(round(as.numeric(ci.auc(type_val, as.numeric(res$val_predictor)))[2],2))
  } else {
    df[i,]$confint_l <- suppressMessages(1-round(as.numeric(ci.auc(type_val, as.numeric(res$val_predictor)))[1],2))
    df[i,]$confint_hi <- suppressMessages(1-round(as.numeric(ci.auc(type_val, as.numeric(res$val_predictor)))[3],2))
    df[i,]$auc <- suppressMessages(1-round(as.numeric(ci.auc(type_val, as.numeric(res$val_predictor)))[2],2))
  }
  # index cpgs
  index_coef <- coef(res$fit.cv, s = "lambda.min")
  names <- rownames(index_coef)
  index_coef <- as.numeric(index_coef)
  names(index_coef) <- names
  index_coef <- index_coef[index_coef!=0][-1]
  df[i,]$n_index <- length(index_coef)
  df[i,]$index_cpgs <- paste(names(index_coef), collapse = ",")
  df[i,]$mqtl_cpgs_index <- sum(names(index_coef) %in% mQTL)
  
  # cleanup and progress
  setTxtProgressBar(pB, i)
  gc()
}

save(df, file = "1-analysis-pipeline/3-output/permutation.Rdata")