# Index training in Discovery set -----------------------
# EGAD00010002079 and EGAD00010002081

# Get Loci
mqtl <- readxl::read_xlsx("1-analysis-pipeline/0-data/10549_2021_6185_MOESM1_ESM.xlsx", sheet = 6, skip = 3) |>  # get mQTL CpGs
  janitor::clean_names() |>
  dplyr::pull(cp_g_site)

# load data
load("~/Dropbox/data/3c/beta_merged.Rdata")
load("1-analysis-pipeline/0-data/data.Rdata")
tr <- data |> 
  dplyr::filter(grepl("training", dataset))

val <- data |> 
  dplyr::filter(dataset == "discovery set, internal validation")

# get intersecting 704 CpGs and training samples
mqtl <- intersect(mqtl, rownames(beta_merged))
beta_tr <- beta_merged[mqtl,tr$basename] # subset
beta_val <- beta_merged[mqtl,val$basename] 
identical(colnames(beta_tr), tr$basename) # check
identical(colnames(beta_val), val$basename) # check
rm(beta_merged);gc()


# Train using 2/3 split and get out of bag AUCs-----------------
source("0-source-code/el-classifier.R")
res_ridge <- el_classifier(beta_tr,
                           as.factor(tr$type),
                           beta_val,
                           as.factor(val$type),
                           alpha = 0,
                           lambda = "unrestricted")
res_ridge$oob_auc # 0.66 - calibration looks ok (slope - 0.994, intercept - 0.054)

res_lasso <- el_classifier(beta_tr,
                           as.factor(tr$type),
                           beta_val,
                           as.factor(val$type),
                           alpha = 1,
                           lambda = "unrestricted")
res_lasso$oob_auc # 0.637

# Finalise on entire dataset-----------------
beta <- cbind(beta_tr, beta_val)
pheno <- rbind(tr, val)
type <- pheno$type
res <- el_classifier(beta,
                     type,
                     beta,
                     type,
                     alpha = 1,
                     lambda = "unrestricted")
res$oob_auc # 0.68

# save coefs
index_coef <- coef(res$fit.cv, s = "lambda.min")
names <- rownames(index_coef)
index_coef <- as.numeric(index_coef)
names(index_coef) <- names
index_coef <- index_coef[index_coef!=0] # get non-zero coefs
saveRDS(index_coef, file="1-analysis-pipeline/2-output/index_coef.Rds")

# compute index in discovery for scaling of the index -------------
index <- index_coef
ind <- na.omit(match(names(index), rownames(beta)))
b <- beta[ind,]
ind <- na.omit(match(rownames(b), names(index)))
w <- index[ind] 
if(!identical(names(w), rownames(b))){
  stop('***names mismatch***')
}
scale <- colSums(b*w, na.rm = TRUE)

# quick plot: check directionality of index
pheno$scale <- scale
pheno |> 
  ggplot(aes(x = type,
             y = scale)) +
  geom_boxplot() # "inverted"
saveRDS(scale, file='1-analysis-pipeline/2-output/scale.Rds')

