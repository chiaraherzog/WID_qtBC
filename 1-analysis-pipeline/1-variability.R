# Step 1: Compute variability of CpGs in Breast and matched blood, cervical, and buccal samples from same individuals.
# Raw data are deposited in EGA (access controlled), under accession numbers EGAS00001005055, EGAS00001005070, and EGAS00001005626.
# Raw data were preprocessed and normalised using the eutopsQC quality control pipeline.

library(dplyr)
library(here)

# 1.1 Breast data variability ------------------------
# Datasets EGAD00010002074, EGAD00010002075, EGAD00010002076

load("~/Dropbox/data/tissue-at-risk/beta_merged.Rdata") # load methylation beta matrix (after eutopsQC)
load("~/Documents/Work/data.nosync/GEO/GSE133985/beta.Rdata")
intersect <- intersect(rownames(beta_merged), rownames(beta))

load("1-analysis-pipeline/0-data/data.Rdata") # load pheno file

dat <- data |> # grab samples of interest
  dplyr::filter(dataset == "breast variability set") |> 
  droplevels()

table(dat$type) # normal, normal adjacent and BRCA1/2 mut carrier samples
beta <- beta_merged
beta <- beta[,match(dat$basename, colnames(beta))] # subset beta if necessary
identical(dat$basename, colnames(beta)) # check identical

sd_breast <- matrixStats::rowSds(as.matrix(beta)) # get cpg-level standard deviation
names(sd_breast) <- rownames(beta) # append names
thresholds <- sapply(c(0.99, 0.98, 0.95, 0.9, 0.85, 0.8), function(i){quantile(sd_breast, i)}) # select top pgs
top_cpgs <- list() # generate list variable
for(i in 1:length(thresholds)){
  top_cpgs[[i]] <-  names(sd_breast)[sd_breast>thresholds[i]]
}

sd_breast <- data.frame(sd_breast = sd_breast) |> # create a flag for top cpgs
  tibble::rownames_to_column("cg") |> 
  dplyr::filter(cg %in% top_cpgs[[6]]) |> 
  dplyr::mutate(top1 = case_when(cg %in% top_cpgs[[1]] ~ "y", TRUE ~ ""),
                top2 = case_when(cg %in% top_cpgs[[2]] ~ "y", TRUE ~ ""),
                top5 = case_when(cg %in% top_cpgs[[3]] ~ "y", TRUE ~ ""),
                top10 = case_when(cg %in% top_cpgs[[4]] ~ "y", TRUE ~ ""),
                top15 = case_when(cg %in% top_cpgs[[5]] ~ "y", TRUE ~ ""),
                top20 = case_when(cg %in% top_cpgs[[6]] ~ "y", TRUE ~ ""))

save(sd_breast, file = "1-analysis-pipeline/1-output/sd_breast.Rdata") # save outpout
rm(beta, beta_merged, dat, thresholds, top_cpgs, i);gc() # clean

# 1.2 Matched surrogate sample variability -----------------------
# Datasets EGAD00010002231, EGAD00010002232, EGAD00010002233
load("~/Dropbox/data/brca-ds1/beta_merged.Rdata") # load methylation beta matrix (after eutopsQC)
beta <- beta_merged

## 1.2.1 Cervical variability ----------
cervical <- data |> # select cervical samples
  dplyr::filter(dataset == "matched surrogate samples" & type == "cervical") |> 
  droplevels()
beta_cervical <- beta[,cervical$basename] # subset beta
identical(cervical$basename, colnames(beta_cervical)) # sanity check

sd_cervical <- data.frame(matrix(nrow = nrow(beta_cervical),
                                 ncol = 6)) # generate df
colnames(sd_cervical) <- c("cg", "sd_all", "sd_0_25", "sd_25_50", "sd_50_75", "sd_75_100") # name columns
sd_cervical$cg <- rownames(beta_cervical) # name rows
sd_cervical$sd_all <- matrixStats::rowSds(as.matrix(beta_cervical)) # variability for ic subgroups
sd_cervical$sd_0_25 <- matrixStats::rowSds(as.matrix(beta_cervical[,cervical$ic <= 0.25]))
sd_cervical$sd_25_50 <- matrixStats::rowSds(as.matrix(beta_cervical[,(cervical$ic > 0.25 & cervical$ic <= 0.50)]))
sd_cervical$sd_50_75 <- matrixStats::rowSds(as.matrix(beta_cervical[,(cervical$ic > 0.5 & cervical$ic <= 0.75)]))
sd_cervical$sd_75_100 <- matrixStats::rowSds(as.matrix(beta_cervical[,(cervical$ic >= 0.75)]))

## 1.2.2 buccal variability ----------
buccal <- data |> # select buccal samples
  dplyr::filter(dataset == "matched surrogate samples" & type == "buccal") |> 
  droplevels()
beta_buccal <- beta[,buccal$basename] # subset beta
identical(buccal$basename, colnames(beta_buccal)) # sanity check

sd_buccal <- data.frame(matrix(nrow = nrow(beta_buccal),
                                 ncol = 6)) # generate df
colnames(sd_buccal) <- c("cg", "sd_all", "sd_0_25", "sd_25_50", "sd_50_75", "sd_75_100") # name columns
sd_buccal$cg <- rownames(beta_buccal) # name rows
sd_buccal$sd_all <- matrixStats::rowSds(as.matrix(beta_buccal)) # variability for ic subgroups
sd_buccal$sd_0_25 <- matrixStats::rowSds(as.matrix(beta_buccal[,buccal$ic < 0.25]))
sd_buccal$sd_25_50 <- matrixStats::rowSds(as.matrix(beta_buccal[,(buccal$ic >= 0.25 & buccal$ic <= 0.50)]))
sd_buccal$sd_50_75 <- matrixStats::rowSds(as.matrix(beta_buccal[,(buccal$ic > 0.5 & buccal$ic < 0.75)]))
sd_buccal$sd_75_100 <- matrixStats::rowSds(as.matrix(beta_buccal[,(buccal$ic >= 0.75)]))

## 1.2.3 blood variability (myeloid v lymphoid) ----------
blood <- data |> # select blood samples
  dplyr::filter(dataset == "matched surrogate samples" & type == "blood") |> 
  droplevels()
beta_blood <- beta[,blood$basename] # subset beta
identical(blood$basename, colnames(beta_blood)) # sanity check

sd_blood <- data.frame(matrix(nrow = nrow(beta_blood),
                                 ncol = 4)) # generate df
colnames(sd_blood) <- c("cg", "sd_all", "sd_0_50", "sd_50_100") # name columns
sd_blood$cg <- rownames(beta_blood) # name rows
sd_blood$sd_all <- matrixStats::rowSds(as.matrix(beta_blood)) # variability for ic subgroups
sd_blood$sd_0_50 <- matrixStats::rowSds(as.matrix(beta_blood[,(blood$hepidish_Eosino+blood$hepidish_Neutro+blood$hepidish_Mono) <= 0.5]))
sd_blood$sd_50_100 <- matrixStats::rowSds(as.matrix(beta_blood[,(blood$hepidish_Eosino+blood$hepidish_Neutro+blood$hepidish_Mono) > 0.5]))

sd_tissues <- list(cervical = sd_cervical,
                   buccal = sd_buccal,
                   blood = sd_blood)

# 1.2.4 Compute median values across breast variability for plotting -----------------
# 1.2.4.1 Cervical ------
groups <- 6
groupnames <- colnames(sd_breast)[3:8]

icgroups <- 5
icgroupnames <- colnames(sd_cervical)[2:6]
df <- as.data.frame(matrix(nrow = groups*icgroups,
                           ncol = 5))
colnames(df) <- c("set", "icgroup", "median", "iqr_lower", "iqr_higher")

df$set <- rep(groupnames, each = icgroups)
df$icgroup <- rep(icgroupnames, groups)

for (i in groupnames){
  tmp <- sd_breast |> 
    dplyr::filter(get(i) == "y") |> 
    pull(cg)
  
  for(j in icgroupnames){
    tmp1 <- sd_cervical |> 
      dplyr::select(cg, any_of(j)) |> 
      dplyr::filter(cg %in% tmp) |> 
      pull(j)
    
    # median
    df[df$set==i & df$icgroup==j,]$median <- median(tmp1)
    df[df$set==i & df$icgroup==j,]$iqr_lower <- unname(quantile(tmp1, 0.25))
    df[df$set==i & df$icgroup==j,]$iqr_higher <- unname(quantile(tmp1, 0.75))
  }
}

cervical <- df
cervical$tissue <- "cervical"

# 1.2.4.2. Buccal ------
groups <- 6
groupnames <- colnames(sd_breast)[3:8]

icgroups <- 5
icgroupnames <- colnames(sd_buccal)[2:6]
df <- as.data.frame(matrix(nrow = groups*icgroups,
                           ncol = 5))
colnames(df) <- c("set", "icgroup", "median", "iqr_lower", "iqr_higher")

df$set <- rep(groupnames, each = icgroups)
df$icgroup <- rep(icgroupnames, groups)

for (i in groupnames){
  tmp <- sd_breast |> 
    dplyr::filter(get(i) == "y") |> 
    pull(cg)
  
  for(j in icgroupnames){
    tmp1 <- sd_buccal |> 
      dplyr::select(cg, any_of(j)) |> 
      dplyr::filter(cg %in% tmp) |> 
      pull(j)
    
    # median
    df[df$set==i & df$icgroup==j,]$median <- median(tmp1)
    df[df$set==i & df$icgroup==j,]$iqr_lower <- unname(quantile(tmp1, 0.25))
    df[df$set==i & df$icgroup==j,]$iqr_higher <- unname(quantile(tmp1, 0.75))
  }
}

buccal <- df
buccal$tissue <- "buccal"

# 1.2.4.3. Blood ------
groups <- 6
groupnames <- colnames(sd_breast)[3:8]

icgroups <- 3
icgroupnames <- colnames(sd_blood)[2:4]
df <- as.data.frame(matrix(nrow = groups*icgroups,
                           ncol = 5))
colnames(df) <- c("set", "icgroup", "median", "iqr_lower", "iqr_higher")

df$set <- rep(groupnames, each = icgroups)
df$icgroup <- rep(icgroupnames, groups)

for (i in groupnames){
  tmp <- sd_breast |> 
    dplyr::filter(get(i) == "y") |> 
    pull(cg)
  
  for(j in icgroupnames){
    tmp1 <- sd_blood |> 
      dplyr::select(cg, any_of(j)) |> 
      dplyr::filter(cg %in% tmp) |> 
      pull(j)
    
    # median
    df[df$set==i & df$icgroup==j,]$median <- median(tmp1)
    df[df$set==i & df$icgroup==j,]$iqr_lower <- unname(quantile(tmp1, 0.25))
    df[df$set==i & df$icgroup==j,]$iqr_higher <- unname(quantile(tmp1, 0.75))
  }
}

blood <- df
blood$tissue <- "blood"

# 1.2.4.4. Merge
variability <- rbind(cervical, buccal, blood)

# Save object for plot in Fig 1
save(variability, file = "1-analysis-pipeline/1-output/variability_breast_matched.Rdata")

# save sd variabilities for plotting (Fig 1)
save(sd_tissues, file = "1-analysis-pipeline/1-output/sd_tissues.Rdata")
