############################################################
# Fig. 7 ###################################################
# Trait variance explained by PGS in weighted subsamples ###
############################################################

library(data.table)
library(ggplot2)
library(reshape2)
library(ggnewscale)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)

transformEAtoEduYears <- function(EA_vector, to="years"){
  
  ea_names <- EA_vector
  ea_names[!(ea_names %in% 1:8)] <- NA
  
  if(to=="years"){
    
    ea_names[ea_names == 8] <- 22
    ea_names[ea_names == 7] <- 20
    ea_names[ea_names == 6] <- 18
    ea_names[ea_names == 5] <- 15
    ea_names[ea_names == 4] <- 13
    ea_names[ea_names %in% c(2, 3)] <- 10
    ea_names[ea_names == 1] <- 7
    ea_names[ea_names == 0] <- 1
    
  }else if(to=="binary"){
    
    ea_names[ea_names < 6] <- 0
    ea_names[ea_names >= 6] <- 1
    
  }else if(to=="normal"){
    
    ea_names <- as.numeric(ea_names)
    ea_names <- qnorm((rank(ea_names+runif(length(ea_names), -0.1, 0.1),na.last="keep")-0.5)
                      /sum(!is.na(ea_names)))
    
  }
  
  ea_names <- as.numeric(ea_names)
  return(ea_names)
  
}

transformSREAtoEduYears <- function(EA_vector, to="years"){
  
  ea <- strsplit(as.character(EA_vector), "(;| )+")
  ea_names <- unlist(lapply(ea, `[`, 1))
  ea_names[!(ea_names %in% 1:8)] <- NA
  
  if(to=="years"){
    
    ea_names[ea_names == 8] <- 22
    ea_names[ea_names %in% c(6, 7)] <- 20
    ea_names[ea_names == 5] <- 15
    ea_names[ea_names == 4] <- 13
    ea_names[ea_names == 3] <- 10
    ea_names[ea_names == 2] <- 7
    ea_names[ea_names == 1] <- 1
    
  }else if(to=="binary"){
    
    ea_names[ea_names < 6] <- 0
    ea_names[ea_names >= 6] <- 1
    
  }else if(to=="normal"){
    
    ea_names <- as.numeric(ea_names)
    ea_names <- qnorm((rank(ea_names+runif(length(ea_names), -0.1, 0.1),na.last="keep")-0.5)
                      /sum(!is.na(ea_names)))
    
  }
  
  ea_names <- as.numeric(ea_names)
  return(ea_names)
  
}

getValue <- function(trait_vector, measurement="first"){
  
  trait <- strsplit(as.character(trait_vector), "(;| )+")
  trait[trait == ""] <- NA
  if(measurement == "first"){
    values <- unlist(lapply(trait, `[`, 1))
  }else if(measurement == "last"){
    values <- unlist(lapply(trait, function(x) tail(x, 1)))
  }
  
  values <- as.numeric(values)
  return(values)
  
}


getOccStat <- function(x) {
  # Extract all numbers from the string
  occupation <- str_extract_all(x, "\\d+")[[1]]
  # Check if we have any numbers; if not, return NA
  if (length(occupation) == 0) return(NA)
  # Extract the first digit of each number and convert to numeric
  categories <- as.numeric(substr(occupation, 1, 1))
  # Check if we have only 0, which means work in armed forces; if yes, return NA
  if (all(categories == 0)) return(NA)
  categories <- categories[which(categories != 0)]
  # Calculate and return the mean of the first digits
  OS <- mean(categories, na.rm = TRUE)
  # Inverse the scale
  OS <- 10 - OS
  
  return(OS)
  
}


getIndepInd <- function(cutoff = 15){
  
  if(cutoff == 15){
    
    ind_01 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/ps_list_unrel.tsv", header=F)
    ind_02 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/s_list_unrel.tsv", header=F)
    ind_11 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p1ps_list_unrel.tsv", header=F)
    ind_12 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p1s_list_unrel.tsv", header=F)
    ind_21 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p2ps_list_unrel.tsv", header=F)
    ind_22 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p2s_list_unrel.tsv", header=F)
    
  }else if(cutoff == 10){
    
    ind_01 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/ps_list_unrel_10.tsv", header=F)
    ind_02 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/s_list_unrel_10.tsv", header=F)
    ind_11 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p1ps_list_unrel_10.tsv", header=F)
    ind_12 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p1s_list_unrel_10.tsv", header=F)
    ind_21 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p2ps_list_unrel_10.tsv", header=F)
    ind_22 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/gcta/data/p2s_list_unrel_10.tsv", header=F)
    
  } else {
    
    
    stop("ERROR: not defined cutoff")
    
  }
  
  ind_01 <- ind_01$V1
  ind_02 <- ind_02$V1
  ind_11 <- ind_11$V1
  ind_12 <- ind_12$V1
  ind_21 <- ind_21$V1
  ind_22 <- ind_22$V1
  
  return(list(ind_01, ind_02, ind_11, ind_12, ind_21, ind_22))
  
}


matchCohorts <- function(ref_ind, biased_ind, n, replacement, max_ratio=1, seed = 1){
  
  # Take samples from biased sample matching with the reference one
  # Let's say "we are accounting for bias moving towards reference distribution"
  set.seed(seed = seed)
  
  # Match EduYears
  s1 <- as.data.table(table(ebb_test[vkood %in% ref_ind, .(EduYears)]))
  s2 <- as.data.table(table(ebb_test[vkood %in% biased_ind, .(EduYears)]))
  s <- merge(s1, s2, by = "V1")
  s[, ratio := N.y/N.x]
  ratio <- s[, min(ratio, max_ratio)]
  s1[, N := floor(N*ratio)]
  
  r2 <- c()
  for(i in 1:n){
    print(i)
    ea <- data.frame()
    for(l in 1:nrow(s1)){
      ea_i <- ebb_test[vkood %in% biased_ind & EduYears==as.numeric(s1[l,1]), ]
      ea_i <- ea_i[sample(x = 1:nrow(ea_i), size = as.numeric(s1[l,2]), replace = replacement), ]
      ea <- rbind(ea, ea_i)
    }
    
    r2 <- c(r2, ea[, cor(EduYears, PRS)^2])
    
  }
  r2_EA <- r2
  
  # match EduYears and Sex
  s1 <- as.data.table(table(ebb_test[vkood %in% ref_ind, .(EduYears, Sex)]))
  s2 <- as.data.table(table(ebb_test[vkood %in% biased_ind, .(EduYears, Sex)]))
  s <- merge(s1, s2, by = c("EduYears", "Sex"))
  s[, ratio := N.y/N.x]
  ratio <- s[, min(ratio, 1)]
  s1[, N := floor(N*ratio)]
  
  r2 <- c()
  for(i in 1:n){
    print(i+100)
    ea <- data.frame()
    for(l in 1:nrow(s1)){
      ea_i <- ebb_test[vkood %in% biased_ind & EduYears==as.numeric(s1[l,1]) & Sex==as.numeric(s1[l,2]), ]
      ea_i <- ea_i[sample(x = 1:nrow(ea_i), size = as.numeric(s1[l,3]), replace = replacement), ]
      ea <- rbind(ea, ea_i)
    }
    
    r2 <- c(r2, ea[, cor(EduYears, PRS)^2])
    
  }
  r2_EA_Sex <- r2
  
  
  return(list(r2_EA=r2_EA, r2_EA_Sex=r2_EA_Sex))
  
}

estbb_filtered <- fread("/gpfs/helios/home/kuznetsi/EBB_project/phenotypes/EstBB_filtered.tsv")
ebb <- fread("/gpfs/helios/home/kuznetsi/EBB_project/phenotypes/query1.tsv")
ebb <- ebb[, c("Person skood", "PersonLocation birthParishName", "PersonLocation residencyParishName", 
               "CONCATSTR(BMIAssembled ageAtBmi)", "CONCATSTR(BMIAssembled bmi)", "CONCATSTR(BMIAssembled height)")]
colnames(ebb) <- c("skood", "ParishBirth", "ParishRes", "Age_meas", "BMI", "Height")
ebb <- merge(estbb_filtered, ebb, by="skood")
ebb <- ebb[Nat=="Eestlane", ]
ebb[, Age2 := Age^2]
ebb[, SxA := Age*Sex]

# Height & BMI
ebb_hb <- ebb[, c("skood", "Age_meas", "BMI", "Height")]
ebb_hb[, Age_first := getValue(trait_vector = Age_meas, measurement = "first")]
ebb_hb[, BMI_first := getValue(trait_vector = BMI, measurement = "first")]
ebb_hb[, Height_first := getValue(trait_vector = Height, measurement = "first")]

ebb_hb <- ebb_hb[!(is.na(Age_first) | is.na(BMI_first) | is.na(Height_first)), ]

ebb_hb[, Age_last := getValue(trait_vector = Age_meas, measurement = "last")]
ebb_hb[, BMI_last := getValue(trait_vector = BMI, measurement = "last")]
ebb_hb[, Height_last := getValue(trait_vector = Height, measurement = "last")]

ebb_hb <- ebb_hb[BMI_first != 0 & BMI_last != 0 & Height_first != 0 & Height_last != 0]
ebb_hb[, BMI_first := log10(BMI_first)]
ebb_hb[, BMI_last := log10(BMI_last)]
ebb_hb[, Height_first := log10(Height_first)]
ebb_hb[, Height_last := log10(Height_last)]

ebb_hb <- ebb_hb[abs(BMI_first-mean(BMI_first))<4*sd(BMI_first) &
                   abs(BMI_last-mean(BMI_last))<4*sd(BMI_last) &
                   abs(Height_first-mean(Height_first))<4*sd(Height_first) &
                   abs(Height_last-mean(Height_last))<4*sd(Height_last),]

ebb <- merge(ebb, ebb_hb, by="skood", all.x = TRUE)


# EA
ebb2 <- fread("/gpfs/helios/home/kuznetsi/EBB_project/phenotypes/query_bmi_edu.tsv")
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA_portrait")
ebb <- merge(ebb, ebb2, by="skood")
ebb[, EduYears := transformEAtoEduYears(EA_portrait, to="years")]
ebb[is.na(EduYears), EduYears := NA]


# OS
ebb2 <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/data/query_OccupationStatus.tsv")
ebb2 <- ebb2[, c("Person skood", "CONCATSTR(Work currentOccupation code)", "CONCATSTR(Work mainOccupation code)")]
colnames(ebb2) <- c("skood", "curOcc", "mainOcc")
ebb2[, curOcc := sapply(curOcc, getOccStat)]
ebb2[, mainOcc := sapply(mainOcc, getOccStat)]
ebb2[, OS := rowMeans(.SD, na.rm = T), .SDcols = c("curOcc", "mainOcc")]
ebb2[is.na(OS), OS := NA]
ebb <- merge(ebb, ebb2, by="skood")


# PCA
pca_est <- fread("/gpfs/helios/home/kuznetsi/EBB_project/data_filtering/pca/pcs_EstBB_estonian")
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
ebb <- merge(ebb, pca_est, by="vkood")


# PGS
prs <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/data/EA4_excl_23andMe_EGCUT/PRS.chrALL.sscore")
colnames(prs) <- c("vkood", "PRS_EA")
ebb <- merge(ebb, prs)

prs <- fread("~/EBB_project/PRSs/pan-UKB/50/50.chrALL.sscore")
colnames(prs) <- c("vkood", "PRS_Height")
ebb <- merge(ebb, prs)

prs <- fread("~/EBB_project/PRSs/pan-UKB/21001/21001.chrALL.sscore")
colnames(prs) <- c("vkood", "PRS_BMI")
ebb <- merge(ebb, prs)


ebb <- ebb[AgeAtAgr > 25,]


ind_list <- getIndepInd(cutoff = 15)
ind_01 <- ind_list[[1]]
ind_02 <- ind_list[[2]]
ind_11 <- ind_list[[3]]
ind_12 <- ind_list[[4]]
ind_21 <- ind_list[[5]]
ind_22 <- ind_list[[6]]



# Census
census <- fread("/gpfs/helios/home/kuznetsi/EA_heritability/data/Cens_edu_combined_RLV302_20240619-103825.csv")
census <- census[Year != 2000,]
census[, Unknown := NULL]
census[Sex=="Male", Sex := "1"]
census[Sex=="Female", Sex := "2"]


# Matching EA categories with census
ebb[, EA_categ_tmp := EA_portrait]
ebb[EA_categ_tmp <= 3, EA_categ := 1]
ebb[EA_categ_tmp > 3 & EA_categ_tmp < 6, EA_categ := 2]
ebb[EA_categ_tmp >= 6 & EA_categ_tmp < 9, EA_categ := 3]
ebb <- ebb[!is.na(EA_categ),] 
ebb <- ebb[vkood %in% c(ind_11, ind_12, ind_21, ind_22),]


# Age categories as in Census
list_phases <- list(c(ind_11, ind_12), c(ind_21, ind_22))
for(i in 1:2){
  
  if(i==1){ Y<-2011 } else{ Y<-2021 }
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 25 & Y - YoB <= 29, Age_group := "25-29"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 30 & Y - YoB <= 34, Age_group := "30-34"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 35 & Y - YoB <= 39, Age_group := "35-39"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 40 & Y - YoB <= 44, Age_group := "40-44"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 45 & Y - YoB <= 49, Age_group := "45-49"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 50 & Y - YoB <= 54, Age_group := "50-54"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 55 & Y - YoB <= 59, Age_group := "55-59"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 60 & Y - YoB <= 64, Age_group := "60-64"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 65 & Y - YoB <= 69, Age_group := "65-69"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 70 & Y - YoB <= 74, Age_group := "70-74"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 75 & Y - YoB <= 79, Age_group := "75-79"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 80 & Y - YoB <= 84, Age_group := "80-84"]
  ebb[vkood %in% list_phases[[i]] & Y - YoB >= 85, Age_group := "85 and older"]
  
  
}

ebb <- ebb[!is.na(Age_group),]


getWeights <- function(ebb, ind, Y, census, correct = "all"){
  
  ebb <- ebb[data.table(vkood = ind), on = "vkood", nomatch = 0]
  
  ebb_stat <- dcast(
    ebb[, .N, by = c("Sex", "EA_categ", "Age_group")],
    Age_group + Sex ~ EA_categ,
    value.var = "N"
  )
  
  census_restr <- census[Year==Y & Age %in% unique(ebb$Age_group),]
  
  a <- as.matrix(ebb_stat[, 3:5])
  b <- as.matrix(census_restr[, 4:6])
  
  
  if(correct == "all"){
    
    a <- a/sum(a)
    b <- b/sum(b)
    
  }else if(correct == "bygroup"){
    
    a <- a/rowSums(a)
    b <- b/rowSums(b)
    
  }else if(correct == "none"){
    
    a <- a/a
    b <- b/b
    
  }
  
  weight <- b/a
  weight[!is.finite(weight)] <- 0
  weight <- data.table(census_restr[, 1:3], weight)
  colnames(weight)[4:6] <- c("1", "2", "3")
  
  weight <- melt(weight, id.vars = c("Year", "Age", "Sex"), variable.name = "EA", value.name = "N")
  colnames(weight)[5] <- "wt"
  weight <- as.data.table(weight)[, .(Age, Sex, EA, wt)]
  weight[, Sex := as.integer(Sex)]
  weight[, EA := as.double(EA)]
  
  ebb <- merge(ebb, weight, by.x = c("Age_group", "Sex", "EA_categ"), by.y = c("Age", "Sex", "EA"))
  
  
  lm_0 <- lm(paste0("EduYears ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("EduYears ~ PRS_EA + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_EA <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  
  lm_0 <- lm(paste0("OS ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("OS ~ PRS_EA + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_OS <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  
  lm_0 <- lm(paste0("Height_first ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("Height_first ~ PRS_Height + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_Height <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  
  lm_0 <- lm(paste0("BMI_first ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("BMI_first ~ PRS_BMI + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_BMI <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  r2_inc <- c(r2_inc_EA, r2_inc_OS, r2_inc_Height, r2_inc_BMI)
  
  return(r2_inc)
  
}


bootstrapR2 <- function(ebb, ind, Y, census, correct = "all", B = 100, seed = 123) {
  
  set.seed(seed)
  # Initialize a B x 4 matrix to store R² for each trait
  r2_boot <- matrix(NA_real_, nrow = B, ncol = 4)
  
  r2_orig <- getWeights(ebb, ind, Y, census, correct)
  
  for (b in 1:B) {
    # Sample with replacement from indices
    ind_boot <- sample(ind, length(ind), replace = TRUE)
    
    # Run the weighted R2 estimation
    r2_boot[b, ] <- getWeights(ebb, ind_boot, Y, census, correct)
  }
  
  # Compute CI
  ci_lower <- apply(r2_boot, 2, quantile, probs = 0.025)
  ci_upper <- apply(r2_boot, 2, quantile, probs = 0.975)
  r2_mean <- apply(r2_boot, 2, mean)
  
  return(list(
    r2_orig = r2_orig,
    r2_distribution = r2_boot,
    r2_mean = r2_mean,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

compareR2 <- function(r2_res, trait_n = 1){
  
  traits <- c("EA", "OS", "Height", "BMI")
  trait <- traits[trait_n]
  
  bootstrap <- data.table(r2_res$r_11$r2_distribution[,trait_n],
                          r2_res$r12$r2_distribution[,trait_n],
                          r2_res$r_21$r2_distribution[,trait_n],
                          r2_res$r_22$r2_distribution[,trait_n])
  colnames(bootstrap) <- c("p1ps", "p1s", "p2ps", "p2s")
  
  p1_s_ps <- bootstrap[, sum(p1s < p1ps)]
  p2_s_ps <- bootstrap[, sum(p2s < p2ps)]
  
  p1p2_s <- bootstrap[, sum(p1s < p2s)]
  p1p2_ps <- bootstrap[, sum(p1ps < p2ps)]
  
  comparison_r2 <- c(p1_s_ps, p2_s_ps, p1p2_s, p1p2_ps)
  pval_r2 <- round(sapply(comparison_r2, function(x) min(x/1000, 1 - x/1000)*2), 3)
  names_r2 <- c("p1_s_ps", "p2_s_ps", "p1p2_s", "p1p2_ps")
  names(pval_r2) <- names_r2
  
  return(pval_r2)
  
}


###############################
# With bootstrap CI
###############################


r_11 <- bootstrapR2(ebb = ebb, ind = ind_11, Y = 2011, census = census, correct = "bygroup", B = 1000)
r_12 <- bootstrapR2(ebb = ebb, ind = ind_12, Y = 2011, census = census, correct = "bygroup", B = 1000)
r_21 <- bootstrapR2(ebb = ebb, ind = ind_21, Y = 2021, census = census, correct = "bygroup", B = 1000)
r_22 <- bootstrapR2(ebb = ebb, ind = ind_22, Y = 2021, census = census, correct = "bygroup", B = 1000)

r2_bygroup <- data.table(weights = "by group", # byEA_in_sex_age_groups
                         trait = rep(c("EA", "OS", "Height", "BMI")),
                         Era = rep(c("p1ps", "p1s", "p2ps", "p2s"), each = 4),
                         r2 = c(r_11$r2_orig, r_12$r2_orig, r_21$r2_orig, r_22$r2_orig),
                         ci_low = c(r_11$ci_lower, r_12$ci_lower, r_21$ci_lower, r_22$ci_lower),
                         ci_up = c(r_11$ci_upper, r_12$ci_upper, r_21$ci_upper, r_22$ci_upper))

r2_list_bygroup <- list(r_11 = r_11, r12 = r_12, r_21 = r_21, r_22 = r_22)


r_11 <- bootstrapR2(ebb = ebb, ind = ind_11, Y = 2011, census = census, correct = "all", B = 1000)
r_12 <- bootstrapR2(ebb = ebb, ind = ind_12, Y = 2011, census = census, correct = "all", B = 1000)
r_21 <- bootstrapR2(ebb = ebb, ind = ind_21, Y = 2021, census = census, correct = "all", B = 1000)
r_22 <- bootstrapR2(ebb = ebb, ind = ind_22, Y = 2021, census = census, correct = "all", B = 1000)

r2_all <- data.table(weights = "overall", # byEA_sex_age
                     trait = rep(c("EA", "OS", "Height", "BMI")),
                     Era = rep(c("p1ps", "p1s", "p2ps", "p2s"), each = 4),
                     r2 = c(r_11$r2_orig, r_12$r2_orig, r_21$r2_orig, r_22$r2_orig),
                     ci_low = c(r_11$ci_lower, r_12$ci_lower, r_21$ci_lower, r_22$ci_lower),
                     ci_up = c(r_11$ci_upper, r_12$ci_upper, r_21$ci_upper, r_22$ci_upper))

r2_list_all <- list(r_11 = r_11, r12 = r_12, r_21 = r_21, r_22 = r_22)


r_11 <- bootstrapR2(ebb = ebb, ind = ind_11, Y = 2011, census = census, correct = "none", B = 1000)
r_12 <- bootstrapR2(ebb = ebb, ind = ind_12, Y = 2011, census = census, correct = "none", B = 1000)
r_21 <- bootstrapR2(ebb = ebb, ind = ind_21, Y = 2021, census = census, correct = "none", B = 1000)
r_22 <- bootstrapR2(ebb = ebb, ind = ind_22, Y = 2021, census = census, correct = "none", B = 1000)

r2_unadj <- data.table(weights = "original",
                       trait = rep(c("EA", "OS", "Height", "BMI")),
                       Era = rep(c("p1ps", "p1s", "p2ps", "p2s"), each = 4),
                       r2 = c(r_11$r2_orig, r_12$r2_orig, r_21$r2_orig, r_22$r2_orig),
                       ci_low = c(r_11$ci_lower, r_12$ci_lower, r_21$ci_lower, r_22$ci_lower),
                       ci_up = c(r_11$ci_upper, r_12$ci_upper, r_21$ci_upper, r_22$ci_upper))

r2_list_unadj <- list(r_11 = r_11, r12 = r_12, r_21 = r_21, r_22 = r_22)


r2 <- rbind(r2_unadj, r2_bygroup, r2_all)
r2 <- as.data.table(r2)

r2_list <- list(unadj = r2_list_unadj, bygroup = r2_list_bygroup, all = r2_list_all)

saveRDS(r2_list, file = "~/EA_heritability/results/r2_15_weighted_ci_alltraits.RDS")
write.table(r2,  "~/EA_heritability/results/r2_15_weighted_ci_alltraits.tsv", col.names = T, row.names = F, quote = F, sep = "\t")



# calculate p-values
r2_list <- readRDS("~/EA_heritability/results/r2_15_weighted_ci_alltraits.RDS")
pval_res <- data.frame()
for(i in 1:4){
  pval <- compareR2(r2_res = r2_list$unadj, trait_n = i)
  pval_res <- rbind(pval_res, c("No weights", pval))
  pval <- compareR2(r2_res = r2_list$bygroup, trait_n = i)
  pval_res <- rbind(pval_res, c("Inside groups", pval))
  pval <- compareR2(r2_res = r2_list$all, trait_n = i)
  pval_res <- rbind(pval_res, c("Overall", pval))
}
pval_res <- cbind(rep(c("EA", "OS", "Height", "BMI"), each = 3), pval_res)
colnames(pval_res) <- c("Trait", "weighting", "p1_s_ps", "p2_s_ps", "p1p2_s", "p1p2_ps")
write.table(pval_res, "~/EA_heritability/figures/paper/r2_15_weighted_ci_alltraits_pval.tsv", 
            col.names = T, row.names = F, quote = F, sep = "\t")


# make plots
r2 <- fread("~/EA_heritability/results/r2_15_weighted_ci_alltraits.tsv")
r2[, Era := factor(rep(rep(c("phase1\npost-soviet", "phase1\nsoviet", "phase2\npost-soviet", "phase2\nsoviet"), each = 4), 3),
                   levels = c("phase1\nsoviet", "phase1\npost-soviet", "phase2\nsoviet", "phase2\npost-soviet"))]
r2[, weights := factor(weights, levels = c("original", "by group", "overall"))]


# grid.newpage()
# grid.draw(leg)
# plot(leg)

# Make palette
# Red and blue shades (darker to lighter or vice versa)
red_shades  <- c("#FF9999", "#FF3333")
blue_shades <- c("#ADD8E6", "#5094CD")
r2[, fill_group := interaction(Era, weights, sep = "_")]
# Assign names matching fill_group levels
custom_colors <- c(
  # Red: phase1 and phase2 soviet
  "phase1\nsoviet_original"       = red_shades[1],
  "phase1\nsoviet_by group"       = red_shades[1],
  "phase1\nsoviet_overall"        = red_shades[1],
  
  "phase2\nsoviet_original"       = red_shades[2],
  "phase2\nsoviet_by group"       = red_shades[2],
  "phase2\nsoviet_overall"        = red_shades[2],
  
  # Blue: post-soviet
  "phase1\npost-soviet_original"  = blue_shades[1],
  "phase1\npost-soviet_by group"  = blue_shades[1],
  "phase1\npost-soviet_overall"   = blue_shades[1],
  
  "phase2\npost-soviet_original"  = blue_shades[2],
  "phase2\npost-soviet_by group"  = blue_shades[2],
  "phase2\npost-soviet_overall"   = blue_shades[2]
)


pl_list <- list()
for(tr in c("EA", "OS", "Height", "BMI")){
  
  pl <- ggplot(r2[trait == tr, ], aes(x = weights, y = r2, fill = fill_group)) + 
    geom_bar(stat='identity', position = position_dodge(width = 0.8), width=.6, color = "black") +
    geom_errorbar(aes(ymin =  ci_low, ymax = ci_up), width = 0.3,
                  position = position_dodge(width = 0.8, preserve = "single")) +
    theme_bw() + theme(text = element_text(size=10),
                       panel.grid.major.x = element_blank(),
                       axis.title.x=element_blank(),
                       legend.position="none") +
    # scale_fill_manual(values = c("#FFE699", "#FDAE6B", "#DA70D6"), labels = c("No weights", "Inside\nsex-age groups", "In overall sample")) +
    scale_fill_manual(values = custom_colors) +
    scale_x_discrete(labels = c("No weighting", "Weigting inside\nsex-age groups", "Weigting in\noverall sample")) + labs(y = bquote(R^2), fill = "Weighting") + ggtitle(tr)
  
  pl_list <- append(pl_list, list(pl))
  
}


empty_plot <- ggplot() + theme_void()
pl_list <- append(pl_list, list(empty_plot))



# Build a legend
# legend_labels <- c("No weights", "Inside\nsex-age groups", "In overall\nsample")
legend_labels <- c("wave 1\nSoviet", "wave 1\npost-Soviet", "wave 2\nSoviet", "wave 2\npost-Soviet")
# Build a legend data frame: two rectangles per row (red + blue)
# legend_df <- data.frame(
#   weight_label = rep(legend_labels, each = 2),
#   color_type = rep(c("red", "blue"), times = 3),
#   xmin = rep(c(0, 0.3), times = 3),
#   xmax = rep(c(0.3, 0.6), times = 3),
#   ymin = rep(2:0+0.1, each = 2),
#   ymax = rep(3:1-0.1, each = 2),
#   fill = c(rbind(red_shades, blue_shades))
# )
legend_df <- data.frame(
  weight_label = legend_labels,
  xmin = rep(0, times = 4),
  xmax = rep(0.6, times = 4),
  ymin = 3:0+0.1,
  ymax = 4:1-0.1,
  fill = c(rbind(red_shades, blue_shades))
)
# We'll add a new data frame for the labels, positioned just right of rectangles
label_df <- data.frame(
  weight_label = legend_labels,
  x = 0.7,  # just right of rectangles
  y = 3:0 + 0.5  # centered vertically in each row
)
legend_title <- data.frame(x = 0, y = 4.2, label = "Wave and era")
leg <- ggplot() +
  geom_rect(data = legend_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "black") +
  geom_text(data = label_df,
            aes(x = x, y = y, label = weight_label),
            hjust = 0, size=3, lineheight = 0.8) +
  geom_text(data = legend_title,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 0, size = 3.5) +  # Title text on top-left
  scale_fill_identity() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +  # wider space for text
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
  theme_void() +
  theme(plot.margin = margin(5, 20, 5, 5))  # bigger right margin

pl_list <- append(pl_list, list(leg))


pdf("~/EA_heritability/figures/paper/figure7.pdf", width = 5.5, height = 6)

print(
  grid.arrange(
    grobs = pl_list,
    layout_matrix = matrix(c(rep(1:4, each=2), 5,5,5,6,6,5,5,5), ncol = 2, byrow = F),
    widths = c(1, .35)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.02, y = 0.73, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.02, y = 0.23, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()

