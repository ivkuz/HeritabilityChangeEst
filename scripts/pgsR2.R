##########################################

# All the calculations with R2 ###########
# 1. Prepare the data for the analysis ###
# 2. Main R2 analysis ####################
# 3. Weighting ###########################
# 4. Matching distributions ##############
# 5. Original sample #####################

##########################################


library(data.table)
library(reshape2)
library(stringr)

# Transfer personal data (register) EA from categories to Years of Education
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


# Transfer answerset (self-reported) EA from categories to Years of Education
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


# get the first or the last value from a variable with several measurements
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


# extract the Ocupational Status values
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


# download the lists of unrelated individuals
getIndepInd <- function(cutoff = 15){
  
  if(cutoff == 15){
    
    ind_01 <- fread("~/EA_heritability/gcta/data/ps_list_unrel.tsv", header=F)
    ind_02 <- fread("~/EA_heritability/gcta/data/s_list_unrel.tsv", header=F)
    ind_11 <- fread("~/EA_heritability/gcta/data/p1ps_list_unrel.tsv", header=F)
    ind_12 <- fread("~/EA_heritability/gcta/data/p1s_list_unrel.tsv", header=F)
    ind_21 <- fread("~/EA_heritability/gcta/data/p2ps_list_unrel.tsv", header=F)
    ind_22 <- fread("~/EA_heritability/gcta/data/p2s_list_unrel.tsv", header=F)
    
  }else if(cutoff == 10){
    
    ind_01 <- fread("~/EA_heritability/gcta/data/ps_list_unrel_10.tsv", header=F)
    ind_02 <- fread("~/EA_heritability/gcta/data/s_list_unrel_10.tsv", header=F)
    ind_11 <- fread("~/EA_heritability/gcta/data/p1ps_list_unrel_10.tsv", header=F)
    ind_12 <- fread("~/EA_heritability/gcta/data/p1s_list_unrel_10.tsv", header=F)
    ind_21 <- fread("~/EA_heritability/gcta/data/p2ps_list_unrel_10.tsv", header=F)
    ind_22 <- fread("~/EA_heritability/gcta/data/p2s_list_unrel_10.tsv", header=F)
    
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


#######################################
# Prepare the data for the analysis ###
#######################################

# Upload main data
estbb_filtered <- fread("~/EBB_project/phenotypes/EstBB_filtered.tsv")
ebb <- fread("~/EBB_project/phenotypes/query1.tsv")
ebb <- ebb[, c("Person skood", "PersonLocation birthParishName", "PersonLocation residencyParishName", 
               "CONCATSTR(BMIAssembled ageAtBmi)", "CONCATSTR(BMIAssembled bmi)", "CONCATSTR(BMIAssembled height)")]
colnames(ebb) <- c("skood", "ParishBirth", "ParishRes", "Age_meas", "BMI", "Height")
ebb <- merge(estbb_filtered, ebb, by="skood")
ebb <- ebb[Nat=="Eestlane", ]

# Process Height & BMI
# Extract the first measurements
ebb_hb <- ebb[, c("skood", "Age_meas", "BMI", "Height")]
ebb_hb[, Age_first := getValue(trait_vector = Age_meas, measurement = "first")]
ebb_hb[, BMI_first := getValue(trait_vector = BMI, measurement = "first")]
ebb_hb[, Height_first := getValue(trait_vector = Height, measurement = "first")]
# Exclude empty values
ebb_hb <- ebb_hb[!(is.na(Age_first) | is.na(BMI_first) | is.na(Height_first)), ]

# Extract the last measurements
ebb_hb[, Age_last := getValue(trait_vector = Age_meas, measurement = "last")]
ebb_hb[, BMI_last := getValue(trait_vector = BMI, measurement = "last")]
ebb_hb[, Height_last := getValue(trait_vector = Height, measurement = "last")]

# Exclude zero values
ebb_hb <- ebb_hb[BMI_first != 0 & BMI_last != 0 & Height_first != 0 & Height_last != 0]
# Transform to the logarithmic scale
ebb_hb[, BMI_first := log10(BMI_first)]
ebb_hb[, BMI_last := log10(BMI_last)]
ebb_hb[, Height_first := log10(Height_first)]
ebb_hb[, Height_last := log10(Height_last)]

# Exclude outliers 
ebb_hb <- ebb_hb[abs(BMI_first-mean(BMI_first))<4*sd(BMI_first) &
                   abs(BMI_last-mean(BMI_last))<4*sd(BMI_last) &
                   abs(Height_first-mean(Height_first))<4*sd(Height_first) &
                   abs(Height_last-mean(Height_last))<4*sd(Height_last),]

ebb <- merge(ebb, ebb_hb, by="skood", all.x = TRUE)


# Upload and process EA
ebb2 <- fread("~/EBB_project/phenotypes/query_bmi_edu.tsv")
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA_portrait")
ebb <- merge(ebb, ebb2, by="skood")
ebb[, EduYears := transformEAtoEduYears(EA_portrait, to="years")]
ebb[is.na(EduYears), EduYears := NA]


# Upload and process OS
ebb2 <- fread("~/EA_heritability/data/query_OccupationStatus.tsv")
ebb2 <- ebb2[, c("Person skood", "CONCATSTR(Work currentOccupation code)", "CONCATSTR(Work mainOccupation code)")]
colnames(ebb2) <- c("skood", "curOcc", "mainOcc")
ebb2[, curOcc := sapply(curOcc, getOccStat)]
ebb2[, mainOcc := sapply(mainOcc, getOccStat)]
ebb2[, OS := rowMeans(.SD, na.rm = T), .SDcols = c("curOcc", "mainOcc")]
ebb2[is.na(OS), OS := NA]
ebb <- merge(ebb, ebb2, by="skood")


# Upload PCA
pca_est <- fread("~/EBB_project/data_filtering/pca/pcs_EstBB_estonian")
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
ebb <- merge(ebb, pca_est, by="vkood")


# Upload polygenic scores (PGS)
prs <- fread("~/EA_heritability/data/EA4_excl_23andMe_EGCUT/PRS.chrALL.sscore")
colnames(prs) <- c("vkood", "PRS_EA")
ebb <- merge(ebb, prs)

prs <- fread("~/EBB_project/PRSs/pan-UKB/50/50.chrALL.sscore")
colnames(prs) <- c("vkood", "PRS_Height")
ebb <- merge(ebb, prs)

prs <- fread("~/EBB_project/PRSs/pan-UKB/21001/21001.chrALL.sscore")
colnames(prs) <- c("vkood", "PRS_BMI")
ebb <- merge(ebb, prs)

# Upload the list of unrelated individuals for each group on Era and Participation Wave
ind_list <- getIndepInd(cutoff = 15)
ind_01 <- ind_list[[1]]
ind_02 <- ind_list[[2]]
ind_11 <- ind_list[[3]]
ind_12 <- ind_list[[4]]
ind_21 <- ind_list[[5]]
ind_22 <- ind_list[[6]]



####################################
###### Main R2 analysis ############
####################################

# Main function for R2 analysis
analyzeR2 <- function(ebb_test, pheno, age, prs, bootstrap_n = 1000){
  
  # Format analysed data
  colnames(ebb_test)[which(colnames(ebb_test)==pheno)] <- "Trait"
  colnames(ebb_test)[which(colnames(ebb_test)==age)] <- "Age"
  colnames(ebb_test)[which(colnames(ebb_test)==prs)] <- "PRS"
  
  ebb_test[, Age2 := Age^2]
  ebb_test[, SxA := Age*Sex]
  
  ebb_test <- ebb_test[!is.na(Trait), ]
  
  # Generate pre-adjusted trait
  lm_prs2 <- lm(paste0("Trait ~ Sex + Age + I(Age^2) + I(Sex*Age) + ", paste("PC", 1:40, sep = "", collapse = " + ")),
                data = ebb_test)
  ebb_test$Trait_adj2 <- ebb_test$Trait - predict(object = lm_prs2, newdata = ebb_test)
  ebb_test[, Trait_adj2 := (Trait_adj2-mean(Trait_adj2))/sd(Trait_adj2)]
  
  # Create empty vectors for summary output values
  rsq <- c()
  rsq_top95 <- c()
  rsq_bottom95 <- c()
  rsq_inc <- c()
  rsq_inc_top95 <- c()
  rsq_inc_bottom95 <- c()
  N <- c()
  dt <- data.table()
  
  rsq_all <- c() # for with main results
  n <- c() # for the table with main results
  
  # Bootstrap for each group
  for(ind in list(ind_02, ind_01, ind_12, ind_11, ind_22, ind_21)){
    
    set.seed(100)
    
    # Create empty vectors for R2 values
    rsq_bootstrap <- c()
    rsq_inc_bootstrap <- c()
    
    # Select individuals from the current group
    ebb_test_bootstrap <- ebb_test[vkood %in% ind]
    
    # Bootstrap
    for(i in 1:bootstrap_n){
      
      # Select subsample
      ebb_test_bootstrap_i <- ebb_test_bootstrap[sample(1:nrow(ebb_test_bootstrap),
                                                        nrow(ebb_test_bootstrap), replace = T), ]
      
      # R2 with pre-adjusted phenotype 
      lm_xx <- lm("Trait_adj2 ~ PRS", data = ebb_test_bootstrap_i)
      rsq_bootstrap <- c(rsq_bootstrap, summary(lm_xx)$r.squared)
      
      # Incremental R2
      lm_0 <- lm(paste0("Trait ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")),
                 data = ebb_test_bootstrap_i)
      lm_0prs <- lm(paste0("Trait ~ PRS + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")),
                    data = ebb_test_bootstrap_i)
      r2_inc <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
      rsq_inc_bootstrap <- c(rsq_inc_bootstrap, r2_inc)
      
    }
    
    # Extract R2 confidence interval
    rsq_top95 <- c(rsq_top95, apply(rsq_bootstrap, 2, quantile, probs = 0.975))
    rsq_bottom95 <- c(rsq_bottom95, apply(rsq_bootstrap, 2, quantile, probs = 0.025))
    rsq_inc_top95 <- c(rsq_inc_top95, sort(rsq_inc_bootstrap)[975])
    rsq_inc_bottom95 <- c(rsq_inc_bottom95, sort(rsq_inc_bootstrap)[26])
    dt <- cbind(dt, rsq_bootstrap, rsq_inc_bootstrap)
    
    # Original R2 for pre-adjusted phenotype 
    lm_xx <- lm("Trait_adj2 ~ PRS", data = ebb_test_bootstrap)
    rsq <- c(rsq, summary(lm_xx)$r.squared)
    rsq_all <- c(rsq_all, summary(lm_xx)$r.squared) # for the table with main results
    
    # Original incremental R2
    lm_0 <- lm(paste0("Trait ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
               data = ebb_test_bootstrap)
    lm_0prs <- lm(paste0("Trait ~ PRS + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                  data = ebb_test_bootstrap)
    r2_inc <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
    rsq_inc <- c(rsq_inc, r2_inc)
    rsq_all <- c(rsq_all, r2_inc) # for the table with main results
    
    N <- c(N, nrow(ebb_test_bootstrap))
    n <- c(n, rep(nrow(ebb_test_bootstrap), 2)) # for the table with main results
    
  }
  
  cohort <- c("s_r2", "s_r2_inc", "ps_r2", "ps_r2_inc",
              "s_r2_p1", "s_r2_inc_p1", "ps_r2_p1", "ps_r2_inc_p1", 
              "s_r2_p2", "s_r2_inc_p2", "ps_r2_p2", "ps_r2_inc_p2") # for the table with main results
  r2 <- rep(c("r2", "r2_inc"), 6) # for the table with main results
  
  # Save bootstrap results
  colnames(dt) <- cohort
  write.table(dt, paste0("~/EA_heritability/results/r2_adj_", pheno, "_1000.tsv"),
              row.names = F, quote = F, sep = "\t")
  
  ################################################
  # save the table with main results #
  name <- c("S", "S", "PS", "PS",
            "p1S", "p1S", "p1PS", "p1PS",
            "p2S", "p2S", "p2PS", "p2PS")
  r2_res <- rbind(name, r2, round(rsq_all, 8), n)
  colnames(r2_res) <- cohort
  write.table(r2_res, paste0("~/EA_heritability/results/r2_adj_log_", pheno, ".tsv"),
              row.names = F, quote = F, sep = "\t")
  
  
}



# Make a copy of the dataset
ebb_test <- ebb

# We analyse individuals older than 25 years at the time of recruitment
ebb_test <- ebb_test[AgeAtAgr > 25,]

ebb_test[, cohort := ifelse(vkood %in% ind_11, "p1ps", ifelse(vkood %in% ind_12, "p1s", ifelse(vkood %in% ind_21, "p2ps", ifelse(vkood %in% ind_22, "p2s", NA))))]
ebb_test[, cohort := factor(cohort, levels = c("p1s", "p1ps", "p2s", "p2ps"))]


ages <- c("Age", "Age", "Age_first", "Age_first")
phenotypes <- c("EA", "OS", "Height_first", "BMI_first")
prses <- c("PRS_EA", "PRS_EA", "PRS_Height", "PRS_BMI")
for(i in 1:4){
  
  pheno <- phenotypes[i]
  age <- ages[i]
  prs <- prses[i]
  
  analyzeR2(ebb_test, pheno=pheno, age = ages[i], prs = prses[i], bootstrap_n = 1000)
  
}

###############################
# Weighting ###################
###############################

# Functions for Weighting #####

# Calculate post-stratification weights based on sex, age, EA and apply to EA, OS, height, and BMI
getWeights <- function(ebb, ind, Y, census, correct = "all"){
  
  # extract (possibly) bootstrap subsample with replacemen (if needed)
  ebb <- ebb[data.table(vkood = ind), on = "vkood", nomatch = 0]
  
  # Count individuals in categories
  ebb_stat <- dcast(
    ebb[, .N, by = c("Sex", "EA_categ", "Age_group")],
    Age_group + Sex ~ EA_categ,
    value.var = "N"
  )
  
  # extract needed subset from census
  census_restr <- census[Year==Y & Age %in% unique(ebb$Age_group),]
  
  # get matrixed with the numbers of the individuals in EstBB and Census
  a <- as.matrix(ebb_stat[, 3:5])
  b <- as.matrix(census_restr[, 4:6])
  
  # Calculate relative proportions of the categories
  # Correct for joint distribution of sex, age and EA
  if(correct == "all"){
    
    a <- a/sum(a)
    b <- b/sum(b)
    
    # Correct for EA distribution in each sex-age category
  }else if(correct == "bygroup"){
    
    a <- a/rowSums(a)
    b <- b/rowSums(b)
    
    # No correction (weights = 1)
  }else if(correct == "none"){
    
    a <- a/a
    b <- b/b
    
  }
  
  # Make weights
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
  
  # Calculate incremental R2
  # EA
  lm_0 <- lm(paste0("EduYears ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("EduYears ~ PRS_EA + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_EA <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  # OS
  lm_0 <- lm(paste0("OS ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("OS ~ PRS_EA + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_OS <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  # Height
  lm_0 <- lm(paste0("Height_first ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("Height_first ~ PRS_Height + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_Height <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  # BMI
  lm_0 <- lm(paste0("BMI_first ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
             data = ebb, weights = wt)
  lm_0prs <- lm(paste0("BMI_first ~ PRS_BMI + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                data = ebb, weights = wt)
  r2_inc_BMI <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
  
  r2_inc <- c(r2_inc_EA, r2_inc_OS, r2_inc_Height, r2_inc_BMI)
  
  return(r2_inc)
  
}


# Bootstrap sampling procedure for R2
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

# Make a copy of the dataset
ebb_test <- ebb

# We analyse individuals older than 25 years at the time of recruitment
ebb_test <- ebb_test[AgeAtAgr > 25,]

# Upload Census data (2011 and 2021)
census <- fread("~/EA_heritability/data/Cens_edu_combined_RLV302_20240619-103825.csv")
census <- census[Year != 2000,]
census[, Unknown := NULL]
census[Sex=="Male", Sex := "1"]
census[Sex=="Female", Sex := "2"]


# Matching EA categories with census
ebb_test[, EA_categ_tmp := EA_portrait]
ebb_test[EA_categ_tmp <= 3, EA_categ := 1]
ebb_test[EA_categ_tmp > 3 & EA_categ_tmp < 6, EA_categ := 2]
ebb_test[EA_categ_tmp >= 6 & EA_categ_tmp < 9, EA_categ := 3]
ebb_test <- ebb_test[!is.na(EA_categ),] 
ebb_test <- ebb_test[vkood %in% c(ind_11, ind_12, ind_21, ind_22),]


# Age categories as in Census
list_phases <- list(c(ind_11, ind_12), c(ind_21, ind_22))
for(i in 1:2){
  
  if(i==1){ Y<-2011 } else{ Y<-2021 }
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 25 & Y - YoB <= 29, Age_group := "25-29"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 30 & Y - YoB <= 34, Age_group := "30-34"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 35 & Y - YoB <= 39, Age_group := "35-39"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 40 & Y - YoB <= 44, Age_group := "40-44"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 45 & Y - YoB <= 49, Age_group := "45-49"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 50 & Y - YoB <= 54, Age_group := "50-54"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 55 & Y - YoB <= 59, Age_group := "55-59"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 60 & Y - YoB <= 64, Age_group := "60-64"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 65 & Y - YoB <= 69, Age_group := "65-69"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 70 & Y - YoB <= 74, Age_group := "70-74"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 75 & Y - YoB <= 79, Age_group := "75-79"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 80 & Y - YoB <= 84, Age_group := "80-84"]
  ebb_test[vkood %in% list_phases[[i]] & Y - YoB >= 85, Age_group := "85 and older"]
  
}

ebb_test <- ebb_test[!is.na(Age_group),]


r_11 <- bootstrapR2(ebb_test = ebb_test, ind = ind_11, Y = 2011, census = census, correct = "bygroup", B = 1000)
r_12 <- bootstrapR2(ebb_test = ebb_test, ind = ind_12, Y = 2011, census = census, correct = "bygroup", B = 1000)
r_21 <- bootstrapR2(ebb_test = ebb_test, ind = ind_21, Y = 2021, census = census, correct = "bygroup", B = 1000)
r_22 <- bootstrapR2(ebb_test = ebb_test, ind = ind_22, Y = 2021, census = census, correct = "bygroup", B = 1000)

r2_bygroup <- data.table(weights = "by group", # byEA_in_sex_age_groups
                         trait = rep(c("EA", "OS", "Height", "BMI")),
                         Era = rep(c("p1ps", "p1s", "p2ps", "p2s"), each = 4),
                         r2 = c(r_11$r2_orig, r_12$r2_orig, r_21$r2_orig, r_22$r2_orig),
                         ci_low = c(r_11$ci_lower, r_12$ci_lower, r_21$ci_lower, r_22$ci_lower),
                         ci_up = c(r_11$ci_upper, r_12$ci_upper, r_21$ci_upper, r_22$ci_upper))

r2_list_bygroup <- list(r_11 = r_11, r12 = r_12, r_21 = r_21, r_22 = r_22)


r_11 <- bootstrapR2(ebb_test = ebb_test, ind = ind_11, Y = 2011, census = census, correct = "all", B = 1000)
r_12 <- bootstrapR2(ebb_test = ebb_test, ind = ind_12, Y = 2011, census = census, correct = "all", B = 1000)
r_21 <- bootstrapR2(ebb_test = ebb_test, ind = ind_21, Y = 2021, census = census, correct = "all", B = 1000)
r_22 <- bootstrapR2(ebb_test = ebb_test, ind = ind_22, Y = 2021, census = census, correct = "all", B = 1000)

r2_all <- data.table(weights = "overall", # byEA_sex_age
                     trait = rep(c("EA", "OS", "Height", "BMI")),
                     Era = rep(c("p1ps", "p1s", "p2ps", "p2s"), each = 4),
                     r2 = c(r_11$r2_orig, r_12$r2_orig, r_21$r2_orig, r_22$r2_orig),
                     ci_low = c(r_11$ci_lower, r_12$ci_lower, r_21$ci_lower, r_22$ci_lower),
                     ci_up = c(r_11$ci_upper, r_12$ci_upper, r_21$ci_upper, r_22$ci_upper))

r2_list_all <- list(r_11 = r_11, r12 = r_12, r_21 = r_21, r_22 = r_22)


r_11 <- bootstrapR2(ebb_test = ebb_test, ind = ind_11, Y = 2011, census = census, correct = "none", B = 1000)
r_12 <- bootstrapR2(ebb_test = ebb_test, ind = ind_12, Y = 2011, census = census, correct = "none", B = 1000)
r_21 <- bootstrapR2(ebb_test = ebb_test, ind = ind_21, Y = 2021, census = census, correct = "none", B = 1000)
r_22 <- bootstrapR2(ebb_test = ebb_test, ind = ind_22, Y = 2021, census = census, correct = "none", B = 1000)

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


######################################
# Matching distributions #############
######################################

matchCohorts <- function(ref_ind, biased_ind, n, replacement, max_ratio=1, seed = 1){
  
  # The matching function
  # Take samples from biased sample matching with the reference one
  # Let's say "we are accounting for bias moving towards reference distribution"
  set.seed(seed = seed)
  
  # Match EduYears
  s1 <- as.data.table(table(ebb_test[vkood %in% ref_ind, .(EduYears)]))
  s2 <- as.data.table(table(ebb_test[vkood %in% biased_ind, .(EduYears)]))
  s <- merge(s1, s2, by = "EduYears")
  s[, ratio := N.y/N.x]
  ratio <- s[, min(ratio, max_ratio)]
  s1[, N := floor(N*ratio)]
  
  r2 <- c()
  for(i in 1:n){
    ea <- data.frame()
    for(l in 1:nrow(s1)){
      ea_i <- ebb_test[vkood %in% biased_ind & EduYears==as.numeric(s1[l,1]), ]
      ea_i <- ea_i[sample(x = 1:nrow(ea_i), size = as.numeric(s1[l,2]), replace = replacement), ]
      ea <- rbind(ea, ea_i)
    }
    
    r2 <- c(r2, ea[, cor(EduYears_adj2, PRS)^2])
    
  }
  r2_EA <- r2
  
  print("First part done")
  
  # match EduYears and Sex
  s1 <- as.data.table(table(ebb_test[vkood %in% ref_ind, .(EduYears, Sex)]))
  s2 <- as.data.table(table(ebb_test[vkood %in% biased_ind, .(EduYears, Sex)]))
  s <- merge(s1, s2, by = c("EduYears", "Sex"))
  s[, ratio := N.y/N.x]
  ratio <- s[, min(ratio, 1)]
  s1[, N := floor(N*ratio)]
  
  r2 <- c()
  for(i in 1:n){
    ea <- data.frame()
    for(l in 1:nrow(s1)){
      ea_i <- ebb_test[vkood %in% biased_ind & EduYears==as.numeric(s1[l,1]) & Sex==as.numeric(s1[l,2]), ]
      ea_i <- ea_i[sample(x = 1:nrow(ea_i), size = as.numeric(s1[l,3]), replace = replacement), ]
      ea <- rbind(ea, ea_i)
    }
    
    r2 <- c(r2, ea[, cor(EduYears_adj2, PRS)^2])
    
  }
  r2_EA_Sex <- r2
  
  
  return(list(r2_EA=r2_EA, r2_EA_Sex=r2_EA_Sex))
  
}



# Make a copy of the dataset
ebb_test <- ebb

# We analyse individuals older than 25 years at the time of recruitment
ebb_test <- ebb_test[AgeAtAgr > 25,]

# Pre-adjust the EA phenotype
lm_EA <- lm(paste0("EduYears ~ Sex + Age + I(Age^2) + I(Sex*Age) + ", paste("PC", 1:40, sep = "", collapse = " + ")),
            data = ebb_test)
ebb_test$EduYears_adj2 <- ebb_test$EduYears - predict(object = lm_EA, newdata = ebb_test)
ebb_test[, EduYears_adj2 := scale(EduYears_adj2)]



r2_list <- list()
for(i in 3:6){
  for(j in 3:6){
    print(c(i, j))
    r2_matched <- matchCohorts(ref_ind = ind_list[[i]], biased_ind = ind_list[[j]], n = 1000, replacement=T, max_ratio=10)
    r2_ref <- ebb_test[vkood %in% ind_list[[i]], cor(PRS, EduYears_adj2)^2]
    r2_biased <- ebb_test[vkood %in% ind_list[[j]], cor(PRS, EduYears_adj2)^2]
    r2_list_tmp <- list(r2_matched, rw_ref = r2_ref, r2_biased = r2_biased)
    r2_list <- append(r2_list, list(r2_list_tmp))
  }
}

saveRDS(r2_list, "~/EA_heritability/scripts/paper/r2_matched.RDS")




######################################
# Original sample ####################
######################################

# Main function for R2 analysis
analyzeR2_origSample <- function(ebb_test, pheno, age, prs, analysis_title, pcs = 40, bootstrap_n = 1000){
  
  # Format analysed data
  colnames(ebb_test)[which(colnames(ebb_test)==pheno)] <- "Trait"
  colnames(ebb_test)[which(colnames(ebb_test)==age)] <- "Age"
  colnames(ebb_test)[which(colnames(ebb_test)==prs)] <- "PRS"
  
  ebb_test[, Age2 := Age^2]
  ebb_test[, SxA := Age*Sex]
  
  ebb_test <- ebb_test[!is.na(Trait), ]
  
  # Generate pre-adjusted trait
  lm_prs2 <- lm(paste0("Trait ~ Sex + Age + I(Age^2) + I(Sex*Age) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
                data = ebb_test)
  ebb_test$Trait_adj2 <- ebb_test$Trait - predict(object = lm_prs2, newdata = ebb_test)
  ebb_test[, Trait_adj2 := (Trait_adj2-mean(Trait_adj2))/sd(Trait_adj2)]
  
  # Create empty vectors for summary output values
  rsq <- c()
  rsq_top95 <- c()
  rsq_bottom95 <- c()
  rsq_inc <- c()
  rsq_inc_top95 <- c()
  rsq_inc_bottom95 <- c()
  N <- c()
  dt <- data.table()
  
  rsq_all <- c() # for with main results
  n <- c() # for the table with main results
  
  # Bootstrap for each group
  for(ind in list(ind_02, ind_01, ind_12, ind_11, ind_22, ind_21)){
    
    set.seed(100)
    
    # Create empty vectors for R2 values
    rsq_bootstrap <- c()
    rsq_inc_bootstrap <- c()
    
    # Select individuals from the current group
    ebb_test_bootstrap <- ebb_test[vkood %in% ind]
    
    # Bootstrap
    for(i in 1:bootstrap_n){
      
      # Select subsample
      ebb_test_bootstrap_i <- ebb_test_bootstrap[sample(1:nrow(ebb_test_bootstrap),
                                                        nrow(ebb_test_bootstrap), replace = T), ]
      
      # R2 with pre-adjusted phenotype 
      lm_xx <- lm("Trait_adj2 ~ PRS", data = ebb_test_bootstrap_i)
      rsq_bootstrap <- c(rsq_bootstrap, summary(lm_xx)$r.squared)
      
      # Incremental R2
      lm_0 <- lm(paste0("Trait ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
                 data = ebb_test_bootstrap_i)
      lm_0prs <- lm(paste0("Trait ~ PRS + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
                    data = ebb_test_bootstrap_i)
      r2_inc <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
      rsq_inc_bootstrap <- c(rsq_inc_bootstrap, r2_inc)
      
    }
    
    # Extract R2 confidence interval
    rsq_top95 <- c(rsq_top95, apply(rsq_bootstrap, 2, quantile, probs = 0.975))
    rsq_bottom95 <- c(rsq_bottom95, apply(rsq_bootstrap, 2, quantile, probs = 0.025))
    rsq_inc_top95 <- c(rsq_inc_top95, sort(rsq_inc_bootstrap)[975])
    rsq_inc_bottom95 <- c(rsq_inc_bottom95, sort(rsq_inc_bootstrap)[26])
    dt <- cbind(dt, rsq_bootstrap, rsq_inc_bootstrap)
    
    # Original R2 for pre-adjusted phenotype 
    lm_xx <- lm("Trait_adj2 ~ PRS", data = ebb_test_bootstrap)
    rsq <- c(rsq, summary(lm_xx)$r.squared)
    rsq_all <- c(rsq_all, summary(lm_xx)$r.squared) # for the table with main results
    
    # Original incremental R2
    lm_0 <- lm(paste0("Trait ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")), 
               data = ebb_test_bootstrap)
    lm_0prs <- lm(paste0("Trait ~ PRS + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")), 
                  data = ebb_test_bootstrap)
    r2_inc <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
    rsq_inc <- c(rsq_inc, r2_inc)
    rsq_all <- c(rsq_all, r2_inc) # for the table with main results
    
    N <- c(N, nrow(ebb_test_bootstrap))
    n <- c(n, rep(nrow(ebb_test_bootstrap), 2)) # for the table with main results
    
  }
  
  colnames(dt) <- c("all_r2", "all_r2_inc", 
                    "s_15_r2", "s_15_r2_inc", "ps_15_r2", "ps_15_r2_inc", 
                    "s_10_r2", "s_10_r2_inc", "ps_10_r2", "ps_10_r2_inc",
                    "all_est_r2", "all_est_r2_inc", 
                    "s_15_est_r2", "s_15_est_r2_inc", "ps_15_est_r2", "ps_15_est_r2_inc", 
                    "s_10_est_r2", "s_10_est_r2_inc", "ps_10_est_r2", "ps_10_est_r2_inc") # for the table with main results
  r2 <- rep(c("r2", "r2_inc"), 10) # for the table with main results
  
  # Save bootstrap results
  colnames(dt) <- cohort
  write.table(dt, paste0("~/EA_heritability/results/r2_adj_", pheno, "_origSample.tsv"),
              row.names = F, quote = F, sep = "\t")
  
  ################################################
  # save the table with main results #
  name <- c("All", "All", "S15", "S15", "PS15", "PS15", "S10", "S10", "PS10", "PS10",
            "All", "All", "S15", "S15", "PS15", "PS15", "S10", "S10", "PS10", "PS10")
  r2_res <- rbind(name, r2, round(rsq_all, 10), n)
  colnames(r2_res) <- cohort
  write.table(r2_res, paste0("~/EA_heritability/results/r2_adj_log_", pheno, "_origSample.tsv"),
              row.names = F, quote = F, sep = "\t")
  
  
}


# Make a copy of the dataset
ebb_test <- ebb_test

# Filter for the original sample
origSample <- read_excel("~/EA_heritability/data/origSamples/2024-06-18_KRimfeld/I16_2024-06-18.xlsx")
ebb_test <- ebb_test[vkood %in% origSample$Vcode1, ]
ebb_test <- ebb_test[YoA-YoB >= 25, ]

# Years of Education from questionnaire 
ebb_test[, EA_a_order := transformSREAtoEduYears(EA_answerset, to="order")]
ebb_test[, EA_p_years := EduYears]

# Define the original groups
ind_all <- ebb_test[, vkood]
ind15_s <- ebb_test[YoB<1976, vkood]
ind15_ps <- ebb_test[YoB>=1976, vkood]
ind10_s <- ebb_test[YoB<1981, vkood]
ind10_ps <- ebb_test[YoB>=1981, vkood]
ind_all_est <- ebb_test[Nat == "Eestlane", vkood]
ind15_s_est <- ebb_test[Nat == "Eestlane" & YoB<1976, vkood]
ind15_ps_est <- ebb_test[Nat == "Eestlane" & YoB>=1976, vkood]
ind10_s_est <- ebb_test[Nat == "Eestlane" & YoB<1981, vkood]
ind10_ps_est <- ebb_test[Nat == "Eestlane" & YoB>=1981, vkood]

ind_list <- list(ind_all, ind15_s, ind15_ps, ind10_s, ind10_ps, ind_all_est, ind15_s_est, ind15_ps_est, ind10_s_est, ind10_ps_est)

phenotypes <- c("EA_p_years", "EA_a_order", "OS")
for(i in 1:4){
  
  pheno <- phenotypes[i]
  analyzeR2_origSample(ebb_test, pheno=pheno, age = "Age", prs = "PRS_EA", pc_n = 40, bootstrap_n = 1000)
  
}




