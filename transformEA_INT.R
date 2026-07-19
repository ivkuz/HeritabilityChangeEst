library(data.table)
library(readr)
library(sjmisc)

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



# Upload main data
estbb_filtered <- fread("~/EBB_project/phenotypes/EstBB_filtered.tsv")
ebb <- fread("~/EBB_project/phenotypes/query1.tsv")
ebb <- ebb[, c("Person skood", "PersonLocation birthParishName", "PersonLocation residencyParishName", 
               "CONCATSTR(BMIAssembled ageAtBmi)", "CONCATSTR(BMIAssembled bmi)", "CONCATSTR(BMIAssembled height)")]
colnames(ebb) <- c("skood", "ParishBirth", "ParishRes", "Age_meas", "BMI", "Height")
ebb <- merge(estbb_filtered, ebb, by="skood")
ebb <- ebb[Nat=="Eestlane", ]



# Upload and process EA
ebb2 <- fread("~/EBB_project/phenotypes/query_bmi_edu.tsv")
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA_portrait")
ebb <- merge(ebb, ebb2, by="skood")
ebb[, EduYears := transformEAtoEduYears(EA_portrait, to="years")]
ebb[is.na(EduYears), EduYears := NA]


# Upload PCA
pca_est <- fread("~/EBB_project/data_filtering/pca/pcs_EstBB_estonian")
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
ebb <- merge(ebb, pca_est, by="vkood")


cohorts <- read_lines("~/EA_heritability/gcta/data/file_list_cohorts.tsv")

for(cohort in cohorts){
  
  if(str_contains(cohort, "_ldak_list_all_15.tsv")){
    pref <- "all"
    suff <- sub(".tsv", "", sub("_ldak_list_all", "", cohort))

  }else if(str_contains(cohort, "ldak_list_all_")){
    pref <- "all"
    suff <- sub(".tsv", "", sub("ldak_list_all_", "", cohort))
    
  }else if(str_contains(cohort, "ldak_kin0.05_") &
           str_contains(cohort, "_15.keep")){
    pref <- "kin0.05"
    suff <- sub(".keep", "", sub("ldak_kin0.05_", "", cohort))
    
  }else if(str_contains(cohort, "ldak_kin0.05_noAgeAtAgr_")){
    pref <- "kin0.05"
    suff <- sub(".tsv.keep", "", sub("ldak_kin0.05_", "", cohort))
  }
  
  ind <- fread(paste0("~/EA_heritability/gcta/data/", cohort))
  
  ebb_sub <- ebb[vkood %in% ind$V2, ]
  
  # Select individuals from the current group
  ebb_sub_m <- ebb_sub[Sex == 1, ]
  ebb_sub_f <- ebb_sub[Sex == 2, ]
  
  # Pre-adjust trait by sex
  # Male
  lm_m <- lm(paste0("EduYears ~ Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")),
                  data = ebb_sub_m)
  ebb_sub_m$EduYears_adj2 <- ebb_sub_m$EduYears - predict(object = lm_m, newdata = ebb_sub_m)
  ebb_sub_m[, EduYears_INT := qnorm((rank(EduYears_adj2)-0.5)/nrow(ebb_sub_m))]
  # Female
  lm_f <- lm(paste0("EduYears ~ Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")),
             data = ebb_sub_f)
  ebb_sub_f$EduYears_adj2 <- ebb_sub_f$EduYears - predict(object = lm_f, newdata = ebb_sub_f)
  ebb_sub_f[, EduYears_INT := qnorm((rank(EduYears_adj2)-0.5)/nrow(ebb_sub_f))]
  
  ebb_sub <- rbind(ebb_sub_m, ebb_sub_f)
  
  write.table(data.table(0, ebb_sub[, .(vkood, EduYears_INT)]),
              paste0("~/EA_heritability/gcta/data/phenoLDAK_EA_INT_", pref, "_", suff, ".tsv"),
              row.names = F, col.names = F, quote = F, sep = "\t")
  
}
