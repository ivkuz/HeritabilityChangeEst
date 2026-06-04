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

settlement <- fread("~/EA_heritability/data/query_settletype.tsv")
colnames(settlement) <- c("skood", "settlCode", "settlName")
ebb <- merge(ebb, settlement, by="skood")

# We analyse individuals older than 25 years at the time of recruitment
ebb <- ebb[AgeAtAgr > 25,]

covar <- fread("~/EA_heritability/gcta/data/covarLDAK.tsv", header = F)
covar[, V1 := V2]
covar <- merge(covar, ebb[, .(vkood, Age_last)], by.x = "V1", by.y = "vkood")
covar[, V1 := 0]
write.table(covar, "~/EA_heritability/gcta/data/covarBV_GCTA.tsv", col.names = F, row.names = F, quote = F, sep = "\t")

unrelated <- fread("/gpfs/space/GI/GV/EGCUT_data/genotype_data/GSA_arrays/4_latest_freeze/additional_files/king/unrelated_set/unrelated_degree_2unrelated.txt", header = F)
unrelated[, V1:=0]
write.table(unrelated, "~/EA_heritability/gcta/data/unrelatedBV_GCTA.tsv", col.names = F, row.names = F, quote = F, sep = "\t")

# pheno <- ebb[vkood %in% unrelated$V2, .(0, vkood, EduYears, OS, Height_last, BMI_last)]
pheno <- ebb[, .(0, vkood, EduYears, OS, Height_last, BMI_last)]

selected_ind <- pheno[vkood %in% unrelated$V2 & complete.cases(pheno), .(V1, vkood)]  
write.table(selected_ind, "~/EA_heritability/gcta/data/unrelated_complete_BV_GCTA.tsv", col.names = F, row.names = F, quote = F, sep = "\t")


ind_s_15 <- ebb[YoB < 1991 - 15, vkood]
ind_ps_15 <- ebb[YoB >= 1991 - 15, vkood]
ind_s_10 <- ebb[YoB < 1991 - 10, vkood]
ind_ps_10 <- ebb[YoB >= 1991 - 10, vkood]


# pheno file for GCTA
pheno[vkood %in% ind_s_15, `:=`(EduYears_s_15=EduYears, OS_s_15=OS, Height_s_15=Height_last, BMI_s_15=BMI_last)]
pheno[vkood %in% ind_s_10, `:=`(EduYears_s_10=EduYears, OS_s_10=OS, Height_s_10=Height_last, BMI_s_10=BMI_last)]
pheno[vkood %in% ind_ps_15, `:=`(EduYears_ps_15=EduYears, OS_ps_15=OS, Height_ps_15=Height_last, BMI_ps_15=BMI_last)]
pheno[vkood %in% ind_ps_10, `:=`(EduYears_ps_10=EduYears, OS_ps_10=OS, Height_ps_10=Height_last, BMI_ps_10=BMI_last)]

pheno[, `:=`(EduYears=NULL, OS=NULL, Height_last=NULL, BMI_last=NULL)]

write.table(pheno, "~/EA_heritability/gcta/data/phenoBV_GCTA.tsv", col.names = F, row.names = F, quote = F, sep = "\t")


t<-1
n<-45850
mem <- (t+4)*n*n*8/1023^3
mem



