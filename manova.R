#####################################################################
# MANOVA analysis: compares jointly EA and OS averages between eras #
#####################################################################

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


# Upload the list of unrelated individuals for each group on Era and Participation Wave
ind_list <- getIndepInd(cutoff = 15)
ind_01 <- ind_list[[1]]
ind_02 <- ind_list[[2]]
ind_11 <- ind_list[[3]]
ind_12 <- ind_list[[4]]
ind_21 <- ind_list[[5]]
ind_22 <- ind_list[[6]]


# Make a copy of the dataset
ebb_test <- ebb

# We analyse individuals older than 25 years at the time of recruitment
ebb_test <- ebb_test[AgeAtAgr > 25,]

ebb_test[, cohort := ifelse(vkood %in% ind_01, "ps", ifelse(vkood %in% ind_02, "s", NA))]
ebb_test[, cohort := factor(cohort, levels = c("s", "ps"))]
ebb_test[, cohort_ph := ifelse(vkood %in% ind_11, "p1ps", ifelse(vkood %in% ind_12, "p1s", ifelse(vkood %in% ind_21, "p2ps", ifelse(vkood %in% ind_22, "p2s", NA))))]

hist(ebb_test[cohort == "s", EduYears])
hist(ebb_test[cohort == "ps", EduYears])
hist(ebb_test[cohort == "s", OS])
hist(ebb_test[cohort == "ps", OS])

library(ggpubr)
ggqqplot(ebb_test[cohort == "s", EduYears])
ggqqplot(ebb_test[cohort == "ps", EduYears])
ggqqplot(ebb_test[cohort == "s", OS])
ggqqplot(ebb_test[cohort == "ps", OS])

# shapiro.test(ebb_test[cohort == "s", EduYears])
# shapiro.test(ebb_test[cohort == "ps", EduYears])
# shapiro.test(ebb_test[cohort == "s", OS])
# shapiro.test(ebb_test[cohort == "ps", OS])


# Both phases
res.man <- manova(cbind(EduYears, OS) ~ cohort, data = ebb_test[!is.na(cohort),])
summary(res.man)
summary.aov(res.man)

# Phase 1
res.man_p1 <- manova(cbind(EduYears, OS) ~ cohort_ph, data = ebb_test[cohort_ph %in% c("p1s", "p1ps"),])
summary(res.man_p1)
summary.aov(res.man_p1)

# Phase 2
res.man_p2 <- manova(cbind(EduYears, OS) ~ cohort_ph, data = ebb_test[cohort_ph %in% c("p2s", "p2ps"),])
summary(res.man_p2)
summary.aov(res.man_p2)

# library(vegan)
# df <- ebb_test[!is.na(ebb_test$cohort) &!is.na(ebb_test$EduYears) & !is.na(ebb_test$OS), ]
# set.seed(42)
# df_sub <- df[sample(nrow(df), 1000), ]
# 
# res.man.perm <- adonis2(
#   df[, c("EduYears", "OS")] ~ cohort,
#   data = df_sub,
#   method = "euclidean",
#   permutations = 999,
#   by = "margin"
# )
# 


res.an.ea <- aov(EduYears ~ cohort, data = ebb_test)
summary(res.an.ea)
summary.aov(res.an.ea)

res.an.os <- aov(OS ~ cohort, data = ebb_test)
summary(res.an.os)
summary.aov(res.an.os)


library(ggplot2)
ggplot(ebb_test[!is.na(cohort),], aes(fill = cohort, x=EduYears)) + geom_histogram(aes(y = after_stat(density)), position = position_dodge())
ggplot(ebb_test[!is.na(cohort),], aes(fill = cohort, x=OS)) + geom_histogram(aes(y = after_stat(density)), position = position_dodge())
