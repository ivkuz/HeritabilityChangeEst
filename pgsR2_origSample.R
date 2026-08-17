# https://www.sciencedirect.com/science/article/pii/S0264999399000346

library(data.table)
library(ggplot2)
library(reshape2)
library(readxl)
library(stringr)


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


filter_relatives <- function(rel, vcodes, k=100){
  rel <- rel[ID1 %in% vcodes & ID2 %in% vcodes,]
  if(nrow(rel) == 0){
    return(vcodes)
  }
  u <- length(unique(c(rel$ID1, rel$ID2)))
  exclude <- c()
  i=0
  n=2
  while(n>1){
    i = i + 1
    print(n)
    rel_tab <- as.data.table(table(c(rel$ID1, rel$ID2)))
    rel_tab <- rel_tab[order(N, decreasing=T),]
    n <- rel_tab[1, N]
    hub <- rel_tab[N == n, V1][1:min(length(rel_tab[N == n, V1]), k)]
    exclude <- c(exclude, hub)
    rel <- rel[!(ID1 %in% hub) & !(ID2 %in% hub),]
  }
  exclude <- c(exclude, rel$ID1)
  non_relatives <- vcodes[!(vcodes %in% exclude)]
  return(non_relatives)
}

getOccStat <- function(x) {
  if (is.na(x)) return(NA)
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


analyzeR2 <- function(ebb_test, ind_list, pheno){
  
  pc_n = 40
  ebb_test <- ebb_test[!is.na(EduYears), ]
  lm_prs2 <- lm(paste0("EduYears ~ Sex + Age + I(Age^2) + I(Sex*Age) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
                # lm_prs2 <- lm("EduYears ~ Sex + Age + I(Age^2) + I(Sex*Age)",
                data = ebb_test)
  ebb_test$EduYears_adj2 <- ebb_test$EduYears - predict(object = lm_prs2, newdata = ebb_test)
  ebb_test[, EduYears_adj2 := (EduYears_adj2-mean(EduYears_adj2))/sd(EduYears_adj2)]
  
  rsq <- c()
  rsq_top95 <- c()
  rsq_bottom95 <- c()
  rsq_inc <- c()
  rsq_inc_top95 <- c()
  rsq_inc_bottom95 <- c()
  N <- c()
  dt <- data.table()
  
  rsq_all <- c() # added to calculate a table with main results
  n <- c() # added to calculate a table with main results
  
  j=0
  for(ind in ind_list){
    
    j <- j+1
    
    set.seed(100)
    rsq_bootstrap <- c()
    rsq_inc_bootstrap <- c()
    ebb_test_bootstrap <- ebb_test[vkood %in% ind]
    
    for(i in 1:1000){
      
      print(c(j, i))
      ebb_test_bootstrap_i <- ebb_test_bootstrap[sample(1:nrow(ebb_test_bootstrap), 
                                                        nrow(ebb_test_bootstrap), replace = T), ]
      lm_xx <- lm("EduYears_adj2 ~ PGI_EA", data = ebb_test_bootstrap_i)
      rsq_bootstrap <- c(rsq_bootstrap, summary(lm_xx)$r.squared)
      
      lm_0 <- lm(paste0("EduYears ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
                 data = ebb_test_bootstrap_i)
      lm_0prs <- lm(paste0("EduYears ~ PGI_EA + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
                    data = ebb_test_bootstrap_i)
      r2_inc <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
      rsq_inc_bootstrap <- c(rsq_inc_bootstrap, r2_inc)
    }
    
    rsq_top95 <- c(rsq_top95, quantile(rsq_bootstrap, probs = 0.975))
    rsq_bottom95 <- c(rsq_bottom95, quantile(rsq_bootstrap, probs = 0.025))
    rsq_inc_top95 <- c(rsq_inc_top95, quantile(rsq_inc_bootstrap, probs = 0.975))
    rsq_inc_bottom95 <- c(rsq_inc_bottom95, quantile(rsq_inc_bootstrap, probs = 0.025))
    
    dt <- cbind(dt, rsq_bootstrap, rsq_inc_bootstrap)
    
    lm_xx <- lm("EduYears_adj2 ~ PGI_EA", data = ebb_test_bootstrap)
    rsq <- c(rsq, summary(lm_xx)$r.squared)
    rsq_all <- c(rsq_all, summary(lm_xx)$r.squared) # added to calculate a table with main results
    
    lm_0 <- lm(paste0("EduYears ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
               data = ebb_test_bootstrap)
    lm_0prs <- lm(paste0("EduYears ~ PGI_EA + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:pc_n, sep = "", collapse = " + ")),
                  data = ebb_test_bootstrap)
    r2_inc <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
    rsq_inc <- c(rsq_inc, r2_inc)
    rsq_all <- c(rsq_all, r2_inc) # added to calculate a table with main results
    
    N <- c(N, nrow(ebb_test_bootstrap))
    n <- c(n, rep(nrow(ebb_test_bootstrap), 2)) # added to calculate a table with main results
    
  }
  
  colnames(dt) <- c("all_r2", "all_r2_inc", 
                    "s_15_r2", "s_15_r2_inc", "ps_15_r2", "ps_15_r2_inc", 
                    "s_10_r2", "s_10_r2_inc", "ps_10_r2", "ps_10_r2_inc",
                    "all_est_r2", "all_est_r2_inc", 
                    "s_15_est_r2", "s_15_est_r2_inc", "ps_15_est_r2", "ps_15_est_r2_inc", 
                    "s_10_est_r2", "s_10_est_r2_inc", "ps_10_est_r2", "ps_10_est_r2_inc")
  
  write.table(dt, paste0("~/EA_heritability/results/revision/r2_adj_", pheno, "_origSample_newPGS_1000.tsv"),
              row.names = F, quote = F, sep = "\t")
  
  cohort <- c("all_r2", "all_r2_inc", 
              "s_15_r2", "s_15_r2_inc", "ps_15_r2", "ps_15_r2_inc", 
              "s_10_r2", "s_10_r2_inc", "ps_10_r2", "ps_10_r2_inc",
              "all_est_r2", "all_est_r2_inc", 
              "s_15_est_r2", "s_15_est_r2_inc", "ps_15_est_r2", "ps_15_est_r2_inc", 
              "s_10_est_r2", "s_10_est_r2_inc", "ps_10_est_r2", "ps_10_est_r2_inc") # added to calculate a table with main results
  r2 <- rep(c("r2", "r2_inc"), 10) # added to calculate a table with main results
  
  
  name <- c("soviet 15", "post-soviet 15", "soviet 10", "post-soviet 10", "soviet 25", "transition")
  
  ################################################
  # added to calculate a table with main results #
  name <- c("All", "All", "S15", "S15", "PS15", "PS15", "S10", "S10", "PS10", "PS10",
            "All", "All", "S15", "S15", "PS15", "PS15", "S10", "S10", "PS10", "PS10")
  r2_res <- rbind(name, r2, round(rsq_all, 10), n)
  colnames(r2_res) <- cohort
  write.table(r2_res, paste0("~/EA_heritability/results/revision/r2_adj_", pheno, "_origSample_newPGS.tsv"),
              row.names = F, quote = F, sep = "\t")
  ################################################
  
}




codes <- fread("~/EBB_project/phenotypes/svcodes.tsv")
ebb <- merge(fread("~/EBB_project/phenotypes/query1.tsv"), #, encoding = "UTF-8"
             fread("~/EBB_project/phenotypes/query2.tsv"), #, encoding = "UTF-8"
             by = "Person skood")
ebb <- merge(ebb, fread("~/EA_heritability/data/query_OccupationStatus.tsv"), all = TRUE)
ebb <- merge(codes, ebb, by="Person skood")
ebb[, lastAge := 2022 - `Person birthYear`]


selected_columns <- c("Sample vkood", "Person skood",
                      "Person birthYear", "Person gender code",
                      "PersonLocation birthCounty name",
                      "PersonLocation residencyCounty name",
                      "PersonPortrait nationality name",
                      "lastAge", "Person ageAtAgreement",
                      "CONCATSTR(Education highestEducationLevel code)", 
                      "CONCATSTR(Work currentOccupation code)", "CONCATSTR(Work mainOccupation code)")

ebb <- ebb[, ..selected_columns]
colnames(ebb) <- c("vkood", "skood", "YoB", "Sex", "PoB", "PoR", "Nat", "Age", "AgeAtAgr", "EA_answerset", "curOcc", "mainOcc")

ebb[, Age2 := Age^2]
ebb[, SxA := Age*Sex]
ebb[, YoA := YoB+AgeAtAgr]

pca <- fread("~/EA_heritability/data/EstBB_additional_files_2022-09-14.eigenvec")
colnames(pca)[1] <- "vkood"
ebb <- merge(ebb, pca, by="vkood")


ebb2 <- fread("~/EBB_project/phenotypes/query_bmi_edu.tsv")
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA_portrait")
ebb <- merge(ebb, ebb2, by="skood")

origSample <- read_excel("~/EA_heritability/data/origSamples/2024-06-18_KRimfeld/I16_2024-06-18_15506GD_VK.xlsx")
ebb <- ebb[vkood %in% origSample$Vcode1, ]
ebb <- ebb[YoA-YoB >= 25, ]

ebb[, EA_a_order := transformSREAtoEduYears(EA_answerset, to="order")]
ebb[, EA_a_years := transformSREAtoEduYears(EA_answerset, to="years")]
ebb[, EA_a_binary := transformSREAtoEduYears(EA_answerset, to="binary")]
ebb[, EA_a_norm := transformSREAtoEduYears(EA_answerset, to="normal")]
ebb[, EA_p_order := transformEAtoEduYears(EA_portrait, to="order")]
ebb[, EA_p_years := transformEAtoEduYears(EA_portrait, to="years")]
ebb[, EA_p_binary := transformEAtoEduYears(EA_portrait, to="binary")]
ebb[, EA_p_norm := transformEAtoEduYears(EA_portrait, to="normal")]

ebb[, curOcc := sapply(curOcc, getOccStat)]
ebb[, mainOcc := sapply(mainOcc, getOccStat)]
ebb[, OS := rowMeans(.SD, na.rm = T), .SDcols = c("curOcc", "mainOcc")]
ebb[is.na(OS), OS := NA]



pgi_ea <-  fread("/gpfs/space/GI/GV/Projects/PGI_repository_v2/SSGAC_PGI_Repository_v2_EstBB.txt", 
                 select = c("IID", "PGI_EA"))
colnames(pgi_ea)[1] <- "vkood"
ebb <- merge(ebb, pgi_ea, by="vkood")



ind_all <- ebb[, vkood]
ind15_s <- ebb[YoB<1976, vkood]
ind15_ps <- ebb[YoB>=1976, vkood]
ind10_s <- ebb[YoB<1981, vkood]
ind10_ps <- ebb[YoB>=1981, vkood]
ind_all_est <- ebb[Nat == "Eestlane", vkood]
ind15_s_est <- ebb[Nat == "Eestlane" & YoB<1976, vkood]
ind15_ps_est <- ebb[Nat == "Eestlane" & YoB>=1976, vkood]
ind10_s_est <- ebb[Nat == "Eestlane" & YoB<1981, vkood]
ind10_ps_est <- ebb[Nat == "Eestlane" & YoB>=1981, vkood]



ind_list <- list(ind_all, ind15_s, ind15_ps, ind10_s, ind10_ps, ind_all_est, ind15_s_est, ind15_ps_est, ind10_s_est, ind10_ps_est) # , ind_tr_s, ind_tr_tr

ebb_test <- ebb

phenotypes <- c("EA_p_years", "EA_a_order", "OS")
for(pheno in phenotypes){
  
  colnames(ebb_test)[which(colnames(ebb_test)==pheno)] <- "EduYears"
  analyzeR2(ebb_test, ind_list=ind_list, pheno=pheno)
  colnames(ebb_test)[which(colnames(ebb_test)=="EduYears")] <- pheno
  
}

