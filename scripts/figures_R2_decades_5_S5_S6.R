
###################################################################
# Plot Fig. 5, S5 and S6 ##########################################
# Variance explained by PGS in the birth cohorts by half-decade ### 
###################################################################


library(data.table)
library(ggplot2)
library(reshape2)
library(psychometric)
library(stringr)
library(grid)
library(gridExtra)


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


r2Decades <- function(ebb_test, trait, age, prs, title, xl = "", ylims = NULL, step = NULL){
  
  # Format the data
  colnames(ebb_test)[which(colnames(ebb_test)==trait)] <- "Trait"
  colnames(ebb_test)[which(colnames(ebb_test)==age)] <- "Age"
  colnames(ebb_test)[which(colnames(ebb_test)==prs)] <- "PRS"
  
  ebb_test <- ebb_test[!is.na(Trait), ]
  
  # list for result plots
  pl_list <- list()
  
  # k is the bin size
  for(k in c(5,10)){
    
    # for each bin
    r2_decades <- data.frame()
    for(i in seq(1930, 1990, k)){
      
      # Calculate incremental R2 for wave 1 bins
      ebb_test_bin <- ebb_test[(cohort=="p1s" | cohort=="p1ps") & YoB > i & YoB <= i+k, ]
      n1 <- nrow(ebb_test_bin)
      
      if(n1 > 10){
        lm_0 <- lm(paste0("Trait ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                   data = ebb_test_bin)
        lm_0prs <- lm(paste0("Trait ~ PRS + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                      data = ebb_test_bin)
        c1 <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
        r2_decades <- rbind(r2_decades, c(i, paste0(i+1, "-", i+k), c1, n1, "phase1"))
        
      }
      
      # Calculate incremental R2 for wave 2 bins
      ebb_test_bin <- ebb_test[(cohort=="p2s" | cohort=="p2ps") & YoB > i & YoB <= i+k, ]
      n2 <- nrow(ebb_test_bin)
      
      if(n2 > 10){
        lm_0 <- lm(paste0("Trait ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                   data = ebb_test_bin)
        lm_0prs <- lm(paste0("Trait ~ PRS + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                      data = ebb_test_bin)
        c2 <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
        r2_decades <- rbind(r2_decades, c(i, paste0(i+1, "-", i+k), c2, n2, "phase2"))
        
      }
      
      # Calculate incremental R2 for both wave bins
      ebb_test_bin <- ebb_test[indep_set == TRUE & YoB > i & YoB <= i+k, ]
      n3 <- nrow(ebb_test_bin)
      
      if(n3 > 10){
        lm_0 <- lm(paste0("Trait ~ Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                   data = ebb_test_bin)
        lm_0prs <- lm(paste0("Trait ~ PRS + Sex + Age + Sex*Age + I(Age^2) + ", paste("PC", 1:40, sep = "", collapse = " + ")), 
                      data = ebb_test_bin)
        c3 <- summary(lm_0prs)$r.squared - summary(lm_0)$r.squared
        r2_decades <- rbind(r2_decades, c(i, paste0(i+1, "-", i+k), c3, n3, "all"))
        
      }
      
    }
    colnames(r2_decades) <- c("Birth", "YoB", "R2", "N", "phase")
    
    r2_decades <- data.table(r2_decades)
    r2_decades[, R2 := as.numeric(R2)]
    r2_decades[, N := as.numeric(N)]
    # Calculate R2 95% CIs using Fisher’s Z-transformation
    r2_decades[, FT_0.025 := mapply(function(r, n) CIr(r = r, n = n, level = 0.95)[1]^2, sqrt(R2), N)]
    r2_decades[, FT_0.975 := mapply(function(r, n) CIr(r = r, n = n, level = 0.95)[2]^2, sqrt(R2), N)]
    

    # Make plot with merged waves (phases)
    pl1 <- ggplot(r2_decades[N>=100 & phase == "all",], aes(x=Birth, y=R2, colour = Birth)) + 
      geom_point() + 
      geom_errorbar(aes(ymin=FT_0.025, 
                        ymax=FT_0.975), 
                    width=.2, linewidth = 0.3) +
      theme_bw() + ggtitle(title) + ylab(bquote(R^2)) +
      theme(text = element_text(size = 10),
            title = element_text(size = 7.5),
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            # axis.title.x=element_blank(),
            legend.position = "none")  + xlab(xl)
    if(!is.null(ylims)){pl1 <- pl1 + scale_y_continuous(breaks = seq(0, 2, by = step), limits = ylims)}
    
    if(k == 10){
      pl1 <- pl1 +
        scale_color_manual(values = c(rep(rgb(254, 0, 0, maxColorValue = 254), 4), "black",
                                      rep(rgb(80, 148, 205, maxColorValue = 254), 3))) +
        scale_x_discrete(breaks=r2_decades$Birth, labels=r2_decades$YoB)
    } else if( k == 5){
      pl1 <- pl1 +
        scale_color_manual(values = c(rep(rgb(254, 0, 0, maxColorValue = 254), 9),
                                      rep(rgb(80, 148, 205, maxColorValue = 254), 4))) +
        scale_x_discrete(labels = function(x) {
          # Match factor levels with labels from YoB, show every second one
          idx <- match(x, r2_decades$Birth)
          lbls <- r2_decades$YoB[idx]
          ifelse(seq_along(lbls) %% 2 == 1, lbls, "")
        })
    }
    
    # Make plot with separate waves (phases)
    pl2 <- ggplot(r2_decades[N>=100 & phase != "all",], aes(x=Birth, y=R2, group=phase, color=phase)) + 
      geom_point(position = position_dodge(w = 0.5)) + 
      geom_errorbar(aes(ymin=FT_0.025, 
                        ymax=FT_0.975), 
                    width=.2, linewidth = 0.3,
                    position = position_dodge(w = 0.5)) +
      theme_bw() + ggtitle(title) + ylab(bquote(R^2)) + 
      theme(text = element_text(size = 10),
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            # axis.title.x=element_blank(),
            legend.position = "none") +
      scale_color_manual(values=c("#ad8d49", "#38785e"))  + xlab(xl)
    
    if(k == 10){
      pl2 <- pl2 +
        scale_x_discrete(breaks=r2_decades$Birth, labels=r2_decades$YoB)
    } else if( k == 5){
      pl2 <- pl2 +
        scale_x_discrete(labels = function(x) {
          # Match factor levels with labels from YoB, show every second one
          idx <- match(x, r2_decades$Birth)
          lbls <- r2_decades$YoB[idx]
          ifelse(seq_along(lbls) %% 2 == 1, lbls, "")
        })
    }

    
    pl_list <- append(pl_list, list(ggplotGrob(pl1), ggplotGrob(pl2)))
    # pl_list <- append(pl_list, list(pl1, pl2))
    
  }
  
  
  return(pl_list)
  
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

# We analyse individuals older than 25 years at the time of recruitment
ebb_test <- ebb[AgeAtAgr > 25,]

# Keep unrelated subsample
ind_all <- fread("~/EBB_project/data_filtering/non_relatives_est_perfect.tsv")
ebb_test[vkood %in% ind_all$vkood, indep_set := TRUE]


# 15 years
ind_01 <- fread("~/EA_heritability/gcta/data/ps_list_unrel.tsv", header=F)
ind_02 <- fread("~/EA_heritability/gcta/data/s_list_unrel.tsv", header=F)
ind_11 <- fread("~/EA_heritability/gcta/data/p1ps_list_unrel.tsv", header=F)
ind_12 <- fread("~/EA_heritability/gcta/data/p1s_list_unrel.tsv", header=F)
ind_21 <- fread("~/EA_heritability/gcta/data/p2ps_list_unrel.tsv", header=F)
ind_22 <- fread("~/EA_heritability/gcta/data/p2s_list_unrel.tsv", header=F)

ind_01 <- ind_01$V1
ind_02 <- ind_02$V1
ind_11 <- ind_11$V1
ind_12 <- ind_12$V1
ind_21 <- ind_21$V1
ind_22 <- ind_22$V1

# Assign a group by era and wave for future sepatation
ebb_test[, cohort := ifelse(vkood %in% ind_11, "p1ps", ifelse(vkood %in% ind_12, "p1s", ifelse(vkood %in% ind_21, "p2ps", ifelse(vkood %in% ind_22, "p2s", NA))))]
ebb_test[, cohort := factor(cohort, levels = c("p1s", "p1ps", "p2s", "p2ps"))]


# Calculate R2 values
pl1 <- r2Decades(ebb_test = ebb_test, trait = "Height_first", age = "Age_first", prs = "PRS_Height", title = "Height")
pl2 <- r2Decades(ebb_test = ebb_test, trait = "BMI_first", age = "Age_first", prs = "PRS_BMI", title = "BMI")
pl3 <- r2Decades(ebb_test = ebb_test, trait = "EduYears", age = "Age", prs = "PRS_EA", title = "EA")
pl4 <- r2Decades(ebb_test = ebb_test, trait = "OS", age = "Age", prs = "PRS_EA", title = "OS")



# Make plots
# Including Supplementary Fig. 6
pdf("~/EA_heritability/figures/paper/r2_ages.pdf", width=5, height=5)

for(i in 1:4){

  # plot_list <- list(ggplotGrob(pl3[[i]]), ggplotGrob(pl4[[i]]),
  #                   ggplotGrob(pl1[[i]]), ggplotGrob(pl2[[i]]))
  plot_list <- list(pl3[[i]], pl4[[i]],
                    pl1[[i]], pl2[[i]])

  print(
    grid.arrange(
      grobs = plot_list,
      layout_matrix = matrix(1:4, ncol = 2, byrow = T)
    )
  )

  grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
  grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
  grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
  grid.text("d", x = 0.52, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))

}

dev.off()

# Make plots for Height and BMI separately
i <- 1
# plot_list <- list(ggplotGrob(pl1[[i]]), ggplotGrob(pl2[[i]]))
plot_list <- list(pl1[[i]], pl2[[i]])


# Supplementary figure 5
pdf("~/EA_heritability/figures/paper/r2_ages_height_bmi.pdf", width=5, height=2.5)

print(
  grid.arrange(
    grobs = plot_list,
    layout_matrix = matrix(1:2, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.96, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.96, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


# For figure 5
pl3 <- r2Decades(ebb_test = ebb_test, trait = "EduYears", age = "Age", prs = "PRS_EA", title = "", xl = "Decade of birth", ylims = c(0.023, 0.125), step=0.03)
pl4 <- r2Decades(ebb_test = ebb_test, trait = "OS", age = "Age", prs = "PRS_EA", title = "", xl = "Decade of birth", ylims = c(0.015, 0.115), step=0.03)

plot_list <- list(pl3[[1]], pl4[[1]])

saveRDS(plot_list, "~/EA_heritability/figures/paper/files_for_figures/fig5cf.rds")
