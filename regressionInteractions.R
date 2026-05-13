
####################################
# Linear regression - GxE andlysis #
####################################

library(data.table)
library(ggplot2)
library(reshape2)
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
ebb[, EA_p_binary := transformEAtoEduYears(EA_portrait, to="binary")]
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

pca_est <- fread("~/EBB_project/data_filtering/pca/pcs_EstBB_estonian")
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
ebb <- merge(ebb, pca_est, by="vkood")

prs <- fread("~/EA_heritability/data/EA4_excl_23andMe_EGCUT/PRS.chrALL.sscore")

colnames(prs) <- c("vkood", "PRS")
ebb <- merge(ebb, prs)
ebb[, PRS := scale(PRS)]


ebb <- ebb[AgeAtAgr > 25,]



# # 15 years
# ind_01 <- fread("~/EA_heritability/gcta/data/ps_list_unrel.tsv", header=F)
# ind_02 <- fread("~/EA_heritability/gcta/data/s_list_unrel.tsv", header=F)
# ind_11 <- fread("~/EA_heritability/gcta/data/p1ps_list_unrel.tsv", header=F)
# ind_12 <- fread("~/EA_heritability/gcta/data/p1s_list_unrel.tsv", header=F)
# ind_21 <- fread("~/EA_heritability/gcta/data/p2ps_list_unrel.tsv", header=F)
# ind_22 <- fread("~/EA_heritability/gcta/data/p2s_list_unrel.tsv", header=F)

# 10 years
ind_01 <- fread("~/EA_heritability/gcta/data/ps_list_unrel_10.tsv", header=F)
ind_02 <- fread("~/EA_heritability/gcta/data/s_list_unrel_10.tsv", header=F)
ind_11 <- fread("~/EA_heritability/gcta/data/p1ps_list_unrel_10.tsv", header=F)
ind_12 <- fread("~/EA_heritability/gcta/data/p1s_list_unrel_10.tsv", header=F)
ind_21 <- fread("~/EA_heritability/gcta/data/p2ps_list_unrel_10.tsv", header=F)
ind_22 <- fread("~/EA_heritability/gcta/data/p2s_list_unrel_10.tsv", header=F)


ind_01 <- ind_01$V1
ind_02 <- ind_02$V1
ind_11 <- ind_11$V1
ind_12 <- ind_12$V1
ind_21 <- ind_21$V1
ind_22 <- ind_22$V1

ind_all <- fread("~/EBB_project/data_filtering/non_relatives_est_perfect.tsv")

ebb_test <- ebb
ebb_test[, EA_p_years := EduYears]
ebb_test <- ebb_test[vkood %in% ind_all$vkood, ]


ebb_test[, cohort := ifelse(vkood %in% ind_11, "p1ps", ifelse(vkood %in% ind_12, "p1s", ifelse(vkood %in% ind_21, "p2ps", ifelse(vkood %in% ind_22, "p2s", NA))))]
ebb_test[, cohort := factor(cohort, levels = c("p1s", "p1ps", "p2s", "p2ps"))]


#########################
# regression interactions
ebb_test[cohort %in% c("p1ps", "p1s"), phase := 0]
ebb_test[cohort %in% c("p2ps", "p2s"), phase := 1]
ebb_test[cohort %in% c("p1s", "p2s"), era := 0]
ebb_test[cohort %in% c("p1ps", "p2ps"), era := 1]
ebb_test[, YoB := scale(YoB)]
ebb_test[, Sex := Sex-1]



# ---- Helper to make formula string ----
make_formula <- function(outcome, pcs = 40, include_phase = TRUE, include_PRS_YoB = TRUE) {
  base_terms <- c(
    "PRS", "phase", "era", "Sex", "YoB",
    "I(Sex*YoB)", "I(YoB^2)", "I(phase*YoB)", "I(era*YoB)",
    "I(Sex*phase)", "I(Sex*era)", "I(phase*era)",
    "I(PRS*phase)", "I(PRS*era)", "I(PRS*Sex)", "I(PRS*YoB)"
  )
  
  # optionally drop all phase interactions
  if (!include_phase)
    base_terms <- setdiff(base_terms, grep("phase", base_terms, value = TRUE))
  
  # optionally drop PRS*YoB interaction
  if (!include_PRS_YoB)
    base_terms <- setdiff(base_terms, "I(PRS*YoB)")
  
  pcs_str <- paste0("PC", 1:pcs, collapse = " + ")
  paste(outcome, "~", paste(c(base_terms, pcs_str), collapse = " + "))
}

# ---- Fit model safely ----
fit_model <- function(formula_str, data, family = NULL) {
  if (is.null(family)) lm(as.formula(formula_str), data = data)
  else glm(as.formula(formula_str), data = data, family = family)
}

# ---- Extract clean coefficient table ----
extract_coefs <- function(model, pcs = 40) {
  tab <- as.data.table(summary(model)$coefficients)
  tab <- tab[2:(.N - pcs)] # drop intercept and PCs
  setnames(tab, c("beta", "se", "t", "p"))
  tab
}

# ---- Combine models (all, phase1, phase2) ----
run_all_models <- function(outcome, data, family = NULL, pcs = 40, 
                           drop_na_col = outcome, include_PRS_YoB = TRUE) {
  d_all <- data[!is.na(era) & !is.na(phase) & !is.na(get(drop_na_col))]
  d_p1 <- d_all[phase == 0]
  d_p2 <- d_all[phase == 1]
  
  m_all <- fit_model(make_formula(outcome, pcs, TRUE, include_PRS_YoB), d_all, family)
  m_p1  <- fit_model(make_formula(outcome, pcs, FALSE, include_PRS_YoB), d_p1, family)
  m_p2  <- fit_model(make_formula(outcome, pcs, FALSE, include_PRS_YoB), d_p2, family)
  
  list(
    all = extract_coefs(m_all, pcs),
    p1  = extract_coefs(m_p1, pcs),
    p2  = extract_coefs(m_p2, pcs)
  )
}

# ---- Format table for forest plot ----
format_forest_table <- function(tab_list, is_logistic = FALSE, include_PRS_YoB = TRUE) {
  # adjust row selection depending on whether PRS*YoB is present
  n_terms <- if (include_PRS_YoB) c(4, 3, 3) else c(3, 2, 2)
  
  tab_all <- rbind(
    tail(tab_list$all, n_terms[1]),
    tail(tab_list$p1, n_terms[2]),
    tail(tab_list$p2, n_terms[3])
  )
  
  # lways define all possible levels (to keep colors consistent)
  full_levels <- c("Sex", "YoB", "Era", "Wave")
  
  # Define which terms are present
  if (include_PRS_YoB) {
    interact_vals <- c("Wave", "Era", "Sex", "YoB", "Era", "Sex", "YoB", "Era", "Sex", "YoB")
    wave_vals <- c(rep("Waves 1 & 2", 4), rep("Wave 1", 3), rep("Wave 2", 3))
  } else {
    interact_vals <- c("Wave", "Era", "Sex", "Era", "Sex", "Era", "Sex")
    wave_vals <- c(rep("Waves 1 & 2", 3), rep("Wave 1", 2), rep("Wave 2", 2))
  }
  
  tab_all[, `:=`(
    interact = factor(interact_vals, levels = full_levels),  # fixed
    wave = factor(wave_vals, levels = c("Waves 1 & 2", "Wave 1", "Wave 2")),
    sign = p < 0.05
  )]
  
  # Compute CI or OR
  if (is_logistic) {
    tab_all[, `:=`(
      OR = exp(beta),
      lower = exp(beta - 1.96 * se),
      upper = exp(beta + 1.96 * se)
    )]
  } else {
    tab_all[, `:=`(
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se
    )]
  }
  
  tab_all[]
}
# ---- Generic forest plot ----
make_forest_plot <- function(tab, title, x_label = NULL) {
  is_logistic <- "OR" %in% names(tab)
  
  xvar <- if (is_logistic) "OR" else "beta"
  xminvar <- "lower"
  xmaxvar <- "upper"
  
  if (is.null(x_label))
    x_label <- if (is_logistic) "Odds Ratio (95% CI)" else "Beta (95% CI)"
  
  ggplot(tab, aes(y = interact, color = interact, shape = sign)) +
    geom_vline(xintercept = if (is_logistic) 1 else 0,
               linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = .data[[xminvar]], xmax = .data[[xmaxvar]]),
                   height = 0.25, linewidth = 0.8) +
    geom_point(aes(x = .data[[xvar]]), size = 3, fill = "white") +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 19)) +
    facet_wrap(~ wave, scales = "free_y", ncol = 1, strip.position = "top") +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(size=12),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    labs(
      x = x_label,
      y = "Interaction factor",
      title = title
    ) + scale_color_manual(values = c(
      "Sex" = "#D55E00",      # vermilion
      "YoB" = "#E69F00",      # orange
      "Era" = "#56B4E9",      # sky blue
      "Wave" = "#009E73"      # bluish green
    ))
}


# Linear model for years of education
lm_years <- run_all_models("EA_p_years", ebb_test)
tab_lm_all <- format_forest_table(lm_years)
lm_plot <- make_forest_plot(tab_lm_all, "EA (years of education)", "Beta (95% CI)")

lm_years <- run_all_models("EA_p_years", ebb_test, include_PRS_YoB = FALSE)
tab_lm_all <- format_forest_table(lm_years, include_PRS_YoB = FALSE)
lm_plot_noYoB <- make_forest_plot(tab_lm_all, "", "Beta (95% CI)")

# Linear model for occupational status
lm_os <- run_all_models("OS", ebb_test)
tab_lm_OS_all <- format_forest_table(lm_os)
lm_OS_plot <- make_forest_plot(tab_lm_OS_all, "Occupational status", "Beta (95% CI)")

lm_os <- run_all_models("OS", ebb_test, include_PRS_YoB = FALSE)
tab_lm_OS_all <- format_forest_table(lm_os, include_PRS_YoB = FALSE)
lm_OS_plot_noYoB <- make_forest_plot(tab_lm_OS_all, "", "Beta (95% CI)")

plot_lm_list <- list(lm_plot, lm_plot_noYoB, lm_OS_plot, lm_OS_plot_noYoB)

# Logistic model for binary education
glm_binary <- run_all_models("EA_p_binary", ebb_test, family = binomial(link = "logit"))
tab_logm_all <- format_forest_table(glm_binary, is_logistic = TRUE)
logm_plot <- make_forest_plot(tab_logm_all, "EA (university degree)", "Odds Ratio (95% CI)")

glm_binary <- run_all_models("EA_p_binary", ebb_test, 
                             family = binomial(link = "logit"), include_PRS_YoB = FALSE)
tab_logm_all <- format_forest_table(glm_binary, is_logistic = TRUE, include_PRS_YoB = FALSE)
logm_plot_noYoB <- make_forest_plot(tab_logm_all, "", "Odds Ratio (95% CI)")

# Logistic model for binary OS
ebb_test[!is.na(OS), OS_binary := ifelse(OS < 7, 0, 1)]
glm_binary <- run_all_models("OS_binary", ebb_test, family = binomial(link = "logit"))
tab_logm_all <- format_forest_table(glm_binary, is_logistic = TRUE)
logm_OS_plot <- make_forest_plot(tab_logm_all, "Occupational status (3/4 vs. 1/2)", "Odds Ratio (95% CI)")

glm_binary <- run_all_models("OS_binary", ebb_test, 
                             family = binomial(link = "logit"), include_PRS_YoB = FALSE)
tab_logm_all <- format_forest_table(glm_binary, is_logistic = TRUE, include_PRS_YoB = FALSE)
logm_OS_plot_noYoB <- make_forest_plot(tab_logm_all, "", "Odds Ratio (95% CI)")

plot_logm_list <- list(logm_plot, logm_plot_noYoB, logm_OS_plot, logm_OS_plot_noYoB)


pdf("~/EA_heritability/figures/regres_interact_10.pdf", width=7, height=7)

print(
  grid.arrange(
    grobs = plot_lm_list,
    layout_matrix = matrix(1:4, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.52, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))

print(
  grid.arrange(
    grobs = plot_logm_list,
    layout_matrix = matrix(1:4, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.52, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()

