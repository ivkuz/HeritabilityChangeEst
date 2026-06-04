
################################################################
# Supplementary Fig. 1 (GCTA and binary EA) and 2 (cutoff 10) ##
# Heritability in the post-Soviet and Soviet groups and ########
# the groups additionally divided by the wave of participation #
# and by sex and settlement type (rural-town-city) #############
################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


# Plotting h2 bar plot with error bars
plotH2 <- function(reml_res, names, errors = "CI", title = "", lim = NULL, step = 0.1){
  
  n <- length(names)
  if(errors %in% c("SE", "se")){
    i <-  1
  } else if(errors %in% c("CI", "ci")){
    i <- 1.96
  }
  
  n <- length(names)
  tab <- reml_res[name %in% names, ]
  tab[, name := factor(name, levels = names)]
  
  pl <- ggplot(tab, aes(x=name, y=h2, fill = name)) +
    geom_bar(stat = "identity", width=.55, color="black") + 
    geom_errorbar(aes(ymin=(h2-i*se), 
                      ymax=(h2+i*se)), 
                  width=.2, linewidth = 0.3) +
    theme_bw() + theme(text = element_text(size=12),
                       title = element_text(size=10),
                       panel.grid.major.x = element_blank(),
                       axis.title.x=element_blank(),
                       legend.position = "none") +
    ggtitle(title) + # scale_y_continuous(breaks = seq(-1, 4, by = 0.1)) +
    ylab(bquote(h^2))
  
  if(!is.null(lim)){
    pl <- pl + scale_y_continuous(breaks = c(seq(0, -2, by = -step), seq(0, 3, by = step)), limits = lim)
  } else{
    pl <- pl + scale_y_continuous(breaks = c(seq(0, -2, by = -step), seq(0, 3, by = step)))
  }
  
  if(n == 2){
    pl <- pl + scale_fill_manual(values = c(rgb(254, 0, 0, maxColorValue = 254), 
                                            rgb(80, 148, 205, maxColorValue = 254))) +
      # scale_x_discrete(labels = c("S", "PS"))
      scale_x_discrete(labels = c("Soviet", "post-Soviet"))
  }else if(n == 4){
    pl <- pl + scale_fill_manual(values = rep(c(rgb(254, 0, 0, maxColorValue = 254), 
                                                rgb(80, 148, 205, maxColorValue = 254)), 2)) +
      # scale_x_discrete(labels = c("S-w1", "PS-w1", "S-w2", "PS-w2"))
      scale_x_discrete(labels = c("wave 1\nSoviet", "wave 1\npost-Soviet", "wave 2\nSoviet", "wave 2\npost-Soviet"))
  }
  
  return(pl)
}

# Extract h2 ans se from GCTA-GREML output files where h2 were stored
getH2table <- function(rem_files, prefix = "~/EA_heritability/gcta/results/", suffix = ".hsq"){
  
  reml_res <- data.frame()
  for( res_file in reml_files){
    file_path <- paste0(prefix, res_file, suffix)
    if (file.exists(file_path)){
      res = readLines(file_path)
      l <- unlist(strsplit(res[grepl("V(G)/Vp", res, fixed = TRUE)], "\t"))
      h2 <- l[2]
      se <- l[3]
      reml_res <- rbind(reml_res, c(res_file, h2, se))
      
    } else{
      print(paste0(res_file, " doesn't exist"))
    }
  }
  reml_res <- as.data.table(reml_res)
  colnames(reml_res) <- c("name", "h2", "se")
  reml_res[, h2 := as.numeric(h2)]
  reml_res[, se := as.numeric(se)]
  
}

# Calculate p-values for h2 differences
compareH2 <- function(reml_res, names){
  
  # Select names of output files where h2 were stored
  n <- length(names)
  tab <- reml_res[name %in% names, ]
  tab <- tab[order(match(name, names))]
  
  # compare S and PS if n = 2 (no waves)
  if(n == 2){
    p <- pnorm(q = -abs(tab[1, h2] - tab[2, h2]), 
               mean = 0, 
               sd = sqrt(tab[1, se]^2 + tab[2, se]^2))*2

    p_df <- data.frame(p)
    colnames(p_df) = c("s_ps")
    
  # compare S and PS by waves and waves themselves if n = 4 (P and PS in waves 1 and 2)
  }else if(n == 4){
    
    p1_s_ps <- pnorm(q = -abs(tab[1, h2] - tab[2, h2]), 
               mean = 0, 
               sd = sqrt(tab[1, se]^2 + tab[2, se]^2))*2
    p2_s_ps <- pnorm(q = -abs(tab[3, h2] - tab[4, h2]), 
                     mean = 0, 
                     sd = sqrt(tab[3, se]^2 + tab[4, se]^2))*2
    p1p2_s <- pnorm(q = -abs(tab[1, h2] - tab[3, h2]), 
                     mean = 0, 
                     sd = sqrt(tab[1, se]^2 + tab[3, se]^2))*2
    p1p2_ps <- pnorm(q = -abs(tab[2, h2] - tab[4, h2]), 
                     mean = 0, 
                     sd = sqrt(tab[2, se]^2 + tab[4, se]^2))*2

    p_df <- data.frame(p1_s_ps, p2_s_ps, p1p2_s, p1p2_ps)
    colnames(p_df) = c("p1_s_ps", "p2_s_ps", "p1p2_s", "p1p2_ps")
    
  }
  
  return(p_df)
  
}

########
# LDAK #
########

# list of files with LDAK results
reml_files = readLines("~/EA_heritability/gcta/data/ldak_out_files_all.txt")

# Get REML results from the files
reml_res <- data.frame()
for( reml_file in reml_files){
  file_path <- paste0("~/EA_heritability/gcta/results/", reml_file)
  if (file.exists(file_path)){
    res = readLines(file_path)
    l <- unlist(strsplit(res[grepl("Her_K1", res)], " "))
    h2 <- l[2]
    se <- l[3]
    reml_res <- rbind(reml_res, c(reml_file, h2, se))
    
  } else{
    print(paste0(res_file, " doesn't exist"))
  }
}

colnames(reml_res) <- c("name", "h2", "se")
reml_res <- as.data.table(reml_res)
reml_res[, h2 := as.numeric(h2)]
reml_res[, se := as.numeric(se)]


# #############
# # Cutoff 15 #
# #############
# 
# # Cutoff 15, both waves together #
# # EA
# reml_names <- c("ldak_s_15.reml", "ldak_ps_15.reml")
# pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "EA", lim = c(0, 0.4), step = 0.1)
# df12 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # OS
# reml_names <- c("ldak_s_ldak_list_OS_unrel_15_OS.reml", "ldak_ps_ldak_list_OS_unrel_15_OS.reml")
# pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS", lim = c(-0.1, 0.27), step = 0.1)
# df14 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Height
# reml_names <- c("ldak_s_ldak_list_unrel_logHeight.reml", "ldak_ps_ldak_list_unrel_logHeight.reml")
# pl15 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Height", lim = c(0, 0.85), step = 0.2)
# df15 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # BMI
# reml_names <- c("ldak_s_ldak_list_unrel_logBMI.reml", "ldak_ps_ldak_list_unrel_logBMI.reml")
# pl16 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "BMI", lim = c(0, 0.5), step = 0.1)
# df16 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Cutoff 15, by wave #
# # EA
# reml_names <- c("ldak_p1s_15.reml", "ldak_p1ps_15.reml", "ldak_p2s_15.reml", "ldak_p2ps_15.reml")
# pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.4), step = 0.1)
# df22 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # OS
# reml_names <- c("ldak_p1s_ldak_list_OS_unrel_15_OS.reml", "ldak_p1ps_ldak_list_OS_unrel_15_OS.reml",
#                 "ldak_p2s_ldak_list_OS_unrel_15_OS.reml", "ldak_p2ps_ldak_list_OS_unrel_15_OS.reml")
# pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.1, 0.27), step = 0.1)
# df24 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Height
# reml_names <- c("ldak_p1s_ldak_list_unrel_logHeight.reml", "ldak_p1ps_ldak_list_unrel_logHeight.reml",
#                 "ldak_p2s_ldak_list_unrel_logHeight.reml", "ldak_p2ps_ldak_list_unrel_logHeight.reml")
# pl25 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.85), step = 0.2)
# df25 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # BMI
# reml_names <- c("ldak_p1s_ldak_list_unrel_logBMI.reml", "ldak_p1ps_ldak_list_unrel_logBMI.reml",
#                 "ldak_p2s_ldak_list_unrel_logBMI.reml", "ldak_p2ps_ldak_list_unrel_logBMI.reml")
# pl26 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.5), step = 0.1)
# df26 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# 
# # cutoff 15, collect p-value results in a table
# df1 <- list(df12, df14, df15, df16)
# df2 <- list(df22, df24, df25, df26)
# df_p <- data.frame()
# for(i in 1:4){
#   df_p_tmp <- cbind(df1[[i]], df2[[i]])
#   df_p <- rbind(df_p, df_p_tmp)
# }
# df_p <- as.data.table(format(df_p, scientific = TRUE, digits = 2))
# df_p[, Trait := c("EA", "OS", "Height", "BMI")]
# df_p[, cutoff := 15]
# df_p_cutoff15 <- df_p
# 
# 
# # cutoff 15, collect plots
# plots <- list(pl12, pl22, pl14, pl24, pl15, pl25, pl16, pl26)
# 
# # Make the plot
# pdf("~/EA_heritability/figures/paper/h2_main.pdf", width=5.5, height=6)
# 
# print(
#   grid.arrange(
#     grobs = plots,
#     layout_matrix = matrix(1:8, ncol = 2, byrow = T),
#     widths = c(1, 1.5)
#   )
# )
# 
# grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("b", x = 0.42, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("c", x = 0.02, y = 0.73, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("d", x = 0.42, y = 0.73, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("e", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("f", x = 0.42, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("g", x = 0.02, y = 0.23, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("h", x = 0.42, y = 0.23, gp = gpar(fontsize=14, fontface = "bold"))
# 
# dev.off()





########################
# Supplementary Fig. 2 #
# Cutoff 10 ############
########################

# Cutoff 10, both waves together #
# EA
reml_names <- c("ldak_kin0.05_s_10.reml", "ldak_ps_10.reml")
pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "EA", step = 0.1) #, lim = c(-0.2, 0.9)
df12 <- compareH2(reml_res = reml_res, names = reml_names)

# OS
reml_names <- c("ldak_OS_kin0.05_s_10.reml", "ldak_OS_kin0.05_ps_10.reml")
pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS", step = 0.05) #, lim = c(-0.05, 1.1)
df14 <- compareH2(reml_res = reml_res, names = reml_names)

# Height
reml_names <- c("ldak_Height_kin0.05_s_10.reml", "ldak_Height_kin0.05_ps_10.reml")
pl15 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Height", step = 0.2) #, lim = c(-0.05, 1.1)
df15 <- compareH2(reml_res = reml_res, names = reml_names)

# BMI
reml_names <- c("ldak_BMI_kin0.05_s_10.reml", "ldak_BMI_kin0.05_ps_10.reml")
pl16 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "BMI", step = 0.1) #, lim = c(-0.05, 1.15)
df16 <- compareH2(reml_res = reml_res, names = reml_names)


# Cutoff 10, by wave #
# EA
reml_names <- c("ldak_kin0.05_p1s_10.reml", "ldak_kin0.05_p1ps_10.reml", "ldak_kin0.05_p2s_10.reml", "ldak_kin0.05_p2ps_10.reml")
pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.3) #, lim = c(-0.2, 0.9)
df22 <- compareH2(reml_res = reml_res, names = reml_names)

# OS
reml_names <- c("ldak_OS_kin0.05_p1s_10.reml", "ldak_OS_kin0.05_p1ps_10.reml", "ldak_OS_kin0.05_p2s_10.reml", "ldak_OS_kin0.05_p2ps_10.reml")
pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.3) #, lim = c(-0.05, 1.1)
df24 <- compareH2(reml_res = reml_res, names = reml_names)

# Height
reml_names <- c("ldak_Height_kin0.05_p1s_10.reml", "ldak_Height_kin0.05_p1ps_10.reml", "ldak_Height_kin0.05_p2s_10.reml", "ldak_Height_kin0.05_p2ps_10.reml")
pl25 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.3) #, lim = c(0, 1.15)
df25 <- compareH2(reml_res = reml_res, names = reml_names)

# BMI
reml_names <- c("ldak_BMI_kin0.05_p1s_10.reml", "ldak_BMI_kin0.05_p1ps_10.reml", "ldak_BMI_kin0.05_p2s_10.reml", "ldak_BMI_kin0.05_p2ps_10.reml")
pl26 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.3) #, lim = c(-0.05, 1.15)
df26 <- compareH2(reml_res = reml_res, names = reml_names)


# cutoff 10, collect p-value results in a table
df1 <- list(df12, df14, df15, df16)
df2 <- list(df22, df24, df25, df26)
df_p <- data.frame()
for(i in 1:4){
  df_p_tmp <- cbind(df1[[i]], df2[[i]])
  df_p <- rbind(df_p, df_p_tmp)
}
df_p <- as.data.table(format(df_p, scientific = TRUE, digits = 2))
df_p[, Trait := c("EA", "OS", "Height", "BMI")]
df_p[, cutoff := 10]
df_p_cutoff10 <- df_p

# Save p-value table for both cutoffs
# df_p <- rbind(df_p_cutoff15, df_p_cutoff10)
write.table(df_p_cutoff10, "~/EA_heritability/figures/paper/revision/h2_main_pval_10.tsv",
            row.names = F, quote = F, sep = "\t")


# cutoff 10, collect plots
plots <- list(pl12, pl22, pl14, pl24, pl15, pl25, pl16, pl26)

# Make the plot 
pdf("~/EA_heritability/figures/paper/revision/h2_main_10.pdf", width=5.5, height=6)

print(
  grid.arrange(
    grobs = plots,
    layout_matrix = matrix(1:8, ncol = 2, byrow = T),
    widths = c(1, 1.5)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.42, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.73, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.42, y = 0.73, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("e", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("f", x = 0.42, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("g", x = 0.02, y = 0.23, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("h", x = 0.42, y = 0.23, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()





########################
# Supplementary Fig. 1 #
########################

########
# GCTA #
########

# no waves
reml_files <- c("s", "ps")
reml_res_gcta <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = "_kin0.05.hsq")
pl11 <- plotH2(reml_res = reml_res_gcta, names = reml_res_gcta$name, errors = "CI", title = "EA, GCTA", step = 0.1) #, lim = c(0, 0.4)
df11 <- compareH2(reml_res = reml_res_gcta, names = reml_res_gcta$name)

# by wave
reml_files <- c("p1s", "p1ps", "p2s", "p2ps")
reml_res_gcta <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = "_kin0.05.hsq")
pl21 <- plotH2(reml_res = reml_res_gcta, names = reml_res_gcta$name, errors = "CI", title = "", step = 0.1) #, lim = c(0, 0.4)
df21 <- compareH2(reml_res = reml_res_gcta, names = reml_res_gcta$name)



# binary EA, cutoff 15
reml_names <- c("ldak_EA_binary_kin0.05_s_15.reml.liab", "ldak_EA_binary_kin0.05_ps_15.reml.liab")
pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Univ. degree, 15 yo", step = 0.1) #, lim = c(-0.04, 0.5)
df12 <- compareH2(reml_res = reml_res, names = reml_names)

# binary EA, cutoff 15, phases (waves)
reml_names <- c("ldak_EA_binary_kin0.05_p1s_15.reml.liab", "ldak_EA_binary_kin0.05_p1ps_15.reml.liab",
                "ldak_EA_binary_kin0.05_p2s_15.reml.liab", "ldak_EA_binary_kin0.05_p2ps_15.reml.liab")
pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.2) #, lim = c(-0.04, 0.5)
df22 <- compareH2(reml_res = reml_res, names = reml_names)

# binary EA, cutoff 10
reml_names <- c("ldak_EA_binary_kin0.05_s_10.reml.liab", "ldak_EA_binary_kin0.05_ps_10.reml.liab")
pl13 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Univ. degree, 10 yo", step = 0.1) #, lim = c(-0.4, 1.4)
df13 <- compareH2(reml_res = reml_res, names = reml_names)

# binary EA, cutoff 10, phases (waves)
reml_names <- c("ldak_EA_binary_kin0.05_p1s_10.reml.liab", "ldak_EA_binary_kin0.05_p1ps_10.reml.liab",
                "ldak_EA_binary_kin0.05_p2s_10.reml.liab", "ldak_EA_binary_kin0.05_p2ps_10.reml.liab")
pl23 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.5) #, lim = c(-0.35, 1.35)
df23 <- compareH2(reml_res = reml_res, names = reml_names)


# p-value table
df1 <- list(df11, df12, df13)
df2 <- list(df21, df22, df23)
df_p <- data.frame()
for(i in 1:3){
  df_p_tmp <- cbind(df1[[i]], df2[[i]])
  df_p <- rbind(df_p, df_p_tmp)
}
df_p <- as.data.table(format(df_p, scientific = TRUE, digits = 2))
df_p[, Trait := c("EA, GCTA", "EA, bin, 15", "EA, bin, 10")]
write.table(df_p, "~/EA_heritability/figures/paper/revision/h2_suppl_pval.tsv",
            row.names = F, quote = F, sep = "\t")

plots <- list(pl11, pl21, pl12, pl22, pl13, pl23)

# Make the plot
pdf("~/EA_heritability/figures/paper/revision/h2_suppl.pdf", width=5.5, height=4.5)
print(
  grid.arrange(
    grobs = plots,
    layout_matrix = matrix(1:6, ncol = 2, byrow = T),
    widths = c(1, 1.5)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.42, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.667, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.42, y = 0.667, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("e", x = 0.02, y = 0.333, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("f", x = 0.42, y = 0.333, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()




# By sex

# Male cutoff 15
reml_names <- c("ldak_kin0.05_s_male_15.reml", "ldak_kin0.05_ps_male_15.reml")
# pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "EA")
pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Men, cutoff 15", lim = c(0, 0.32), step = 0.1)
df12 <- compareH2(reml_res = reml_res, names = reml_names)

# Female cutoff 15
reml_names <- c("ldak_kin0.05_s_female_15.reml", "ldak_kin0.05_ps_female_15.reml")
# pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS")
pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Women, cutoff 15", lim = c(0, 0.32), step = 0.1)
df14 <- compareH2(reml_res = reml_res, names = reml_names)

# Male cutoff 10
reml_names <- c("ldak_kin0.05_s_male_10.reml", "ldak_kin0.05_ps_male_10.reml")
# pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "EA")
pl15 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Men, cutoff 10", lim = c(0, 0.32), step = 0.1)
df15 <- compareH2(reml_res = reml_res, names = reml_names)

# Female cutoff 10
reml_names <- c("ldak_kin0.05_s_female_10.reml", "ldak_kin0.05_ps_female_10.reml")
# pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS")
pl16 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Women, cutoff 10", lim = c(0, 0.32), step = 0.1)
df16 <- compareH2(reml_res = reml_res, names = reml_names)


# # Phases Male cutoff 15
# reml_names <- c("ldak_s_p1male_15.reml", "ldak_ps_p1male_15.reml", "ldak_s_p2_male_15.reml", "ldak_ps_p2_male_15.reml")
# pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.2)
# # pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.4), step = 0.1)
# df22 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Phases Female cutoff 15
# reml_names <- c("ldak_s_p1female_15.reml", "ldak_ps_p1female_15.reml", "ldak_s_p2_female_15.reml", "ldak_ps_p2_female_15.reml")
# pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.2)
# # pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.1, 0.27), step = 0.1)
# df24 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Phases Male cutoff 10
# reml_names <- c("ldak_s_p1male_10.reml", "ldak_ps_p1male_10.reml", "ldak_s_p2_male_10.reml", "ldak_ps_p2_male_10.reml")
# pl25 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.5)
# # pl25 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.85), step = 0.2)
# df25 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Phases Female cutoff 10
# reml_names <- c("ldak_s_p1female_10.reml", "ldak_ps_p1female_10.reml", "ldak_s_p2_female_10.reml", "ldak_ps_p2_female_10.reml")
# pl26 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.5)
# # pl26 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.5), step = 0.1)
# df26 <- compareH2(reml_res = reml_res, names = reml_names)


# table p-values
df1 <- list(df12, df14, df15, df16)
# df2 <- list(df22, df24, df25, df26)
df_p <- data.frame()
for(i in 1:4){
  # df_p_tmp <- cbind(df1[[i]], df2[[i]])
  df_p_tmp <- cbind(df1[[i]])
  df_p <- rbind(df_p, df1[[i]])
}
df_p <- as.data.table(format(df_p, scientific = TRUE, digits = 2))
df_p[, Sex := c("M", "W", "M", "W")]
df_p[, cutoff := c(15, 15, 10, 10)]
write.table(df_p, "~/EA_heritability/figures/paper/revision/h2_bysex_pval.tsv",
            row.names = F, quote = F, sep = "\t")


plots1 <- list(pl12, pl14, pl15, pl16)
# plots2_15 <- list(pl12, pl22, pl14, pl24)
# plots2_10 <- list(pl15, pl25, pl16, pl26)

pdf("/gpfs/helios/home/kuznetsi/EA_heritability/figures/paper/revision/h2_bysex.pdf", width=4.4, height=3)

print(
  grid.arrange(
    grobs = plots1,
    layout_matrix = matrix(1:4, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.52, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()

# pdf("/gpfs/helios/home/kuznetsi/EA_heritability/figures/paper/revision/h2_bysex_all_15_10.pdf", width=5.5, height=3)
# 
# print(
#   grid.arrange(
#     grobs = plots2_15,
#     layout_matrix = matrix(1:4, ncol = 2, byrow = T),
#     widths = c(1, 1.5)
#   )
# )
# 
# grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("b", x = 0.42, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("d", x = 0.42, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# 
# 
# print(
#   grid.arrange(
#     grobs = plots2_10,
#     layout_matrix = matrix(1:4, ncol = 2, byrow = T),
#     widths = c(1, 1.5)
#   )
# )
# 
# 
# grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("b", x = 0.42, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("d", x = 0.42, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# 
# dev.off()



# By settlement type

# Rural cutoff 15
reml_names <- c("ldak_kin0.05_s_rural_15.reml", "ldak_kin0.05_ps_rural_15.reml")
# pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "EA")
pl11 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Rural, cutoff 15", lim = c(0, 0.31), step = 0.1)
df11 <- compareH2(reml_res = reml_res, names = reml_names)

# Town cutoff 15
reml_names <- c("ldak_kin0.05_s_town_15.reml", "ldak_kin0.05_ps_town_15.reml")
# pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS")
pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Town, cutoff 15", lim = c(0, 0.31), step = 0.1)
df12 <- compareH2(reml_res = reml_res, names = reml_names)

# City cutoff 15
reml_names <- c("ldak_kin0.05_s_city_15.reml", "ldak_kin0.05_ps_city_15.reml")
# pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS")
pl13 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "City, cutoff 15", lim = c(0, 0.31), step = 0.1)
df13 <- compareH2(reml_res = reml_res, names = reml_names)


# Rural cutoff 10
reml_names <- c("ldak_kin0.05_s_rural_10.reml", "ldak_kin0.05_ps_rural_10.reml")
# pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "EA")
pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Rural, cutoff 10", lim = c(0, 0.39), step = 0.1)
df14 <- compareH2(reml_res = reml_res, names = reml_names)

# Town cutoff 10
reml_names <- c("ldak_kin0.05_s_town_10.reml", "ldak_kin0.05_ps_town_10.reml")
# pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS")
pl15 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Town, cutoff 10", lim = c(0, 0.39), step = 0.1)
df15 <- compareH2(reml_res = reml_res, names = reml_names)

# City cutoff 10
reml_names <- c("ldak_kin0.05_s_city_10.reml", "ldak_kin0.05_ps_city_10.reml")
# pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS")
pl16 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "City, cutoff 10", lim = c(0, 0.39), step = 0.1)
df16 <- compareH2(reml_res = reml_res, names = reml_names)


# # Phases Rural cutoff 15
# reml_names <- c("ldak_s_p1rural_15.reml", "ldak_ps_p1rural_15.reml", "ldak_s_p2rural_15.reml", "ldak_ps_p2rural_15.reml")
# pl21 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.2)
# # pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.45), step = 0.1)
# df21 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Phases Town cutoff 15
# reml_names <- c("ldak_s_p1town_15.reml", "ldak_ps_p1town_15.reml", "ldak_s_p2town_15.reml", "ldak_ps_p2town_15.reml")
# pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.5)
# # pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.1, 0.27), step = 0.1)
# df22 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Phases City cutoff 15
# reml_names <- c("ldak_s_p1city_15.reml", "ldak_ps_p1city_15.reml", "ldak_s_p2city_15.reml", "ldak_ps_p2city_15.reml")
# pl23 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", step = 0.2)
# # pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.1, 0.27), step = 0.1)
# df23 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# 
# # Phases Rural cutoff 10
# reml_names <- c("ldak_s_p1rural_10.reml", "ldak_ps_p1rural_10.reml", "ldak_s_p2rural_10.reml", "ldak_ps_p2rural_10.reml")
# pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "")
# # pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 0.45), step = 0.1)
# df24 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Phases Town cutoff 10
# reml_names <- c("ldak_s_p1town_10.reml", "ldak_ps_p1town_10.reml", "ldak_s_p2town_10.reml", "ldak_ps_p2town_10.reml")
# pl25 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "")
# # pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.1, 0.27), step = 0.1)
# df25 <- compareH2(reml_res = reml_res, names = reml_names)
# 
# # Phases City cutoff 10
# reml_names <- c("ldak_s_p1town_10.reml", "ldak_ps_p1town_10.reml", "ldak_s_p2town_10.reml", "ldak_ps_p2town_10.reml")
# pl26 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "")
# # pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.1, 0.27), step = 0.1)
# df26 <- compareH2(reml_res = reml_res, names = reml_names)


# table p-values
df1 <- list(df11, df12, df13)
# df2 <- list(df21, df22, df23)
df2 <- list(df14, df15, df16)
df_p <- data.frame()
for(i in 1:3){
  df_p_tmp <- cbind(df1[[i]], df2[[i]])
  df_p <- rbind(df_p, df_p_tmp)
}
df_p <- as.data.table(format(df_p, scientific = TRUE, digits = 2))
colnames(df_p) <- c("s_ps_15", "s_ps_10")
df_p[, Settlement := c("Rural", "Town", "Urban")]
write.table(df_p, "~/EA_heritability/figures/paper/revision/h2_rural-urban_pval.tsv",
            row.names = F, quote = F, sep = "\t")


# # cutoff 15
# df1 <- list(df11, df12, df13)
# df2 <- list(df21, df22, df23)
# df_p <- data.frame()
# for(i in 1:3){
#   df_p_tmp <- cbind(df1[[i]], df2[[i]])
#   df_p <- rbind(df_p, df_p_tmp)
# }
# df_p <- as.data.table(format(df_p, scientific = TRUE, digits = 2))
# df_p[, Settelment := c("Rural", "Town", "City")]
# df_p[, cutoff := c(15, 15, 15)]
# write.table(df_p, "~/EA_heritability/figures/paper/revision/h2_rural-urban_pval.tsv",
#             row.names = F, quote = F, sep = "\t")


# plots <- list(pl12, pl22, pl14, pl24, pl15, pl25, pl16, pl26)
# plots <- list(pl11, pl21, pl12, pl22, pl13, pl23)
plots <- list(pl11, pl14, pl12, pl15, pl13, pl16)

# pdf("/gpfs/helios/home/kuznetsi/EA_heritability/figures/BGA/h2_main.pdf", width=5.5, height=10)
pdf("/gpfs/helios/home/kuznetsi/EA_heritability/figures/paper/revision/h2_rural-urban.pdf", width=4.4, height=4.5)

print(
  grid.arrange(
    grobs = plots,
    layout_matrix = matrix(1:6, ncol = 2, byrow = T),
    widths = c(1, 1)
  )
)

grid.text("a", x = 0.02, y = 0.97, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.42, y = 0.97, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.64, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.42, y = 0.64, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("e", x = 0.02, y = 0.30, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("f", x = 0.42, y = 0.30, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()





