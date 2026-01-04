
################################################################
# Supplementary Fig. 1 (GCTA and binary EA) and 2 (cutoff 10) ##
# Heritability in the post-Soviet and Soviet groups and ########
# the groups additionally divided by the wave of participation #
################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


# Plotting h2 bar plot with error bars
plotH2 <- function(reml_res, names, errors = "CI", title = "", lim = NULL, step = 0.1){
  
  # Choose whether plot SE or 95% CI error bars
  if(errors %in% c("SE", "se")){
    i <-  1
  } else if(errors %in% c("CI", "ci")){
    i <- 1.96
  }
  
  # Select names of output files where h2 were stored
  n <- length(names)
  tab <- reml_res[name %in% names, ]
  tab[, name := factor(name, levels = names)]
  
  # Main plot
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
    ggtitle(title) + ylab(bquote(h^2))
  
  # Additional parameters: unified y axis
  if(!is.null(lim)){
    pl <- pl + scale_y_continuous(breaks = c(seq(0, -2, by = -step), seq(0, 2, by = step)), limits = lim)
  } else{
    pl <- pl + scale_y_continuous(breaks = c(seq(0, -2, by = -step), seq(0, 2, by = step)))
  }
  
  # Additional parameters: colors and x axis labels
  if(n == 2){
    pl <- pl + scale_fill_manual(values = c(rgb(254, 0, 0, maxColorValue = 254), 
                                            rgb(80, 148, 205, maxColorValue = 254))) +
      scale_x_discrete(labels = c("S", "PS"))
  }else if(n == 4){
    pl <- pl + scale_fill_manual(values = rep(c(rgb(254, 0, 0, maxColorValue = 254), 
                                                rgb(80, 148, 205, maxColorValue = 254)), 2)) +
      scale_x_discrete(labels = c("S-w1", "PS-w1", "S-w2", "PS-w2"))
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
reml_files = readLines("~/EA_heritability/gcta/data/ldak_out_files2.txt")

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
reml_names <- c("ldak_s_15.reml", "ldak_ps_10.reml")
pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "EA", lim = c(-0.2, 0.9), step = 0.3)
df12 <- compareH2(reml_res = reml_res, names = reml_names)

# OS
reml_names <- c("ldak_s_ldak_list_OS_unrel_10_OS.reml", "ldak_ps_ldak_list_OS_unrel_10_OS.reml")
pl14 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "OS", lim = c(-0.05, 1.1), step = 0.3)
df14 <- compareH2(reml_res = reml_res, names = reml_names)

# Height
reml_names <- c("ldak_s_ldak_list_unrel_10_logHeight.reml", "ldak_ps_ldak_list_unrel_10_logHeight.reml")
pl15 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Height", lim = c(-0.05, 1.1), step = 0.3)
df15 <- compareH2(reml_res = reml_res, names = reml_names)

# BMI
reml_names <- c("ldak_s_ldak_list_unrel_10_logBMI.reml", "ldak_ps_ldak_list_unrel_10_logBMI.reml")
pl16 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "BMI", lim = c(-0.05, 1.15), step = 0.3)
df16 <- compareH2(reml_res = reml_res, names = reml_names)


# Cutoff 10, by wave #
# EA
reml_names <- c("ldak_p1s_10.reml", "ldak_p1ps_10.reml", "ldak_p2s_10.reml", "ldak_p2ps_10.reml")
pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.2, 0.9), step = 0.3)
df22 <- compareH2(reml_res = reml_res, names = reml_names)

# OS
reml_names <- c("ldak_p1s_ldak_list_OS_unrel_10_OS.reml", "ldak_p1ps_ldak_list_OS_unrel_10_OS.reml",
                "ldak_p2s_ldak_list_OS_unrel_10_OS.reml", "ldak_p2ps_ldak_list_OS_unrel_10_OS.reml")
pl24 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.05, 1.1), step = 0.3)
df24 <- compareH2(reml_res = reml_res, names = reml_names)

# Height
reml_names <- c("ldak_p1s_ldak_list_unrel_10_logHeight.reml", "ldak_p1ps_ldak_list_unrel_10_logHeight.reml",
                "ldak_p2s_ldak_list_unrel_10_logHeight.reml", "ldak_p2ps_ldak_list_unrel_10_logHeight.reml")
pl25 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(0, 1.15), step = 0.3)
df25 <- compareH2(reml_res = reml_res, names = reml_names)

# BMI
reml_names <- c("ldak_p1s_ldak_list_unrel_10_logBMI.reml", "ldak_p1ps_ldak_list_unrel_10_logBMI.reml",
                "ldak_p2s_ldak_list_unrel_10_logBMI.reml", "ldak_p2ps_ldak_list_unrel_10_logBMI.reml")
pl26 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.05, 1.15), step = 0.3)
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
df_p <- rbind(df_p_cutoff15, df_p_cutoff10)
write.table(df_p, "~/EA_heritability/figures/paper/h2_main_pval.tsv",
            row.names = F, quote = F, sep = "\t")


# cutoff 10, collect plots
plots <- list(pl12, pl22, pl14, pl24, pl15, pl25, pl16, pl26)

# Make the plot 
pdf("~/EA_heritability/figures/paper/h2_main_10.pdf", width=5.5, height=6)

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
reml_res <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = ".hsq")
pl11 <- plotH2(reml_res = reml_res, names = reml_res$name, errors = "CI", title = "EA, GCTA", lim = c(0, 0.4), step = 0.1)
df11 <- compareH2(reml_res = reml_res, names = reml_res$name)

# by wave
reml_files <- c("p1s", "p1ps", "p2s", "p2ps")
reml_res <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = ".hsq")
pl21 <- plotH2(reml_res = reml_res, names = reml_res$name, errors = "CI", title = "", lim = c(0, 0.4), step = 0.1)
df21 <- compareH2(reml_res = reml_res, names = reml_res$name)



# binary EA, cutoff 15
reml_names <- c("ldak_s_ldak_list_unrel_EA_binary.reml.liab", "ldak_ps_ldak_list_unrel_EA_binary.reml.liab")
pl12 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Univ. degree, 15 yo", lim = c(-0.04, 0.5), step = 0.1)
df12 <- compareH2(reml_res = reml_res, names = reml_names)

# binary EA, cutoff 15, phases (waves)
reml_names <- c("ldak_p1s_ldak_list_unrel_EA_binary.reml.liab", "ldak_p1ps_ldak_list_unrel_EA_binary.reml.liab",
                "ldak_p2s_ldak_list_unrel_EA_binary.reml.liab", "ldak_p2ps_ldak_list_unrel_EA_binary.reml.liab")
pl22 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.04, 0.5), step = 0.1)
df22 <- compareH2(reml_res = reml_res, names = reml_names)

# binary EA, cutoff 10
reml_names <- c("ldak_s_ldak_list_unrel_10_EA_binary.reml.liab", "ldak_ps_ldak_list_unrel_10_EA_binary.reml.liab")
pl13 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Univ. degree, 10 yo", lim = c(-0.4, 1.4), step = 0.3)
df13 <- compareH2(reml_res = reml_res, names = reml_names)

# binary EA, cutoff 10, phases (waves)
reml_names <- c("ldak_p1s_ldak_list_unrel_10_EA_binary.reml.liab", "ldak_p1ps_ldak_list_unrel_10_EA_binary.reml.liab",
                "ldak_p2s_ldak_list_unrel_10_EA_binary.reml.liab", "ldak_p2ps_ldak_list_unrel_10_EA_binary.reml.liab")
pl23 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "", lim = c(-0.35, 1.35), step = 0.3)
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
write.table(df_p, "~/EA_heritability/figures/paper/h2_suppl_pval.tsv",
            row.names = F, quote = F, sep = "\t")

plots <- list(pl11, pl21, pl12, pl22, pl13, pl23)

# Make the plot
pdf("~/EA_heritability/figures/paper/h2_suppl.pdf", width=5.5, height=4.5)
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

