################################################################
# Plot Fig. 2abc, 3ab, 4 #######################################
# Heritability in the post-Soviet and Soviet groups and ########
# the groups additionally divided by the wave of participation #
################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggsignif)

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
    theme_bw() + theme(text = element_text(size=10),
                       title = element_text(size=7.5),
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
      scale_x_discrete(labels = c("Soviet", "post-Soviet"))
  }else if(n == 4){
    pl <- pl + scale_fill_manual(values = rep(c(rgb(254, 0, 0, maxColorValue = 254), 
                                                rgb(80, 148, 205, maxColorValue = 254)), 2)) +
      scale_x_discrete(labels = c("wave 1\nSoviet", "wave 1\npost-Soviet", 
                                  "wave 2\nSoviet", "wave 2\npost-Soviet"))
  }
  
  return(pl)
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


# For figure 2

# EA, 15
reml_names <- c("ldak_s_15.reml", "ldak_ps_15.reml")
pl1 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Educational Attainment", lim = c(0, 0.32))
df1 <- compareH2(reml_res = reml_res, names = reml_names)
pl2a <- pl1 + geom_signif(comparisons=list(c("ldak_s_15.reml", "ldak_ps_15.reml")),
                          annotations = paste0("p=", format(df1[1, 1], scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.295, tip_length = c(0.3, 0.9),
                          step_increase = 0.1)

# OS, 15
reml_names <- c("ldak_s_ldak_list_OS_unrel_15_OS.reml", "ldak_ps_ldak_list_OS_unrel_15_OS.reml")
pl2 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Occupational Status", lim = c(0, 0.2))
df2 <- compareH2(reml_res = reml_res, names = reml_names)
pl2b <- pl2 + geom_signif(comparisons=list(c("ldak_s_ldak_list_OS_unrel_15_OS.reml", "ldak_ps_ldak_list_OS_unrel_15_OS.reml")),
                          annotations = paste0("p=", format(df2[1, 1], scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.185, tip_length = c(0.4, 0.25),
                          step_increase = 0.1)


# binary EA, cutoff 15
reml_names <- c("ldak_s_ldak_list_unrel_EA_binary.reml.liab", "ldak_ps_ldak_list_unrel_EA_binary.reml.liab")
pl3 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "University degree", lim = c(0, 0.38))
df3 <- compareH2(reml_res = reml_res, names = reml_names)

pl2c <- pl3 + geom_signif(comparisons=list(c("ldak_s_ldak_list_unrel_EA_binary.reml.liab", "ldak_ps_ldak_list_unrel_EA_binary.reml.liab")),
                          annotations = paste0("p=", format(df3[1, 1], scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.35, tip_length = c(0.07, 0.47),
                          step_increase = 0.1)


saveRDS(list(pl2a, pl2b, pl2c), "~/EA_heritability/figures/paper/files_for_figures/fig2abc.rds")




# For figure 3

# Height
reml_names <- c("ldak_s_ldak_list_unrel_logHeight.reml", "ldak_ps_ldak_list_unrel_logHeight.reml")
pl1 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Height", lim = c(0, 0.7))
df1 <- compareH2(reml_res = reml_res, names = reml_names)
pl3a <- pl1 + geom_signif(comparisons=list(c("ldak_s_ldak_list_unrel_logHeight.reml", "ldak_ps_ldak_list_unrel_logHeight.reml")),
                          annotations = paste0("p=", format(df1[1, 1], scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.65, tip_length = c(0.1, 0.7),
                          step_increase = 0.1)

# BMI
reml_names <- c("ldak_s_ldak_list_unrel_logBMI.reml", "ldak_ps_ldak_list_unrel_logBMI.reml")
pl2 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "BMI", lim = c(0, 0.4))
df2 <- compareH2(reml_res = reml_res, names = reml_names)
pl3b <- pl2 + geom_signif(comparisons=list(c("ldak_s_ldak_list_unrel_logBMI.reml", "ldak_ps_ldak_list_unrel_logBMI.reml")),
                          annotations = paste0("p=", format(df2[1, 1], scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.375, tip_length = c(0.15, 0.75),
                          step_increase = 0.1)

saveRDS(list(pl3a, pl3b), "~/EA_heritability/figures/paper/files_for_figures/fig3ab.rds")



# For figure 4

# Cutoff 15, by wave #
# EA
reml_names <- c("ldak_p1s_15.reml", "ldak_p1ps_15.reml", "ldak_p2s_15.reml", "ldak_p2ps_15.reml")
pl1 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Educational Attainment")
df1 <- compareH2(reml_res = reml_res, names = reml_names)
pl4a <- pl1 + geom_signif(comparisons=list(c("ldak_p2s_15.reml", "ldak_p2ps_15.reml")),
                          annotations = paste0("p=", format(df1[1, 2], scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.3, tip_length = c(0.1, 0.2),
                          step_increase = 0.1)


# OS
reml_names <- c("ldak_p1s_ldak_list_OS_unrel_15_OS.reml", "ldak_p1ps_ldak_list_OS_unrel_15_OS.reml", 
                "ldak_p2s_ldak_list_OS_unrel_15_OS.reml", "ldak_p2ps_ldak_list_OS_unrel_15_OS.reml")
pl2 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Occupational Status", step = 0.1)
df2 <- compareH2(reml_res = reml_res, names = reml_names)
pl4b <- pl2


# Height
reml_names <- c("ldak_p1s_ldak_list_unrel_logHeight.reml", "ldak_p1ps_ldak_list_unrel_logHeight.reml",
                "ldak_p2s_ldak_list_unrel_logHeight.reml", "ldak_p2ps_ldak_list_unrel_logHeight.reml")
pl3 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "Height", lim = c(0, 0.92), step = 0.2)
df3 <- compareH2(reml_res = reml_res, names = reml_names)
pl4c <- pl3 + geom_signif(comparisons=list(c("ldak_p2s_ldak_list_unrel_logHeight.reml", "ldak_p2ps_ldak_list_unrel_logHeight.reml"),
                                           c("ldak_p1s_ldak_list_unrel_logHeight.reml", "ldak_p2s_ldak_list_unrel_logHeight.reml")),
                          annotations = paste0("p=", format(c(df3[1, 2], df3[1, 3]), scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.75, tip_length = c(0.25, 0.35, 0.8, 0.05),
                          step_increase = 0.25)

# BMI
reml_names <- c("ldak_p1s_ldak_list_unrel_logBMI.reml", "ldak_p1ps_ldak_list_unrel_logBMI.reml",
                "ldak_p2s_ldak_list_unrel_logBMI.reml", "ldak_p2ps_ldak_list_unrel_logBMI.reml")
pl4 <- plotH2(reml_res = reml_res, names = reml_names, errors = "CI", title = "BMI")
df4 <- compareH2(reml_res = reml_res, names = reml_names)
pl4d <- pl4 + geom_signif(comparisons=list(c("ldak_p2s_ldak_list_unrel_logBMI.reml", "ldak_p2ps_ldak_list_unrel_logBMI.reml")),
                          annotations = paste0("p=", format(df4[1, 2], scientific = TRUE, digits = 2)),
                          textsize=3, size=0.5, vjust = -0.3,
                          y_position = 0.4, tip_length = c(0.1, 0.22),
                          step_increase = 0.1)


# cutoff 10, collect plots
plots <- list(pl4a, pl4b, pl4c, pl4d)

# Make the plot 
pdf("~/EA_heritability/figures/paper/figure4.pdf", width=5.5, height=5)

print(
  grid.arrange(
    grobs = plots,
    layout_matrix = matrix(1:4, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.52, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()
