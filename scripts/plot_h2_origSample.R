
############################################################################
# SF 1 #####################################################################
# Heritability in the post-Soviet and Soviet groups in the original sample #
############################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

# Plotting h2 bar plot with error bars
plotH2 <- function(reml_res, names, errors = "CI", title = ""){
  
  if(errors %in% c("SE", "se")){
    i <-  1
  } else if(errors %in% c("CI", "ci")){
    i <- 1.96
  }
  
  tab <- reml_res[name %in% names, ]
  tab[, name := factor(name, levels = names)]

  pl <- ggplot(tab, aes(x=name, y=h2)) +
    geom_bar(fill=c(rgb(38, 65, 140, maxColorValue = 254), rep(c(rgb(254, 0, 0, maxColorValue = 254), 
                                                                 rgb(80, 148, 205, maxColorValue = 254)), 2)),
             stat = "identity", width=.55, color="black") + 
    geom_errorbar(aes(ymin=(h2-i*se), 
                      ymax=(h2+i*se)), 
                  width=.2, linewidth = 0.3) + 
    theme_bw() + theme(text = element_text(size=12),
                       title = element_text(size=10),
                       panel.grid.major.x = element_blank(),
                       axis.title.x=element_blank()) +
    ggtitle(title) + scale_y_continuous(breaks = seq(-1, 4, by = 0.4)) +
    scale_x_discrete(labels = c("All", "S10", "PS10", "S15", "PS15")) +
    ylab(bquote(h^2))
  
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
  
  tab <- reml_res[name %in% names, ]
  tab <- tab[order(match(name, names))]
  
  cutoff10 <- pnorm(q = -abs(tab[2, h2] - tab[3, h2]), 
                    mean = 0, 
                    sd = sqrt(tab[2, se]^2 + tab[3, se]^2))*2
  cutoff15 <- pnorm(q = -abs(tab[3, h2] - tab[4, h2]), 
                    mean = 0, 
                    sd = sqrt(tab[3, se]^2 + tab[4, se]^2))*2
  
  p_df <- data.frame(cutoff10, cutoff15)
  colnames(p_df) = c("cutoff10", "cutoff15")
  
  return(p_df)
  
}


# Upload GCTA h2 results
reml_files <- c("all", "s10", "ps10", "s15", "ps15")

reml_res <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = "_norm_origSample_unrelEst.hsq")
pl11 <- plotH2(reml_res = reml_res, names = reml_res$name, errors = "CI", title = "GCTA, INT EA, Estonians")
df11 <- compareH2(reml_res = reml_res, names = reml_res$name)

reml_res <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = "_origSample_unrelEst.hsq")
pl12 <- plotH2(reml_res = reml_res, names = reml_res$name, errors = "CI", title = "GCTA, EA, Estonians")
df12 <- compareH2(reml_res = reml_res, names = reml_res$name)


reml_res <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = "_norm_origSample_unrel05.hsq")
pl21 <- plotH2(reml_res = reml_res, names = reml_res$name, errors = "CI", title = "GCTA, INT EA")
df21 <- compareH2(reml_res = reml_res, names = reml_res$name)

reml_res <- getH2table(reml_files, prefix = "~/EA_heritability/gcta/results/", suffix = "_origSample_unrel05.hsq")
pl22 <- plotH2(reml_res = reml_res, names = reml_res$name, errors = "CI", title = "GCTA, EA")
df22 <- compareH2(reml_res = reml_res, names = reml_res$name)


# Upload LDAK h2 results
reml_res <- data.frame()
reml_files <- c("all", "s10", "ps10", "s15", "ps15", "all_OS", "s10_OS", "ps10_OS", "s15_OS", "ps15_OS")
for( res_file in reml_files){
  file_path <- paste0("~/EA_heritability/gcta/results/ldak_", res_file, "_originalSample.reml")
  res = readLines(file_path)
  l <- unlist(strsplit(res[grepl("Her_K1", res)], " "))
  h2 <- l[2]
  se <- l[3]
  reml_res <- rbind(reml_res, c(res_file, h2, se))
}

colnames(reml_res) <- c("name", "h2", "se")
reml_res <- as.data.table(reml_res)
reml_res[, h2 := as.numeric(h2)]
reml_res[, se := as.numeric(se)]


names <- c("all", "s10", "ps10", "s15", "ps15")
pl32 <- plotH2(reml_res = reml_res, names = names, errors = "CI", title = "LDAK, EA, Estonians")
df32 <- compareH2(reml_res = reml_res, names = reml_res$name)

names <- c("all_OS", "s10_OS", "ps10_OS", "s15_OS", "ps15_OS")
pl31 <- plotH2(reml_res = reml_res, names = names, errors = "CI", title = "LDAK, OS, Estonians")
df31 <- compareH2(reml_res = reml_res, names = reml_res$name)


# Make and save p-value table
dfs <- list(df11, df12, df21, df22, df32, df31)
df_p <- data.frame()
for(df_p_tmp in dfs){
  df_p <- rbind(df_p, df_p_tmp)
}
df_p <- as.data.table(format(df_p, scientific = TRUE, digits = 2))
df_p[, Trait := c("GCTA, INT EA, Estonians", "GCTA, EA, Estonians",
                  "GCTA, INT EA", "GCTA, EA",
                  "LDAK, EA, Estonians", "LDAK, OS, Estonians")]
write.table(df_p, "~/EA_heritability/figures/paper/h2_origSample_pval.tsv",
            row.names = F, quote = F, sep = "\t")



# Make the plot
plots <- list(pl11, pl12, pl21, pl22, pl31, pl32)

pdf("~/EA_heritability/figures/paper/h2_origSample.pdf", width=5.5, height=6)
print(
  grid.arrange(
    grobs = plots,
    layout_matrix = matrix(1:6, ncol = 2, byrow = T)
  )
)


grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.6466, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.52, y = 0.6466, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("e", x = 0.02, y = 0.3133, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("f", x = 0.52, y = 0.3133, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()

