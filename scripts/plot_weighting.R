
############################################################
# Fig. 6 ###################################################
# Trait variance explained by PGS in weighted subsamples ###
############################################################



library(data.table)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)


# Calculate p-values for R2 differences
compareR2 <- function(r2_res, trait_n = 1){
  
  traits <- c("EA", "OS", "Height", "BMI")
  trait <- traits[trait_n]
  
  bootstrap <- data.table(r2_res$r_11$r2_distribution[,trait_n],
                          r2_res$r12$r2_distribution[,trait_n],
                          r2_res$r_21$r2_distribution[,trait_n],
                          r2_res$r_22$r2_distribution[,trait_n])
  colnames(bootstrap) <- c("p1ps", "p1s", "p2ps", "p2s")
  
  p1_s_ps <- bootstrap[, sum(p1s < p1ps)]
  p2_s_ps <- bootstrap[, sum(p2s < p2ps)]
  
  p1p2_s <- bootstrap[, sum(p1s < p2s)]
  p1p2_ps <- bootstrap[, sum(p1ps < p2ps)]
  
  comparison_r2 <- c(p1_s_ps, p2_s_ps, p1p2_s, p1p2_ps)
  pval_r2 <- round(sapply(comparison_r2, function(x) min(x/1000, 1 - x/1000)*2), 3)
  names_r2 <- c("p1_s_ps", "p2_s_ps", "p1p2_s", "p1p2_ps")
  names(pval_r2) <- names_r2
  
  return(pval_r2)
  
}


# calculate p-values
r2_list <- readRDS("~/EA_heritability/results/r2_15_weighted_ci_alltraits.RDS")
pval_res <- data.frame()
for(i in 1:4){
  pval <- compareR2(r2_res = r2_list$unadj, trait_n = i)
  pval_res <- rbind(pval_res, c("No weights", pval))
  pval <- compareR2(r2_res = r2_list$bygroup, trait_n = i)
  pval_res <- rbind(pval_res, c("Inside groups", pval))
  pval <- compareR2(r2_res = r2_list$all, trait_n = i)
  pval_res <- rbind(pval_res, c("Overall", pval))
}
pval_res <- cbind(rep(c("EA", "OS", "Height", "BMI"), each = 3), pval_res)
colnames(pval_res) <- c("Trait", "weighting", "p1_s_ps", "p2_s_ps", "p1p2_s", "p1p2_ps")
write.table(pval_res, "~/EA_heritability/figures/paper/r2_15_weighted_ci_alltraits_pval.tsv", 
            col.names = T, row.names = F, quote = F, sep = "\t")


# make plots
r2 <- fread("~/EA_heritability/results/r2_15_weighted_ci_alltraits.tsv")
r2[, Era := factor(rep(rep(c("phase1\npost-soviet", "phase1\nsoviet", "phase2\npost-soviet", "phase2\nsoviet"), each = 4), 3),
                   levels = c("phase1\nsoviet", "phase1\npost-soviet", "phase2\nsoviet", "phase2\npost-soviet"))]
r2[, weights := factor(weights, levels = c("original", "by group", "overall"))]

# Make main plots
pl_list <- list()
for(tr in c("EA", "OS", "Height", "BMI")){
  
  pl <- ggplot(r2[trait == tr, ], aes(x = Era, y = r2, fill = weights)) + 
    geom_bar(stat='identity', position = position_dodge(width = 0.8), width=.6, color = "black") +
    geom_errorbar(aes(ymin =  ci_low, ymax = ci_up), width = 0.3,
                  position = position_dodge(width = 0.8, preserve = "single")) +
    theme_bw() + theme(text = element_text(size=10),
                       panel.grid.major.x = element_blank(),
                       axis.title.x=element_blank(),
                       legend.position="none") +
    scale_fill_manual(values = c("#FFE699", "#FDAE6B", "#DA70D6"), labels = c("No weights", "Inside\nsex-age groups", "In overall sample")) +
    scale_x_discrete(labels = c("S-w1", "PS-w1", "S-w2", "PS-w2")) + labs(y = bquote(R^2), fill = "Weighting") + ggtitle(tr)
  
  pl_list <- append(pl_list, list(pl))
  
}

# Empty plots for the grid plot
empty_plot <- ggplot() + theme_void()
pl_list <- append(pl_list, list(empty_plot))

# Extract legend
pl0 <- ggplot(r2, aes(x = Era, y = r2, fill = weights)) + 
  geom_bar(stat='identity', position = position_dodge(width = 0.8), width=.6, color = "black") +
  scale_fill_manual(values = c("#FFE699", "#FDAE6B", "#DA70D6"), 
                    labels = c("No weights", "Inside\nsex-age groups", "In overall sample")) +
  theme(text = element_text(size=10))
leg <- get_legend(pl0)


pl_grobs <- lapply(pl_list, ggplotGrob)
pl_grobs <- append(pl_grobs, list(leg_plot))


pdf("~/EA_heritability/figures/paper/r2_15_weighted_ci_alltraits.pdf", width = 5.5, height = 6)

print(
  grid.arrange(
    grobs = pl_grobs,
    layout_matrix = matrix(c(1:4, 6,5,5,6), ncol = 2, byrow = F),
    widths = c(1, .3)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.02, y = 0.73, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.02, y = 0.23, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()

