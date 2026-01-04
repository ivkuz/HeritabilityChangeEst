
############################################################
# Fig. 7 ###################################################
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


# Make palette
# Red and blue shades (darker to lighter or vice versa)
red_shades  <- c("#FF9999", "#FF3333")
blue_shades <- c("#ADD8E6", "#5094CD")
r2[, fill_group := interaction(Era, weights, sep = "_")]
# Assign names matching fill_group levels
custom_colors <- c(
  # Red: phase1 and phase2 soviet
  "phase1\nsoviet_original"       = red_shades[1],
  "phase1\nsoviet_by group"       = red_shades[1],
  "phase1\nsoviet_overall"        = red_shades[1],
  
  "phase2\nsoviet_original"       = red_shades[2],
  "phase2\nsoviet_by group"       = red_shades[2],
  "phase2\nsoviet_overall"        = red_shades[2],
  
  # Blue: post-soviet
  "phase1\npost-soviet_original"  = blue_shades[1],
  "phase1\npost-soviet_by group"  = blue_shades[1],
  "phase1\npost-soviet_overall"   = blue_shades[1],
  
  "phase2\npost-soviet_original"  = blue_shades[2],
  "phase2\npost-soviet_by group"  = blue_shades[2],
  "phase2\npost-soviet_overall"   = blue_shades[2]
)


pl_list <- list()
for(tr in c("EA", "OS", "Height", "BMI")){
  
  pl <- ggplot(r2[trait == tr, ], aes(x = weights, y = r2, fill = fill_group)) + 
    geom_bar(stat='identity', position = position_dodge(width = 0.8), width=.6, color = "black") +
    geom_errorbar(aes(ymin =  ci_low, ymax = ci_up), width = 0.3,
                  position = position_dodge(width = 0.8, preserve = "single")) +
    theme_bw() + theme(text = element_text(size=10),
                       panel.grid.major.x = element_blank(),
                       axis.title.x=element_blank(),
                       legend.position="none") +
    scale_fill_manual(values = custom_colors) +
    scale_x_discrete(labels = c("No weighting", "Weigting inside\nsex-age groups", "Weigting in\noverall sample")) + labs(y = bquote(R^2), fill = "Weighting") + ggtitle(tr)
  
  pl_list <- append(pl_list, list(pl))
  
}


empty_plot <- ggplot() + theme_void()
pl_list <- append(pl_list, list(empty_plot))



# Build a legend
legend_labels <- c("wave 1\nSoviet", "wave 1\npost-Soviet", "wave 2\nSoviet", "wave 2\npost-Soviet")
legend_df <- data.frame(
  weight_label = legend_labels,
  xmin = rep(0, times = 4),
  xmax = rep(0.6, times = 4),
  ymin = 3:0+0.1,
  ymax = 4:1-0.1,
  fill = c(rbind(red_shades, blue_shades))
)
# add a new data frame for the labels, positioned just right of rectangles
label_df <- data.frame(
  weight_label = legend_labels,
  x = 0.7,  # just right of rectangles
  y = 3:0 + 0.5  # centered vertically in each row
)
legend_title <- data.frame(x = 0, y = 4.2, label = "Wave and era")
leg <- ggplot() +
  geom_rect(data = legend_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "black") +
  geom_text(data = label_df,
            aes(x = x, y = y, label = weight_label),
            hjust = 0, size=3, lineheight = 0.8) +
  geom_text(data = legend_title,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 0, size = 3.5) +  # Title text on top-left
  scale_fill_identity() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 3)) +  # wider space for text
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
  theme_void() +
  theme(plot.margin = margin(5, 20, 5, 5))  # bigger right margin

pl_list <- append(pl_list, list(leg))


pdf("~/EA_heritability/figures/paper/figure7.pdf", width = 5.5, height = 6)

print(
  grid.arrange(
    grobs = pl_list,
    layout_matrix = matrix(c(rep(1:4, each=2), 5,5,5,6,6,5,5,5), ncol = 2, byrow = F),
    widths = c(1, .35)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.02, y = 0.73, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.02, y = 0.23, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()

