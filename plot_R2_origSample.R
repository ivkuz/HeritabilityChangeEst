
###########################################################################################
# SF 3 ####################################################################################
# Variance explained by PGS in the post-Soviet and Soviet groups in the original sample ###
###########################################################################################


library(data.table)
library(psychometric)
library(ggplot2)
library(grid)
library(gridExtra)

# Plot R2 bars with error bars
plotR2 <- function(r2_res, names, r2_var = "r2_inc", title = "", lim = NA, step = 0.05){
  
  n <- length(names)
  tab <- r2_res[name %in% names & var == r2_var, ]
  tab[, name := factor(name, levels = names)]
  
  pl <- ggplot(tab, aes(x=name, y=r2, fill = name)) +
    geom_bar(stat = "identity", width=.55, color="black") + 
    geom_errorbar(aes(ymin=bstr_0.025, 
                      ymax=bstr_0.975), 
                  width=.2, linewidth = 0.3) +
    theme_bw() + theme(text = element_text(size=12),
                       title = element_text(size=10),
                       panel.grid.major.x = element_blank(),
                       axis.title.x=element_blank(),
                       legend.position = "none") +
    ggtitle(title) + ylab(bquote(R^2)) + 
    scale_fill_manual(values = c(rgb(38, 65, 140, maxColorValue = 254), rep(c(rgb(254, 0, 0, maxColorValue = 254), 
                                                                                  rgb(80, 148, 205, maxColorValue = 254)), 2))) +
    scale_x_discrete(labels = c("All", "S10", "PS10", "S15", "PS15"))
  
  if(!is.na(lim)){
    pl <- pl + scale_y_continuous(breaks = seq(0, 0.5, by = step), limits = c(0, lim))
  } else{
    pl <- pl + scale_y_continuous(breaks = seq(0, 0.5, by = step))
  }
  
  return(pl)
}


# Extract R2 table with CIs from the input with bootstrap results
makeR2 <- function(r2_res, bootstrap, r2_var = "r2_inc", title = "", lim = NA, step = 0.05){
  
  bootstrap_result <- bootstrap[, lapply(.SD, function(x) list(
    median = median(x, na.rm = TRUE),
    quantile_0025 = quantile(x, 0.025, na.rm = TRUE),
    quantile_0975 = quantile(x, 0.975, na.rm = TRUE)
  ))]
  
  r2_res <- rbind(r2_res, bootstrap_result)
  
  r2_res <- as.data.table(t(r2_res))
  colnames(r2_res) <- c("name", "var", "r2", "N", "bstr_0.5", "bstr_0.025", "bstr_0.975")
  r2_res[, r2 := as.numeric(r2)]
  r2_res[, N := as.numeric(N)]
  r2_res[, bstr_0.5 := as.double(bstr_0.5)]
  r2_res[, bstr_0.025 := as.double(bstr_0.025)]
  r2_res[, bstr_0.975 := as.double(bstr_0.975)]
  r2_res[, FT_0.025 := mapply(function(r, n) CIr(r = r, n = n, level = 0.95)[1]^2, sqrt(r2), N)]
  r2_res[, FT_0.975 := mapply(function(r, n) CIr(r = r, n = n, level = 0.95)[2]^2, sqrt(r2), N)]
  
  pl1 <- plotR2(r2_res = r2_res, names = c("All", "S10", "PS10", "S15", "PS15"), r2_var = r2_var, title = title, lim = lim, step = step)

  pl_list <- list(pl1)
  return(pl_list)
  
}

# Calculate bootstrap p-values for the difference between R2
compareR2 <- function(r2_res, bootstrap){
  
  s_ps_10 <- bootstrap[, sum(s_10_r2 < ps_10_r2)]
  s_ps_10_inc <- bootstrap[, sum(s_10_r2_inc < ps_10_r2_inc)]
  s_ps_15 <- bootstrap[, sum(s_15_r2 < ps_15_r2)]
  s_ps_15_inc <- bootstrap[, sum(s_15_r2_inc < ps_15_r2_inc)]

  s_ps_10_est <- bootstrap[, sum(s_10_est_r2 < ps_10_est_r2)]
  s_ps_10_est_inc <- bootstrap[, sum(s_10_est_r2_inc < ps_10_est_r2_inc)]
  s_ps_15_est <- bootstrap[, sum(s_15_est_r2 < ps_15_est_r2)]
  s_ps_15_est_inc <- bootstrap[, sum(s_15_est_r2_inc < ps_15_est_r2_inc)]
  
  comparison_r2 <- c(s_ps_10, s_ps_15, s_ps_10_est, s_ps_15_est,
                     s_ps_10_inc, s_ps_15_inc, s_ps_10_est_inc, s_ps_15_est_inc)
  pval_r2 <- round(sapply(comparison_r2, function(x) min(x/1000, 1 - x/1000)*2), 3)
  names_r2 <- c("s_ps_10", "s_ps_15", "s_ps_10_est", "s_ps_15_est",
                "s_ps_10_inc", "s_ps_15_inc", "s_ps_10_est_inc", "s_ps_15_est_inc")
  res <- rbind(c("r2_inc", pval_r2[5:8]), c("r2", pval_r2[1:4]))
  colnames(res) <- c("r2", names_r2[1:4])
  
  return(res)
  
}


r2_res_v <- c("r2_adj_EA_p_years_origSample.tsv", "r2_adj_EA_a_order_origSample.tsv", "r2_adj_OS_origSample.tsv")
bootstrap_v <- c("r2_adj_EA_p_years_origSample_1000.tsv", "r2_adj_EA_a_order_origSample_1000.tsv", "r2_adj_OS_origSample_1000.tsv")
titles <- c("EA years", "EA categories", "OS")

plot_list_inc <- list()
plot_list <- list()
pval_dt <- data.table()
for(i in 1:3){

  # 1:10 are with all nationalities; 11:20 are with self-reported Estonians
  r2_res <- fread(paste0("~/EA_heritability/results/", r2_res_v[i]))
  bootstrap <- fread(paste0("~/EA_heritability/results/", bootstrap_v[i]))
  plot_list_inc <- append(plot_list_inc, makeR2(r2_res = r2_res[, 1:10], bootstrap = bootstrap[, 1:10], 
                                                r2_var = "r2_inc", title = titles[i]))
  plot_list <- append(plot_list, makeR2(r2_res = r2_res[, 1:10], bootstrap = bootstrap[, 1:10], 
                                        r2_var = "r2", title = titles[i]))
  plot_list_inc <- append(plot_list_inc, makeR2(r2_res = r2_res[, 11:20], bootstrap = bootstrap[, 11:20], 
                                                r2_var = "r2_inc", title = paste0(titles[i], ", Estonians")))
  plot_list <- append(plot_list, makeR2(r2_res = r2_res[, 11:20], bootstrap = bootstrap[, 11:20], 
                                        r2_var = "r2", title = paste0(titles[i], ", Estonians"))) #, 
  pval_dt <- rbind(pval_dt, compareR2(r2_res = r2_res, bootstrap = bootstrap))
  
}

pval_dt <- as.data.table(pval_dt)
pval_dt[, trait := rep(titles, each = 2)]
write.table(pval_dt, "~/EA_heritability/figures/paper/r2_origSample_pval.tsv",
            row.names = F, quote = F, sep = "\t")


pdf("~/EA_heritability/figures/paper/r2_origSample_1incr_2adj.pdf", width=5.5, height=6)

# Plot incremental R2
print(
  grid.arrange(
    grobs = plot_list_inc,
    layout_matrix = matrix(1:6, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.6466, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.52, y = 0.6466, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("e", x = 0.02, y = 0.3133, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("f", x = 0.52, y = 0.3133, gp = gpar(fontsize=14, fontface = "bold"))


# Plot R2 for pre-adjusted traits
print(
  grid.arrange(
    grobs = plot_list,
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




