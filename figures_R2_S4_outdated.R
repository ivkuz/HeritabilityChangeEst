
##########################################################################
# Plot Fig. S4 ###########################################################
# Trait variance explained by PGS in the post-Soviet and Soviet groups ###
##########################################################################


library(data.table)
library(psychometric)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggsignif)


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
    ggtitle(title) + ylab(bquote(R^2))
  
  if(!is.na(lim)){
    pl <- pl + scale_y_continuous(breaks = seq(0, 0.5, by = step), limits = c(0, lim))
  } else{
    pl <- pl + scale_y_continuous(breaks = seq(0, 0.5, by = step))
  }
  
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
  
  pl1 <- plotR2(r2_res = r2_res, names = c("S", "PS"), r2_var = r2_var, title = title, lim = lim, step = step)
  pl2 <- plotR2(r2_res = r2_res, names = c("p1S", "p1PS", "p2S", "p2PS"), r2_var = r2_var, title = "", lim = lim, step = step)
  
  pl_list <- list(pl1, pl2)
  return(pl_list)
  
}

# Calculate bootstrap p-values for the difference between R2
compareR2 <- function(r2_res, bootstrap){
  
  s_ps <- bootstrap[, sum(s_r2 < ps_r2)]
  s_ps_inc <- bootstrap[, sum(s_r2_inc < ps_r2_inc)]
  
  p1_s_ps <- bootstrap[, sum(s_r2_p1 < ps_r2_p1)]
  p1_s_ps_inc <- bootstrap[, sum(s_r2_inc_p1 < ps_r2_inc_p1)]
  p2_s_ps <- bootstrap[, sum(s_r2_p2 < ps_r2_p2)]
  p2_s_ps_inc <- bootstrap[, sum(s_r2_inc_p2 < ps_r2_inc_p2)]
  
  p1p2_s <- bootstrap[, sum(s_r2_p1 < s_r2_p2)]
  p1p2_s_inc <- bootstrap[, sum(s_r2_inc_p1 < s_r2_inc_p2)]
  p1p2_ps <- bootstrap[, sum(ps_r2_p1 < ps_r2_p2)]
  p1p2_ps_inc <- bootstrap[, sum(ps_r2_inc_p1 < ps_r2_inc_p2)]
  
  comparison_r2 <- c(s_ps, p1_s_ps, p2_s_ps, p1p2_s, p1p2_ps,
                     s_ps_inc, p1_s_ps_inc, p2_s_ps_inc, p1p2_s_inc, p1p2_ps_inc)
  pval_r2 <- round(sapply(comparison_r2, function(x) min(x/1000, 1 - x/1000)*2), 3)
  names_r2 <- c("s_ps", "p1_s_ps", "p2_s_ps", "p1p2_s", "p1p2_ps",
                "s_ps_inc", "p1_s_ps_inc", "p2_s_ps_inc", "p1p2_s_inc", "p1p2_ps_inc")
  res <- rbind(c("r2_inc", pval_r2[6:10]), c("r2", pval_r2[1:5]))
  colnames(res) <- c("r2", names_r2[1:5])
  
  return(res)
  
}


# files with R2 values and bootstrap results
r2_res_v <- c("r2_adj_EA_p_years_15y.tsv", "r2_adj_OS_15y.tsv", 
              "r2_adj_log_Height_first_ukbb.tsv", "r2_adj_log_BMI_first_ukbb.tsv")
bootstrap_v <- c("r2_adj_EA_p_years_1000.tsv", "r2_adj_OS_OS_15y_1000.tsv", 
                 "r2_adj_log_Height_first_1000_ukbb.tsv", "r2_adj_log_BMI_first_1000_ukbb.tsv")

# parameters for plotting
titles <- c("EA", "OS", "Height", "BMI") 
lim_list_inc <- c(0.105, 0.10, 0.15, 0.13)
lim_list <- c(0.11, 0.10, 0.32, 0.13)
steps = c(0.05, 0.05, 0.1, 0.05)

# manage data
plot_list_inc <- list()
plot_list <- list()
pval_dt <- data.table()
for(i in 1:4){
  
  r2_res <- fread(paste0("~/EA_heritability/results/", r2_res_v[i]))
  bootstrap <- fread(paste0("~/EA_heritability/results/", bootstrap_v[i]))
  plot_list_inc <- append(plot_list_inc, makeR2(r2_res = r2_res, bootstrap = bootstrap, 
                                                r2_var = "r2_inc", title = titles[i], 
                                                lim = lim_list_inc[i], step = 0.05))
  plot_list <- append(plot_list, makeR2(r2_res = r2_res, bootstrap = bootstrap, 
                                        r2_var = "r2", title = titles[i], 
                                        lim = lim_list[i], step = steps[i]))
  pval_dt <- rbind(pval_dt, compareR2(r2_res = r2_res, bootstrap = bootstrap))
  
}

# save p-value table
pval_dt <- as.data.table(pval_dt)
pval_dt[, trait := rep(titles, each = 2)]
write.table(pval_dt, "~/EA_heritability/results/r2_pval.tsv",
            row.names = F, quote = F, sep = "\t")


# Plot Supplementary Fig. 4

pdf("~/EA_heritability/figures/paper/r2_main.pdf", width=5.5, height=6)

# # incremental R2
# print(
#   grid.arrange(
#     grobs = plot_list_inc,
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

# R2 for pre-adjusted traits
print(
  grid.arrange(
    grobs = plot_list,
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




