
####################################################
# Fig. 2d, 3cd  and SF 3 ###########################
# Heritability across different birth year cohorts #
####################################################

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

# Upload LDAK h2 results for decade bins
plots <- list()
traits <- c("EA", "EA_binary_liab", "OS", "logHeight", "logBMI")
trait_names <- c("Educational Attainment", "University degree", "Occupational Status", "Height", "BMI")
limits <- c(0.32, NA, NA, 0.7, 0.4)
traits_for_table <- c("EA", "EA_binary", "OS", "Height", "BMI")
results_tab <- data.frame()
for(i in 1:5){
  
  trait <- traits[i]
  age_res <- fread(paste0("~/EA_heritability/gcta/results/ldak_age_h2_", trait, ".tsv"))
  colnames(age_res) <- c("YoB", "h2", "se", "Size", "Mega_Intensity", "Int_SD")
  age_res$YoB <- factor(c("1986-1995", "1981-1990", "1976-1985", "1971-1980", "1966-1975", 
                          "1961-1970", "1956-1965", "1951-1960", "1946-1955", "1941-1950", "1991-2000"), 
                        levels = c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                                   "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"))
  
  age_res[YoB %in% c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", "1966-1975"), 
          Era := "Soviet"]
  age_res[YoB == "1971-1980", Era := "Transition"]
  age_res[YoB %in% c("1976-1985", "1981-1990", "1986-1995", "1991-2000"), 
          Era := "post-Soviet"]
  age_res$Era <- factor(age_res$Era, levels = c("Soviet", "Transition", "post-Soviet"))
  
  age_res_tab <- age_res[, Trait := traits_for_table[i]]
  age_res_tab <- age_res_tab[order(YoB), .(Trait, YoB, Era, h2, se)]
  results_tab <- rbind(results_tab, age_res_tab)
  
  if(trait == "OS"){
    age_res[YoB == "1991-2000", h2 := NA]
  }
  
  # # Make the LDAK h2 plots
  # pl <- ggplot(age_res,
  #              aes(x = YoB, y = h2)) +
  #   geom_point(aes(shape = YoB, color = YoB),
  #              show.legend = FALSE) +
  #   geom_errorbar(aes(ymin = (h2 - 1.96 * se),
  #                     ymax = (h2 + 1.96 * se),
  #                     linewidth = YoB, color = YoB),
  #                 width = 0.2) + # , linewidth = 0.5
  #   theme_bw() +
  #   theme(text = element_text(size = 10),
  #         panel.grid.major.x = element_blank(),
  #         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  #         axis.title.x=element_blank(),
  #         legend.position = "none") +
  #   scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
  #   scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
  #   scale_color_manual(values = c(rep(rgb(254, 0, 0, maxColorValue = 254), 6), rgb(38, 65, 140, maxColorValue = 254),
  #                                 rep(rgb(80, 148, 205, maxColorValue = 254), 4))) +
  #   ylab(bquote(h^2)) + ggtitle(trait_names[i]) + # xlab("Year of Birth") +
  #   scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, ""))
  
  # Make the LDAK h2 plots
  pl <- ggplot(age_res,
               aes(x = YoB, y = h2, color = Era)) +
    geom_point(aes(shape = YoB)) +
    geom_errorbar(aes(ymin = (h2 - 1.96 * se),
                      ymax = (h2 + 1.96 * se),
                      linewidth = YoB),
                  width = 0.2) + # , linewidth = 0.5
    theme_bw() +
    theme(text = element_text(size = 10),
          title = element_text(size=8),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
    scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
    scale_color_manual(values = c(rgb(254, 0, 0, maxColorValue = 254), rgb(38, 65, 140, maxColorValue = 254),
                                  rgb(80, 148, 205, maxColorValue = 254))) +
    ylab(bquote(h^2)) + xlab("Decade of birth") + ggtitle(trait_names[i]) + # xlab("Year of Birth") +
    scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
    guides(linewidth = "none", 
           shape = "none")
  
  if(!is.na(limits[i])){
    pl <- pl + ylim(c(0, limits[i]))
  }
  
  
  plots <- append(plots, list(pl))
  
}

# for the main figure
plots1 <- list(plots[[1]], plots[[3]], plots[[4]], plots[[5]])


saveRDS(plots[1], "~/EA_heritability/figures/paper/files_for_figures/fig2d.rds")
saveRDS(c(plots[4], plots[5]), "~/EA_heritability/figures/paper/files_for_figures/fig3cd.rds")

write.table(results_tab, "~/EA_heritability/figures/paper/revision/h2_estimates_YoBbins_alltraits_unrel_2nd_degree.tsv",
            row.names = F, quote = F, sep = "\t")


# Upload GCTA h2 results for decade bins
age_res <- fread("~/EA_heritability/gcta/results/age_h2.tsv")
colnames(age_res) <- c("YoB", "h2", "se")
age_res$YoB <- factor(c("1986-1995", "1941-1950", "1991-2000", "1981-1990", "1976-1985", 
                        "1971-1980", "1966-1975", "1961-1970", "1956-1965", "1951-1960", "1946-1955"), 
                      levels = c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                                 "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"))

age_res <- age_res[order(YoB), ]
write.table(age_res, "~/EA_heritability/figures/paper/revision/h2_estimates_YoBbins_EA_gcta_unrel_2nd_degree.tsv",
            row.names = F, quote = F, sep = "\t")

# Make the GCTA h2 plot
pl1 <- ggplot(age_res, 
              aes(x = YoB, y = h2)) + 
  geom_point(aes(shape = YoB, color = YoB)) + 
  geom_errorbar(aes(ymin = (h2 - 1.96 * se), 
                    ymax = (h2 + 1.96 * se),
                    linewidth = YoB, color = YoB),  
                width = 0.2) +
  theme_bw() + 
  theme(text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x=element_blank(),
        legend.position = "none") +
  scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
  scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
  scale_color_manual(values = c(rep(rgb(254, 0, 0, maxColorValue = 254), 6), rgb(38, 65, 140, maxColorValue = 254),
                                rep(rgb(80, 148, 205, maxColorValue = 254), 4))) +
  ylab(bquote(h^2)) + ggtitle("EA, GCTA") + # xlab("Year of Birth") +
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, ""))

# for the supplementary figure
plots2 <- list(pl1, plots[[2]])


# # Figure 2
# pdf("~/EA_heritability/figures/paper/h2_ages3.pdf", width=7/3*2, height = 5)
# 
# print(
#   grid.arrange(
#     grobs = plots1,
#     layout_matrix = matrix(1:4, ncol = 2, byrow = T)
#   )
# )
# 
# grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("b", x = 0.52, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# grid.text("d", x = 0.52, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
# 
# dev.off()


# SF 3
pdf("~/EA_heritability/figures/paper/h2_ages3_suppl.pdf", width=7/3*2, height = 2.5)

print(
  grid.arrange(
    grobs = plots2,
    layout_matrix = matrix(1:2, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.96, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.96, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()





############################################################
# Relatedness threshold 0.05 ###############################
############################################################

# Upload LDAK h2 results for decade bins
plots <- list()
traits <- c("EA", "EA_binary_liab", "OS", "Height", "BMI")
trait_names <- c("Educational Attainment", "University degree", "Occupational Status", "Height", "BMI")
limits <- c(0.32, NA, NA, 0.73, 0.4)
traits_for_table <- c("EA", "EA_binary", "OS", "Height", "BMI")
results_tab <- data.frame()
for(i in 1:5){
  
  trait <- traits[i]
  # age_res <- fread(paste0("~/EA_heritability/gcta/results/ldak_age_h2_kin0.05_", trait, ".tsv"))
  age_res <- fread(paste0("~/EA_heritability/gcta/results/ldak_age_h2_kin0.05_noAgeAtAgr_", trait, ".tsv"))
  colnames(age_res) <- c("YoB", "h2", "se", "Size", "Mega_Intensity", "Int_SD")
  age_res$YoB <- factor(c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                          "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"), 
                        levels = c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                                   "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"))
  
  age_res[YoB %in% c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", "1966-1975"), 
          Era := "Soviet"]
  age_res[YoB == "1971-1980", Era := "Transition"]
  age_res[YoB %in% c("1976-1985", "1981-1990", "1986-1995", "1991-2000"), 
          Era := "post-Soviet"]
  age_res$Era <- factor(age_res$Era, levels = c("Soviet", "Transition", "post-Soviet"))
  
  age_res_tab <- age_res[, Trait := traits_for_table[i]]
  age_res_tab <- age_res_tab[order(YoB), .(Trait, YoB, Era, h2, se)]
  results_tab <- rbind(results_tab, age_res_tab)
  
  if(trait == "OS"){
    age_res[YoB == "1991-2000", h2 := NA]
  }
  # age_res <- age_res[YoB != "1991-2000", ]
  
  
  # Make the LDAK h2 plots
  pl <- ggplot(age_res,
               aes(x = YoB, y = h2, color = Era)) +
    geom_point(aes(shape = YoB)) +
    geom_errorbar(aes(ymin = (h2 - 1.96 * se),
                      ymax = (h2 + 1.96 * se),
                      linewidth = YoB),
                  width = 0.2) + # , linewidth = 0.5
    theme_bw() +
    theme(text = element_text(size = 10),
          title = element_text(size=8),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
    scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
    scale_color_manual(values = c(rgb(254, 0, 0, maxColorValue = 254), rgb(38, 65, 140, maxColorValue = 254),
                                  rgb(80, 148, 205, maxColorValue = 254))) +
    ylab(bquote(h^2)) + xlab("Decade of birth")  + # xlab("Year of Birth") +
    scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
    guides(linewidth = "none", 
           shape = "none")
  
  if(!is.na(limits[i])){
    pl <- pl + ylim(c(0, limits[i]))
  }
  if(trait == "EA_binary_liab"){
    pl <- pl + theme(legend.position = "none")
  }
  if(trait %in% c("EA", "EA_binary_liab")){
    pl <- pl + ggtitle(trait_names[i])
  } else{
    pl <- pl + ggtitle("")
  }
  
  
  plots <- append(plots, list(pl))
  
}

# for the main figure
plots1 <- list(plots[[1]], plots[[3]], plots[[4]], plots[[5]])


# saveRDS(plots[1], "~/EA_heritability/figures/paper/revision/files_for_figures/fig2d.rds")
# saveRDS(c(plots[4], plots[5]), "~/EA_heritability/figures/paper/revision/files_for_figures/fig3cd.rds")
saveRDS(plots[1], "~/EA_heritability/figures/paper/revision/files_for_figures/fig2d_noAgeAtAgr.rds")
saveRDS(c(plots[4], plots[5]), "~/EA_heritability/figures/paper/revision/files_for_figures/fig3cd_noAgeAtAgr.rds")

# write.table(results_tab, "~/EA_heritability/figures/paper/revision/h2_estimates_YoBbins_alltraits_kin0.05.tsv",
#             row.names = F, quote = F, sep = "\t")
write.table(results_tab, "~/EA_heritability/figures/paper/revision/h2_estimates_YoBbins_alltraits_kin0.05_noAgeAtAgr.tsv",
            row.names = F, quote = F, sep = "\t")


# Upload GCTA h2 results for decade bins
age_res <- fread("~/EA_heritability/gcta/results/age_h2_kin0.05.tsv")
colnames(age_res) <- c("YoB", "h2", "se")
age_res$YoB <- factor(c("1941-1950", "1991-2000", "1986-1995", "1981-1990", "1976-1985", 
                        "1971-1980", "1966-1975", "1961-1970", "1956-1965", "1951-1960", "1946-1955"), 
                      levels = c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                                 "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"))

age_res <- age_res[order(YoB), ]
write.table(age_res, "~/EA_heritability/figures/paper/revision/h2_estimates_YoBbins_EA_gcta_kin0.05.tsv",
            row.names = F, quote = F, sep = "\t")

# age_res <- age_res[YoB != "1991-2000", ]

# Make the GCTA h2 plot
pl1 <- ggplot(age_res, 
              aes(x = YoB, y = h2)) + 
  geom_point(aes(shape = YoB, color = YoB)) + 
  geom_errorbar(aes(ymin = (h2 - 1.96 * se), 
                    ymax = (h2 + 1.96 * se),
                    linewidth = YoB, color = YoB),  
                width = 0.2) +
  theme_bw() + 
  theme(text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x=element_blank(),
        legend.position = "none") +
  scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
  scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
  scale_color_manual(values = c(rep(rgb(254, 0, 0, maxColorValue = 254), 6), rgb(38, 65, 140, maxColorValue = 254),
                                rep(rgb(80, 148, 205, maxColorValue = 254), 4))) +
  ylab(bquote(h^2)) + ggtitle("EA, GCTA") + # xlab("Year of Birth") +
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, ""))

# for the supplementary figure
plots2 <- list(pl1, plots[[2]])



# SF 3
pdf("~/EA_heritability/figures/paper/revision/h2_ages3_suppl_noAgeAtAgr.pdf", width=7/3*2, height = 2.5)

print(
  grid.arrange(
    grobs = plots2,
    layout_matrix = matrix(1:2, ncol = 2, byrow = T)
  )
)

grid.text("a", x = 0.02, y = 0.96, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.96, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


############################################################
# INT EA ###################################################
############################################################

trait <- "EA"
trait_name <- "Educational Attainment"

# Upload LDAK h2 results for decade bins
age_res <- fread("~/EA_heritability/gcta/results/ldak_age_h2_kin0.05_noAgeAtAgr_EA_INT.tsv")
colnames(age_res) <- c("YoB", "h2", "se", "Size", "Mega_Intensity", "Int_SD")
age_res$YoB <- factor(c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                        "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"), 
                      levels = c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                                 "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"))

age_res[YoB %in% c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", "1966-1975"), 
        Era := "Soviet"]
age_res[YoB == "1971-1980", Era := "Transition"]
age_res[YoB %in% c("1976-1985", "1981-1990", "1986-1995", "1991-2000"), 
        Era := "post-Soviet"]
age_res$Era <- factor(age_res$Era, levels = c("Soviet", "Transition", "post-Soviet"))

age_res_tab <- age_res[order(YoB), .(YoB, Era, h2, se)]
write.table(age_res_tab, "~/EA_heritability/figures/paper/revision/h2_estimates_YoBbins_INT_EA_kin0.05_noAgeAtAgr.tsv",
            row.names = F, quote = F, sep = "\t")

# Make the LDAK h2 plots
pl <- ggplot(age_res,
             aes(x = YoB, y = h2, color = Era)) +
  geom_point(aes(shape = YoB)) +
  geom_errorbar(aes(ymin = (h2 - 1.96 * se),
                    ymax = (h2 + 1.96 * se),
                    linewidth = YoB),
                width = 0.2) + # , linewidth = 0.5
  theme_bw() +
  theme(text = element_text(size = 10),
        title = element_text(size=8),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
  scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
  scale_color_manual(values = c(rgb(254, 0, 0, maxColorValue = 254), rgb(38, 65, 140, maxColorValue = 254),
                                rgb(80, 148, 205, maxColorValue = 254))) +
  ylab(bquote(h^2)) + xlab("Decade of birth")  + # xlab("Year of Birth") +
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
  guides(linewidth = "none", 
         shape = "none")

pl

pdf("~/EA_heritability/figures/paper/revision/h2_ages_kin0.05_INT.pdf", width=5.5, height=2.5)

plot(pl)

dev.off()




# Increasing bins around 1976
# Upload LDAK h2 results for decade bins
bin_res <- fread("~/EA_heritability/gcta/results/ldak_bin_h2_kin0.05_EA.tsv")
colnames(bin_res) <- c("name", "h2", "se", "Size", "Mega_Intensity", "Int_SD")
bin_res[, bin_size := factor(c(rep(c(2,4,6,8), 2), 10, 10, rep(c(12,14,16,18,20), 2)), levels = seq(2, 20, 2))]
bin_res[, Era := factor(c(rep(c("post-Soviet", "Soviet"), each = 4), "Soviet", "post-Soviet",
                          rep(c("post-Soviet", "Soviet"), each = 5)), levels = c("Soviet", "post-Soviet"))]


# Make the LDAK h2 plots
pl <- ggplot(bin_res,
             aes(x = bin_size, y = h2, color = Era)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = (h2 - 1.96 * se),
                    ymax = (h2 + 1.96 * se)),
                width = 0.2, position = position_dodge(width = 0.5)) + # , linewidth = 0.5
  theme_bw() +
  theme(text = element_text(size = 10),
        title = element_text(size=8),
        panel.grid.major.x = element_blank()) +
  scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
  scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
  scale_color_manual(values = c(rgb(254, 0, 0, maxColorValue = 254),
                                rgb(80, 148, 205, maxColorValue = 254))) +
  ylab(bquote(h^2)) + xlab("Birth-year bin size (years since 1976)")  + 
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
  guides(linewidth = "none", shape = "none")


# calculate the difference between post-Soviet and Soviet h2
bin_res <- bin_res[order(bin_size,Era),]
bin_res <- bin_res[, .(bin_size, Era, h2, se)]
delta_h2 <- c()
se <- c()
for(i in 1:10){
  delta_h2 <- c(delta_h2, bin_res[2*i, h2] - bin_res[2*i-1, h2])
  se <- c(se, sqrt(bin_res[2*i, se^2] + bin_res[2*i-1, se^2]))
}
bin_res_dif <- data.table(bin_size=factor(seq(2, 20, 2)), delta_h2, se)
bin_res_dif[, p := pnorm(-delta_h2, 0, se)]
bin_res_dif[, `-log(p-value)` := -log10(p)]

library(viridis)
pl_dif <- ggplot(bin_res_dif,
             aes(x = bin_size, y = delta_h2, color = `-log(p-value)`)) +
  geom_point() +
  geom_errorbar(aes(ymin = (delta_h2 - 1.96 * se),
                    ymax = (delta_h2 + 1.96 * se)),
                width = 0.2, position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(text = element_text(size = 10),
        title = element_text(size=8),
        panel.grid.major.x = element_blank()) +
  scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
  scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
  ylab(bquote(Delta~h^2)) + xlab("Birth-year bin size (years since 1976)")  + 
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
  scale_colour_viridis(
    option = "viridis",
    direction = 1
  ) +
  guides(linewidth = "none", shape = "none")



plots <- list(pl, pl_dif)

pdf("~/EA_heritability/figures/paper/revision/h2_bandwidths.pdf", width=5.5, height = 5)

print(
  grid.arrange(
    grobs = plots,
    layout_matrix = matrix(c(1,1,2,3), ncol = 2, byrow = T),
    widths = c(1, 0.043)
  )
)

grid.text("a", x = 0.02, y = 0.96, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.02, y = 0.46, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()

bin_res[, `:=`(h2 = round(h2, digits = 3), se = round(se, digits = 3))]
bin_res_dif[, `:=`(delta_h2 = round(delta_h2, digits = 3), se = round(se, digits = 3), 
                   p = round(p, digits = 3), `-log(p-value)` = round(`-log(p-value)`, digits = 3))]
write.table(bin_res, "~/EA_heritability/figures/paper/revision/h2_bandwidths.tsv", row.names = F, quote = F, sep = "\t")
write.table(bin_res_dif, "~/EA_heritability/figures/paper/revision/h2_bandwidths_dif.tsv", row.names = F, quote = F, sep = "\t")
