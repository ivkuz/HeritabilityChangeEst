
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


# Upload GCTA h2 results for decade bins
age_res <- fread("~/EA_heritability/gcta/results/age_h2.tsv")
colnames(age_res) <- c("YoB", "h2", "se")
age_res$YoB <- factor(c("1986-1995", "1941-1950", "1991-2000", "1981-1990", "1976-1985", 
                        "1971-1980", "1966-1975", "1961-1970", "1956-1965", "1951-1960", "1946-1955"), 
                      levels = c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", 
                                 "1966-1975", "1971-1980", "1976-1985", "1981-1990", "1986-1995", "1991-2000"))

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
