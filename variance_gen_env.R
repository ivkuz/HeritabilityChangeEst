library(data.table)
library(readr)
library(ggplot2)
library(grid)
library(gridExtra)

#########################
# kin 0.05 noAgeatAgr ###
#########################

traits <- c("kin0.05_noAgeAtAgr", "OS_kin0.05_noAgeAtAgr", "Height_kin0.05_noAgeAtAgr", "BMI_kin0.05_noAgeAtAgr")
names <- c("EA", "OS", "Height", "BMI")
plots <- list()
results <- data.frame()

for(j in 1:4){
  
  trait <- traits[j]
  if(trait == "kin0.05_noAgeAtAgr"){
    trait_name <- "Educational Attainment"
  }else if(trait == "OS_kin0.05_noAgeAtAgr"){
    trait_name <- "Occupational Status"
  }else {
    trait_name <- sub("_kin0.05_noAgeAtAgr", "", trait)
  }

  res <- data.table()
  years <- c()
  for(i in seq(1940, 1990, 5)){
    
    res_tmp <- fread(paste0("~/EA_heritability/gcta/results/ldak_", trait, "_", i, ".vars"), header = T)
    res_tmp[, Component := c("G", "E")]
    year <- paste(i+1, i+10, sep = "-")
    res <- rbind(res, data.table(res_tmp, YoB = year))
    
    years <- c(years, year)
    
  }
  
  res[, YoB := factor(YoB, levels = years)]
  res[, Component := factor(Component, levels = c("G", "E"))]
  
  res[YoB %in% c("1941-1950", "1946-1955", "1951-1960", "1956-1965","1961-1970", "1966-1975"), 
      Era := "Soviet"]
  res[YoB == "1971-1980", Era := "Transition"]
  res[YoB %in% c("1976-1985", "1981-1990", "1986-1995", "1991-2000"), 
      Era := "post-Soviet"]
  res$Era <- factor(res$Era, levels = c("Soviet", "Transition", "post-Soviet"))
  res$Trait <- names[j]
  res <- res[, .(Trait, Era, YoB, Component, Variance, SE)]
  
  results <- rbind(results, res)
  
  pl <- ggplot(res,
               aes(x = YoB, y = Variance, color = Era)) +
    geom_point(aes(shape = Component)) +
    geom_errorbar(aes(ymin = (Variance - 1.96 * SE),
                      ymax = (Variance + 1.96 * SE),
                      linewidth = YoB),
                  width = 0.2) + # , linewidth = 0.5
    theme_bw() +
    theme(text = element_text(size = 10),
          title = element_text(size=8),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    # scale_shape_manual(values = c(rep(c(19, 1), 5), 19)) +
    scale_linewidth_manual(values = c(rep(c(0.5, 0.2), 5), 0.5)) +
    scale_color_manual(values = c(rgb(254, 0, 0, maxColorValue = 254), rgb(38, 65, 140, maxColorValue = 254),
                                  rgb(80, 148, 205, maxColorValue = 254))) +
    ylab("Variance") + xlab("Decade of birth") + ggtitle(trait_name) + 
    scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")) +
    guides(linewidth = "none")
  
  plots <- append(plots, list(pl))
  
}

write.table(results, "~/EA_heritability/figures/paper/revision/var_GE_ages_kin0.05_noAgeatAgr.tsv", 
            row.names = F, quote = F, sep = "\t")

pdf("~/EA_heritability/figures/paper/revision/var_GE_ages_kin0.05_noAgeatAgr.pdf", width=7, height = 6)

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

