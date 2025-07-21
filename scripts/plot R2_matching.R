
#####################################################################
# Plot Fig. 5 and SF 9 ##############################################
# Variance of EA explained by PGSEA in subsamples from the groups ###
# (defined by wave participation and era) selected based on #########
# the distribution of EA in another subcohort #######################
#####################################################################


library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)


r2_list <- readRDS("~/EA_heritability/scripts/paper/r2_matched.RDS")

plot_list_EA <- list()
plot_list_EA_Sex <- list()
for(l in r2_list){

  plot_list_EA <- c(plot_list_EA,
                    list(ggplot(data.table(r2=l[[1]]$r2_EA), aes(x=r2)) + 
                           geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 0.002, fill = "#888") +
                           geom_vline(xintercept=as.numeric(l[[2]]), color="red", linetype = "dashed", linewidth = 0.5) + # linewidth = 1
                           geom_vline(xintercept=as.numeric(l[[3]]), color="black", linetype = "dashed", linewidth = 0.5) +
                           theme_bw() +
                           theme(axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 text = element_text(size = 10)) +
                         xlim(0.05, 0.13) + ylim(0, 0.33)
                    )
  )
  plot_list_EA_Sex <- c(plot_list_EA_Sex,
                        list(ggplot(data.table(r2=l[[1]]$r2_EA_Sex), aes(x=r2)) + 
                               geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 0.002, fill = "#888") +
                               geom_vline(xintercept=as.numeric(l[[2]]), color="red", linetype = "dashed", linewidth = 0.5) +
                               geom_vline(xintercept=as.numeric(l[[3]]), color="black", linetype = "dashed", linewidth = 0.5) +
                               theme_bw() +
                               theme(axis.title.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     text = element_text(size = 10)) +
                             xlim(0.04, 0.15) + ylim(0, 0.32)
                        )
  )
}


pdf("~/EA_heritability/figures/paper/matching_replacement2.pdf", width=6, height = 6)

print(
  grid.arrange(
    grobs = plot_list_EA,
    layout_matrix = matrix(c(6,5,8,7,2,1,4,3,14,13,16,15,10,9,12,11), ncol=4, byrow = F)
  )
)
print(
  grid.arrange(
    grobs = plot_list_EA_Sex,
    layout_matrix = matrix(c(6,5,8,7,2,1,4,3,14,13,16,15,10,9,12,11), ncol=4, byrow = F)
  )
)

dev.off()



