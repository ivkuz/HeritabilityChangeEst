library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggsignif)


# fig2abc <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig2abc.rds")
fig2abc <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig2abc_pedigree.rds")
# fig2d <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig2d.rds")
# fig2d <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig2d_noAgeAtAgr.rds")
fig2d <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig2d_pedigree_noAgeAtAgr.rds")

fig2 <- c(fig2abc, fig2d)

# pdf("~/EA_heritability/figures/paper/revision/figure2.pdf", width=5.5, height=5)
# pdf("~/EA_heritability/figures/paper/revision/figure2_noAgeAtAgr2d.pdf", width=5.5, height=5)
pdf("~/EA_heritability/figures/paper/revision/figure2_pedigree_noAgeAtAgr2d.pdf", width=5.5, height=5)

print(
  grid.arrange(
    grobs = fig2,
    layout_matrix = matrix(c(1:3, 4,4,4), nrow = 2, byrow = T),
    heights = c(1, 1.2)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.35, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.69, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.02, y = 0.53, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()



# fig3ab <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig3ab.rds")
fig3ab <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig3ab_pedigree.rds")
# fig3cd <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig3cd.rds")
# fig3cd <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig3cd_noAgeAtAgr.rds")
fig3cd <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig3cd_pedigree_noAgeAtAgr.rds")

fig3 <- c(fig3ab, fig3cd)

# pdf("~/EA_heritability/figures/paper/revision/figure3.pdf", width=5.5, height=6)
# pdf("~/EA_heritability/figures/paper/revision/figure3_noAgeAtAgr3cd.pdf", width=5.5, height=6)
pdf("~/EA_heritability/figures/paper/revision/figure3_peedigree_noAgeAtAgr3cd.pdf", width=5.5, height=6)

print(
  grid.arrange(
    grobs = fig3,
    layout_matrix = matrix(c(1,3,NA,3,2,4,NA,4), ncol = 2, byrow = T),
    widths = c(1, 2), heights = c(1, 0.2, 1, 0.2)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.35, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.02, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.35, y = 0.48, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()



# fig5abde <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5abde_weakPGS.rds")
fig5abde <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5abde.rds")
fig5abde <- list(ggplotGrob(fig5abde[[1]]), ggplotGrob(fig5abde[[2]]), ggplotGrob(fig5abde[[3]]), ggplotGrob(fig5abde[[4]]))
# fig5cf <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5cf_weakPGS.rds")
# fig5cf <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5cf_weakPGS_noAgeAtAgr.rds")
fig5cf <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5cf.rds")
# fig5cf <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5cf_noAgeAtAgr.rds")

fig5 <- c(fig5abde, fig5cf)

# pdf("~/EA_heritability/figures/paper/revision/figure5_weakPGS.pdf", width=7.5, height=5.25)
# pdf("~/EA_heritability/figures/paper/revision/figure5_weakPGS_noAgeAtAgr.pdf", width=7.5, height=5.25)
pdf("~/EA_heritability/figures/paper/revision/figure5.pdf", width=7.5, height=5.25)
# pdf("~/EA_heritability/figures/paper/revision/figure5_noAgeAtAgr.pdf", width=7.5, height=5.25)

print(
  grid.arrange(
    grobs = fig5,
    layout_matrix = matrix(c(1,2,5, 3,4,6, NA,NA,6), nrow = 3, byrow = T),
    widths = c(1, 1.55, 1.55), heights = c(1, 1, 0.04)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.26, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.64, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.02, y = 0.5, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("e", x = 0.26, y = 0.5, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("f", x = 0.64, y = 0.5, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


# INT

fig5abde <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5_INT_abde.rds")
fig5abde <- list(ggplotGrob(fig5abde[[1]]), ggplotGrob(fig5abde[[2]]), ggplotGrob(fig5abde[[3]]), ggplotGrob(fig5abde[[4]]))
fig5cf <-  readRDS("~/EA_heritability/figures/paper/revision/files_for_figures/fig5_INT_cf.rds")

fig5 <- c(fig5abde, fig5cf)

pdf("~/EA_heritability/figures/paper/revision/figure5_INT.pdf", width=7.5, height=5.25)

print(
  grid.arrange(
    grobs = fig5,
    layout_matrix = matrix(c(1,2,5, 3,4,6, NA,NA,6), nrow = 3, byrow = T),
    widths = c(1, 1.55, 1.55), heights = c(1, 1, 0.04)
  )
)

grid.text("a", x = 0.02, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("b", x = 0.26, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("c", x = 0.64, y = 0.98, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("d", x = 0.02, y = 0.5, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("e", x = 0.26, y = 0.5, gp = gpar(fontsize=14, fontface = "bold"))
grid.text("f", x = 0.64, y = 0.5, gp = gpar(fontsize=14, fontface = "bold"))

dev.off()


