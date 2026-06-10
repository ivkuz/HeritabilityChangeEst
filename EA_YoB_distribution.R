library(data.table)
library(ggplot2)


ebb <- fread("~/EBB_project/phenotypes/EstBB_filtered.tsv")
ebb <- ebb[Nat=="Eestlane", ]

# Upload and process EA
ebb2 <- fread("~/EBB_project/phenotypes/query_bmi_edu.tsv")
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA_portrait")
ebb <- merge(ebb, ebb2, by="skood")
ebb <- ebb[!is.na(EA_portrait),]

ebb <- ebb[AgeAtAgr > 25,]

ea_freq <- ebb[, .N, by = .(YoB, EA_portrait)][, prop := N / sum(N), by = YoB]

pdf("~/EA_heritability/figures/paper/revision/YoB_threshold.pdf")
ggplot(ea_freq[YoB > 1925 & YoB <= 1995, ], aes(x=YoB, y = N, color = as.character(EA_portrait), group = EA_portrait)) + 
  geom_line() + theme_bw() + guides(color=guide_legend(title="EA"))
ggplot(ea_freq[YoB > 1925 & YoB <= 1995, ], aes(x=YoB, y = prop, color = as.character(EA_portrait), group = EA_portrait)) + 
  geom_line() + theme_bw() + guides(color=guide_legend(title="EA"))

ggplot(ea_freq[YoB > 1925 & YoB <= 1995 & EA_portrait == 8, ], aes(x=YoB, y = N, color = as.character(EA_portrait), group = EA_portrait)) + 
  geom_line() + theme_bw() + guides(color=guide_legend(title="EA"))
ggplot(ea_freq[YoB > 1925 & YoB <= 1995 & EA_portrait == 8, ], aes(x=YoB, y = prop, color = as.character(EA_portrait), group = EA_portrait)) + 
  geom_line() + theme_bw() + guides(color=guide_legend(title="EA"))

dev.off()

# ggplot(ea_freq[YoB > 1925 & YoB <= 1990, ], aes(x=YoB, y = N, color = as.character(EA_portrait), group = EA_portrait)) + geom_line()
# ggplot(ea_freq[YoB > 1925 & YoB <= 1990, ], aes(x=YoB, y = prop, color = as.character(EA_portrait), group = EA_portrait)) + geom_line()
# 
# ggplot(ea_freq[YoB > 1925 & YoB <= 1985, ], aes(x=YoB, y = N, color = as.character(EA_portrait), group = EA_portrait)) + geom_line()
# ggplot(ea_freq[YoB > 1925 & YoB <= 1985, ], aes(x=YoB, y = prop, color = as.character(EA_portrait), group = EA_portrait)) + geom_line()
# 
# ebb[YoB <= 1985, min(AgeAtAgr)]
# 
# ea_freq2 <- ebb[, .N, by = .(AgeAtAgr, EA_portrait)][, prop := N / sum(N), by = AgeAtAgr]
# 
# ggplot(ea_freq2[AgeAtAgr < 70 & AgeAtAgr > 25, ], aes(x=AgeAtAgr, y = N, color = as.character(EA_portrait), group = EA_portrait)) + geom_line()
# ggplot(ea_freq2[AgeAtAgr < 70 & AgeAtAgr > 25, ], aes(x=AgeAtAgr, y = prop, color = as.character(EA_portrait), group = EA_portrait)) + geom_line()
