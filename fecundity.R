library(data.table)

estbb_filtered <- fread("~/EBB_project/phenotypes/EstBB_filtered.tsv")
ebb <- fread("~/EBB_project/phenotypes/query1.tsv")
ebb <- ebb[, c("Person skood", "PersonLocation birthParishName", "PersonLocation residencyParishName", 
               "CONCATSTR(BMIAssembled ageAtBmi)", "CONCATSTR(BMIAssembled bmi)", "CONCATSTR(BMIAssembled height)", "CONCATSTR(Education highestEducationLevel code)")]
colnames(ebb) <- c("skood", "ParishBirth", "ParishRes", "Age_meas", "BMI", "Height", "EA_answerset")
ebb <- merge(estbb_filtered, ebb, by="skood")
ebb <- ebb[Nat=="Eestlane", ]

pca_est <- fread("~/EBB_project/data_filtering/pca/pcs_EstBB_estonian")
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
ebb <- merge(ebb, pca_est, by="vkood")


# Upload polygenic scores (PGS)
pgi_ea <-  fread("/gpfs/space/GI/GV/Projects/PGI_repository_v2/SSGAC_PGI_Repository_v2_EstBB.txt", 
                 select = c("IID", "PGI_EA"))
colnames(pgi_ea)[1] <- "vkood"
ebb <- merge(ebb, pgi_ea)

# Maximum unrelated set, KING degree 1 and 2 excluded
unrelatedSet <- fread("/gpfs/space/GI/GV/EGCUT_data/genotype_data/GSA_arrays/4_latest_freeze/additional_files/king/unrelated_set/unrelated_degree_2unrelated.txt",
                      header = F, col.names = c("vkood", "vkood2"))

# all unrelated
ebb_unrel <- ebb[vkood %in% unrelatedSet$vkood, ]

lm_prs <- lm(paste0("PGI_EA ~ Sex + ", paste("PC", 1:40, sep = "", collapse = " + ")),
              data = ebb_unrel)
ebb_unrel$PGI_EA_adj <- ebb_unrel$PGI_EA - predict(object = lm_prs, newdata = ebb_unrel)
ebb_unrel[, PGI_EA_adj := scale(PGI_EA_adj)]


meanPGI <- ebb_unrel[
  , .(meanPGI = mean(PGI_EA_adj),
      sdPGI = sd(PGI_EA_adj),
      N = .N),
  by = YoB
]

meanPGI[, sePGI := sdPGI / sqrt(N)]
setorder(meanPGI, YoB)

# Keep YoBxSex with >= 10 individuals
meanPGI <- meanPGI[YoB >= 1920 & YoB <= 2002, ]

meanPGI[, rollMeanPGI := frollmean(meanPGI, n = 5, align = "center")]
meanPGI[, rollN := frollsum(N,   n = 5, align = "center")]
meanPGI[, rollSEPGI := 1/sqrt(rollN)]
write.table(meanPGI, "~/EA_heritability/figures/paper/revision/PGI_EA_YoB.tsv",
            row.names = F, quote = F, sep = "\t")


pdf("~/EA_heritability/figures/paper/revision/PGI_EA_YoB.pdf", width=5.5, height=3)

ggplot(meanPGI, aes(x = YoB, rollMeanPGI)) + 
  geom_pointrange(aes(ymin = rollMeanPGI - rollSEPGI, ymax = rollMeanPGI + rollSEPGI), size = 0.15, linewidth = 0.3) + 
  theme_bw() + ylab(expression("5-year sliding window" ~ bar(PGS[EA]))) + xlab("Year of birth") +
  theme(
    text = element_text(size = 10),
    title = element_text(size=8)
  )


dev.off()
