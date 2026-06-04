
###############################################################
# Make files with lists of individuals (related or unrelated) #
###############################################################


library(data.table)

filter_relatives <- function(rel, vcodes, k=100){
  rel <- rel[ID1 %in% vcodes & ID2 %in% vcodes,]
  if(nrow(rel) == 0){
    return(vcodes)
  }
  u <- length(unique(c(rel$ID1, rel$ID2)))
  exclude <- c()
  i=0
  n=2
  while(n>1){
    i = i + 1
    # print(n)
    rel_tab <- as.data.table(table(c(rel$ID1, rel$ID2)))
    rel_tab <- rel_tab[order(N, decreasing=T),]
    n <- rel_tab[1, N]
    hub <- rel_tab[N == n, V1][1:min(length(rel_tab[N == n, V1]), k)]
    exclude <- c(exclude, hub)
    rel <- rel[!(ID1 %in% hub) & !(ID2 %in% hub),]
  }
  exclude <- c(exclude, rel$ID1)
  non_relatives <- vcodes[!(vcodes %in% exclude)]
  return(non_relatives)
}


estbb_filtered <- fread("~/EBB_project/phenotypes/EstBB_filtered.tsv")
ebb <- fread("~/EBB_project/phenotypes/query1.tsv")
ebb <- ebb[, c("Person skood", "PersonLocation birthParishName", "PersonLocation residencyParishName", 
               "CONCATSTR(Education highestEducationLevel code)")]
colnames(ebb) <- c("skood", "ParishBirth", "ParishRes", "EA_answerset")
ebb <- merge(estbb_filtered, ebb, by="skood")
ebb <- ebb[Nat=="Eestlane", ]
ebb[, Age2 := Age^2]
ebb[, SxA := Age*Sex]

ebb2 <- fread("~/EBB_project/phenotypes/query_bmi_edu.tsv")
ebb2 <- ebb2[, c("Person skood", "PersonPortrait lastEducation code")]
colnames(ebb2) <- c("skood", "EA_portrait")
ebb <- merge(ebb, ebb2, by="skood")

pca_est <- fread("~/EBB_project/data_filtering/pca/pcs_EstBB_estonian")
pca_est <- pca_est[,c("IID", paste0("PC", 1:100))]
colnames(pca_est)[1] <- "vkood"
pca_rus <- fread("~/EBB_project/data_filtering/pca/pcs_EstBB_russian")
pca_rus <- pca_rus[,c("IID", paste0("PC", 1:100))]
colnames(pca_rus)[1] <- "vkood"
pca <- rbind(pca_est, pca_rus)
ebb <- merge(ebb, pca, by="vkood")

settlement <- fread("~/EA_heritability/data/query_settletype.tsv")
colnames(settlement) <- c("skood", "settlCode", "settlName")
ebb <- merge(ebb, settlement, by="skood")




king <- "~/EBB_project/GWAS/relatedness/king.kin0"
rel <- fread(king)

ebb_test <- ebb[AgeAtAgr > 25 & Nat =="Eestlane",]


# Unfiltered
phases <- rep(c("", "p1", "p2"), each = 3)
eras <- rep(c("all", "ps", "s"), 3)
for(age in c(10, 15)){
  
  ind_00  <- data.table(0, ebb_test[,                                vkood]) # moved to fle without age manually 
  ind_01  <- data.table(0, ebb_test[YoB >= 1991 - age,               vkood])
  ind_02  <- data.table(0, ebb_test[YoB < 1991 - age,                vkood])
  ind_10  <- data.table(0, ebb_test[                    YoA < 2016,  vkood]) # moved to fle without age manually 
  ind_11  <- data.table(0, ebb_test[YoB >= 1991 - age & YoA < 2016,  vkood])
  ind_12  <- data.table(0, ebb_test[YoB < 1991 - age &  YoA < 2016,  vkood])
  ind_20  <- data.table(0, ebb_test[                    YoA >= 2016, vkood]) # moved to fle without age manually 
  ind_21  <- data.table(0, ebb_test[YoB >= 1991 - age & YoA >= 2016, vkood])
  ind_22  <- data.table(0, ebb_test[YoB < 1991 - age &  YoA >= 2016, vkood])
  
  ind_list <- list(ind_00, ind_01, ind_02, ind_10, ind_11, ind_12, ind_20, ind_21, ind_22)
  
  for(i in 1:9){
    
    write.table(ind_list[[i]], paste0("~/EA_heritability/gcta/data/", phases[i], eras[i], "_ldak_list_all_", age, ".tsv"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    
  }
  
}

# Age bins unfiltered
for(year in seq(1940, 1990, 5)){
  
  ind <- data.table(0, ebb_test[YoB > year & YoB <= year + 10, vkood])
  
  write.table(ind, paste0("~/EA_heritability/gcta/data/ldak_list_all_", year, ".tsv"),
              row.names = F, col.names = F, quote = F, sep = "\t")
  
}



# Unfiltered by sex
sexes <- c("male", "female")
for(age in c(10, 15)){
  for(s in c(1, 2)){

    ind_01  <- data.table(0, ebb_test[Sex==s & YoB >= 1991 - age, vkood])
    ind_02  <- data.table(0, ebb_test[Sex==s & YoB < 1991 - age, vkood])

    write.table(ind_01, paste0("~/EA_heritability/gcta/data/", "ps", "_ldak_list_all_", sexes[s], "_", age, ".tsv"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(ind_02, paste0("~/EA_heritability/gcta/data/", "s", "_ldak_list_all_", sexes[s], "_", age, ".tsv"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    
  }
}


# Unfiltered by settlement type
settl <- c("city", "town", "rural")
settl_codes <- c("L", "V", "M")
for(age in c(10, 15)){
  for(s in 1:3){
    
    ind_01  <- data.table(0, ebb_test[settlCode==settl_codes[s] & YoB >= 1991 - age, vkood])
    ind_02  <- data.table(0, ebb_test[settlCode==settl_codes[s] & YoB < 1991 - age, vkood])
    
    write.table(ind_01, paste0("~/EA_heritability/gcta/data/", "ps", "_ldak_list_all_", settl[s], "_", age, ".tsv"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(ind_02, paste0("~/EA_heritability/gcta/data/", "s", "_ldak_list_all_", settl[s], "_", age, ".tsv"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    
  }
}



# # Sex
# print("Sex 0")
# ind_0_male   <- filter_relatives(rel = rel, vcodes = ebb_test[Sex == 1, vkood], k = 1) # 51455
# ind_0_female <- filter_relatives(rel = rel, vcodes = ebb_test[Sex == 2, vkood], k = 1) # 96366
# print("Sex 1")
# ind_1_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoA <  2016 & Sex == 1, vkood], k = 1) # 10338
# ind_1_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoA <  2016 & Sex == 2, vkood], k = 1) # 20129 5
# print("Sex 2")
# ind_2_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 & Sex == 1, vkood], k = 1) # 41117 5
# ind_2_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 & Sex == 2, vkood], k = 1) # 76237 10
# 
# # type of settlement
# print("Settlement 0")
# ind_0_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[settlCode=="M", vkood], k = 1) # 5
# ind_0_town   <- filter_relatives(rel = rel, vcodes = ebb_test[settlCode=="V", vkood], k = 1)
# ind_0_city   <- filter_relatives(rel = rel, vcodes = ebb_test[settlCode=="L", vkood], k = 1) # 10
# print("Settlement 1")
# ind_1_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoA < 2016 & settlCode=="M", vkood], k = 1)
# ind_1_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoA < 2016 & settlCode=="V", vkood], k = 1)
# ind_1_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoA < 2016 & settlCode=="L", vkood], k = 1)
# print("Settlement 2")
# ind_2_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 & settlCode=="M", vkood], k = 1) # 5
# ind_2_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 & settlCode=="V", vkood], k = 1)
# ind_2_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 & settlCode=="L", vkood], k = 1) # 10
# 
# # place of birth
# print("PoB 0")
# ind_0_BO <- filter_relatives(rel = rel, vcodes = ebb_test[!(PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
# ind_0_BC <- filter_relatives(rel = rel, vcodes = ebb_test[ (PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
# print("PoB 1")
# ind_1_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoA <  2016 & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1) 
# ind_1_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoA <  2016 &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1) 
# print("PoB 2")
# ind_2_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
# ind_2_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
# 
# # place of residence
# print("PoR 0")
# ind_0_RO <- filter_relatives(rel = rel, vcodes = ebb_test[!(PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
# ind_0_RC <- filter_relatives(rel = rel, vcodes = ebb_test[ (PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
# print("PoR 1")
# ind_1_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoA <  2016 & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1) 
# ind_1_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoA <  2016 &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1) 
# print("PoR 2")
# ind_2_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
# ind_2_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoA >= 2016 &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
# 
# 
# 
# l <- list(ind_0_male, ind_0_female, ind_1_male, ind_1_female, ind_2_male, ind_2_female,
#           ind_0_rural, ind_0_town, ind_0_city, ind_1_rural, ind_1_town, ind_1_city,  
#           ind_2_rural, ind_2_town, ind_2_city,
#           ind_0_BO, ind_0_BC, ind_1_BO, ind_1_BC, ind_2_BO, ind_2_BC,
#           ind_0_RO, ind_0_RC, ind_1_RO, ind_1_RC, ind_2_RO, ind_2_RC)
# name <- c("male", "female", "p1male", "p1female", "p2_male", "p2_female",
#           "rural", "town", "city", "p1rural", "p1town", "p1city",  
#           "p2rural", "p2town", "p2city",
#           "BO", "BC", "p1BO", "p1BC", "p2BO", "p2BC",
#           "RO", "RC", "p1RO", "p1RC", "p2RO", "p2RC")


# for(age in c(10, 15)){
#   # Sex
#   print("Sex 0")
#   ind_0_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & Sex == 1, vkood], k = 1) # 51455
#   ind_0_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & Sex == 2, vkood], k = 1) # 96366
#   print("Sex 1")
#   ind_1_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA <  2016 & Sex == 1, vkood], k = 1) # 10338
#   ind_1_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA <  2016 & Sex == 2, vkood], k = 1) # 20129 5
#   print("Sex 2")
#   ind_2_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 & Sex == 1, vkood], k = 1) # 41117 5
#   ind_2_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 & Sex == 2, vkood], k = 1) # 76237 10
#   
#   # type of settlement
#   print("Settlement 0")
#   ind_0_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & settlCode=="M", vkood], k = 1) # 5
#   ind_0_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & settlCode=="V", vkood], k = 1)
#   ind_0_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & settlCode=="L", vkood], k = 1) # 10
#   print("Settlement 1")
#   ind_1_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA < 2016 & settlCode=="M", vkood], k = 1)
#   ind_1_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA < 2016 & settlCode=="V", vkood], k = 1)
#   ind_1_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA < 2016 & settlCode=="L", vkood], k = 1)
#   print("Settlement 2")
#   ind_2_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 & settlCode=="M", vkood], k = 1) # 5
#   ind_2_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 & settlCode=="V", vkood], k = 1)
#   ind_2_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 & settlCode=="L", vkood], k = 1) # 10
#   
#   # place of birth
#   print("PoB 0")
#   ind_0_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
#   ind_0_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
#   print("PoB 1")
#   ind_1_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA <  2016 & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1) 
#   ind_1_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA <  2016 &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1) 
#   print("PoB 2")
#   ind_2_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
#   ind_2_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
#   
#   # place of residence
#   print("PoR 0")
#   ind_0_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
#   ind_0_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
#   print("PoR 1")
#   ind_1_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA <  2016 & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1) 
#   ind_1_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA <  2016 &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1) 
#   print("PoR 2")
#   ind_2_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
#   ind_2_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB >= 1991 - age & YoA >= 2016 &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10 
#   
#   
#   l <- list(ind_0_male, ind_0_female, ind_1_male, ind_1_female, ind_2_male, ind_2_female,
#             ind_0_rural, ind_0_town, ind_0_city, ind_1_rural, ind_1_town, ind_1_city,  
#             ind_2_rural, ind_2_town, ind_2_city,
#             ind_0_BO, ind_0_BC, ind_1_BO, ind_1_BC, ind_2_BO, ind_2_BC,
#             ind_0_RO, ind_0_RC, ind_1_RO, ind_1_RC, ind_2_RO, ind_2_RC)
#   name <- c("male", "female", "p1male", "p1female", "p2_male", "p2_female",
#             "rural", "town", "city", "p1rural", "p1town", "p1city",  
#             "p2rural", "p2town", "p2city",
#             "BO", "BC", "p1BO", "p1BC", "p2BO", "p2BC",
#             "RO", "RC", "p1RO", "p1RC", "p2RO", "p2RC")
#   
#   for(i in 1:length(l)){
#     ind <- l[[i]]
#     n <- name[i]
#     write.table(data.frame(ind, ind),
#                 paste0("~/EA_heritability/gcta/data/ps_", n, "_list_unrel_", age, ".tsv"),
#                 quote = F, sep = "\t", col.names = F, row.names = F)
#   }
#   
# }
# 


for(age in c(10, 15)){
  # Sex
  print("Sex 0")
  ind_0_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & Sex == 1, vkood], k = 1) # 51455
  ind_0_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & Sex == 2, vkood], k = 1) # 96366
  print("Sex 1")
  ind_1_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA <  2016 & Sex == 1, vkood], k = 1) # 10338
  ind_1_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA <  2016 & Sex == 2, vkood], k = 1) # 20129 5
  print("Sex 2")
  ind_2_male   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 & Sex == 1, vkood], k = 1) # 41117 5
  ind_2_female <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 & Sex == 2, vkood], k = 1) # 76237 10

  # type of settlement
  print("Settlement 0")
  ind_0_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & settlCode=="M", vkood], k = 1) # 5
  ind_0_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & settlCode=="V", vkood], k = 1)
  ind_0_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & settlCode=="L", vkood], k = 1) # 10
  print("Settlement 1")
  ind_1_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA < 2016 & settlCode=="M", vkood], k = 1)
  ind_1_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA < 2016 & settlCode=="V", vkood], k = 1)
  ind_1_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA < 2016 & settlCode=="L", vkood], k = 1)
  print("Settlement 2")
  ind_2_rural  <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 & settlCode=="M", vkood], k = 1) # 5
  ind_2_town   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 & settlCode=="V", vkood], k = 1)
  ind_2_city   <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 & settlCode=="L", vkood], k = 1) # 10

  # place of birth
  print("PoB 0")
  ind_0_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
  ind_0_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
  print("PoB 1")
  ind_1_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA <  2016 & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1)
  ind_1_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA <  2016 &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1)
  print("PoB 2")
  ind_2_BO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 & !(PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10
  ind_2_BC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 &  (PoB %in% c("Harju", "Tartu")), vkood], k = 1) # 10

  # place of residence
  print("PoR 0")
  ind_0_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10
  ind_0_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10
  print("PoR 1")
  ind_1_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA <  2016 & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1)
  ind_1_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA <  2016 &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1)
  print("PoR 2")
  ind_2_RO <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 & !(PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10
  ind_2_RC <- filter_relatives(rel = rel, vcodes = ebb_test[YoB < 1991 - age & YoA >= 2016 &  (PoR %in% c("Harju", "Tartu")), vkood], k = 1) # 10


  l <- list(ind_0_male, ind_0_female, ind_1_male, ind_1_female, ind_2_male, ind_2_female,
            ind_0_rural, ind_0_town, ind_0_city, ind_1_rural, ind_1_town, ind_1_city,
            ind_2_rural, ind_2_town, ind_2_city,
            ind_0_BO, ind_0_BC, ind_1_BO, ind_1_BC, ind_2_BO, ind_2_BC,
            ind_0_RO, ind_0_RC, ind_1_RO, ind_1_RC, ind_2_RO, ind_2_RC)
  name <- c("male", "female", "p1male", "p1female", "p2_male", "p2_female",
            "rural", "town", "city", "p1rural", "p1town", "p1city",
            "p2rural", "p2town", "p2city",
            "BO", "BC", "p1BO", "p1BC", "p2BO", "p2BC",
            "RO", "RC", "p1RO", "p1RC", "p2RO", "p2RC")

  for(i in 1:length(l)){
    ind <- l[[i]]
    n <- name[i]
    write.table(data.frame(ind, ind),
                paste0("~/EA_heritability/gcta/data/s_", n, "_list_unrel_", age, ".tsv"),
                quote = F, sep = "\t", col.names = F, row.names = F)
  }

}




inputfile <- "~/EA_heritability/gcta/data/ind_list_files.txt"
outfile <- "~/EA_heritability/gcta/data/ldak_out_files.txt"
covarfile <- "~/EA_heritability/gcta/data/ldak_covar_files.txt"
name <- c("male", "female", "p1male", "p1female", "p2_male", "p2_female",
          "rural", "town", "city", "p1rural", "p1town", "p1city",
          "p2rural", "p2town", "p2city",
          "BO", "BC", "p1BO", "p1BC", "p2BO", "p2BC",
          "RO", "RC", "p1RO", "p1RC", "p2RO", "p2RC")
name_sex <- c("male", "female", "p1male", "p1female", "p2_male", "p2_female")
for(era in c("s", "ps")){
  for(age in c(10, 15)){
    for(n in name){
      # write(paste0("~/EA_heritability/gcta/data/", era, "_", n, "_list_unrel_", age, ".tsv"), file=inputfile, append=TRUE)
      # write(paste0("~/EA_heritability/gcta/results/ldak_", era, "_", n, "_", age), file=outfile, append=TRUE)
      if(n %in% name_sex){
        write("~/EA_heritability/gcta/data/covarLDAK_nosex.tsv", file=covarfile, append=TRUE)
      } else{
        write("~/EA_heritability/gcta/data/covarLDAK.tsv", file=covarfile, append=TRUE)
      }
      # writeLines(paste0("~/EA_heritability/gcta/data/", era, "_", n, "_list_unrel_", age, ".tsv"), inputfiles)
      # writeLines(paste0("~/EA_heritability/gcta/results/ldak_", era, "_", n, "_", age), outfiles)

      # ind <- fread(paste0("~/EA_heritability/gcta/data/", era, "_", n, "_list_unrel_", age, ".tsv"), header = F)
      # write.table(data.frame(0, ind$V1),
      #             paste0("~/EA_heritability/gcta/data/", era, "_", n, "_list_unrel_", age, ".tsv"),
      #             quote = F, sep = "\t", col.names = F, row.names = F)
    }

  }
}

for(phase in c("", "p1", "p2")){
  for(era in c("s", "ps")){
    for(age in c("_10", "")){
      # write(paste0("~/EA_heritability/gcta/data/", phase, era, "_ldak_list_unrel", age, ".tsv"), file=inputfile, append=TRUE)
      if(age==""){
        age <- "_15"
      }
      # write(paste0("~/EA_heritability/gcta/results/ldak_", phase, era, age), file=outfile, append=TRUE)
      write("~/EA_heritability/gcta/data/covarLDAK.tsv", file=covarfile, append=TRUE)
      

      # ind <- fread(paste0("~/EA_heritability/gcta/data/", phase, era, "_list_unrel", age, ".tsv"), header = F)
        # write.table(data.frame(0, ind$V1),
        #             paste0("~/EA_heritability/gcta/data/", phase, era, "_ldak_list_unrel", age, ".tsv"),
        #             quote = F, sep = "\t", col.names = F, row.names = F)
    }
  }
}




head(ebb)
