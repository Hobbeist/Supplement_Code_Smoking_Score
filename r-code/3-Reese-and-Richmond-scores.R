#==========================================================
# Creating the Reese and the Richmond scores: 
# Example Raine Study (done in the same way for NFBC1986 and NFBC1966)
#
# 1. Reese
# 2. Richmond 568 CpGs
# 3. Richmond 19 CpGs
# 
#==========================================================

# 1. Reese Score

# Load the Raine study CpG data for Reese, excluding missing values
Raine_reese_raw <- na.omit(readRDS("/path/to/raine/data.rds"))

# Loop to apply the Reese betas. Weights taken from the  supplemental material Table S1 of 
# Reese SE, Zhao S, Wu MC, Joubert BR, Parr CL, Haberg SE, et al. 2017. 
# DNA methylation score as a biomarker in newborns for sustained maternal smoking during pregnancy. Environmental health perspectives 125:760-766.
y   <- NULL
SMK <- c()

for (i in 1:length(Raine_reese_raw[,1])) {
  
  y <- -1.667*Raine_reese_raw[i,"cg00709966"] - 0.191*Raine_reese_raw[i,"cg02256631"] + 2.706*Raine_reese_raw[i,"cg02482603"] + 1.786*Raine_reese_raw[i,"cg04103532"] + 14.027*Raine_reese_raw[i,"cg04180046"] +   
    2.318*Raine_reese_raw[i,"cg04506190"] + 6.210*Raine_reese_raw[i,"cg05549655"] -10.909*Raine_reese_raw[i,"cg05575921"] +1.142*Raine_reese_raw[i,"cg08698721"] - 6.330*Raine_reese_raw[i,"cg09743950"] -
    4.963*Raine_reese_raw[i,"cg10799846"] - 0.370*Raine_reese_raw[i,"cg11864574"] + 3.847*Raine_reese_raw[i,"cg12186702"] + 1.514*Raine_reese_raw[i,"cg13834112"] - 0.963*Raine_reese_raw[i,"cg13893782"] -
    6.304*Raine_reese_raw[i,"cg14179389"] + 6.361*Raine_reese_raw[i,"cg14351425"] + 5.050*Raine_reese_raw[i,"cg14633298"] + 2.286*Raine_reese_raw[i,"cg14743346"] - 2.912*Raine_reese_raw[i,"cg17397069"] +
    5.245*Raine_reese_raw[i,"cg19381766"] - 0.773*Raine_reese_raw[i,"cg22154659"] - 0.254*Raine_reese_raw[i,"cg22802102"] + 0.011*Raine_reese_raw[i,"cg23304605"] - 3.903*Raine_reese_raw[i,"cg25189904"] -
    46.991*Raine_reese_raw[i,"cg25949550"] -0.246*Raine_reese_raw[i,"cg26764244"] + 0.836*Raine_reese_raw[i,"cg27291468"]
  SMK <- c(SMK,y)
  
}


#==========================================================
# Richmond Score: 568 CpGs
#==========================================================
# load Joubert CpGs (taken from the Supplement 4, Excel table S3 of 
# Joubert BR, Felix JF, Yousefi P, Bakulski KM, Just AC, Breton C, et al. 2016. 
# DNA methylation in newborns and maternal smoking in pregnancy: Genome-wide consortium meta-analysis. American journal of human genetics 98:680-696.)
#==========================================================
# The data used is the named "Meta-Analysis of sustained smoking and newborn methylation adjusted for cell type" in the supplement
# load the Joubert data
joubert.cpg           <- read_csv("data/Joubert-cpgs.csv")
bonferroni.threshold  <- 1.07613e-07

# Select the CpGs that Richmond used (Bonferroni corrected); this will be the 568 CpGs
richmond.cpg <- joubert.cpg %>%
  filter(P_2 < bonferroni.threshold) %>%
  select(CpG,Coef_2) %>%
  rename(beta = Coef_2)

# Save all the CpG names so it can be selected from NFBC
saveRDS(richmond.cpg$CpG, file="data/richmond-cpgs_normal.rds")

# This data contains the CpGs necessary to create the Richmond score
# The data only contains the 995 Caucasian Raine Study participants

richmond.data <- read_csv("data/richmond-score-data-raine.csv")

# Calculate the Richmond score

SCORE     <- NULL
score     <- NULL
cpg       <- NULL

# In case of unavailability of some CpGs, first check which CpGs are available
availCpGs <- names(richmond.data)[names(richmond.data) %in% richmond.cpg$CpG]

for (i in availCpGs) {
  
  CPG <- as.numeric(richmond.cpg %>% filter(CpG %in% i) %>% 
                      select(beta))
  score <- cbind(score, as.numeric(unlist((CPG * richmond.data[,i]))))
  cpg <- c(cpg, i)
  
}

score                             <- as.data.frame(score)
names(score)                      <- cpg
richmond.data$richmond_568  <- rowSums(score) 


#==========================================================
# Richmond score: 19 CpGs
#==========================================================
# (taken from the Supplement 4, Excel table S3 of 
# Joubert BR, Felix JF, Yousefi P, Bakulski KM, Just AC, Breton C, et al. 2016. 
# DNA methylation in newborns and maternal smoking in pregnancy: 
# Genome-wide consortium meta-analysis. American journal of human genetics 98:680-696.)
#==========================================================
# The data used is the named "Meta-Analysis of sustained smoking and methylation in older children" in the supplement

# load Joubert CpGs
joubert.cpg           <- read_csv("data/Joubert-cpgs-older-kids.csv")
bonferroni.threshold  <- 1.07613e-07

# Select the CpGs that Richmond used (Bonferroni corrected); this will be the 19 CpGs.
richmond.cpg <- joubert.cpg %>%
  filter(P_old < bonferroni.threshold) %>%
  select(CpG,Coef_old) %>%
  rename(beta = Coef_old)

# Load Richmond score data, created on Nimbus. This data
# contains the CpGs necessary to create the Richmond score, although some of them are missing.
# The data only contains the 995 Caucasian Raine Study participants

richmond.data <- read_csv("data/richmond-score-data-raine.csv")


# Calculate the Richmond score

SCORE     <- NULL
score     <- NULL
cpg       <- NULL
availCpGs <- names(richmond.data)[names(richmond.data) %in% richmond.cpg$CpG]

for (i in availCpGs) {
  
  CPG <- as.numeric(richmond.cpg %>% filter(CpG %in% i) %>% 
                      select(beta))
  score <- cbind(score, as.numeric(unlist((CPG * richmond.data[,i]))))
  cpg <- c(cpg, i)
  
}

score                             <- as.data.frame(score)
names(score)                      <- cpg
richmond.data$richmond_19  <- rowSums(score) 

#=========================================================================
# END
#=========================================================================




