#MODIFIED DEAP LINEAR MIXED EFFECT MODEL SCRIPT
#Max M Owens, PhD - University of Vermont (Updated Dec 2020)

#set this to your working directory
setwd("/home/max/Documents/linear_mixed_model_abcd/")

#load libraries you need
library(gamm4)
library(MuMIn)

#load in mixed model functions
mixed_model <- dget("mixed_model.R")
mixed_model_nocov <- dget("mixed_model_nocov.R")
matrixify <- dget("matrixify.R")

#read lists of split half subs
half1_subs <- readLines("subj_list_half1.txt")
half2_subs <- readLines("subj_list_half2.txt")

#read in data
data <- readRDS(paste0("/home/max/Documents/linear_mixed_model_abcd/nda2.0.1.Rds"))
data <- data[which(data$eventname == "baseline_year_1_arm_1"), ]

#variables to keep in dataset
ivs <- readLines("/home/max/Documents/linear_mixed_model_abcd/rsvarnames_rds.txt")
dvs <- readLines("/home/max/Documents/linear_mixed_model_abcd/rs_networks.txt")
covs <- c("race.4level", "sex", "high.educ.bl", "household.income.bl", "age", "iqc_rsfmri_all_mean_motion")
re_covs <- c("rel_family_id", "mri_info_device.serial.number")
proc_vars <- c("src_subject_id", "eventname", "fsqc_qc", "iqc_rsfmri_good_ser", "rsfmri_var_ntpoints")
allvars <- c(proc_vars, ivs, dvs, re_covs, covs)
data <- data[c(allvars)]

#remove subjects with missing data
missing_dvs <- data[!complete.cases(data[c(dvs)]), ]
data <- data[complete.cases(data[c(dvs)]), ]
missing_ivs <- data[!complete.cases(data[c(ivs)]), ]
data <- data[complete.cases(data[c(ivs)]), ]
missing_covs <- data[!complete.cases(data[c(covs)]), ]
data <- data[complete.cases(data[c(covs)]), ]
missing_recovs <- data[!complete.cases(data[c(re_covs)]), ]
data <- data[complete.cases(data[c(re_covs)]), ]
data <- data[complete.cases(data), ]

#remove phillips
philips_exclude <- readLines("philips.txt")
data <- data[!(data$src_subject_id %in% philips_exclude), ]

######movementcensor#######

data_censor <- data[data$rsfmri_var_ntpoints < 375, ]
data <- data[data$rsfmri_var_ntpoints >= 375, ]
data_notenoughscans <- data[data$iqc_rsfmri_good_ser <= 1, ]
data <- data[data$iqc_rsfmri_good_ser > 1, ]
data_motionperTR <- data[data$iqc_rsfmri_all_mean_motion >= 0.9, ]
data <- data[data$iqc_rsfmri_all_mean_motion < 0.9, ]
fsqc_reject <- data[data$fsqc_qc != "accept", ]
data <- data[data$fsqc_qc == "accept", ]

#set categorical variables as factors
factor_ls <- c("rel_family_id", "ksads_14_853_p", "ksads_16_897_p", 
               "ksads_16_898_p", "ksads_1_840_p", "ksads_22_969_p")
data[factor_ls] <- lapply(data[factor_ls], as.factor)

#make split half datasets
data1 <- data[(data$src_subject_id %in% half1_subs),]
data2 <- data[(data$src_subject_id %in% half2_subs),]

#set DV, IVs, and COVs
dv <- "rsfmri_cor_network.gordon_default_network.gordon_dorsalattn"

#run mixed model function
stat_list <-(lapply(ivs, mixed_model, covs=covs, y=dv, data=data))
stat_matrix <- matrixify(stat_list)

#run mixed model function on split half
stat_list1 <-(lapply(ivs, mixed_model, covs=covs, y=dv, data=data1))
stat_matrix1 <- matrixify(stat_list1)

stat_list2 <-(lapply(ivs, mixed_model, covs=covs, y=dv, data=data2))
stat_matrix2 <- matrixify(stat_list2)

#write stat matrix to csv
write.csv(stat_matrix, "anticorrelation_results.csv")
write.csv(stat_matrix1, "anticorrelation_results_half1.csv")
write.csv(stat_matrix2, "anticorrelation_results_half2.csv")
