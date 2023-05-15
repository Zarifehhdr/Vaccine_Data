# Load the data from the RDS file
install.packages("readRDS")
library(readRDS)
data_extended_old <- readRDS("~/Downloads/VaccineData/extendedOld_norm_batchCorrectedFromYoung_withResponse_eset.rds")
data_young <- readRDS("~/Downloads/VaccineData/young_norm_withResponse_eset.rds")
data_old <- readRDS("~/Downloads/VaccineData/old_norm_batchCorrectedFromYoung_withResponse_eset.rds")
#data_all < -readRDS("~/Downloads/VaccineData/all_norm_withResponse_eset.rds")

# View the data
head(my_data)