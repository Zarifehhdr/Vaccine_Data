######### Transcriptomic and Immune response data Subjects age 18-50, 
#####cross-study normalized and batch corrected expression
# Load the RDS file
data_young <- readRDS("/ix/djishnu/Zarifeh/ML_MM/young_norm_withResponse_eset.rds")
str(data_young)
young_X <- data_young@assayData$exprs  #Transcriptomic data of young cohort
dim(young_X) # 10086  2689
young_Y <- data_young@phenoData@data$ImmResp_postVax_value_MFC #Immune Response of young cohort

######### Transcriptomic and Immune response data Subjects age 60-90
# Load the RDS file
data_old <- readRDS("/ix/djishnu/Zarifeh/ML_MM/old_norm_batchCorrectedFromYoung_withResponse_eset.rds")
str(data_old)
old_X <- data_old@assayData$exprs  #Transcriptomic data of old cohort
dim(old_X) # 14016   859
old_Y <- data_old@phenoData@data$ImmResp_postVax_value_MFC #Immune Response of old cohort

######### Transcriptomic and Immune response data Subjects age 50-90
# Load the RDS file
data_extended_old <- readRDS("/ix/djishnu/Zarifeh/ML_MM/extendedOld_norm_batchCorrectedFromYoung_withResponse_eset.rds")
str(data_extended_old)
extended_old_X <- data_extended_old@assayData$exprs  #Transcriptomic data of extended old cohort
dim(extended_old_X) # 14016   938
extended_old_Y <- data_extended_old@phenoData@data$ImmResp_postVax_value_MFC #Immune Response of extended old cohort

######### Transcriptomic and Immune response data of all participants
# Load the RDS file
data_all <- readRDS("/ix/djishnu/Zarifeh/ML_MM/all_norm_withResponse_eset.rds")
str(data_all)
all_X <- data_all@assayData$exprs  #Transcriptomic data of extended all cohort
dim(all_X) # 10086   4104
all_Y <- data_all@phenoData@data$ImmResp_postVax_value_MFC #Immune Response of all

###########################First start to work on young and old cohorts
#######################Young cohort
################## Removing the baseline from data 
# Get the column names from the exprs matrix
col_names <- colnames(young_X)
# Identify the indices of columns to remove that correspond to day 0
indices_to_remove <- grep("_0_Days_", col_names)
# Remove the columns corresponding to day 0 from the expression and immune response data
young_X <- young_X[, -indices_to_remove]
dim(young_X) # 10086  1999
young_Y <- young_Y[-indices_to_remove]
length(young_Y) #1999
#############
# Load the SummarizedExperiment package
library(SummarizedExperiment)
# Get the study accessions from the phenoData metadata
study_accessions <- data_young@phenoData@data[["study_accession"]]
# Count the number of occurrences of each study accession
accession_counts <- table(study_accessions)
# Print the accession counts
print(accession_counts)

#study_accessions for young 
#SDY1119 SDY1260 SDY1264 SDY1276 SDY1289 SDY1294 SDY1325 SDY1328 SDY1364 SDY1370 SDY1529  SDY180  SDY212  SDY224 
#67      90      87     816      84     109       4      50      36      48     180     156      28      55 
#SDY269  SDY270  SDY400  SDY404  SDY520   SDY56   SDY61   SDY63  SDY640   SDY80  SDY984 
#163      83      60      64      51      30      27      42      44     251      64 
# Get the study accessions from the phenoData metadata
study_accessions <- data_old@phenoData@data[["study_accession"]]
# Count the number of occurrences of each study accession
accession_counts <- table(study_accessions)
# Print the accession counts
print(accession_counts)

#study_accessions for old
#SDY1119 SDY1325 SDY1328  SDY212  SDY400  SDY404  SDY520   SDY56   SDY63  SDY640   SDY80  SDY984 
#66       4     270      60       60      92      43        S118      30      35       5      76 
study_accessions <- data_extended_old@phenoData@data[["study_accession"]]
# Count the number of occurrences of each study accession
accession_counts <- table(study_accessions)
Print the accession counts
#study_accessions of extended old
#SDY1119 SDY1325 SDY1328  SDY212  SDY400  SDY404  SDY520   SDY56   SDY63  SDY640   SDY80  SDY984 
#110      14     270      60      60      92      43       118      30      35      30      76 
############################
# two flu datasets: SDY1119 and SDY400, SDY212, SDY400, SDY404, SDY56, SDY520, SDY63, SDY640, SDY80
#Zoster : SDY984
# Meningococcal:SDY1325 
# Hepatit B: SDY1328

# Loop over each unique study accession
######################### Providing the X and Y data for the selected Studies
sub_young_X <- t(young_X[, grep("\\.1119_|\\.400_|\\.984_|\\.1325_|\\.1328_|\\.212_|\\.404_|\\.56_|\\.520_|\\.63_|\\.80_", colnames(young_X))])
sub_young_Y <- young_Y[grep("\\.1119_|\\.400_|\\.984_|\\.1325_|\\.1328_|\\.212_|\\.404_|\\.56_|\\.520_|\\.63_|\\.80_", colnames(young_X))]
# Making the Ysubset a matrix
sub_young_Y <- matrix(sub_young_Y, ncol = 1, dimnames = list(rownames(sub_young_X), "names"))

# Dimension of the subset matrix
dim(sub_young_X) #503 10086
dim(sub_young_Y) #503 1
#### Zero Filtering :::::: No genes to remove
mt_cols1 <- grep("^MT-",colnames(sub_young_X),value = TRUE) # There is no mitochondrial gene
mt_cols2 <- grep("^RP[LS]",colnames(sub_young_X),value = TRUE) # There are 58 ribosomal genes 
# Removing the ribosomal genes
sub_young_X <- sub_young_X[, !colnames(sub_young_X) %in% mt_cols2] #503 10028
write.csv(sub_young_X, file = "/ix/djishnu/Zarifeh/ML_MM/X_young.csv", row.names = TRUE)
write.csv(sub_young_Y, file = "/ix/djishnu/Zarifeh/ML_MM/Y_young.csv", row.names = TRUE)

#######################Old cohort
################## Removing the baseline from data 
# Get the column names from the exprs matrix
col_names <- colnames(old_X)
# Identify the indices of columns to remove that correspond to day 0
indices_to_remove <- grep("_0_Days_", col_names)
# Remove the columns corresponding to day 0 from the expression and immune response data
old_X <- old_X[, -indices_to_remove]
dim(old_X) # 14016   524
old_Y <- old_Y[-indices_to_remove]
length(old_Y) #524
#############
######################### Providing the X and Y data for the selected Studies
sub_old_X <- t(old_X[, grep("\\.1119_|\\.400_|\\.984_|\\.1325_|\\.1328_|\\.212_|\\.404_|\\.56_|\\.520_|\\.63_|\\.80_", colnames(old_X))])
sub_old_Y <- old_Y[grep("\\.1119_|\\.400_|\\.984_|\\.1325_|\\.1328_|\\.212_|\\.404_|\\.56_|\\.520_|\\.63_|\\.80_", colnames(old_X))]
# Making the Ysubset a matrix
sub_old_Y <- matrix(sub_old_Y, ncol = 1, dimnames = list(rownames(sub_old_X), "names"))

# Dimension of the subset matrix
dim(sub_old_X) #498 14016
dim(sub_old_Y) #498 1
#### Zero Filtering :::::: No genes to remove

mt_cols1 <- grep("^MT-",colnames(sub_old_X),value = TRUE) # There is no mitochondrial gene
mt_cols2 <- grep("^RP[LS]",colnames(sub_old_X),value = TRUE) # There are 74 ribosomal genes 
# Removing the ribosomal genes
sub_old_X <- sub_old_X[, !colnames(sub_old_X) %in% mt_cols2] #498 13943
write.csv(sub_old_X, file = "/ix/djishnu/Zarifeh/ML_MM/X_old.csv", row.names = TRUE)
write.csv(sub_old_Y, file = "/ix/djishnu/Zarifeh/ML_MM/Y_old.csv", row.names = TRUE)

################### Concat the young and old data for selected studies
### Some shared genes are available

###############################
# Read the csv files as data frames
X1 <- read.csv("/ix/djishnu/Zarifeh/ML_MM/X_young.csv", row.names = 1, header = TRUE) #503 10028 
X2 <- read.csv("/ix/djishnu/Zarifeh/ML_MM/X_old.csv", row.names = 1, header = TRUE)   #498 13943

# Find the common column names between X1 and X2
common_cols <- intersect(colnames(X1), colnames(X2))   # 10028 shared genes

# Remove common columns from X1 and X2
X1_new <- X1[, setdiff(colnames(X1), common_cols)]
X2_new <- X2[, setdiff(colnames(X2), common_cols)]

# Merging not shared columns
X = merge(X1_new, X2_new, all=TRUE, by=0)
dim(X)   #1001 3916
# Set the Row.names column as row names
rownames(X) <- X$Row.names
X$Row.names <- NULL

for (col in 1:length(common_cols)) {
  X <- cbind(X, 0)
}
dim(X)   # 1001 13943

colnames(X)[(ncol(X)-length(common_cols)+1):ncol(X)] <- paste0(common_cols)

# Inserting shared genes
v=rbind(X1[,common_cols], X2[,common_cols])
v
# get the row indices of a that match the row names of v
indices <- match(rownames(X), rownames(v))

# reorder the rows of v
v_sorted <- v[indices, ]
X[,common_cols] <- v_sorted

X[is.na(X)] <- 0
#####Writing CSV
write.csv(X, file = "/ix/djishnu/Zarifeh//ML_MM/X_young_old_concat.csv", row.names = TRUE, col.names = TRUE)

#####################################################
Y1 <- sub_young_Y
Y2 <- sub_old_Y
# Combine Y matrices  
Y_young_old_concat <- rbind(Y1, Y2)

# View the resulting matrix
dim(Y_young_old_concat)

write.csv(Y_young_old_concat, file = "/ix/djishnu/Zarifeh//ML_MM/Y_young_old_concat.csv", row.names = TRUE)



