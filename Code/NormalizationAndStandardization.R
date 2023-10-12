###############################################################################
##### Load Data ###############################################################
###############################################################################
# Provide a path to Metabolon .xlsx file. 
metabolon_path <- "../data/UNAZ-0501-22VW_ DATA TABLES.xlsx" 

#  Load raw peak data
peak_data <- read.xlsx(metabolon_path, sheet = "Peak Area Data") 


###############################################################################
###### Standardization and Normalization ######################################
###############################################################################

# 1 Median standardization
med_std <- median_standardization(peak_data = peak_data)

# 2 Minimum Value imputation
impute <- min_val_impute(med_std)

# 3 Log Transformation
log_trans_met <- log_transformation(impute)