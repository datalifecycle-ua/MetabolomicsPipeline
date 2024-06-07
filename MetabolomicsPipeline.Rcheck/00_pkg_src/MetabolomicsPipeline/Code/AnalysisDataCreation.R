# Load required packages for analysis
source("Code/Setup.R")

################################################################################
### Load Data ##################################################################
################################################################################

# Load Metabolon data from the Excel file
dat <- loadMetPipeFromMetabolon(metabolon_path = "data/UNAZ-0501-22VW_ DATA TABLES.xlsx" )

# load additional metadata Set as NULL is not applicable.
meta_data_additional <- read.xlsx("data/AdditionalVars.xlsx")


# 2. Merge additional vars to the meta data
if(nrow(meta_data_additional)>0){
  
  meta_data <- meta_data_additional %>% 
    left_join(dat@meta,"PARENT_SAMPLE_NAME")
  
  # Update meta data slot
  dat@meta <- meta_data
}


################################################################################
### Enter Metadata Variables For Analysis ######################################
################################################################################
metadata_variables <- c("PARENT_SAMPLE_NAME", 
                        "GROUP_NAME", 
                        "TIME1", 
                        "Gender") 


################################################################################
### Add Analysis Data To MetPipe Object ########################################
################################################################################
# 2. Create analysis data
dat@analysis <- dat@meta %>% 
  select(all_of(metadata_variables)) %>% 
  left_join(dat@standardized_peak,"PARENT_SAMPLE_NAME") 


################################################################################
### Create Table 1 #############################################################
################################################################################
# 5. Create table 1
tbl1 <- table1(~ TIME1 + GROUP_NAME| Gender 
               , data = dat@analysis) 

# 6. Display table 1
tbl1


# Save data
saveRDS(dat, file = "data/dat.Rds")

