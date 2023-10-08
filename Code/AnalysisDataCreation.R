################################################################################
### Enter Data Paths ###########################################################
################################################################################

#  Metabolon excel file 
met_excel <- "data/UNAZ-0501-22VW_ DATA TABLES.xlsx"

#  Provide path to additional metadata (if no additional meta data the set = NULL)
additional_metaPath <- "data/AdditionalVars.xlsx" 


################################################################################
### Enter Metadata Variables For Analysis ######################################
################################################################################
metadata_variables <- c("PARENT_SAMPLE_NAME", 
                        "GROUP_NAME", 
                        "TIME1", 
                        "Gender") 


################################################################################
### Load Data ##################################################################
################################################################################

# Create MetPipe Object
dat <- loadMetPipeFromMetabolon(metabolon_path = met_excel)



# 2. Merge additional vars to the meta data
if(!is.null(additional_metaPath)){
  
  # Read in additional meta data 
  meta_data_additional <- read.xlsx(additional_metaPath) 
  
  meta_data <- meta_data_additional %>% 
    left_join(dat@meta,"PARENT_SAMPLE_NAME")
  
  # Update meta data slot
  dat@meta <- meta_data
}


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
# 1. Choose the meta data variable name.
var1 = "GROUP_NAME"

var2 = "TIME1" 

# 2. Choose a single variable to stratify the table by.
stratified_var = "Gender" 

# 3. Assign the variable a label name
label(dat@analysis[,var1]) <- "Treatment Group" 

label(dat@analysis[,var2]) <- "Time" 

# 4. Assign a stratified variable a label name for the table. 
label(dat@analysis[,stratified_var]) = "Gender" 

# 5. Create table 1
tbl1 <- table1(~ TIME1 + GROUP_NAME| Gender 
               , data = dat@analysis) 

# 6. Display table 1
print(tbl1)

################################################################################
### Save MetPipe Object ########################################################
################################################################################

# saveRDS(dat, file = "data/demo.Rds")


