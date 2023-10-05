################################################################################
### Enter Data Paths ###########################################################
################################################################################

#  Metabolon excel file 
met_excel <- "Data/Metabolon/UNAZ-0501-22VW_ DATA TABLES.xlsx"

#  Provide path to additional metadata (if no additional meta data the set = NULL)
additional_metaPath <- "Data/Metabolon/AdditionalVars.xlsx" 

# Provide path for analysis data output
data_output = "Data/Processed/"

# Provide path for table 1 output
table_output <- "Outputs/Tables/"
################################################################################
### Enter Metadata Variables ###################################################
################################################################################
metadata_variables <- c("PARENT_SAMPLE_NAME", 
                        "GROUP_NAME", 
                        "TIME1", 
                        "Gender") 


################################################################################
### Load Data ##################################################################
################################################################################

# Read in sample metadata
meta_data <- read.xlsx(met_excel,sheet = "Sample Meta Data") 


# Read Log transformed data
peak_data_log <- read.xlsx(met_excel,sheet = "Log Transformed Data") 


# 2. Merge additional vars to the meta data
if(!is.null(additional_metaPath)){
  
  # Read in additional meta data 
  meta_data_additional <- read.xlsx(additional_metaPath) 
  
  meta_data <- meta_data_additional %>% 
    left_join(meta_data,"PARENT_SAMPLE_NAME") 
}


################################################################################
### Merge Peak and Metadata ####################################################
################################################################################
# 2. Create analysis data
analysis_data <- meta_data %>% 
  select(all_of(metadata_variables)) %>% 
  left_join(peak_data_log,"PARENT_SAMPLE_NAME") 


################################################################################
### Save Analysis Data #########################################################
################################################################################
# Save analysis data
write.csv(analysis_data,paste0(data_output,"analysis_data.csv"), row.names = F) 


################################################################################
### Create Table 1 #############################################################
################################################################################
# 1. Choose the meta data variable name.
var1 = "GROUP_NAME"

var2 = "TIME1" 

# 2. Choose a single variable to stratify the table by.
stratified_var = "Gender" 

# 3. Assign the variable a label name
label(analysis_data[,var1]) <- "Treatment Group" 

label(analysis_data[,var2]) <- "Time" 

# 4. Assign a stratified variable a label name for the table. 
label(analysis_data[,stratified_var]) = "Gender" 

# 5. Create table 1
tbl1 <- table1(~ TIME1 + GROUP_NAME| Gender 
               , data = analysis_data) 

# 6. Display table 1
print(tbl1)

#7. Save table 1
t1flex(tbl1) %>% #<3>
  save_as_docx(path = paste0(table_output,"table1.docx"))


