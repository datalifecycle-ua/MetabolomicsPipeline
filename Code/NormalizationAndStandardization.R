###############################################################################
##### Enter Info ###############################################################
###############################################################################

# Provide a path to Metabolon .xlsx file. 
metabolon_path <- "Data/Metabolon/UNAZ-0501-22VW_ DATA TABLES.xlsx" 

# Output data path
Outfolder = "Data/Processed/"

###############################################################################
##### Load Data ###############################################################
###############################################################################

#  Load raw peak data
peak_data <- read.xlsx(metabolon_path, sheet = "Peak Area Data") # <2>



###############################################################################
#### Plots Before Standardization and Normalization ###########################
###############################################################################
# Select the first five metabolites for the box plot. 
metabolites <- names(peak_data)[2:6] 


#  Create the boxplot data
plot_data <- peak_data %>% 
  reshape2::melt(id.vars="PARENT_SAMPLE_NAME") %>% 
  filter(variable %in% metabolites) 


# Plot the box plots for each of the five metabolites. 
ggplot(plot_data,aes(x=variable,y=value)) + 
  geom_boxplot() + #<3>
  theme_bw() #<3>


###############################################################################
###### Standardization and Normalization ######################################
###############################################################################

# Median standardization
med_std <- median_standardization(peak_data = peak_data)

# Minimum Value imputation
impute <- min_val_impute(med_std)

# Log Transformation
log_trans_met <- log_transformation(impute)


###############################################################################
#### Plots After Standardization and Normalization ###########################
###############################################################################
metabolites <- names(peak_data)[2:6] 


# Create the boxplot data
plot_data <- log_trans_met %>% 
  reshape2::melt(id.vars="PARENT_SAMPLE_NAME") %>% 
  filter(variable %in% metabolites) 


#  Plot the box plots for each of the five metabolites. 
ggplot(plot_data,aes(x=variable,y=value)) + 
  geom_boxplot() + #<3>
  theme_bw() #<3>

###############################################################################
#### Save Data  ###############################################################
###############################################################################

write.csv(log_trans_met, file = paste0(Outfolder,"Log_transformed_data.csv"), 
          row.names = F) 

