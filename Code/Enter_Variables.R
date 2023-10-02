## ----Enter_Variables----
# 1. Define the variable in the analysis data that we want to stratify 
strata_var <- "Gender" #<1>


# 2. Define the form variable
form <- "GROUP_NAME*TIME1" #<2>



# Enter Non-metabolite variables
non_metabolite <- c("PARENT_SAMPLE_NAME", "GROUP_NAME", #<3>
                                           "TIME1", "Gender") #<3>

metabolites <- names(analysis_data)[!names(analysis_data) %in% non_metabolite] #<4>
