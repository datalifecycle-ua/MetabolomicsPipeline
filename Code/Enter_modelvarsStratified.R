## ----Enter_modelvarsStratified----
# Enter block variable name
blockvar = "TIME1" #<1>

# Enter treatment variable
treatvar = "GROUP_NAME" #<1>


# 1. Name of the variable to stratify by.
stratified_var = "Gender" # <2>

# 3. Find unique values of this variable
uni_vals <- unique(analysis_data[,stratified_var]) #<3>
