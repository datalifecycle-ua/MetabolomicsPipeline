## ----Enter_tableOneLabels----
# 1. Choose the meta data variable name.
var1 = "GROUP_NAME" #<1>

var2 = "TIME1" # <1>


# 2. Choose a single variable to stratify the table by.
stratified_var = "Gender" # <2>


# 3. Assign the variable a label name
label(analysis_data[,var1]) <- "Treatment Group" #<3>

label(analysis_data[,var2]) <- "Time" #<3>


# 4. Assign a stratified variable a label name for the table. 
label(analysis_data[,stratified_var]) = "Gender" #<4>
