## ----readRawData----
# 2. Read in sample metadata
meta_data <- read.xlsx(met_excel,sheet = "Sample Meta Data") # <2>


# 3. Read Log transformed data
peak_data_log <- read.xlsx(met_excel,sheet = "Log Transformed Data") #<3>
