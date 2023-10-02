## ----read_analysisAndChem----
# 1. Analysis Data
analysis_data <- read.csv(analysis_dat_path, check.names = F) #<1>


# 2. Read chemical annotations
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") #<2>
