## ----Function_dataMutation----
# 1. Read in analysis data
analysis_data <- read.csv(analysis, check.names = F) #<1>

# 2. Reach in chemical annotation data
# Read chemical annotation file
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") #<2>


# 3. Prep data for line plots
analysis_data <- analysis_data %>% #<3>
  filter(Gender=="Male") %>% #<3>
  mutate(TIME1 = case_when(TIME1== "PreSymp" ~ 1, #<3>
                           TIME1 =="Onset" ~ 2, #<3>
                       TIME1 == "End" ~ 3)) #<3>
