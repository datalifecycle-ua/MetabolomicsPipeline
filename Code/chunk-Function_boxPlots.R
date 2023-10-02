## ----Function_boxPlots----
# 1. Read in analysis data
analysis_data <- read.csv(analysis, check.names = F) #<1>

# 2. Reach in chemical annotation data
# Read chemical annotation file
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") # <2>


# 3. Run the subpathway_boxplots function
box_1 <- subpathway_boxplots(analysis_data =  analysis_data, #<3>
                chem_data = chem_data, #<3>
                subpathway = subpath, #<3>
                X = TIME1, #<3>
                groupBy = GROUP_NAME, #<3>
                # Additional analysis data filtering options. 
                Gender == "Male", GROUP_NAME != "Control") #<3>

# 4. Customize labels for the boxplots
box_1 + #<4>
  labs(x = "Time", y="Scaled Intensity", color="Treatment Group", #<4>
       title = "Metabolite boxplots for the Gamma-glutamyl Amino Acid Subpathway") #<4>

# 5. Save boxplots. You can change the file type for example change pdf to png. 
ggsave(filename = paste0(plot_save,"/Plots/Metabolite_boxplots_",subpath,".pdf")) #<5>
