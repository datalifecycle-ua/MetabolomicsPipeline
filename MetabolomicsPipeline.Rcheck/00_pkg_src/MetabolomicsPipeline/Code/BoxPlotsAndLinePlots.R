# Load required packages for analysis
source("Code/Setup.R")

# Load Data from AnalysisDataCreation.R
dat <- readRDS("data/dat.Rds")

################################################################################
### BoxPlots ###################################################################
################################################################################

subpathway_boxplots(dat, subpathway = "Lactoyl Amino Acid", block_var=TIME1,
                    treat_var = GROUP_NAME, Gender =="Female")


################################################################################
## Line plots ##################################################################
################################################################################

# Set up data
dat@analysis$TIME1 <- as.numeric(factor(dat@analysis$TIME1,
                                        levels = c("PreSymp","Onset","End")))

# Create line plots 
subpathway_lineplots(dat, subpathway = "Lactoyl Amino Acid",
                     block_var= TIME1, treat_var = GROUP_NAME, Gender=="Female" )
