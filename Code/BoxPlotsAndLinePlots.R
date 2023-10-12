################################################################################
### BoxPlots ###################################################################
################################################################################

subpathway_boxplots(dat, subpathway = "Lactoyl Amino Acid", X=TIME1,
                    groupBy = GROUP_NAME, Gender =="Female")


################################################################################
## Line plots ##################################################################
################################################################################

# Set up data
dat@analysis$TIME1 <- as.numeric(factor(dat@analysis$TIME1,
                                        levels = c("PreSymp","Onset","End")))

# Create line plots 
subpathway_lineplots(dat, subpathway = "Lactoyl Amino Acid",
                     X= TIME1,groupBy = GROUP_NAME, Gender=="Female" )