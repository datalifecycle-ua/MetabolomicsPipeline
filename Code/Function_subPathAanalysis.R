## ----Function_subPathAanalysis----
# 1. Run the subpathway function
path_dat <- subpathway_analysis(analysis_data,#<1>
                            chem_data = chem_data,#<1>
                                  block_var = "TIME1",#<1>
                                  treat_var = "GROUP_NAME") #<1>
