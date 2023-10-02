## ----Function_linePlots----
# 4. Run subpathway_line plots function.
line_plots <- subpathway_lineplots(analysis_data = analysis_data, #<4>
                     chem_data = chem_data, #<4>
                     subpathway = subpath, #<4>
                     X = TIME1, #<4>
                     groupBy = GROUP_NAME) #<4>


# 5. Edit asthetics of plots
line_plots + #<5>
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("PreSymp", "Onset", "End")) +#<5>
  labs(x= "Time",y="Scaled Intensity",color = "Treatment Group")#<5>


# 6. Save plots 
ggsave(filename = paste0(plot_save ,"/LinePlots_",subpath,".pdf")) #<6>
