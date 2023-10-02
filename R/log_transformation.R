log_transformation <- function(peak_data){
  
  # 1. Log transform all of the values
  peak_data_log <- peak_data %>% 
    mutate_if(is.numeric, log) 
  
  # Return log tranformed data
  return(peak_data_log)
}