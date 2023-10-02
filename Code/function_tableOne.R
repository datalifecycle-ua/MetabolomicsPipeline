## ----function_tableOne----
# 1. Creates table1
tbl1 <- table1(~ TIME1 + GROUP_NAME| Gender # <1>
               , data = analysis_data) #<2>
