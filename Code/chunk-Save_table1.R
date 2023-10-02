## ----Save_table1----
# 3. saves table 1 
t1flex(tbl1) %>% #<3>
  save_as_docx(path = paste0("../Outputs/Tables/","table1.docx")) #<3>
