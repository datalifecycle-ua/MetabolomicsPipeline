## ----runPCA----
# 3. Read in analysis data
analysis_data <- read.csv(analysis_data_path, check.names = F) #<3>


# 4. Create PCA data containing only metabolite data
pca_dat <- analysis_data[,!names(analysis_data) %in% non_metabolite] #<4>


# 5. Run PCA of the pca_dat matrix containing only the metabolites.   
res.pca <- PCA(pca_dat, #<5>
               graph = F)#<5>


# 6. Create figure 
fviz_pca_ind(res.pca, #<6>
             label = "none",#<6>
             habillage = as.factor(analysis_data[,meta_var])) #<6>


# 7. Save figure
ggsave(paste0(out,"/PCA.pdf"), width = 10, height = 10)#<7>
