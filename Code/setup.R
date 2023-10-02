## ----setup----
# Set global chunk obptions
knitr::opts_chunk$set(fig.width=12, fig.height=12, warning = F)

# Data 
library(tidyverse)
library(openxlsx)
library(kableExtra)
library(table1)
library(flextable)
library(reshape2)


# Heat Map
library(pheatmap)
library(RColorBrewer)
library(rstatix)
library(plotly)

# PCA
library(factoextra)
library(FactoMineR)


## Combining plots
library(ggpubr)
library(grid)
library(gridExtra)

## Image package
library(magick)

## Pairwise comparisons
library(emmeans)

# Metabolomics Pipeline source files
files.sources = list.files("../R")
sapply(paste0("../R/",files.sources), source)
