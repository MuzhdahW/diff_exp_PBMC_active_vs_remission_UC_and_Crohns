install.packages("dplyr")
library(dplyr)
install.packages("Matrix")
library(Matrix)
install.packages("Seurat")
library(Seurat)
install.packages("patchwork")
library(patchwork)
install.packages("devtools")
install.packages("dplyr")
library(dplyr)
install.packages("tibble")
library(tibble)
#devtools::install_github("immunogenomics/presto")
#install.packages("presto")
#library(presto)
#install.packages("GEOquery")
#library(GEOquery)
# Since the data being utilized is a pre-made RDS object, I am simply going to use CURL to open the file
# and move directly to secondary analysis. 
# geoquery did not exist for the newest version of R, so I decided to just load the file from my local machine
#pbmc_rds <- readRDS(url("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330122/matrix/"))
#Finding differentially expressed features 
data_file <- "diff_exp_of_PBMCs_in_IBD_patients_UC_patients_and_healthy_controls/data/GSE330122_integrated_seurat_object_all12donors.rds"
pbmc_rds <- readRDS(data_file)

# let's clear out the bad reads by removing the NA's. We are looking at the names of the rows in the meta.data table that are complete
pbmc_complete_rows <- rownames(pbmc_rds@meta.data)[complete.cases(pbmc_rds@meta.data)]
# and we are using those complete rows to subset the rds file, but in the final file this will be removing columns, not rows.
# this is because cells are rows in the meta data table and columns in the actual rds
cleaned_pbmc <- subset(pbmc_rds,cells=pbmc_complete_rows)
# now swe split the data based on ulcerative_colitis(UC) or crohns(CD). We will then compare active vs remission
# disease



# 2. Get Crohn's rows by checking for barcodes that contain "CD" and then subsetting the rds file based on those 
pbmc_crohns_rows <- rownames(cleaned_pbmc@meta.data)[grepl("CD", cleaned_pbmc@meta.data$barcode)]
crohns_subset <- subset(cleaned_pbmc, cells = pbmc_crohns_rows)

# 3. Get UC rows the same way
pbmc_uc_rows <- rownames(cleaned_pbmc@meta.data)[grepl("UC", cleaned_pbmc@meta.data$barcode)]
uc_subset <- subset(cleaned_pbmc, cells = pbmc_uc_rows)

# now we see how many levels are in the orig.ident of each subset
print(unique(crohns_subset@meta.data$seurat_clusters))
print(unique(uc_subset@meta.data$seurat_clusters))
# we have twelve cell clusters. These, however, are independent from "active" or "remission" disease.
# in order to compare "active" or "remission" status within a cell cluster ( ie cluster 0 active vs remission or cluster 1 active vs remission)
# I need to create a third column in the metadata table. I'm going to call this column "active disease"
# and list it as "true" or "false". I will then combine that w the factor column, so I can compare
# cluster 0 active vs cluster 0 remission, for example


# just to play around for now I'm going to copy the dataframes. I am going to use dplyr to 
# create the new columns. but I need to use rownames_to_column and column_to_rowname to preserve the 
# barcodes, because the barcodes are the rownames from the metadata table. 

copy_crohns <- crohns_subset
copy_uc <- uc_subset

# using dplyr to add in a column that has " active " or "remission" so I can compare active vs remission
# within a given cluster more easily!We can call that column " status" and have it be " active " or "remission" 
copy_uc@meta.data <- copy_uc@meta.data %>% 
rownames_to_column(var = "cell_barcodes")%>% 
mutate(status = if_else(grepl('Active',barcode),"active","remission")) %>%
column_to_rownames(var = "cell_barcodes") 


copy_crohns@meta.data <- copy_crohns@meta.data %>% 
  rownames_to_column(var = "cell_barcodes")%>% 
  mutate(status = if_else(grepl('Active',barcode),"active","remission")) %>%
  column_to_rownames(var = "cell_barcodes")

