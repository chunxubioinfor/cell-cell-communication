library(rjson)
library(jsonlite)
library(gdata)
library(dplyr)
library(tidyr)
library(Ipaper)
setwd('~/Desktop/CCI_DK/ccc/')
options(stringsAsFactors = FALSE)
# Curation of KEGG pathways
gsea_kegg <- fromJSON("./kegg.json",simplifyDataFrame = TRUE)
kegg_symbols <- names(gsea_kegg)
kegg_ids <- c()
for (i in 1:length(kegg_symbols)){
  kegg_ids <- c(kegg_ids,gsea_kegg[[kegg_symbols[i]]]$exactSource)
}
kegg_symbol_id <- data.frame(kegg_symbols,kegg_ids)
write.csv(kegg_symbol_id,'./tmp.csv')

# Curation of WIKI pathway
gsea_wiki <- fromJSON("./wiki.json",simplifyDataFrame = TRUE)
wiki_symbols <- names(gsea_wiki)
wiki_ids <- c()
for (i in 1:length(wiki_symbols)){
  wiki_ids <- c(wiki_ids,gsea_wiki[[wiki_symbols[i]]]$exactSource)
}
wiki_symbol_id <- data.frame(wiki_symbols,wiki_ids)
write.csv(wiki_symbol_id,'./tmp.csv')


# Then manually curate the KEGG pathway and WIKi pathway
# This step took me around two weeks
# The results were saved in curation_pathway.csv
# Read into the file and reshape the dataset
curation_pathway_kegg <- read.xls('./curation_pathway.xlsx',sheet = 1,na.strings=c("NA","#DIV/0!",''))
curation_pathway_kegg <- curation_pathway_kegg %>% filter(receptor != 'NA') %>%
  select(gsea_symbol,id,receptor)  ## Filter out the invalid rows 
curation_pathway_wiki <- read.xls('./curation_pathway.xlsx',sheet = 2,na.strings=c("NA","#DIV/0!",''))
curation_pathway_wiki <- curation_pathway_wiki %>% filter(receptor != 'NA') %>%
  select(gsea_symbol,id,receptor)
dim(curation_pathway_kegg)
dim(curation_pathway_wiki)
print(paste('The number of valid pathways is',nrow(curation_pathway_kegg)+nrow(curation_pathway_wiki)))

# Reshape the dateset
curation_pathway_kegg <- curation_pathway_kegg %>% as_tibble() %>% separate_rows(receptor,sep = ',')
curation_pathway_wiki <- curation_pathway_wiki %>% as_tibble() %>% separate_rows(receptor,sep = ',')
curation_pathway <- rbind(curation_pathway_kegg,curation_pathway_wiki)

# Filter out the non-receptor (maybe it is caused by the curation mistakes)
# 
cor_matrix_spearman <- readRDS("~/Desktop/CCI_DK/ccc/cor_matrix_spearman.rds")
receptor_ref <- colnames(cor_matrix_spearman)
curation_pathway_filtered <- filter(curation_pathway,receptor %in% receptor_ref)
write.csv(curation_pathway_filtered,'./curation_pathway_filtered.csv')

# Calculation of average ranking
# Define a function to calculate the average ranking given a correlation data.frame and a validation set

ranking_cal <- function(ranking_df,curation_pathway,file_name){
  pb <- txtProgressBar(style=3)
  average_ranking <- c()
  for (i in 1:ncol(ranking_df)){
    setTxtProgressBar(pb, i/ncol(ranking_df))
    receptor_s <- colnames(ranking_df)[i]
    receptor_curation_df <- filter(curation_pathway,receptor == receptor_s)
    pathway_2_receptor <- unique(receptor_curation_df$gsea_symbol)
    ranking <- c()
    for (j in 1:length(pathway_2_receptor)){
      pathway <- pathway_2_receptor[j]
      ranking <- c(ranking,nrow(ranking_df) +1 - rank(ranking_df[,i])[pathway])
    }
    ranking <- mean(ranking)
    average_ranking <- c(average_ranking,ranking)
  }
  print(mean(average_ranking,na.rm = TRUE))
  average_ranking_df <- data.frame(receptor = colnames(ranking_df),average_ranking)
  write.csv(average_ranking_df,file = paste('./',file_name,'.csv',sep = ''))
  p <- ggplot(average_ranking_df, aes(x = average_ranking)) +
    geom_density(color = 'black', fill = 'gray') +
    geom_vline(xintercept = mean(average_ranking,na.rm = TRUE),color = 'red')
  p
  h <- ggplot(average_ranking_df,aes(x = average_ranking)) + 
    geom_histogram(bins = 100) + 
    geom_vline(xintercept = mean(average_ranking,na.rm = TRUE),color = 'red') +
    geom_vline(xintercept = 10,color = 'blue')
  write_fig(h,paste('./',file_name,'.png',sep = ''))
  return(average_ranking)
  
  close(pb)
}

## 1. The valid pathways 
receptor_val <- intersect(receptor_ref,unique(curation_pathway$receptor))
matrix_spearman_kegg_wiki <- cor_matrix_spearman[unique(curation_pathway$gsea_symbol),receptor_val]
ranking_cal(matrix_spearman_kegg_wiki,curation_pathway_filtered,file_name = 'average_ranking_hist_1')

## 2. All KEGG and WIKI pathways
matrix_spearman_kegg_wiki_all <- cor_matrix_spearman[c(kegg_symbols,wiki_symbols),receptor_val]
ranking_cal(matrix_spearman_kegg_wiki_all,curation_pathway_filtered,file_name = 'average_ranking_hist_2')

## 3. All the pathways
matrix_spearman <- cor_matrix_spearman[,receptor_val]
ranking_cal(matrix_spearman,curation_pathway_filtered,file_name = 'average_ranking_hist_3')

## Make a attempt to use median rather than mean
# Define a function to calculate the average ranking given a correlation data.frame and a validation set
ranking_cal <- function(ranking_df,curation_pathway,file_name){
  pb <- txtProgressBar(style=3)
  average_ranking <- c()
  for (i in 1:ncol(ranking_df)){
    setTxtProgressBar(pb, i/ncol(ranking_df))
    receptor_s <- colnames(ranking_df)[i]
    receptor_curation_df <- filter(curation_pathway,receptor == receptor_s)
    pathway_2_receptor <- unique(receptor_curation_df$gsea_symbol)
    ranking <- c()
    for (j in 1:length(pathway_2_receptor)){
      pathway <- pathway_2_receptor[j]
      ranking <- c(ranking,nrow(ranking_df) +1 - rank(ranking_df[,i])[pathway])
    }
    ranking <- median(ranking,na.rm = TRUE)
    average_ranking <- c(average_ranking,ranking)
  }
  print(median(average_ranking,na.rm = TRUE))
  average_ranking_df <- data.frame(receptor = colnames(ranking_df),average_ranking)
  write.csv(average_ranking_df,file = paste('./',file_name,'.csv',sep = ''))
  p <- ggplot(average_ranking_df, aes(x = average_ranking)) +
    geom_density(color = 'black', fill = 'gray') +
    geom_vline(xintercept = mean(average_ranking,na.rm = TRUE),color = 'red')
  p
  h <- ggplot(average_ranking_df,aes(x = average_ranking)) + 
    geom_histogram(bins = 100) + 
    geom_vline(xintercept = median(average_ranking,na.rm = TRUE),color = 'red')
  write_fig(h,paste('./',file_name,'.png',sep = ''))
  return(average_ranking)
  close(pb)
}
ranking_cal(matrix_spearman_kegg_wiki,curation_pathway_filtered,file_name = 'median_ranking_hist_1')
ranking_cal(matrix_spearman_kegg_wiki_all,curation_pathway_filtered,file_name = 'median_ranking_hist_2')
ranking_cal(matrix_spearman,curation_pathway_filtered,file_name = 'median_ranking_hist_3')

# The performance is quite great
# At the pathway point of view
# Define a function to calculate the average ranking given a correlation data.frame and a validation set
ranking_cal <- function(ranking_df,curation_pathway,file_name,focus){
  pb <- txtProgressBar(style=3)
  average_ranking <- c()
  if (focus == 'receptor'){
    for (i in 1:ncol(ranking_df)){
      setTxtProgressBar(pb, i/ncol(ranking_df))
      receptor_s <- colnames(ranking_df)[i]
      receptor_curation_df <- filter(curation_pathway,receptor == receptor_s)
      pathway_2_receptor <- unique(receptor_curation_df$gsea_symbol)
      ranking <- c()
      for (j in 1:length(pathway_2_receptor)){
        pathway <- pathway_2_receptor[j]
        ranking <- c(ranking,nrow(ranking_df) +1 - rank(ranking_df[,i])[pathway])
      }
      ranking <- median(ranking,na.rm = TRUE)
      average_ranking <- c(average_ranking,ranking)
    }
    print(median(average_ranking,na.rm = TRUE))
    average_ranking_df <- data.frame(receptor = colnames(ranking_df),average_ranking)
    write.csv(average_ranking_df,file = paste('./',file_name,'.csv',sep = ''))
    h <- ggplot(average_ranking_df,aes(x = average_ranking)) + 
      geom_histogram(bins = 100) + 
      geom_vline(xintercept = median(average_ranking,na.rm = TRUE),color = 'red')
    write_fig(h,paste('./',file_name,'.png',sep = ''))
    return(average_ranking)
    close(pb)
    }
    else if (focus == 'pathway'){
      for (i in 1:nrow(ranking_df)){
        setTxtProgressBar(pb, i/nrow(ranking_df))
        pathway_s <- rownames(ranking_df)[i]
        pathway_curation_df <- filter(curation_pathway,gsea_symbol == pathway_s)
        receptor_2_pathway <- unique(pathway_curation_df$receptor)
        ranking <- c()
        for (j in 1:length(receptor_2_pathway)){
          receptor <- receptor_2_pathway[j]
          ranking <- c(ranking,ncol(ranking_df)+1 - rank(ranking_df[i,])[receptor])
        }
        ranking <- median(ranking,na.rm = TRUE)
        average_ranking <- c(average_ranking,ranking)
      }
      print(median(average_ranking,na.rm = TRUE))
      average_ranking_df <- data.frame(pathway = rownames(ranking_df),average_ranking)
      write.csv(average_ranking_df,file = paste('./',file_name,'.csv',sep = ''))
      h <- ggplot(average_ranking_df,aes(x = average_ranking)) + 
        geom_histogram(bins = 100) + 
        geom_vline(xintercept = median(average_ranking,na.rm = TRUE),color = 'red')
      write_fig(h,paste('./',file_name,'.png',sep = ''))
      return(average_ranking)
      close(pb)
    }
  }
receptor_val <- intersect(receptor_ref,unique(curation_pathway$receptor))
matrix_spearman_kegg_wiki <- cor_matrix_spearman[unique(curation_pathway$gsea_symbol),receptor_val]
ranking_cal(matrix_spearman_kegg_wiki,curation_pathway_filtered,file_name = 'median_ranking_hist_p_1',focus = 'pathway')
ranking_cal(cor_matrix_spearman[unique(curation_pathway$gsea_symbol),],curation_pathway_filtered,file_name = 'median_ranking_hist_p_2',focus = 'pathway')

# Try applying different algorithms

# For the pathway
# Randomly divide the curation dataset into two groups, training data (2/3) and testing data (1/3)
all_pathway_kegg <- unique(curation_pathway_kegg$gsea_symbol)
all_pathway_wiki <- unique(curation_pathway_wiki$gsea_symbol)
testing_pathway_kegg <- sample(all_pathway_kegg,trunc(length(all_pathway_kegg)/3))
testing_pathway_wiki <- sample(all_pathway_wiki,trunc(length(all_pathway_wiki)/3))
testing_set <- c(testing_pathway_kegg,testing_pathway_wiki)  # testing data,just leave alone
training_pathway_kegg <- setdiff(all_pathway_kegg,testing_pathway_kegg)
training_pathway_wiki <- setdiff(all_pathway_wiki,testing_pathway_wiki)
training_set <- setdiff(unique(curation_pathway$gsea_symbol),testing_set)

# Cross validation 
# Define a function which group the training data into train and validation set for cross-validation step
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}

# Define a function which calculate the performance score in the ranking method
top_score <- function(ranking_df,curation_pathway,pathway,ranking_thr,file_name){
  score_mtx <- c()
  # Extract the top receptors to the certain pathway from the matrix
  top_list <- names(sort(ranking_df[pathway,],decreasing = TRUE)[1:ranking_thr])
  print(top_list)
  # Extract the curated receptors as reference
  receptor_curation <- filter(curation_pathway,gsea_symbol == pathway)$receptor
  print(receptor_curation)
  tp <- length(intersect(top_list,receptor_curation))
  fp <- length(top_list) - length(intersect(top_list,receptor_curation))
  tn <- ncol(ranking_df) - length(receptor_curation) - fp
  sensitivity <- tp/length(receptor_curation)
  specificity <- tn/(ncol(ranking_df) - length(receptor_curation))
  FDR <- fp/length(top_list)
  score_mtx <- c(sensitivity,specificity,FDR)
  return(score_mtx)
}


for (i in 1:length(training_set)){
  pathway <- training_set[i]
  score <- top_score(cor_matrix_spearman,curation_pathway_filtered,pathway,25,'j')
  print(score)
}
ligand_receptor <- colnames(cor_matrix_pearson)

top_score_p <- function(ranking_df,curation_pathway,receptor_a,ranking_thr,file_name){
  score_mtx <- c()
  # Extract the top receptors to the certain pathway from the matrix
  top_list <- names(sort(ranking_df[,receptor_a],decreasing = TRUE)[1:ranking_thr])
  print(top_list)
  # Extract the curated receptors as reference
  pathway_curation <- filter(curation_pathway,receptor == receptor_a)$gsea_symbol
  print(pathway_curation)
  tp <- length(intersect(top_list,pathway_curation))
  fp <- length(top_list) - length(intersect(top_list,pathway_curation))
  tn <- nrow(ranking_df) - length(pathway_curation) - fp
  sensitivity <- tp/length(pathway_curation)
  specificity <- tn/(nrow(ranking_df) - length(pathway_curation))
  FDR <- fp/length(top_list)
  score_mtx <- c(sensitivity,specificity,FDR)
  return(score_mtx)
}

for (i in 1:ncol(matrix_spearman_kegg_wiki_all)){
  receptor_a <- colnames(matrix_spearman_kegg_wiki_all)[i]
  score <- top_score_p(matrix_spearman_kegg_wiki_all,curation_pathway_filtered,receptor_a,25,'j')
  print(score)
}


