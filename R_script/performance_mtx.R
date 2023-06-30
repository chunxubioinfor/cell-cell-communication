library(readxl)
library(dplyr)
library(tidyr)
library(progress)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Based on the Average Ranking result in which the plot is extremely left skewed,
# we decided to develop an algorithm applied onto cor_matrix to determine whether correlated
# Here, we introduce two types of algorithms:
# 1. A more straightforward one -- Iteration of the cut-off value
# 2. A more relative one -- Iteration of top value

## Randomly divide the curation dataset into two groups, testing data (2/3) and validation data (1/3)
curation_pathway <- read_excel('./curation_pathway.xlsx',sheet = "Curation_Pathway_for_R",col_names = TRUE)
curation_pathway <- curation_pathway %>% as_tibble() %>% separate_rows(receptors,sep = ',')  #reshape the dataset
curation_pathway_filtered <- filter(curation_pathway,receptors %in% ligand_receptor)
curation_pathway_filtered$index <- 1:nrow(curation_pathway_filtered)
write.csv(curation_pathway_filtered,'./curation_pathway_filtered.csv')
validation_pathway <- curation_pathway_filtered %>% group_by(source) %>% sample_frac(size = 1/3) #leave alone for final validation
testing_pathway<- anti_join(curation_pathway_filtered,validation_pathway,by = 'index')


## Apply two algorithms to calculate the performance matrix

# Define a function which transfer the correaltion coefficient matrix into correlated or not matrix
correlated_or_not_mtx <- function(cor_matrix, cut_off_value){
  cor_matrix[cor_matrix >= cut_off_value] <- 'TRUE'
  cor_matrix[cor_matrix < cut_off_value] <- 'FALSE'
  return(cor_matrix)
}

correlated_or_not_mtx <- function(cor_matrix, algorithm, value,view){
  if(algorithm == 'cut_off'){
    cut_off_value <- value
    cor_matrix[cor_matrix >= cut_off_value] <- 'TRUE'
    cor_matrix[cor_matrix < cut_off_value] <- 'FALSE'
    return(cor_matrix)
  }
  else if(algorithm == 'top'){
    top_value <- value
    cor_or_not_mtx <- matrix(NA, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix),
                             dimnames = list(rownames(cor_matrix),colnames(cor_matrix)))
    if(view == 'receptor'){
      for (i in 1:ncol(cor_matrix)) {
        col <- cor_matrix[, i]
        sorted_col <- sort(col, decreasing = TRUE)
        threshold <- sorted_col[top_value]
        cor_or_not_mtx[, i] <- col >= threshold
      }
    }
    else if(view == 'pathway'){
      for (i in 1:nrow(cor_matrix)) {
        row <- cor_matrix[i,]
        sorted_row <- sort(row, decreasing = TRUE)
        threshold <- sorted_row[top_value]
        cor_or_not_mtx[i,] <- row >= threshold
      }
    }
    return(cor_or_not_mtx)
  }
}
# Define a function which calculate and output performance matrix
perf_mtx <- function(cor_or_not_mtx,view,curation_pathway){
  if(view == 'receptor'){
    perf_matrix <- matrix(0,nrow = 3,ncol = ncol(cor_or_not_mtx),dimnames = list(c('sensitivity','specificity','FDR'),colnames(cor_or_not_mtx)))
    for(i in 1:ncol(cor_or_not_mtx)){
      correlated_pathway <- rownames(cor_or_not_mtx)[cor_or_not_mtx[,i] == "TRUE"]
      ref_pathway <- filter(curation_pathway,receptors == colnames(cor_or_not_mtx)[i])$gsea_symbol
      tp <- length(intersect(correlated_pathway,ref_pathway))
      fp <- length(correlated_pathway) -tp
      tn <- nrow(cor_or_not_mtx) - length(ref_pathway) - fp
      sensitivity <- tp/length(ref_pathway)
      specificity <- tn/(nrow(cor_or_not_mtx) - length(ref_pathway))
      FDR <- fp/length(correlated_pathway)
      perf_matrix[,i] <- c(sensitivity,specificity,FDR)
    }
    return(perf_matrix)
  }
  else if(view == 'pathway'){
    perf_matrix <- matrix(0,nrow = nrow(cor_or_not_mtx),ncol = 3,dimnames = list(rownames(cor_or_not_mtx),c('sensitivity','specificity','FDR')))
    for(i in 1:nrow(cor_or_not_mtx)){
      correlated_receptor <- colnames(cor_or_not_mtx)[cor_or_not_mtx[i,] == "TRUE"]
      ref_receptor <- filter(curation_pathway,gsea_symbol == rownames(cor_or_not_mtx)[i])$receptors
      tp <- length(intersect(correlated_receptor,ref_receptor))
      fp <- length(correlated_receptor) -tp
      tn <- nrow(cor_or_not_mtx) - length(ref_receptor) - fp
      sensitivity <- tp/length(ref_receptor)
      specificity <- tn/(nrow(cor_or_not_mtx) - length(ref_receptor))
      FDR <- fp/length(correlated_receptor)
      perf_matrix[i,] <- c(sensitivity,specificity,FDR)
    }
    perf_matrix <- t(perf_matrix)
    return(perf_matrix)
  }
}

# Test two functions
t <- correlated_or_not_mtx(cor_matrix_spearman,"cut_off",0.5,'receptor')
a <- perf_mtx(t,'pathway',curation_pathway_filtered)

# Subset the correlation matrix to test two algorithms 
receptors_subset <- unique(curation_pathway_filtered$receptors)
pathways_subset <- unique(curation_pathway_filtered$gsea_symbol)
cor_matrix_subset <- cor_matrix_spearman[pathways_subset,receptors_subset]

# 1. Algorithm I: Iteration of the cut-off value at the view of receptor;
# save all the performance matrix into a list: perf_mtx_ls_co
perf_mtx_ls_co <- list()
pb <- progress_bar$new(total = 15)
for (cut_off in seq(0.2, 0.9, 0.05)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',cut_off,'receptor')
  perf_matrix <- perf_mtx(cor_or_not_mtx,'receptor',curation_pathway_filtered)
  matrix_name <- paste('perf_matrix',cut_off,sep = '_')
  perf_mtx_ls_co[[matrix_name]] <- perf_matrix
  pb$tick()
}


# 2. Algorithm II: Iteration of the top value at the view of receptor
perf_mtx_ls_tp <- list()
pb <- progress_bar$new(total = 4)
for (top in c(1,5,10,25)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'top',top,'receptor')
  perf_matrix <- perf_mtx(cor_or_not_mtx,view = "receptor",curation_pathway_filtered)
  matrix_name <- paste('perf_matrix_top',top,sep = '_')
  perf_mtx_ls_tp[[matrix_name]] <- perf_matrix
  pb$tick()
}

t <- as.data.frame(t(perf_matrix_list[["perf_matrix_0.2"]]))
t <- tibble::rownames_to_column(t,'receptor')

# First meeting
# Performance Matrix Visulaization

## Generate a data frame which contains all the info
mtx_to_df <- function(matrix_list){
  output_df <- data.frame(matrix(nrow = 0,ncol = 5))
  colnames(output_df) <- c('algorithm','receptor','sensitivity','specificity','FDR')
  for(i in 1:length(matrix_list)){
    df <- as.data.frame(t(matrix_list[[i]]))
    df$receptor <- row.names(df)
    row.names(df) <- NULL
    df$algorithm <- names(matrix_list)[i]
    output_df <- rbind(output_df,df)
  }
  return(output_df)
}

perf_df_co <- mtx_to_df(perf_mtx_ls_co)
perf_df_tp <- mtx_to_df(perf_mtx_ls_tp)

## Visualize with a boxplot
perf_df_co <- melt(perf_df_co,id.vars = c('receptor','algorithm'),measure.vars = c('sensitivity','specificity','FDR'),
                  variable.name = 'perf_variable',value.name = 'perf_value')
perf_df_tp <- melt(perf_df_tp,id.vars = c('receptor','algorithm'),measure.vars = c('sensitivity','specificity','FDR'),
                   variable.name = 'perf_variable',value.name = 'perf_value')
perf_df_tp$algorithm <- factor(perf_df_tp$algorithm,levels = c('perf_matrix_top_1','perf_matrix_top_5',
                                                                'perf_matrix_top_10','perf_matrix_top_25'))
viz_boxplot <- function(melted_df,algorithm){
  p <- ggplot(data = melted_df, aes(x = algorithm,y = perf_value,fill=factor(perf_variable)))+
    geom_boxplot()+
    scale_fill_brewer(palette = "Pastel2")+
    scale_x_discrete(name = algorithm)+
    guides(fill=guide_legend(title = NULL))+
    geom_text(stat = "summary", fun = median, vjust = -1, aes(label = round(..y.., 2)),position = position_dodge(0.75))
  print(p)
  ggsave(paste0(algorithm,'.png'),p)
}

viz_boxplot(perf_df_co,"Cut-off Value")
viz_boxplot(perf_df_tp,"Ranking Threshold")

## Then,from the point of view of pathway
perf_mtx_ls_co_pw <- list()
pb <- progress_bar$new(total = 15)
for (cut_off in seq(0.2, 0.9, 0.05)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',cut_off,'pathway')
  perf_matrix <- perf_mtx(cor_or_not_mtx,'pathway',curation_pathway_filtered)
  matrix_name <- paste('perf_matrix',cut_off,sep = '_')
  perf_mtx_ls_co_pw[[matrix_name]] <- perf_matrix
  pb$tick()
}

perf_mtx_ls_tp_pw <- list()
pb <- progress_bar$new(total = 4)
for (top in c(1,5,10,25)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'top',top,'pathway')
  perf_matrix <- perf_mtx(cor_or_not_mtx,view = "pathway",curation_pathway_filtered)
  matrix_name <- paste('perf_matrix_top',top,sep = '_')
  perf_mtx_ls_tp_pw[[matrix_name]] <- perf_matrix
  pb$tick()
}

perf_df_co_pw <- mtx_to_df(perf_mtx_ls_co_pw)
perf_df_tp_pw <- mtx_to_df(perf_mtx_ls_tp_pw)

perf_df_co_pw <- melt(perf_df_co_pw,id.vars = c('receptor','algorithm'),measure.vars = c('sensitivity','specificity','FDR'),
                   variable.name = 'perf_variable',value.name = 'perf_value')
perf_df_tp_pw <- melt(perf_df_tp_pw,id.vars = c('receptor','algorithm'),measure.vars = c('sensitivity','specificity','FDR'),
                   variable.name = 'perf_variable',value.name = 'perf_value')
perf_df_tp_pw$algorithm <- factor(perf_df_tp_pw$algorithm,levels = c('perf_matrix_top_1','perf_matrix_top_5',
                                                               'perf_matrix_top_10','perf_matrix_top_25'))
viz_boxplot(perf_df_co_pw,"Cut-off Value")
viz_boxplot(perf_df_tp_pw,"Ranking Threshold")


## Finally, save the important output incl. all the performance matrix
saveRDS(list(perf_mtx_ls_co,perf_mtx_ls_tp),'./perf_mtx_R.rds')
saveRDS(list(perf_mtx_ls_co_pw,perf_mtx_ls_tp_pw),'./perf_mtx_P.rds')


## To be continued
