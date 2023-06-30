## Provide some examples to demonstrate the biological significance of this method
## Pick critical receptors whose corresponding pathways exhibit significant biological information
library(dplyr)

# load the performance matrix
perf_mtx_R <- readRDS('./perf_mtx_R.rds')
perf_mtx_co_0.5 <- as.data.frame(t(perf_mtx_R[[1]][["perf_matrix_0.5"]]))
cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',0.5,'receptor')

# select some well-performed receptors
well_R <- filter(perf_mtx_co_0.5,sensitivity >= 0.5,specificity >= 0.5,FDR <= 0.7)


# Pay attention to some common signaling pathways


# For one certain critical receptor, define a function to compare the predicted and curated pathways
pred_n_cur <- function(cor_or_not_mtx,receptor){
  pred_P <- rownames(cor_or_not_mtx)[cor_or_not_mtx[,receptor]=='TRUE']
  cur_P <- filter(curation_pathway_filtered,receptors == receptor)$gsea_symbol
  pred_n_cur_ls <- list(predicted = pred_P,curation = cur_P)
  common_P <- intersect(pred_P,cur_P)
  print(paste0('The predicted ',length(pred_P),' correlated pathways:'))
  print(pred_P)
  print(paste0('The curated ',length(cur_P),' correlated pathways:'))
  print(cur_P)
  print('Common Pathways:')
  print(common_P)
}

# R1: TLR4
pred_n_cur(cor_or_not_mtx,'TLR4')
pred_n_cur(cor_or_not_mtx,'FZD2')
pred_n_cur(cor_or_not_mtx,'FZD9')
pred_n_cur(cor_or_not_mtx,'ITGA7')
pred_n_cur(cor_or_not_mtx,'ROR2')
pred_n_cur(cor_or_not_mtx,'CD8A')
pred_n_cur(cor_or_not_mtx,'PDGFRA')
