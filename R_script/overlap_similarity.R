# Prove how our predicted correlated pathways are internally similar

# Importing packages
library(BiocManager)
BiocManager::install("OmnipathR")
BiocManager::install("cogena")
library(OmnipathR)
library(cogena)
library(proxy)
library(ggplot2)
library(progress)
library(colorspace)
library(clusterProfiler)

## Getting a list where each entry is the genes in a pathway 
## and the name of that entry is the name of the pathway
data_folder <- '/Users/hanchunxu/Desktop/cell_cell_comm/data/'
result_folder <- '/Users/hanchunxu/Desktop/cell_cell_comm/result/'
cp_gs_ls <- cogena::gmt2list(paste0(data_folder,'c2.cp.v2022.1.Hs.symbols.gmt'))
h_gs_ls <- cogena::gmt2list(paste0(data_folder,'h.all.v2022.1.Hs.symbols.gmt'))
go_gs_ls <- cogena::gmt2list(paste0(data_folder,'go.slim.v2023.1.Hs.symbols.gmt'))
gsOi <- c(cp_gs_ls,h_gs_ls,go_gs_ls)

## Generating the binary matrix
allGenes <- unique(unlist(gsOi))
gsMat <- do.call(
  cbind,
  lapply(gsOi, function(x) {
    matrix(as.integer(allGenes %in% x))
  })
)
colnames(gsMat) <- names(gsOi)
rownames(gsMat) <- allGenes
dim(gsMat) # 19515 * 3152
saveRDS(gsMat,paste0(result_folder,'gsMat.rds'))

## Calculating the overlap coefficient matrix
jcMat <- as.matrix(dist(t(gsMat), method = 'binary'))
ocMat <- as.matrix(proxy::simil(x = t(gsMat), method="Simpson"))
ocMat[is.na(ocMat)] <- 0

## Define a function which calculate the z score distribution
calc_z_score <- function(ocMat,cor_or_not_mtx){
  z_score_df <- data.frame(receptor = colnames(cor_or_not_mtx),
                           z_score = rep(NA,ncol(cor_or_not_mtx)),
                           algorithm = rep(NA,ncol(cor_or_not_mtx)))
  # Calculate the population mean and sd
  mu <- mean(ocMat,na.rm = TRUE)
  sigma <- sd(ocMat,na.rm = TRUE)
  # For each receptor, calculate the z score
  for(i in 1:ncol(cor_or_not_mtx)){
    correlated_pathways <- rownames(cor_or_not_mtx)[cor_or_not_mtx[,i] == "TRUE"]
    ocMat_subset <- ocMat[correlated_pathways,correlated_pathways]
    score <- mean(ocMat_subset,na.rm = TRUE)
    z_score <- (score-mu)/sigma
    z_score_df$z_score[i] <- z_score
  }
  return(z_score_df)
}

## Check the Z-score of a given receptor
z_score_R <- function(ocMat,cor_or_not_mtx,receptor){
  correlated_pathways <- rownames(cor_or_not_mtx)[cor_or_not_mtx[,receptor] == "TRUE"]
  mu <- mean(ocMat,na.rm = TRUE)
  sigma <- sd(ocMat,na.rm = TRUE)
  ocMat_subset <- ocMat[correlated_pathways,correlated_pathways]
  score <- mean(ocMat_subset,na.rm = TRUE)
  z_score <- (score-mu)/sigma
  print(z_score)
}

eg_R <- c('FZD2','FZD6','FZD9','TLR4','ITGA7','ROR2','CD8A')
for (receptor in eg_R){
  z_score_R(ocMat,cor_or_not_mtx,receptor)
}


## Iterate all the 19 algorithms
## Algorithm I: cut-off value
cor_matrix_subset <- readRDS(paste0(result_folder,'cor_mtx_curation_sub.rds'))
pb <- progress_bar$new(total = 15)
z_score_ls <- list()
for (cut_off in seq(0.2, 0.9, 0.05)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',cut_off,'receptor')
  z_score_df <- calc_z_score(ocMat,cor_or_not_mtx)
  df_name <- as.character(cut_off)
  z_score_df$algorithm <- rep(df_name,nrow(z_score_df))
  z_score_ls[[df_name]] <- z_score_df
  pb$tick()
}
## Visualize
z_cut_off_df <- do.call(rbind,z_score_ls)
z_cut_off_dp <- ggplot(z_cut_off_df,aes(x=z_score,fill=algorithm,alpha=1/10)) +
  geom_density(position="identity")+
  scale_fill_discrete_qualitative(palette = "Set3")+
  guides(alpha = FALSE,fill=guide_legend(title = "cut-off"))+
  theme_minimal()
z_cut_off_dp
ggsave(filename = paste0(result_folder,"z_score_d_cut_off.png"),plot = z_cut_off_dp)
z_cut_off_bp <- ggplot(z_cut_off_df,aes(x=algorithm,y=z_score,fill=algorithm))+
  geom_boxplot()+
  theme_minimal()+
  labs(x="cut-off",y='Z score')+
  theme(axis.title = element_text(size = 15),legend.position="none",
        axis.text = element_text(size = 12))+
  scale_fill_manual(values = colorRampPalette(c('#e6eaf0','#08316b'))(15))+
  geom_text(stat = "summary", fun = median, vjust = -1, 
            aes(label = round(..y.., 2)),position = position_dodge(0.75))
z_cut_off_bp
ggsave(filename = paste0(result_folder,"z_score_boxp_cut_off.png"),plot = z_cut_off_bp)


## Algorithm II: Ranking Threshold
pb <- progress_bar$new(total = 4)
z_score_ls <- list()
for (top in c(1,5,10,25)){
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'top',top,'receptor')
  z_score_df <- calc_z_score(ocMat,cor_or_not_mtx)
  df_name <- as.character(top)
  z_score_df$algorithm <- rep(df_name,nrow(z_score_df))
  z_score_ls[[df_name]] <- z_score_df
  pb$tick()
}
## ## visuzlise
z_ranking_df <- do.call(rbind,z_score_ls)
z_ranking_df$algorithm <- factor(z_ranking_df$algorithm,levels = c('1','5','10','25'))
z_ranking_dp <- ggplot(z_ranking_df,aes(x=z_score,
                                      fill=factor(algorithm,levels = c('1','5','10','25')),
                                      alpha=1/10)) +
  geom_density(position="identity")+
  scale_fill_discrete_qualitative(palette = "Set3") +
  guides(alpha = FALSE,fill=guide_legend(title = "ranking-threshold"))+
  theme_minimal()
z_ranking_dp
ggsave(filename = paste0(result_folder,"z_score_d_ranking_thr.png"),plot = z_ranking_dp)
z_ranking_bp <- ggplot(z_ranking_df,aes(x=algorithm,y=z_score,
                                        fill=algorithm))+
  geom_boxplot()+
  theme_minimal()+
  labs(x="ranking threshold",y='Z score')+
  theme(axis.title = element_text(size = 15),legend.position="none",
        axis.text = element_text(size = 12))+
  scale_fill_manual(values = colorRampPalette(c('#e6eaf0','#08316b'))(4))+
  geom_text(stat = "summary", fun = median, vjust = -1, 
            aes(label = round(..y.., 2)),position = position_dodge(0.75))
z_ranking_bp
ggsave(filename = paste0(result_folder,"z_score_boxp_ranking_thr.png"),plot = z_ranking_bp)

## further determine the final cur-off value at range of (0.6-0.7)
z_score_ls <- list()
nc <- c()
for (cut_off in seq(0.6, 0.7, 0.01)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',cut_off,'receptor')
  z_score_df <- calc_z_score(ocMat,cor_or_not_mtx)
  df_name <- as.character(cut_off)
  z_score_df$algorithm <- rep(df_name,nrow(z_score_df))
  z_score_ls[[df_name]] <- z_score_df
  nc <- c(nc,sum(cor_or_not_mtx == 'TRUE',na.rm = TRUE))
}
z_cut_off_df <- do.call(rbind,z_score_ls)
z_cut_off_bp <- ggplot(z_cut_off_df,aes(x=algorithm,y=z_score,fill=algorithm))+
  geom_boxplot()+
  theme_minimal()+
  labs(x="cut-off",y='Z score')+
  theme(axis.title = element_text(size = 15),legend.position="none",
        axis.text = element_text(size = 12))+
  scale_fill_manual(values = colorRampPalette(c('#e6eaf0','#08316b'))(15))+
  geom_text(stat = "summary", fun = median, vjust = -1, 
            aes(label = round(..y.., 2)),position = position_dodge(0.75))
z_cut_off_bp

perf_df <- data.frame(algorithm = as.character(seq(0.6, 0.7, 0.01)),
                      z_score_median = unlist(lapply(z_score_ls,function(df) median(df$z_score,na.rm = TRUE))),
                      sensitivity = unlist(lapply(perf_mtx_ls_co,function(mtx) median(mtx['sensitivity',],na.rm = TRUE))),
                      specificity = unlist(lapply(perf_mtx_ls_co,function(mtx) median(mtx['specificity',],na.rm = TRUE))),
                      num_of_connections = nc)
cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',0.6,'receptor')
sum(cor_or_not_mtx == 'TRUE',na.rm = TRUE)

cor_or_not_mtx <- as.numeric(cor_or_not_mtx)
curation_pathway_filtered <- read.csv(paste0(data_folder,'curation_pathway_filtered.csv'))
perf_mtx_ls_co <- list()
for (cut_off in seq(0.6, 0.7, 0.01)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',cut_off,'receptor')
  perf_matrix <- perf_mtx(cor_or_not_mtx,'receptor',curation_pathway_filtered)
  matrix_name <- paste('perf_matrix',cut_off,sep = '_')
  perf_mtx_ls_co[[matrix_name]] <- perf_matrix
}
perf_df_co <- mtx_to_df(perf_mtx_ls_co)

