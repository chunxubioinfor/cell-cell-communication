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

## Getting a list where each entry is the genes in a pathway 
## and the name of that entry is the name of the pathway
cp_gs_ls <- cogena::gmt2list('c2.cp.v2022.1.Hs.symbols.gmt')
h_gs_ls <- cogena::gmt2list('h.all.v2022.1.Hs.symbols.gmt')
go_gs_ls <- cogena::gmt2list('go.slim.v2023.1.Hs.symbols.gmt')
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
saveRDS(gsMat,'./gsMat.rds')

## Calculating the overlap coefficient matrix
jcMat <- as.matrix(dist(t(gsMat), method = 'binary'))
ocMat <- as.matrix(proxy::simil(x = t(gsMat), method="Simpson"))

## Define a function which calculate the z score distribution
calc_z_score <- function(ocMat,cor_or_not_mtx){
  z_score_df <- data.frame(receptor = colnames(cor_or_not_mtx),z_score = rep(NA,ncol(cor_or_not_mtx)),
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

## Test
tmp_df <- calc_z_score(ocMat,cor_or_not_mtx)
tmp_p <- ggplot(tmp_df,aes(x=z_score)) +
  geom_density()


## Iterate all the algorithms
pb <- progress_bar$new(total = 19)
z_score_ls <- list()
for (cut_off in seq(0.2, 0.9, 0.05)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'cut_off',cut_off,'receptor')
  z_score_df <- calc_z_score(ocMat,cor_or_not_mtx)
  df_name <- paste('z_score_df',cut_off,sep = '_')
  z_score_df$algorithm <- rep(df_name,nrow(z_score_df))
  z_score_ls[[df_name]] <- z_score_df
  pb$tick()
}
for (top in c(1,5,10,25)){
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_subset,'top',top,'receptor')
  z_score_df <- calc_z_score(ocMat,cor_or_not_mtx)
  df_name <- paste('z_score_df',top,sep = '_')
  z_score_df$algorithm <- rep(df_name,nrow(z_score_df))
  z_score_ls[[df_name]] <- z_score_df
  pb$tick()
}

## visualize the z score
z_score_df <- do.call(rbind,z_score_ls)
z_score_p <- ggplot(z_score_df,aes(x=z_score,fill=algorithm,alpha=1/10)) +
  geom_density(position="identity")+
  scale_fill_discrete_qualitative(palette = "cold")
z_score_p <- ggplot(z_score_df,aes(x=z_score,fill=algorithm,alpha=1/10)) +
  geom_density(position="identity")+
  scale_fill_brewer(palette = "Set3")+
  theme_minimal()


## To be continued 

