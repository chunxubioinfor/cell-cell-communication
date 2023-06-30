## This is a R script run on server ##
## The transposation of the matrix takes too much time##
## So this is an improved version which postpone the transpose to last step##


# install and import packages
library(BiocManager)
library(IsoformSwitchAnalyzeR)
library(GSVA)
library(msigdbr)
library(rhdf5)
library(OmnipathR)
library(tidyverse)
library(snow)
library(WGCNA)
library(clusterProfiler)
options(stringsAsFactors = F) 


# Filter the low preformance archs4 data
start_time <- Sys.time()
setwd('/home/projects/kvs_ccc/')
dir.create('./output')
sink('./output/output.txt')
archs4file= "/home/databases/archs4/v11/human_transcript_v11_tpm.h5"
arch4meta <- 
  data.frame(
    study             = h5read(archs4file, "/meta/samples/series_id"),
    id                = h5read(archs4file, "/meta/samples/geo_accession"),
    organism          = h5read(archs4file, "/meta/samples/organism_ch1"),
    source            = h5read(archs4file, "/meta/samples/source_name_ch1"),
    title             = h5read(archs4file, "/meta/samples/title"),
    characteristics   = h5read(archs4file, "/meta/samples/characteristics_ch1"),
    singleCellProb    = h5read(archs4file, "/meta/samples/singlecellprobability"),
    type              = h5read(archs4file, "/meta/samples/type"),
    library_selection = h5read(archs4file, "/meta/samples/library_selection"),
    readsAlligned     = h5read(archs4file, "/meta/samples/readsaligned"),
    readsTotal        = h5read(archs4file, "/meta/samples/readstotal"),
    submissionDate    = h5read(archs4file, "/meta/samples/submission_date")
  )
arch4metaFilt <- 
  arch4meta %>% 
  as_tibble() %>% 
  group_by(study) %>% 
  mutate(
    allignFrac = readsAlligned / readsTotal,
  ) %>% 
  filter(
    organism == 'Homo sapiens',
    type == 'SRA',
    singleCellProb <= 1/3,
    library_selection == 'cDNA',
    allignFrac > 0.5,
    readsAlligned >= 5e6
  )

print(paste('There are',nrow(arch4metaFilt),'samples available!'))

# Extract the available samples from ARCHS4_transcript_file
sample_id <- arch4metaFilt$id
samples <- h5read(archs4file, "/meta/samples/geo_accession")
print(length(samples))
transcripts <- h5read(archs4file,'meta/transcripts/transcripts')
print(length(transcripts))
sample_locations = which(samples %in% sample_id)
chunk_number <- 17
sample_loc <- split(sample_locations,cut(seq_along(sample_locations),chunk_number,labels = FALSE))
transcript_expression_list <- list()
for (i in 1:17){
  transcript_data <- h5read(archs4file, "data/expression", index=list(sample_loc[[i]], 1:length(transcripts)))
  transposedData <- transposeBigData(transcript_data, blocksize = 20000)
  transcript_expression_list[[i]] <- transposedData
  print(paste('The',i,'is done!'))
}
transcript_expression <- do.call(cbind, transcript_expression_list)
colnames(transcript_expression) <- transcripts
rownames(transcript_expression) <- samples[sample_locations]
print('The generation of transcripts expression matrix has been done!')
print(paste('The dimension of the transcripts expression matrix is',dim(transcript_expression),sep=' '))


# Quantification from transcripts to genes level
## Import GTF file
aSwitchList <- importGTF(pathToGTF='/home/databases/archs4/v11/Homo_sapiens.GRCh38.87.chr_patch_hapl_scaff.gtf.gz')

## Create dataframe with associations
isoGene <- unique(aSwitchList$isoformFeatures[c('isoform_id','gene_name')])
colnames(isoGene) <- c('isoform_id','gene_id')
gene_expression <- IsoformSwitchAnalyzeR::isoformToGeneExp(transcript_expression,
                                                           isoformGeneAnnotation = isoGene,
                                                           quiet = FALSE)
end_time <- Sys.time()
print(end_time - start_time)
print('The generation of genes expression matrix has been done!')
print(paste('The dimension of the gene expression matrix is',dim(gene_expression),sep=' '))
write_rds(gene_expression,'./output/gene_expression_matrix.rds',compress = 'gz')
print('The gene expression matrix has been saved in /home/projects/kvs_ccc/output/gene_expression_matrix.rds.gz')


# Quantification from genes to gene-sets level
## Import gene-sets from MSigDB (a combination of three gene sets -- CP + GO slim + Hallmarks)
cp_gene_sets <- clusterProfiler::read.gmt('./data/c2.cp.v2022.1.Hs.symbols.gmt')
h_gene_sets <- clusterProfiler::read.gmt('./data/h.all.v2022.1.Hs.symbols.gmt')
go_gene_sets_df <- OmnipathR::go_annot_slim(organism = 'human',
                                            slim = 'agr',
                                            aspects = c('C','F','P'),
                                            cache = TRUE)
saveRDS(go_gene_sets_df,'./go_gene_sets_df.rds')
go_gene_sets <- dplyr::select(go_gene_sets_df,go_id,db_object_symbol)
colnames(go_gene_sets) <- c('term','gene')
gene_sets <- rbind(cp_gene_sets,h_gene_sets,go_gene_sets)
gset_list <- split(gene_sets$gene, gene_sets$term)
print('The preparation of gene sets matrix generation is done!')

## Bring the matrix to gene-sets level
start_time_gset <- Sys.time()

gene_expression_matrix <- as.matrix(gene_expression)
gene_set_expression <- gsva(gene_expression_matrix,
                            gset.idx.list = gset_list,
                            method = 'ssgsea',
                            verbose = T,
                            ssgsea.norm = F,
                            parallel.sz = 16)
print('The generation of gene sets expression matrix has been done!')
print(paste('The dimension of the gene expression matrix is',dim(gene_set_expression),sep=' '))
write_rds(gene_set_expression,'./output/gene_set_expression_matrix.rds.gz',compress = 'gz')
print('The gene expression matrix has been saved in /home/projects/kvs_ccc/output/gene_set_expression_matrix.rds.gz')

end_time_gset <- Sys.time()
print(paste('Done! The overall time cost is',end_time_gset - start_time_gset))
sink()
