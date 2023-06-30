# install and import packages
library(BiocManager)
library(IsoformSwitchAnalyzeR)
library(GSVA)
library(GSVAdata)
library(msigdbr)
library(rhdf5)
library(biomaRt)
library(GSA)
library(clusterProfiler)
options(stringsAsFactors = F) 

setwd('/home/projects/kvs_ccc/')
ARCHS4_transcript_file = "/home/databases/archs4/v11/human_transcript_v11_tpm.h5"
human_transcript_v11 <- H5Fopen(ARCHS4_transcript_file)
h5dump(human_transcript_v11,load=FALSE)    # Dump the content of an HDF5 file
samples <- human_transcript_v11$meta$samples$geo_accession
transcripts <- human_transcript_v11$meta$transcripts$transcripts
transcript_expression <- t(h5read(human_transcript_v11, "data/expression"))
H5close()
rownames(transcript_expression) <- transcripts
colnames(transcript_expression) <- samples
write.csv(transcript_expression,'./transcript_expression.csv')
print('The generation of transcripts expression matrix has been done!')
saveRDS(transcript_expression,'./transcript_expression_matrix.rds')

# Quantification from transcripts to genes level
### Import GTF file
aSwitchList <- importGTF(pathToGTF='/home/databases/archs4/v11/Homo_sapiens.GRCh38.87.chr_patch_hapl_scaff.gtf')

### Create dataframe with associations
isoGene <- unique(aSwitchList$isoformFeatures[c('isoform_id','gene_name')])
colnames(isoGene) <- c('isoform_id','gene_id')

### Overwrite with gene names
#isoGene$gene_name <- aSwitchList$isoformFeatures$gene_name[match(
#  isoGene$isoform_id, aSwitchList$isoformFeatures$isoform_id
#)]

gene_expression <- IsoformSwitchAnalyzeR::isoformToGeneExp(transcript_expression,
                                                           isoformGeneAnnotation = isoGene,
                                                           quiet = FALSE)
write.csv(gene_expression,'./gene_expression.csv')
print('The generation of genes expression matrix has been done!')
saveRDS(gene_expression,'./gene_expression_matrix.rds')

# Quantification from genes to gene-sets level
## Import gene-sets from MSigDB (a combination of three gene sets -- CP + GO(BP+MF+CC) + Hallmarks)
cp_gene_sets <- clusterProfiler::read.gmt('./data/c2.cp.v2022.1.Hs.symbols.gmt')
h_gene_sets <- clusterProfiler::read.gmt('./data/h.all.v2022.1.Hs.symbols.gmt')
go_gene_sets <- clusterProfiler::read.gmt('./data/c5.go.v2022.1.Hs.symbols.gmt')
gene_sets <- rbind(cp_gene_sets,h_gene_sets,go_gene_sets)
gset_list <- split(gene_sets$gene, gene_sets$term)

## Bring the matrix to genesets level
gene_expression_matrix <- as.matrix(gene_expression)
gene_set_expression <- gsva(gene_expression_matrix,
                            gset.idx.list = gset_list,
                            method = 'ssgsea',
                            kcdf = 'Poisson',
                            verbose = T)
write.csv(gene_set_expression,'./gset_matrix/gset_expression.csv')
print('The generation of gene sets expression matrix has been done!')
print('The gene sets expression matrix was strored in /home/projects/kvs_ccc//gset_matrix/gset_expression.csv')



