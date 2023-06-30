# write a gmt file containing GO slim gene set

go_slim_df <- OmnipathR::go_annot_slim(organism = 'human',
                                            slim = 'agr',
                                            aspects = c('C','F','P'),
                                            cache = TRUE)
go_slim_df <- dplyr::select(go_slim_df,go_id,db_object_symbol)
colnames(go_slim_df) <- c('GO_id','gene')
go_slim_ls <- split(go_slim_df$gene, go_slim_df$GO_id)

file <- "go.slim.v2023.1.Hs.symbols.gmt"
gs_ls <- go_slim_ls
write.gmt <- function(gs_ls,file){
  sink(file)
  lapply(names(gs_ls), function(i){
    cat( paste(c(i,'NA',gs_ls[[i]]),collapse='\t') )
    cat('\n')
  })
  sink()
}
write.gmt(gs_ls,file)

