#===========================
# function: rGene19()
# Purpose: grab gene IDs and coordinates (Ensembl and/or UCSC) overlapping a query region of the hg19 human genome build.

#=============
# arguments
# input.range : gene coordinate ranges as GRanges object. Seqnames should be of format "chr#"
# searchwindow (T/F) : whether to use a search window up/downstream of inputted range
# windowsize : how far to look up and downstream, if searchwindow==T
# get.tx.coor (T/F) : whether to grab UCSC gene coordinates, from txdb package
# get.en.coor (T/F) : whether to grab ensembl db coordinates, from ensdb package
# returngr (T/F) : whether to return gene info as granges objects
# show.module.versions (T/F) : whether to show and store module versions used in query

rGene19 <- function(input.range,
                    searchwindow=T,
                    windowsize=2e+5,
                    get.tx.coor=T,
                    get.en.coor=T,
                    show.module.versions=T,
                    returngr=T){

packages.list <- c()

if(get.tx.coor==T){
  #library(org.Hs.eg.db)
  require("org.Hs.eg.db")
  packages.list <- c(packages.list,"org.Hs.eg.db")
  #columns(org.Hs.eg.db)
  
  require("TxDb.Hsapiens.UCSC.hg19.knownGene")
  packages.list <- c(packages.list,"TxDb.Hsapiens.UCSC.hg19.knownGene")
}

if(get.en.coor==T){
  require("EnsDb.Hsapiens.v75")
  edb <- EnsDb.Hsapiens.v75; 
  packages.list <- c(packages.list,"EnsDb.Hsapiens.v75")
}

listreturn <- list()

if(searchwindow==T){
  rangeix <- makeGRangesFromDataFrame(data.frame(start=start(input.range)-windowsize,
                                                 end=end(input.range)+windowsize,
                                                 chr=as.character(seqnames(input.range))))
} else{
  rangeix <- makeGRangesFromDataFrame(data.frame(start=start(input.range),
                                                 end=end(input.range),
                                                 chr=as.character(seqnames(input.range))))
}; listreturn <- c(listreturn,rangeix); names(listreturn)[length(listreturn)] <- "query.gr"

if(get.tx.coor==T){
  message("querying ucsc db gene info...")
  gene.txdb <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) # returns large GRanges object
  gid.ez.df <- subsetByOverlaps(gene.txdb,rangeix)
  mcols(gid.ez.df)$genesymbol <- "NA"
  for(z in 1:length(gid.ez.df)){
    try(mcols(gid.ez.df)$genesymbol[z] <- suppressMessages(unique(select(org.Hs.eg.db, gid.ez.df$gene_id[z], c("SYMBOL"), "ENTREZID")$SYMBOL)),
        silent=T)
  }
  
  
  if(!returngr==T){
    gid.ez.df <- data.frame(ez_symbol=gid.ez.df$gene_id,
                            genesymbol=gid.ez.df$genesymbol,
                            start=start(gid.ez.df),
                            end=end(gid.ez.df),
                            chr=seqnames(gid.ez.df),
                            strand=strand(gid.ez.df),
                            stringsAsFactors = FALSE)
    
  }
  
  
  listreturn <- c(listreturn,gid.ez.df); names(listreturn)[length(listreturn)] <- "tx.gr"
}

if(get.en.coor==T){
  message("querying ensembl db gene info...")
  genes.ens <- genes(edb); 
  genes.ens <- genes.ens[genes.ens$gene_biotype=="protein_coding"]; #unique(seqnames(genes.ens))
  
  genes.ens2 <- makeGRangesFromDataFrame(data.frame(start=start(genes.ens),
                                                    end=end(genes.ens),
                                                    seqnames=paste0("chr",seqnames(genes.ens)),
                                                    strand=strand(genes.ens),
                                                    genesymbol=genes.ens$gene_name,
                                                    ens.symbol=genes.ens$gene_id),
                                         keep.extra.columns = T)
  
  gid.ens <- subsetByOverlaps(genes.ens2,rangeix)
  gid.ens.df <- data.frame(ens.symbol=gid.ens$ens.symbol,
                           start=start(gid.ens),
                           end=end(gid.ens),
                           chr=seqnames(gid.ens),
                           genesymbol=gid.ens$genesymbol,
                           strand=strand(gid.ens),
                           stringsAsFactors = FALSE)
  if(returngr==T){
    gid.ens.df <- makeGRangesFromDataFrame(gid.ens.df,keep.extra.columns = T)
  }
  
  listreturn <- c(listreturn,gid.ens.df); names(listreturn)[length(listreturn)] <- "en.gr"
}

if(show.module.versions==T){
  message("Package versions used:")
  
  packages.list.return <- list()
  for(i in 1:length(packages.list)){
    packages.list.return <- c(packages.list.return,as.character(packageVersion(packages.list[i])))
    names(packages.list.return)[length(packages.list.return)] <- packages.list[i]
    
    message(names(packages.list.return)[i],": ",packages.list.return[i])
  }
  listreturn <- c(listreturn,list(packages.list.return)); names(listreturn)[length(listreturn)] <- "Rpackage.versions"
}

return(listreturn)

}
