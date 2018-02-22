# geneMapR
Repo containing info and workhorse scripts to aid gene ID mapping and coordinate retrieval.

Downloadable using the following in R:   
`require(devtools);   install_github("metamaden/methyPre")`

# Mapping a Gene
Sequence -> Codon -> Transcript -> Isoform -> Isoform Cluster -> Gene Coordinate Range

# Gene IDs
Gene identifiers correspond to any number of different molecular data types, including different types of transcripts.

# Genome Builds
The genome build (and affiliate patch number) corresponds to a designation like "hg##" ([UCSC genome build](https://genome.ucsc.edu/FAQ/FAQreleases.html)) or "GRCh##" (Genome Assembly from the NCBI [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)). Genomes can be stored as a scaffold sequence, [annotated](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/) database, or some combination of the two.

# Programmatic Mapping: R
Using R/Bioconductor, id mappings can readily be generated from org.Hs.eg.db package, using the select() interface inherited from AnnotationDbi package (see [org.Hs.rg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) and [Annotation tutorial](https://www.bioconductor.org/help/workflows/annotation/annotation/) Bioconductor pages for more details).

Here is an example building from the [Bioc help archives](https://support.bioconductor.org/p/59549/):
```
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")

library(org.Hs.eg.db)
columns(org.Hs.eg.db)
#[1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"   "EVIDENCEALL" 
#[10] "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"     
#[19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT" 

symbols <- c("ERBB2","TP53")
head(select(org.Hs.eg.db, symbols, c("ENTREZID","GENENAME","REFSEQ"), "ALIAS"))
#  ALIAS ENTREZID                          GENENAME       REFSEQ
#1 ERBB2     2064 erb-b2 receptor tyrosine kinase 2 NM_001005862
#2 ERBB2     2064 erb-b2 receptor tyrosine kinase 2 NM_001289936
#3 ERBB2     2064 erb-b2 receptor tyrosine kinase 2 NM_001289937
#4 ERBB2     2064 erb-b2 receptor tyrosine kinase 2 NM_001289938
#5 ERBB2     2064 erb-b2 receptor tyrosine kinase 2    NM_004448
#6 ERBB2     2064 erb-b2 receptor tyrosine kinase 2 NP_001005862
```
To extract the gene coordinate ranges, use either the UCSC (TxDB/BS.xx) or ENSEMBL (EnsDb.xx) packages. For instance:
```
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
edb
# EnsDb for Ensembl:
# |Backend: SQLite
# |Db type: EnsDb
# |Type of Gene ID: Ensembl Gene ID
# |Supporting package: ensembldb
# |Db created by: ensembldb package from Bioconductor
# |script_version: 0.2.3
# |Creation time: Tue Nov 15 23:35:19 2016
# |ensembl_version: 75
# |ensembl_host: localhost
# |Organism: homo_sapiens
# |genome_build: GRCh37
# |DBSCHEMAVERSION: 1.0
# | No. of genes: 64102.
# | No. of transcripts: 215647.
# |Protein data available.

columns(edb)
# [1] "ENTREZID"            "EXONID"              "EXONIDX"             "EXONSEQEND"          "EXONSEQSTART"        "GENEBIOTYPE"        
# [7] "GENEID"              "GENENAME"            "GENESEQEND"          "GENESEQSTART"        "INTERPROACCESSION"   "ISCIRCULAR"         
# [13] "PROTDOMEND"          "PROTDOMSTART"        "PROTEINDOMAINID"     "PROTEINDOMAINSOURCE" "PROTEINID"           "PROTEINSEQUENCE"    
# [19] "SEQCOORDSYSTEM"      "SEQLENGTH"           "SEQNAME"             "SEQSTRAND"           "SYMBOL"              "TXBIOTYPE"          
# [25] "TXCDSSEQEND"         "TXCDSSEQSTART"       "TXID"                "TXNAME"              "TXSEQEND"            "TXSEQSTART"         
# [31] "UNIPROTDB"           "UNIPROTID"           "UNIPROTMAPPINGTYPE" 

genex <- "TP53"
select(edb, keys=genex, columns=c("TXID", "TXSEQSTART", "TXBIOTYPE","SEQNAME","GENESEQSTART","GENESEQEND"),keytype="SYMBOL")
#SYMBOL            TXID TXSEQSTART            TXBIOTYPE SEQNAME GENESEQSTART GENESEQEND
#1    TP53 ENST00000413465    7565097       protein_coding      17      7565097    7590856
#2    TP53 ENST00000359597    7569404       protein_coding      17      7565097    7590856
#3    TP53 ENST00000504290    7571720      retained_intron      17      7565097    7590856
#4    TP53 ENST00000510385    7571720      retained_intron      17      7565097    7590856
#5    TP53 ENST00000504937    7571720      retained_intron      17      7565097    7590856
#6    TP53 ENST00000269305    7571720       protein_coding      17      7565097    7590856
#7    TP53 ENST00000455263    7571722       protein_coding      17      7565097    7590856
#8    TP53 ENST00000420246    7571722       protein_coding      17      7565097    7590856
#9    TP53 ENST00000445888    7571739       protein_coding      17      7565097    7590856
#10   TP53 ENST00000576024    7572887       protein_coding      17      7565097    7590856
#11   TP53 ENST00000509690    7576853       protein_coding      17      7565097    7590856
#12   TP53 ENST00000514944    7577535       protein_coding      17      7565097    7590856
#13   TP53 ENST00000574684    7577572 processed_transcript      17      7565097    7590856
#14   TP53 ENST00000505014    7577844      retained_intron      17      7565097    7590856
#15   TP53 ENST00000508793    7578434       protein_coding      17      7565097    7590856
#16   TP53 ENST00000604348    7578480       protein_coding      17      7565097    7590856
#17   TP53 ENST00000503591    7578547       protein_coding      17      7565097    7590856

```
NOTE: gene coordinates are not uncontroversial. Even controlling for genome build, different databases can yield varying exact coordinate ranges. These should, however, cluster predictably and according to annotated transcripts and possible variation in transcript annotations, database updates, curation practices, etc. Thus it often pays to check multiple resources. See helpful links below for some options.

To explore this topic a bit more, we can explore the discrepancy between the transcript ranges and gene sequence ranges with a few examples.

```
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

genes(txdb) # returns large GRanges object
#GRanges object with 23056 ranges and 1 metadata column:
#        seqnames                 ranges strand |     gene_id
#           <Rle>              <IRanges>  <Rle> | <character>
#      1    chr19 [ 58858172,  58874214]      - |           1
#     10     chr8 [ 18248755,  18258723]      + |          10
#    100    chr20 [ 43248163,  43280376]      - |         100
#   1000    chr18 [ 25530930,  25757445]      - |        1000
#  10000     chr1 [243651535, 244006886]      - |       10000
#    ...      ...                    ...    ... .         ...
#   9991     chr9 [114979995, 115095944]      - |        9991
#   9992    chr21 [ 35736323,  35743440]      + |        9992
#   9993    chr22 [ 19023795,  19109967]      - |        9993
#   9994     chr6 [ 90539619,  90584155]      + |        9994
#   9997    chr22 [ 50961997,  50964905]      - |        9997
#  -------
#  seqinfo: 93 sequences (1 circular) from hg19 genome


transcriptsBy(txdb) # returns list of genes organized as granges objects listing transcripts
#GRangesList object of length 23459:
#$1 
#GRanges object with 2 ranges and 2 metadata columns:
#      seqnames               ranges strand |     tx_id     tx_name
#         <Rle>            <IRanges>  <Rle> | <integer> <character>
#  [1]    chr19 [58858172, 58864865]      - |     70455  uc002qsd.4
#  [2]    chr19 [58859832, 58874214]      - |     70456  uc002qsf.2
#
#$10 
#GRanges object with 1 range and 2 metadata columns:
#      seqnames               ranges strand | tx_id    tx_name
#  [1]     chr8 [18248755, 18258723]      + | 31944 uc003wyw.1
#
#$100 
#GRanges object with 1 range and 2 metadata columns:
#      seqnames               ranges strand | tx_id    tx_name
#  [1]    chr20 [43248163, 43280376]      - | 72132 uc002xmj.3
#
#...
#<23456 more elements>
#-------
#seqinfo: 93 sequences (1 circular) from hg19 genome

transcriptsBy(txdb)[["1"]] # returns the first granges object in the GRangeslist
#GRanges object with 2 ranges and 2 metadata columns:
#      seqnames               ranges strand |     tx_id     tx_name
#         <Rle>            <IRanges>  <Rle> | <integer> <character>
#  [1]    chr19 [58858172, 58864865]      - |     70455  uc002qsd.4
#  [2]    chr19 [58859832, 58874214]      - |     70456  uc002qsf.2
#  -------
#  seqinfo: 93 sequences (1 circular) from hg19 genome

transcripts(edb) # returns large Granges object
#GRanges object with 215647 ranges and 6 metadata columns:
#                  seqnames               ranges strand |           tx_id                         tx_biotype tx_cds_seq_start tx_cds_seq_end         gene_id         tx_name
                     <Rle>            <IRanges>  <Rle> |     <character>                        <character>        <integer>      <integer>     <character>     <character>
  #ENST00000456328        1       [11869, 14409]      + | ENST00000456328               processed_transcript             <NA>           <NA> ENSG00000223972 ENST00000456328
  #ENST00000515242        1       [11872, 14412]      + | ENST00000515242 transcribed_unprocessed_pseudogene             <NA>           <NA> ENSG00000223972 ENST00000515242
  #ENST00000518655        1       [11874, 14409]      + | ENST00000518655 transcribed_unprocessed_pseudogene             <NA>           <NA> ENSG00000223972 ENST00000518655
  #ENST00000450305        1       [12010, 13670]      + | ENST00000450305 transcribed_unprocessed_pseudogene             <NA>           <NA> ENSG00000223972 ENST00000450305
  #ENST00000438504        1       [14363, 29370]      - | ENST00000438504             unprocessed_pseudogene             <NA>           <NA> ENSG00000227232 ENST00000438504
   #           ...      ...                  ...    ... .             ...                                ...              ...            ...             ...             ...
  #ENST00000420810        Y [28695572, 28695890]      + | ENST00000420810               processed_pseudogene             <NA>           <NA> ENSG00000224240 ENST00000420810
  #ENST00000456738        Y [28732789, 28737748]      - | ENST00000456738             unprocessed_pseudogene             <NA>           <NA> ENSG00000227629 ENST00000456738
  #ENST00000435945        Y [28740998, 28780799]      - | ENST00000435945             unprocessed_pseudogene             <NA>           <NA> ENSG00000237917 ENST00000435945
  #ENST00000435741        Y [28772667, 28773306]      - | ENST00000435741               processed_pseudogene             <NA>           <NA> ENSG00000231514 ENST00000435741
  #ENST00000431853        Y [59001391, 59001635]      + | ENST00000431853               processed_pseudogene             <NA>           <NA> ENSG00000235857 ENST00000431853
  #-------
  #seqinfo: 273 sequences from GRCh37 genome




```

# Citations (package list)
org.Hs.eg.db  
BSgenome.Hsapiens.UCSC.hg19  
TxDb.Hsapiens.UCSC.hg19.knownGene  
AnnotationDbi  
EnsDb.Hsapiens.v75

## helpful links
https://bioconductor.org/packages/3.7/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf  
https://www.bioconductor.org/help/course-materials/2010/EMBL2010/GenomicRanges.pdf  
https://bioconductor.org/packages/3.7/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf  
https://genome.ucsc.edu/  
http://www.genecards.org/  
 
#
