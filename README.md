# geneMapR
Repo containing info and workhorse scripts to aid gene ID mapping and coordinate retrieval.

# Mapping a Gene
Sequence -> Codon -> Transcript -> Isoform -> Isoform Cluster -> Gene Coordinate Range

# Gene IDs
Gene identifiers correspond to any number of different molecular data types, including different types of transcripts.

# Genome Builds
The genome build (and affiliate patch number) corresponds to a designation like "hg##" ([UCSC genome build](https://genome.ucsc.edu/FAQ/FAQreleases.html)) or "GRCh##" (Genome Assembly from the NCBI [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)). Genomes can be stored as a scaffold sequence, [annotated](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/) database, or some combination of the two.

# Programmic Mapping: R
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
head(select(edb, keys=genex, columns=c("TXID", "TXSEQSTART", "TXBIOTYPE","SEQNAME","GENESEQSTART","GENESEQEND"),keytype="SYMBOL"))
#   SYMBOL            TXID TXSEQSTART       TXBIOTYPE SEQNAME GENESEQSTART GENESEQEND
# 1   TP53 ENST00000413465    7565097  protein_coding      17      7565097    7590856
# 2   TP53 ENST00000359597    7569404  protein_coding      17      7565097    7590856
# 3   TP53 ENST00000504290    7571720 retained_intron      17      7565097    7590856
# 4   TP53 ENST00000510385    7571720 retained_intron      17      7565097    7590856
# 5   TP53 ENST00000504937    7571720 retained_intron      17      7565097    7590856
# 6   TP53 ENST00000269305    7571720  protein_coding      17      7565097    7590856
```
NOTE: gene coordinates are not uncontroversial. Even controlling for genome build, different databases can yield varying exact coordinate ranges. These should, however, cluster predictably and according to annotated transcripts and possible variation in transcript annotations, database updates, curation practices, etc. Thus it often pays to check multiple resources. See helpful links below for some options.

# Citations (package list)
org.Hs.eg.db
BSgenome.Hsapiens.UCSC.hg19
TxDb.Hsapiens.UCSC.hg19.knownGene
AnnotationDbi
EnsDb.Hsapiens.v75

## helpful links
https://bioconductor.org/packages/3.7/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
https://www.bioconductor.org/help/course-materials/2010/EMBL2010/GenomicRanges.pdf
https://genome.ucsc.edu/
http://www.genecards.org/

#
