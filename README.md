# geneMappR
Repo containing info and workhorse scripts to aid gene ID mapping and coordinate retrieval.

# Mapping a Gene
Sequence -> Codon -> Transcript -> Isoform -> Isoform Cluster -> Gene Coordinate Range

# Gene IDs
Gene identifiers correspond to any number of different molecular data types, including different types of transcripts.

# Genome Builds
The genome build (and affiliate patch number) corresponds to a designation like "hg##" ([UCSC genome build](https://genome.ucsc.edu/FAQ/FAQreleases.html)) or "GRCh##" (Genome Assembly from the NCBI [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)). Genomes can be stored a a scaffold sequence, [annotated](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/) database, or some combination of the two.

# Programmic Mapping: R
Using R/Bioconductor, id mappings can readily be generated from org.Hs.eg.db package, using the select() interface inherited from AnnotationDbi package (see [org.Hs.rg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) and [Annotation tutorial](https://www.bioconductor.org/help/workflows/annotation/annotation/) Bioconductor pages for more details).

Here is an example building from the [Bioc help archives](https://support.bioconductor.org/p/59549/):
```
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")

library(org.Hs.eg.db)
columns(org.Hs.eg.db)
#[1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     #"EVIDENCEALL" 
#[10] "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"     #   
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

# Citations (package list)
org.Hs.eg.db
BSgenome.Hsapiens.UCSC.hg19
TxDb.Hsapiens.UCSC.hg19.knownGene
AnnotationDbi
EnsDb.Hsapiens.v75
