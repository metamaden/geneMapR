# geneMappR
Repo containing info and workhorse scripts to aid gene ID mapping and coordinate retrieval.

# Mapping a Gene
Sequence -> Codon -> Transcript -> Isoform -> Isoform Cluster -> Gene Coordinate Range

# Gene IDs
Gene identifiers correspond to any number of different molecular data types, including different types of transcripts.

# Genome Builds
The genome build (and affiliate patch number) corresponds to a designation like "hg##" ([UCSC genome build](https://genome.ucsc.edu/FAQ/FAQreleases.html)) or "GRCh##" (Genome Assembly from the NCBI [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)). Genomes can be stored a a scaffold sequence, [annotated](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/) database, or some combination of the two.

# Programmic Mapping
TBD (to discuss programmatic strategies for gene mapping)

# Citations (package list)
org.Hs.eg.db
BSgenome.Hsapiens.UCSC.hg19
TxDb.Hsapiens.UCSC.hg19.knownGene
AnnotationDbi
EnsDb.Hsapiens.v75
