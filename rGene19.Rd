\name{rGene19}
\alias{rGene19}
\title{
Grab hg19 gene info. 
}
\description{
Grab gene IDs and coordinates (Ensembl and/or UCSC) overlapping a query region of the hg19 human genome build.
}
\usage{
rGene19(input.range,searchwindow=T,windowsize=2e+5,get.tx.coor=T,get.en.coor=T,show.module.versions=T,returngr=T)
}

\arguments{
  \item{input.range}{
 Genomic coordinate range to query, formatted as a GRanges R object with "chr#" seqname format. See example. 
}
  \item{searchwindow}{
(T/F) whether to use a search window up/downstream of inputted range
}
  \item{windowsize}{
  How far to look up and downstream, if searchwindow==T
 }
  \item{get.tx.coor}{
 (T/F) Whether to grab UCSC gene coordinates, from TXDB package
}
  \item{get.en.coor}{
 (T/F) Whether to grab Ensembl DB coordinates, from ENSDB package.
}
  \item{returngr}{
  (T/F) Whether to return gene info as granges objects
 }
  \item{show.module.versions}{
 (T/F) Whether to show and store module versions used in query
}
}
\details{
This function makes use of the R GenomicRanges and GenomicFeatures packages, in addition to core database packages built on the AnnotationDB class.
}
\value{
List of specified return object.
  \item{(granges query|Ensembl df/gr|Txdb df/gr|R package versions)} ...
}
\references{
}
\author{
Sean Maden, smaden@fredhutch.org
}
\note{
}
\seealso{
}
\examples{
# ranges object should be GRanges class R object
# easily assemble this as such 
irx <- makeGRangesFromDataFrame(data.frame(start=100000,end=200000,chr="chr1"))
rGene19(irx)

# returns:
#$query.gr
#GRanges object with 1 range and 0 metadata columns:
#      seqnames            ranges strand
#         <Rle>         <IRanges>  <Rle>
#  [1]     chr1 [-100000, 400000]      *
#  -------
#  seqinfo: 1 sequence from an unspecified genome; no seqlengths

#$tx.gr
#GRanges object with 4 ranges and 2 metadata columns:
#            seqnames           ranges strand |     gene_id  genesymbol
#               <Rle>        <IRanges>  <Rle> | <character> <character>
#  100287102     chr1 [ 11874,  14409]      + |   100287102     DDX11L1
#     653635     chr1 [ 14362,  29961]      - |      653635      WASH7P
#     729737     chr1 [134773, 140566]      - |      729737   LOC729737
#      79501     chr1 [ 69091,  70008]      + |       79501       OR4F5
#  -------
#  seqinfo: 93 sequences (1 circular) from hg19 genome

#$en.gr
#GRanges object with 4 ranges and 2 metadata columns:
#            seqnames           ranges strand |     gene_id  genesymbol
#               <Rle>        <IRanges>  <Rle> | <character> <character>
#  100287102     chr1 [ 11874,  14409]      + |   100287102     DDX11L1
#     653635     chr1 [ 14362,  29961]      - |      653635      WASH7P
#     729737     chr1 [134773, 140566]      - |      729737   LOC729737
#      79501     chr1 [ 69091,  70008]      + |       79501       OR4F5
#  -------
#  seqinfo: 93 sequences (1 circular) from hg19 genome

#$Rpackage.versions
#$Rpackage.versions$org.Hs.eg.db
#[1] "3.4.1"
#
#$Rpackage.versions$TxDb.Hsapiens.UCSC.hg19.knownGene
#[1] "3.2.2"
#
#$Rpackage.versions$EnsDb.Hsapiens.v75
#[1] "2.1.0"
}
