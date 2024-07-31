# CloneTrackerXP_exstraction
Cellecta CloneTraker XP barcode exstraction and manipulation used in [Pellecchia et al. 2024](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01327-2)

``` r
#####################
# Global parameters #
#####################
bc14Path = "/path/to/BC14.txt"
bc30Path = "/path/to/BC30.txt"
bc14 = read.delim(bc14Path,stringsAsFactors = F)
bc30 = read.delim(bc30Path,stringsAsFactors = F)
Rcpp::sourceCpp("/path/to/hamming.cpp", showOutput = T)

source("/path/to/lineage_libraries.R")


###################
# MAIN            #
###################

# Example 1: read BC from scData
# read1 contain cell barcode, while read2 can contain cellecta barcode

df = readBC_scRNA(fqin="path/to/SAMPLE_R1.fastq.gz",
                  fqin2="path/to/SAMPLE_R2.fastq.gz",
                  sample="SAMPLE_NAME",max.mismatch=1,
                  condition="WT",
                  outDir = "/path/to/results.rds")

# Example 2: Read BC from Quant-seq Bulk data
# single end sequencing of the amplicon containg the barcode

df = readBC_QSQ(fqin2="path/to/SAMPLE_R1.fastq.gz",
                  sample="SAMPLE_NAME",
                  max.mismatch=1,
                  condition="WT",
                  outDir = "/path/to/results.rds",
                  fqchunks = (10^6)/2)

```
