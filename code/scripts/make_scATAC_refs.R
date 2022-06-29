# Following the data structures vignette at
# https://github.com/timoast/signac/blob/master/vignettes/data_structures.Rmd

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(Signac)
library(GenomicRanges)

# get seqinfo
seqinfo_hg19 <- Seqinfo("hg19")

# convert EnsDb to GRanges
gene.ranges <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# convert to UCSC style
seqlevelsStyle(gene.ranges) <- "UCSC"
genome(gene.ranges) <- "hg19"

annotation_hg19 <- gene.ranges

save(seqinfo_hg19, annotation_hg19, file = here("data/scATAC/references/hg19_genome_info.Rda"))

# promoter coordinates
gene.ranges <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

promoter.coords <- promoters(gene.ranges, upstream = 2500, downstream = 2500)

save(promoter.coords, file = here("data/scATAC/references/promoter_coordinates.Rda"))
