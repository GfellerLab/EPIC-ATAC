## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("GfellerLab/EPIC-ATAC", build_vignettes=TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  # library(EPICATAC) ## If the package isn't loaded (or use EPICATAC::EPIC and so on).
#  out <- EPIC(bulk = PBMC_ATAC_data$counts, reference = BRef_ATAC, ATAC = TRUE, withOtherCells = F)

## ---- eval = FALSE------------------------------------------------------------
#  ?EPICATAC::EPIC
#  ?EPICATAC::EPICATAC.package

## ---- eval = FALSE------------------------------------------------------------
#  # library(EPICATAC) ## If the package isn't loaded (or use EPICATAC::EPIC and so on).
#  counts = read.table("featureCounts.txt", skip = 1, header = T)
#  rownames(counts) = counts$Geneid
#  counts = counts[, c(7:ncol(counts))]
#  normalized_counts <- get_TPMlike_counts(counts)
#  out <- EPIC(bulk = normalized_counts, reference = BRef_ATAC, ATAC = TRUE)

