## ----eval = FALSE-------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("GfellerLab/EPIC-ATAC", build_vignettes=TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  # library(EPICATAC) ## If the package isn't loaded (or use EPICATAC::EPIC_ATAC and so on).
#  out <- EPIC_ATAC(bulk = PBMC_ATAC_data$counts, reference = atacRef_PBMC, ATAC = TRUE, withOtherCells = F)

## ----eval = FALSE-------------------------------------------------------------
#  ?EPICATAC::EPIC_ATAC
#  ?EPICATAC::EPICATAC.package

## ----eval = FALSE-------------------------------------------------------------
#  # library(EPICATAC) ## If the package isn't loaded (or use EPICATAC::EPIC_ATAC and so on).
#  counts = read.table("featureCounts.txt", skip = 1, header = T)
#  rownames(counts) = counts$Geneid
#  counts = counts[, c(7:ncol(counts))]
#  normalized_counts <- get_TPMlike_counts(counts)
#  out <- EPIC_ATAC(bulk = normalized_counts, reference = atacRef_PBMC, ATAC = TRUE)

