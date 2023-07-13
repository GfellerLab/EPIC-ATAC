## ---- eval = FALSE------------------------------------------------------------
#  # library(EPICatac) ## If the package isn't loaded (or use EPICatac::EPIC and so on).
#  out <- EPIC(bulk = PBMC_ATAC_data, reference = BRef_ATAC, ATAC = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  ?EPICatac::EPIC
#  ?EPICatac::EPICatac.package

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("GfellerLab/EPICatac", build_vignettes=TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  # library(EPICatac) ## If the package isn't loaded (or use EPICatac::EPIC and so on).
#  counts = read.table("featureCounts.txt", skip = 1, header = T)
#  rownames(counts) = counts$Geneid
#  counts = counts[, c(7:ncol(counts))]
#  normalized_counts <- get_TPMlike_counts(counts)
#  out <- EPIC(bulk = normalized_counts, reference = BRef_ATAC, ATAC = TRUE)

