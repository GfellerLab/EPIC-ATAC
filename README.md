EPIC-ATAC
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Description

Package implementing the EPIC method to estimate the proportion of
immune, stromal, endothelial and cancer or other cells from bulk
chromatin accessibility data (ATAC-Seq). EPIC-ATAC is based on reference
ATAC-Seq profiles for the main non-malignant cell types and it predicts
the proportion of these cells and of the remaining “other cells” (that
are mostly cancer cells) for which no reference profile is given.

This framework is an extension of EPIC previously developed for bulk
RNA-Seq data deconvolution and described in the publication from *Racle
et al., 2017* available at <https://elifesciences.org/articles/26476>.
We recommend loading the EPIC package for RNA-Seq deconvolution and
EPICATAC for ATAC-Seq deconvolution.

## Installation

To install EPICATAC, run the following lines in R:

``` r
install.packages("devtools")
devtools::install_github("GfellerLab/EPIC-ATAC", build_vignettes=TRUE)
```

## Usage

The main function in this package is `EPIC_ATAC`. It needs as input a
matrix of TPM-like counts for ATAC-Seq data, i.e., counts normalized for
peaks length and sample depth and rescaled so that the counts of each
sample sum to 10^6 (or the TPM (or RPKM) for gene expression data) from
the samples to deconvolve. Note that you can obtain TPM-like transformed
data from counts data using the function `get_TPMlike_counts`. One can
define the reference profiles to use with the `reference` option (for
ATAC-Seq deconvolution, use the `atacRef_PBMC` for PBMC samples and the
`atacRef_TME` for cancer samples).

``` r
# library(EPICATAC) ## If the package isn't loaded (or use EPICATAC::EPIC_ATAC and so on).
out <- EPIC_ATAC(bulk = PBMC_ATAC_data$counts, reference = atacRef_PBMC, ATAC = TRUE)
```

`out` is a list containing the various cell fractions in each samples as
well as some *data.frame* of the goodness of fit.

Various other options are available and are well documented in the help
pages from EPIC_ATAC:

``` r
?EPICATAC::EPIC_ATAC
?EPICATAC::EPICATAC.package
```

## License

EPICATAC can be used freely by academic groups for non-commercial
purposes. The product is provided free of charge, and, therefore, on an
“*as is*” basis, without warranty of any kind. Please read the file
“*LICENSE*” for details.

If you plan to use EPICATAC (version 1.0) in any for-profit application,
you are required to obtain a separate license. Note that EPICATAC
requires the users to install the following packages: GenomicRanges,
rtracklayer, tidyr, GenomeInfoDb, S4Vectors. The user should verify that
their usage is in accordance with the licenses associated to these
packages. To do so, please contact Nadette Bulgin (<nbulgin@lcr.org>) at
the Ludwig Institute for Cancer Research Ltd.

## Contact information

Aurélie Gabriel (<aurelie.gabriel@unil.ch>), Julien Racle
(<julien.racle@unil.ch>), and David Gfeller (<david.gfeller@unil.ch>).

## FAQ

##### What do the “*otherCells*” represent?

- EPIC-ATAC predicts the proportions of the various cell types for which
  we have reference profiles (and corresponding cell-type specific
  markers). But, depending on the bulk sample, it is possible that some
  other cell types are present in the bulk but absent from the reference
  profiles. EPIC-ATAC returns the proportion of these remaining cells
  under the name “*otherCells*”. In the case of tumor samples, most of
  these other cells would certainly correspond to the cancer cells, but
  it could be that there are also some stromal cells or epithelial cells
  for example.

##### I receive an error message “*attempt to set ‘colnames’ on an object with less than two dimensions*”. What can I do?

- This is certainly that some of your data is a vector instead of a
  matrix. Please make sure that your bulk data is in the form of a
  matrix (and also your reference profiles if using custom ones).

##### What is the meaning of the warning message telling that some cell-types have few marker peaks in common with the bulk data?

- This warning means that for some cell types, less than 5 cell-type
  specific marker peaks are in common with the open regions of the bulk
  data given as input for deconvolution. It is likely that these
  cell-types are not present in the bulk samples or in low proportions
  and that no open regions specific to these cell-types were detected
  when performing peak calling at the bulk level. High proportions of
  these cell-types predicted by EPIC-ATAC should be considered with
  caution.

##### What is the meaning of the error message telling that none of the marker peaks are in common with the bulk data?

- In this case, we suggest that you extract raw counts from your data
  for each of our marker peaks. To do that, you can use featureCounts
  and the following command line:

``` bash
featureCounts -F SAF -O --fracOverlap 0.2 -T 1 -p -a markerPeaks.saf -o featureCounts.txt bam_files
```

In the previous command line, the bam_files correspond to your samples
bam files and the file “*markerPeaks.saf*” contains the coordinates
(hg38) of our marker peaks. This file is available in the EPICATAC
package and the file path can be retrieved using:

``` r
system.file("extdata", "markerPeaks.saf", package="EPICATAC")
```

The “*featureCounts.txt*” file can then be passed to the
get_TPMlike_counts function to obtain normalized counts that can be used
for deconvolution.

``` r
# library(EPICATAC) ## If the package isn't loaded (or use EPICATAC::EPIC_ATAC and so on).
counts = read.table("featureCounts.txt", skip = 1, header = T)
rownames(counts) = counts$Geneid
counts = counts[, c(7:ncol(counts))]
normalized_counts <- get_TPMlike_counts(counts)
out <- EPIC_ATAC(bulk = normalized_counts, reference = atacRef_PBMC, ATAC = TRUE)
```

##### I receive a warning message that “*the optimization didn’t fully converge for some samples*”. What does it mean?

- When estimating the cell proportions EPIC-ATAC performs a least square
  regression between the observed gene expression or peaks accessibility
  of the cell-type specific markers and the expression or accessibility
  of these markers predicted based on the estimated proportions and
  reference profiles of the various cell types.

  When such a warning message appears, it means that the optimization
  didn’t manage to fully converge for this regression, for some of the
  samples. You can then check the “*fit.gof\$convergeCode*” (and
  possibly also “*fit.gof\$convergeMessage*”) that is outputted by
  EPIC_ATAC alongside the cell proportions. This will tell you which
  samples had issue with the convergence (a value of 0 means it
  converged ok, while other values are errors/warnings, their meaning
  can be found in the help of “*optim*” (or “*constrOptim*”) function
  from R (from “*stats*” package) which is used during the optimization
  and we simply forward the message it returns).

  The error code that usually comes is a “1” which means that the
  maximum number of iterations has been reached in the optimization.
  This could mean there is an issue with the bulk data that maybe don’t
  completely follow the assumption of equation (1) from our manuscript.
  From our experience, it seems in practice that even when there was
  such a warning message the proportions were predicted well, it is
  maybe that the optimization just wants to be *too precise*, or maybe
  few of the markers didn’t match well but the rest of markers could be
  used to have a good estimate of the proportions.

  If you have some samples that seem to have strange results, it could
  however be useful to check that the issue is not that these samples
  didn’t converge well. To be more conservative you could also remove
  all the samples that didn’t converge well as these are maybe outliers,
  if it is only a small fraction from your original samples. Another
  possibility would be to change the parameters of the optim/constrOptim
  function to allow for more iterations or maybe a weaker tolerance for
  the convergence, this numpber of iterations can be provided using the
  nb_iter parameter in the EPIC_ATAC function.

##### Who should I contact in case of a technical or other issue?

- Aurélie Gabriel (<aurelie.gabriel@unil.ch>). Please provide as much
  details as possible and ideally send also an example input file
  (and/or reference profiles) that is causing the issue.
