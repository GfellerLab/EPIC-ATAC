context("Testing EPIC function")

test_that("Test of bad inputs", {
          a <- matrix(1:20, nrow=5)
          expect_error(EPIC_ATAC(bulk=a, ATAC = F), "There are only 0 signature")
          refStr <- "unknownRef"
          expect_error(EPIC_ATAC(bulk=melanoma_data$counts, reference=refStr, ATAC = F),
                       paste0(refStr, ".* not part of the allowed references"))
          tRef <- BRef; tRef$sigGenes <- NULL
          expect_error(EPIC_ATAC(bulk=melanoma_data$counts, reference=tRef, ATAC = F),
                       "needs to contain .* 'sigGenes'")
          tRef <- BRef;
          tRef$refProfiles.var <- tRef$refProfiles.var[
            nrow(tRef$refProfiles.var):1,]
          expect_error(EPIC_ATAC(bulk=melanoma_data$counts, reference=tRef, ATAC = F),
                       "dimnames of .*refProfiles.*refProfiles.var.*same")
})

test_that("Test for correct result on melanoma data with default input", {
  testFract <- EPIC_ATAC(melanoma_data$counts, ATAC = F, nb_iter = 500)$cellFractions
  expect_equal(testFract, melanoma_data$cellFractions.pred)
})


