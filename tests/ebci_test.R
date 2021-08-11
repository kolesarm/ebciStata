## Test whether Stata output matches R

## Required R packages: haven, testthat, ebci, RStata

## R Stata needs StataVersion and StataPath, usually in .Rprofile
## options("RStata.StataVersion"=17)
## options("RStata.StataPath"="\"C:\\Program Files\\Stata17\\StataSE-64\"")
## StataPath can also be set using RStata::chooseStataBin()

require("testthat")

RStata::stata("ebreg_test.do")

# Load the output from Stata
eb_out1 <- haven::read_dta(file = "output/eb_out1.dta")
eb_out1 <- eb_out1[is.na(eb_out1$EB_df1)==FALSE, ]
eb_out2 <- haven::read_dta(file = "output/eb_out2.dta")
eb_out2 <- eb_out2[is.na(eb_out2$EB_df1)==FALSE, ]
eb_out3 <- haven::read_dta(file = "output/eb_out3.dta")
eb_out3 <- eb_out3[is.na(eb_out3$EB_df1)==FALSE, ]

############
# Load data and test if output are equal
load("../data/cz.rda")

test_that("Compare ebreg1", {
   eb_out1 <- as.matrix(unname(eb_out1))
   Rout <- ebci::ebci(theta25 ~ stayer25, cz[cz$state=="NY", ], se25, 1/se25^2,
                      alpha=0.05, fs_correction = "none", wopt=TRUE)
   REBdf <- as.matrix(unname(Rout$df))
   diff <- eb_out1 - REBdf
   expect_equal(max(abs(diff)), 0, tolerance=1e-4)
})

test_that("Compare ebreg2", {
   eb_out2 <- as.matrix(unname(eb_out2))
   Rout <- ebci::ebci(theta25 ~ stayer25, cz[cz$state=="NY", ], se25, 1/se25^2,
                      alpha=0.1, wopt=TRUE)
   REBdf <- as.matrix(unname(Rout$df))
   diff <- eb_out2 - REBdf
   expect_equal(max(abs(diff)), 0, tolerance=1e-4)
})

test_that("Compare ebreg3", {
   eb_out3 <- as.matrix(unname(eb_out3))
   Rout <- ebci::ebci(theta25 ~ stayer25, cz[cz$state=="FL", ], se25, 1/se25^2,
                      alpha=0.1, kappa = 10, wopt=TRUE)
   REBdf <- as.matrix(unname(Rout$df))
   diff <- eb_out3 - REBdf
   expect_equal(max(abs(diff)), 0, tolerance=1e-4)
})
