## Testing if Stata output checks out with R
require("testthat")
require("ebci")
require("RStata")

## Run stata do file
#chooseStataBin()
options("RStata.StataVersion"=17)
options("RStata.StataPath"="\"C:\\Program Files\\Stata17\\StataSE-64\"")
stata("ebreg_test.do")

# Load the output from Stata
eb_out1 <- haven::read_dta(file = "eb_out1.dta")
eb_out1 <- eb_out1[is.na(eb_out1$EB_df1)==FALSE, ]
eb_out2 <- haven::read_dta(file = "eb_out2.dta")
eb_out2 <- eb_out2[is.na(eb_out2$EB_df1)==FALSE, ]
eb_out3 <- haven::read_dta(file = "eb_out3.dta")
eb_out3 <- eb_out3[is.na(eb_out3$EB_df1)==FALSE, ]

############
# Load data and test if output are equal
load("cz.rda")

test_that("Compare ebreg1", {
   eb_out1 <- as.matrix(unname(eb_out1))
   Rout <- ebci(theta25 ~ stayer25, cz[cz$state=="NY", ],
                se25, 1/se25^2, alpha=0.05, fs_correction = "none", wopt=TRUE)
   REBdf <- as.matrix(unname(Rout$df))
   diff <- eb_out1 - REBdf
   expect_equal(max(abs(diff)), 0, tolerance=1e-4)
})

test_that("Compare ebreg2", {
   eb_out2 <- as.matrix(unname(eb_out2))
   Rout <- ebci(theta25 ~ stayer25, cz[cz$state=="NY", ],
                se25, 1/se25^2, alpha=0.1, wopt=TRUE)
   REBdf <- as.matrix(unname(Rout$df))
   diff <- eb_out2 - REBdf
   expect_equal(max(abs(diff)), 0, tolerance=1e-4)
})

test_that("Compare ebreg3", {
   eb_out3 <- as.matrix(unname(eb_out3))
   Rout <- ebci(theta25 ~ stayer25, cz[cz$state=="FL", ],
                se25, 1/se25^2, alpha=0.1, kappa = 10, wopt=TRUE)
   REBdf <- as.matrix(unname(Rout$df))
   diff <- eb_out3 - REBdf
   expect_equal(max(abs(diff)), 0, tolerance=1e-4)
})
