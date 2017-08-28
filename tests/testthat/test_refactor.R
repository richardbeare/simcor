## Testing whether my refactoring has broken anything
library(simcor)
source("simcor_orig.R")
seed <- 759
set.seed(seed)
noise.iden <- noisecor(diag(15), epsilon = .5, eidim = 2)
set.seed(seed)
noise.iden.orig <- noisecor.orig(diag(15), epsilon = .5, eidim = 2)

set.seed(seed)
s1 <- simcorTop.orig()
set.seed(seed)
s2 <- simcorTop.orig()

set.seed(seed)
b1 <- simcor()
set.seed(seed)
b2 <- simcor.orig()

set.seed(seed)
h1 <- simcor.H.orig()
set.seed(seed)
h2 <- simcor.H()


test_that("noise addition OK", {
  expect_equal(noise.iden, noise.iden.orig)
}
          )

test_that("Toeplitz structure OK", {
  expect_equal(s1, s2)
}
)

test_that("Block structure OK", {
  expect_equal(b1, b2)
}
)

test_that("Hub structure OK", {
  expect_equal(h1, h2)
}
)
