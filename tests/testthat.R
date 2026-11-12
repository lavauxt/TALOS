# testthat entry point.
#
# R CMD check runs this file, which in turn invokes test_check() to load
# every test-*.R file under tests/testthat/. The library() call requires
# the package to be installed; in dev mode (devtools::load_all()) the
# testthat::test_local() workflow sources the helpers directly.

library(testthat)
library(TALOS)

test_check("TALOS")
