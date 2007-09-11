#
# Makefile.am for top level
#

EXTRA_DIST = taucs.pdf ./test_files/full_example_200 ./test_files/sparse_tests external_libraries testfull testsparse tsnnls-libs.in
SUBDIRS = tsnnls 

if HAVE_LAPACK
  TESTS = testfull testsparse tsnnls/f77blas_test tsnnls/genb_test
else
  TESTS = testfull testsparse tsnnls/f77blas_test
endif

bin_SCRIPTS = tsnnls-libs
