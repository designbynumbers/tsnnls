#
# Makefile.am for /tsnnls directory
#

AM_CFLAGS = -g -Wall

SUBDIRS = taucs_basic

pkginclude_HEADERS = tsnnls.h lsqr.h
noinst_HEADERS = tsnnls_blas_wrappers.h taucs_basic/taucs_common.h
nobase_pkginclude_HEADERS = taucs_basic/taucs.h taucs_basic/colamd.h taucs_basic/f2c.h taucs_basic/taucs_config_build.h taucs_basic/taucs_config_tests.h taucs_basic/taucs_private.h
lib_LTLIBRARIES = libtsnnls.la

libtsnnls_la_SOURCES = taucs_basic/taucs_logging.c \
	taucs_basic/taucs_malloc.c \
	taucs_basic/taucs_ccs_order.c \
	taucs_basic/taucs_ccs_ops.c \
	taucs_basic/taucs_vec_base.c \
	taucs_basic/taucs_complex.c \
	taucs_basic/colamd.c \
	taucs_basic/amdbar.c \
	taucs_basic/amdexa.c \
	taucs_basic/amdtru.c \
	taucs_basic/genmmd.c \
	taucs_basic/taucs_timer.c \
	taucs_basic/taucs_sn_llt.c \
	taucs_basic/taucs_ccs_base.c \
	tlsqr.c \
        tsnnls.c \
	lsqr.c

if HAVE_LAPACK 
  check_PROGRAMS = f77blas_test genb_test
  genb_test_SOURCES = genb_test.c
  genb_test_LDFLAGS = -ltsnnls $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
else 
  check_PROGRAMS = tsnnls_test f77blas_test
endif

tsnnls_test_SOURCES = tsnnls_test.c
tsnnls_test_LDFLAGS = -ltsnnls $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

#nodist_PROGRAMS = cblas_test
#cblas_test_SOURCES = cblas_test.c
#cblas_test_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) 

f77blas_test_SOURCES = f77blas_test.c
f77blas_test_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

#if DARWIN
#tsnnls_test_LDFLAGS = -ltsnnls -framework vecLib
#else
#tsnnls_test_LDFLAGS = -L../external_libraries -ltsnnls -llapacklinux -lblaslinux -lf77blas -lcblas -latlas -lg2c -lm
#endif

