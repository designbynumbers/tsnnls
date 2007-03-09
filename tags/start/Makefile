
all: libtsnnls test

CC=gcc
CFLAGS=-O3

TAUCS_OBJS=tsnnls/taucs_basic/taucs_logging.o \
	tsnnls/taucs_basic/taucs_malloc.o \
	tsnnls/taucs_basic/taucs_ccs_order.o \
	tsnnls/taucs_basic/taucs_ccs_ops.o \
	tsnnls/taucs_basic/taucs_vec_base.o \
	tsnnls/taucs_basic/taucs_complex.o \
	tsnnls/taucs_basic/colamd.o \
	tsnnls/taucs_basic/amdbar.o \
	tsnnls/taucs_basic/amdexa.o \
	tsnnls/taucs_basic/amdtru.o \
	tsnnls/taucs_basic/genmmd.o \
	tsnnls/taucs_basic/taucs_timer.o \
	tsnnls/taucs_basic/taucs_sn_llt.o \
	tsnnls/taucs_basic/taucs_ccs_base.o

TSNNLS_OBJS=tsnnls/tlsqr.o \
        tsnnls/tsnnls.o \
	tsnnls/lsqr.o

OBJS=$(TAUCS_OBJS) $(TSNNLS_OBJS)

ifeq ($(shell uname),Darwin)
	LDFLAGS=-framework vecLib
else
	LDFLAGS= -Lexternal_libraries -llapacklinux -lblaslinux -lf77blas -lcblas -latlas -lg2c -lm
endif

clean:
	(cd tsnnls ; make clean)
	rm -f tsnnls_test libtsnnls.a

tsnnls: $(OBJS)
	(cd tsnnls ; make all "CFLAGS=$(CFLAGS)" )

libtsnnls: tsnnls
	ar cru libtsnnls.a $(OBJS)
	ranlib libtsnnls.a

test: $(OBJS)
	$(CC) $(CFLAGS) tsnnls_test.c -o tsnnls_test -L. -ltsnnls $(LDFLAGS)
	./tsnnls_test test_files/full_example_200 0 0
	./tsnnls_test test_files/sparse_tests 0 1

install: 
	cp ./tsnnls/tsnnls.h /usr/local/include/tsnnls.h
	cp ./libtsnnls.a /usr/local/lib/libtsnnls.a



