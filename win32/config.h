/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

#define ADD_PREFIX(x) x

/* Name of daxpy after mangling */
#define DAXPY_F77 ADD_PREFIX(daxpy_)

/* Name of ddot after mangling */
#define DDOT_F77 ADD_PREFIX(ddot_)

/* Name of dgels after mangling */
#define DGELS_F77 ADD_PREFIX(dgels_)

/* Name of dgemm after mangling */
#define DGEMM_F77 ADD_PREFIX(dgemm_)

/* Name of dgemv after mangling */
#define DGEMV_F77 ADD_PREFIX(dgemv_)

/* Name of dnrm2 after mangling */
#define DNRM2_F77 ADD_PREFIX(dnrm2_)

/* Name of dpotrf after mangling */
#define DPOTRF_F77 ADD_PREFIX(dpotrf_)

/* Name of dscal after mangling */
#define DSCAL_F77 ADD_PREFIX(dscal_)

/* Name of dsymv after mangling */
#define DSYMV_F77 ADD_PREFIX(dsymv_)

/* Name of dsyrk after mangling */
#define DSYRK_F77 ADD_PREFIX(dsyrk_)

/* Name of dtrsm after mangling */
#define DTRSM_F77 ADD_PREFIX(dtrsm_)

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the <argtable2.h> header file. */
/* #undef HAVE_ARGTABLE2_H */

/* Define to 1 if you have the <atlas/clapack.h> header file. */
/* #undef HAVE_ATLAS_CLAPACK_H */

/* ATLAS (instead of full) LAPACK */
/* #undef HAVE_ATLAS_LAPACK */

/* Define if you have a BLAS library. */
/* #undef HAVE_BLAS */

/* Define to 1 if you have the <clapack.h> header file. */
/* #undef HAVE_CLAPACK_H */

/* Define to 1 if you have the <ctype.h> header file. */
#define HAVE_CTYPE_H 1

/* Defined if we are in the Apple environment. */
#undef HAVE_DARWIN

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Defined if we have full clapack. */
#undef HAVE_FULL_CLAPACK

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
/* #undef HAVE_LAPACK */

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the `memcpy' function. */
#define HAVE_MEMCPY 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have the `rand' function. */
#define HAVE_RAND 1

/* Define to 1 if you have the `srand' function. */
#define HAVE_SRAND 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/resource.h> header file. */
#undef HAVE_SYS_RESOURCE_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#undef HAVE_SYS_TIME_H

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `time' function. */
#define HAVE_TIME 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#undef HAVE_UNISTD_H

/* Define to 1 if you have the <vecLib/clapack.h> header file. */
#undef HAVE_VECLIB_CLAPACK_H

/* Name of package */
#define PACKAGE "libtsnnls"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "cantarel@math.uga.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "libtsnnls"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "libtsnnls 2.01.2"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "libtsnnls"

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.01.2"

/* The size of `double *', as computed by sizeof. */
#define SIZEOF_DOUBLE_P 4

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long int', as computed by sizeof. */
#define SIZEOF_LONG_INT 4

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "2.01.2"

/* Defined if we have the argtable2 library */
#undef WITH_ARGTABLE2 

/* Define if using the dmalloc debugging malloc package */
/* #undef WITH_DMALLOC */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
