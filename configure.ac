#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl


AC_PREREQ([2.57])
AC_COPYRIGHT([
------------------------------- libtsnnls.a ---------------------------------

                  Copyright (c) 2006 University of Georgia.

 Distributed under terms of the GNU Lesser General Public License version 3
                       or later (see COPYING and COPYING.LESSER).
])

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(libtsnnls, 2.4.5, cantarel@math.uga.edu)
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE

#
# The version number below has a very specific meaning.  It is in the form
# a:b:c.  The first number, a, is the interface version number.  Any time we
# add or remove a function, or change the order of parameters or change the way
# that the library reads or writes files, this number needs to be incremented.
# The second number, b, is the implementation number *of this interface*.  When
# a changes, b goes back to zero.  This is the number that needs to change
# every time we put out a new library distribution.  The third number, c, is
# the interface backward-compatibility count.  Suppose that we simply add a
# needed function to the library.  The new interface still provides everything
# that the old version provided and in the same way.  So the
# backward-compatibility count gets incremented.  Of course it gets reset to 0
# when a new interface is not binary-compatible with the old one.
#
# Due to various broken ideas, we have to have two versions of this number, one
# with colons below and one with periods above.  Please keep them in sync!
#
AC_SUBST(LIB_VERSION, 2:4:5)
AC_CONFIG_SRCDIR([configure])
AC_CONFIG_HEADERS([config.h])

# We need to override the default cflags here. We don't want to compile with -02 by default
CFLAGS=`echo $CFLAGS | sed s/-O2//g`

# Checks for programs.
LT_INIT
AC_PROG_CC([clang gcc])
AM_PROG_CC_C_O

# Checks for libraries
AC_CHECK_LIB([m],[cos])

# Checks for functions
AC_CHECK_FUNCS([memset memcpy rand srand time isnormal isnan])

# Checks for FORTRAN.
#AC_PROG_FC([gfortran g77])		
#AC_FC_LIBRARY_LDFLAGS
#LIBS="$LIBS $FCLIBS -lm"

# Checks for Darwin.
# AC_CANONICAL_HOST
#AC_MSG_CHECKING(how to link with BLAS and LAPACK)

dnl Tsnnls is now assuming that you're building with OpenBlas. This entails
dnl converting from the clapack.h header (giving C prototypes for lapack functions)
dnl by doing whatever f2c did) to the lapacke.h header, which is at least
dnl *maintained* and attempts to provide a consistent C interface, so we have
dnl a shot at keeping everything going.

AC_SEARCH_LIBS([cblas_dgemv],[openblas], ,AC_MSG_ERROR([tsnnls 2.4.5 and higher require OpenBlas]))

#LAPACK_LIBS="-lopenblas"

dnl case "${host_os}" in
dnl   darwin*)
dnl     echo -framework accelerate
dnl     AM_CONDITIONAL([DARWIN],[test 1 = 1])
dnl     AM_CONDITIONAL([HAVE_LAPACK],[test x1 = x1])
dnl     AC_DEFINE([HAVE_DARWIN],[1],[Defined if we are in the Apple environment.])
dnl     # In this case, we need to use -framework accelerate (instead of the older
dnl     # -framework vecLib. Per the Apple docs, -framework is a linker option and
dnl     # not a compiler option, so it goes in LDFLAGS or as a library. Either
dnl     # seems to work.
dnl     LIBS="$LIBS -framework accelerate"

dnl     # Now there's a lot of lazy symbol binding going on in clang, so this
dnl     # actually works in the end even though we don't find the clapack header
dnl     # at include time
    
dnl #   LDFLAGS="-framework accelerate $LDFLAGS"
dnl #   AC_SUBST([LAPACK_LIBS],["-framework accelerate"])
dnl     AC_SUBST([BLAS_LIBS],[" "])
dnl     # $FLIBS was set above by AC_F77_PROG  
dnl     # How can I work around the damaged g77 on Darwin by not linking lcrt2.o?
dnl     # By editing flibs to remove -lcrt2.o, of course!  
dnl     # (I realize that this is a total hack.)

dnl     FLIBS=`echo $FLIBS | sed s/-lcrt2.o//g` 
dnl     ;;
dnl   *)
dnl     echo linux
dnl     AM_CONDITIONAL([DARWIN],[test 0 == 1])
dnl     # Gets ready to compile with 

dnl     ACX_BLAS
dnl     if test "$acx_blas_ok" = no; then 
dnl 	  AC_MSG_ERROR([tsnnls requires an F77BLAS. The ATLAS project
dnl will provide these on any given system if no other is
dnl available. Configure has searched the usual locations and can't find
dnl one. You can use --with-blas=<lib> to point to a blas library in a
dnl nonstandard spot. This test may fail even if you have BLAS if you
dnl don't have a FORTRAN compiler (preferably g77) on your default path.])
dnl     fi

dnl     # Now look for a full (or minimal ATLAS) lapack
dnl     AC_FC_FUNC(dpotrf)
dnl     ACX_LAPACK(AM_CONDITIONAL([HAVE_LAPACK],[test x1 = x1]),AC_CHECK_LIB([lapack],$dpotrf,[atlas_lapack=yes 
dnl        LAPACK_LIBS="-llapack -lcblas -lf77blas -latlas"
dnl        AC_DEFINE_UNQUOTED([HAVE_ATLAS_LAPACK],[1],[ATLAS (instead of full) LAPACK])
dnl        AM_CONDITIONAL([HAVE_LAPACK],[test 1 == 0])
dnl        ],[atlas_lapack=no],[-lcblas -lf77blas -latlas $FLIBS -lm]))
dnl     AC_MSG_CHECKING([for minimal ATLAS lapack])
dnl     echo $atlas_lapack
    
dnl   esac

#LIBS="$LAPACK_LIBS"
 
# We'll need to call F77 functions from C, find out how to mangle their names
dnl AC_FC_FUNC([daxpy],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="daxpy_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DAXPY_F77],[$fname],[Name of daxpy after mangling])
dnl AC_MSG_CHECKING([name of daxpy after mangling])
dnl echo "$fname"

dnl AC_FC_FUNC([dsyrk],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dsyrk_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DSYRK_F77],[$fname],[Name of dsyrk after mangling])
dnl AC_MSG_CHECKING([name of dsyrk after mangling])
dnl echo "$fname"

dnl AC_FC_FUNC([dsymv],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dsymv_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DSYMV_F77],[$fname],[Name of dsymv after mangling])
dnl AC_MSG_CHECKING([name of dsymv after mangling])
dnl echo "$fname"

dnl # dgemm

dnl AC_FC_FUNC([dgemm],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dgemm_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DGEMM_F77],[$fname],[Name of dgemm after mangling])
dnl AC_MSG_CHECKING([name of dgemm after mangling])
dnl echo "$fname"

dnl # dgemv

dnl AC_FC_FUNC([dgemv],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dgemv_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DGEMV_F77],[$fname],[Name of dgemv after mangling])
dnl AC_MSG_CHECKING([name of dgemv after mangling])
dnl echo "$fname"

dnl # dtrsm
 
dnl AC_FC_FUNC([dtrsm],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dtrsm_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DTRSM_F77],[$fname],[Name of dtrsm after mangling])
dnl AC_MSG_CHECKING([name of dtrsm after mangling])
dnl echo "$fname"

dnl # dpotrf

dnl AC_FC_FUNC([dpotrf],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dpotrf_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DPOTRF_F77],[$fname],[Name of dpotrf after mangling])
dnl AC_MSG_CHECKING([name of dpotrf after mangling])
dnl echo "$fname"

dnl # dgels

dnl AC_FC_FUNC([dgels],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dgels_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DGELS_F77],[$fname],[Name of dgels after mangling])
dnl AC_MSG_CHECKING([name of dgels after mangling])
dnl echo "$fname"

dnl # dnrm2

dnl AC_FC_FUNC([dnrm2],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dnrm2_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DNRM2_F77],[$fname],[Name of dnrm2 after mangling])
dnl AC_MSG_CHECKING([name of dnrm2 after mangling])
dnl echo "$fname"

dnl # ddot

dnl AC_FC_FUNC([ddot],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="ddot_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DDOT_F77],[$fname],[Name of ddot after mangling])
dnl AC_MSG_CHECKING([name of ddot after mangling])
dnl echo "$fname"

dnl # dscal

dnl AC_FC_FUNC([dscal],[fname])
dnl if test "$fname" = "unknown"; then   # Assume g77 conventions 
dnl   fname="dscal_"
dnl fi
dnl AC_DEFINE_UNQUOTED([DSCAL_F77],[$fname],[Name of dscal after mangling])
dnl AC_MSG_CHECKING([name of dscal after mangling])
dnl echo "$fname"

# Other checks for libraries.
AC_SEARCH_LIBS([arg_parse],[argtable2],[have_argtable2=true],[have_argtable2=false])

if test "x${have_argtable2}" = xtrue; then
   AC_DEFINE([WITH_ARGTABLE2],[1],[Defined if we have the argtable2 library])
   AM_CONDITIONAL([ARGTABLE],[test 1 == 1])
else
   AM_CONDITIONAL([ARGTABLE],[test 0 == 1])
fi

# Checks for header files.
dnl AC_HEADER_STDC

dnl AC_CHECK_HEADERS([clapack.h vecLib/clapack.h Accelerate/Accelerate.h],AC_DEFINE([HAVE_FULL_CLAPACK],[1],[Defined if we have a full clapack.]),AC_DEFINE([HAVE_FULL_CLAPACK],[],[Defined if we have full clapack.]))
dnl AC_CHECK_HEADERS([atlas/clapack.h])
dnl AC_CHECK_HEADERS([string.h float.h stdlib.h stdio.h sys/time.h sys/resource.h limits.h math.h time.h ctype.h malloc.h sys/malloc.h])

AC_CHECK_HEADERS([cblas.h lapacke.h])
AC_CHECK_HEADERS([argtable2.h])
AC_CHECK_HEADERS([malloc.h])
AC_CHECK_FUNCS([mallinfo])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Since we are using Autoconf 2.59, go through some stupidity to try to 
# find which int type is 32 bits. This should be able to be replaced in 
# Autoconf 2.60 with something like

#AC_TYPE_INT32_T   # Lapack/f77 will expect type integer to be a 32-bit int
#AC_TYPE_UINT32_T  # make sure that we have this type (and the unsigned version

AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])
AC_CHECK_SIZEOF([double *])

# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Enable debugging support for dmalloc, if we configure --with-dmalloc
AM_WITH_DMALLOC

# Output
#AC_CONFIG_FILES([tsnnls-libs])    # Write linker info to file. 
AC_CONFIG_FILES([Makefile
		tsnnls/Makefile
		tsnnls/taucs_basic/Makefile])
AC_OUTPUT
