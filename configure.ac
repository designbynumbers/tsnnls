#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl


AC_PREREQ(2.57)
AC_COPYRIGHT(
------------------------------- libtsnnls.a ---------------------------------

                  Copyright (c) 2006 University of Georgia.

 Distributed under terms of the GNU General Public License version 2.1
                       or later (see COPYING.LIB-2.1).
)

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(libtsnnls, 1.0.1, cantarel@math.uga.edu)
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
AC_SUBST(LIB_VERSION, 1:0:1)
AC_CONFIG_SRCDIR([configure])
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for functions
AC_CHECK_FUNCS([memset memcpy rand srand time])

# Checks for FORTRAN.
AC_PROG_F77		
AC_F77_LIBRARY_LDFLAGS

# Checks for Darwin.
AC_CANONICAL_HOST
AC_MSG_CHECKING(how to link with BLAS and LAPACK)

case "${host_os}" in
  darwin*)
    echo -framework vecLib
    AM_CONDITIONAL([DARWIN],[test 1 = 1])
    AM_CONDITIONAL([HAVE_LAPACK],[test x1 = x1])
    AC_DEFINE([HAVE_DARWIN],[1],[Defined if we are in the Apple environment.])
#    LIBS="$LIBS -framework vecLib"
    AC_SUBST([LAPACK_LIBS],["-framework vecLib"])
    AC_SUBST([BLAS_LIBS],[" "])
    # $FLIBS was set above by AC_F77_PROG  
    # How can I work around the damaged g77 on Darwin by not linking lcrt2.o?
    # By editing flibs to remove -lcrt2.o, of course!  
    # (I realize that this is a total hack.)

    FLIBS=`echo $FLIBS | sed s/-lcrt2.o//g` 
    ;;
  *)
    echo linux
    AM_CONDITIONAL([DARWIN],[test 0 == 1])
    # Gets ready to compile with weird Fortran stuff

    #BLAS_FLAGS="-L/usr/lib64/atlas -L/usr/lib/atlas -I/usr/include/atlas"
    #BLAS_FLAGS="-L/home/cantarel/ATLAS/lib/Linux_UNKNOWNSSE2 -I/home/cantarel/ATLAS/include"
    #AC_MSG_WARN([TEMPORARILY linking manually with ATLAS BLAS in home account. FIX ME!!!!])
    #Should swap back in the commented line once BLAS issues are resolved on george. 
    #LDFLAGS="$LDFLAGS $BLAS_FLAGS"

    ACX_BLAS
    if test "$acx_blas_ok" = no; then 
	  AC_MSG_ERROR([tsnnls requires an F77BLAS. The ATLAS project
will provide these on any given system if no other is
available. Configure has searched the usual locations and can't find
one. You can use --with-blas=<lib> to point to a blas library in a
nonstandard spot. This test may fail even if you have BLAS if you
don't have a FORTRAN compiler (preferably g77) on your default path.])
    fi

    # Now look for a full (or minimal ATLAS) lapack
    AC_F77_FUNC(dpotrf)
    ACX_LAPACK(AM_CONDITIONAL([HAVE_LAPACK],[test 1 == 1]),AC_CHECK_LIB([lapack],$dpotrf,[atlas_lapack=yes 
       LAPACK_LIBS="-llapack -lcblas -lf77blas -latlas"
       AC_DEFINE_UNQUOTED([HAVE_ATLAS_LAPACK],[1],[ATLAS (instead of full) LAPACK])
       ],[atlas_lapack=no],[-lcblas -lf77blas -latlas $FLIBS -lm]))
    AC_MSG_CHECKING([for minimal ATLAS lapack])
    echo $atlas_lapack
    
  esac

# We'll need to call F77 functions from C, find out how to mangle their names
AC_F77_FUNC([daxpy],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="daxpy_"
fi
AC_DEFINE_UNQUOTED([DAXPY_F77],[$fname],[Name of daxpy after mangling])
AC_MSG_CHECKING([name of daxpy after mangling])
echo "$fname"

AC_F77_FUNC([dsyrk],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dsyrk_"
fi
AC_DEFINE_UNQUOTED([DSYRK_F77],[$fname],[Name of dsyrk after mangling])
AC_MSG_CHECKING([name of dsyrk after mangling])
echo "$fname"

AC_F77_FUNC([dsymv],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dsymv_"
fi
AC_DEFINE_UNQUOTED([DSYMV_F77],[$fname],[Name of dsymv after mangling])
AC_MSG_CHECKING([name of dsymv after mangling])
echo "$fname"

# dgemm

AC_F77_FUNC([dgemm],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dgemm_"
fi
AC_DEFINE_UNQUOTED([DGEMM_F77],[$fname],[Name of dgemm after mangling])
AC_MSG_CHECKING([name of dgemm after mangling])
echo "$fname"

# dgemv

AC_F77_FUNC([dgemv],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dgemv_"
fi
AC_DEFINE_UNQUOTED([DGEMV_F77],[$fname],[Name of dgemv after mangling])
AC_MSG_CHECKING([name of dgemv after mangling])
echo "$fname"

# dtrsm
 
AC_F77_FUNC([dtrsm],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dtrsm_"
fi
AC_DEFINE_UNQUOTED([DTRSM_F77],[$fname],[Name of dtrsm after mangling])
AC_MSG_CHECKING([name of dtrsm after mangling])
echo "$fname"

# dpotrf

AC_F77_FUNC([dpotrf],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dpotrf_"
fi
AC_DEFINE_UNQUOTED([DPOTRF_F77],[$fname],[Name of dpotrf after mangling])
AC_MSG_CHECKING([name of dpotrf after mangling])
echo "$fname"

# dgels

AC_F77_FUNC([dgels],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dgels_"
fi
AC_DEFINE_UNQUOTED([DGELS_F77],[$fname],[Name of dgels after mangling])
AC_MSG_CHECKING([name of dgels after mangling])
echo "$fname"

# dnrm2

AC_F77_FUNC([dnrm2],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dnrm2_"
fi
AC_DEFINE_UNQUOTED([DNRM2_F77],[$fname],[Name of dnrm2 after mangling])
AC_MSG_CHECKING([name of dnrm2 after mangling])
echo "$fname"

# ddot

AC_F77_FUNC([ddot],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="ddot_"
fi
AC_DEFINE_UNQUOTED([DDOT_F77],[$fname],[Name of ddot after mangling])
AC_MSG_CHECKING([name of ddot after mangling])
echo "$fname"

# dscal

AC_F77_FUNC([dscal],[fname])
if test "$fname" = "unknown"; then   # Assume g77 conventions 
  fname="dscal_"
fi
AC_DEFINE_UNQUOTED([DSCAL_F77],[$fname],[Name of dscal after mangling])
AC_MSG_CHECKING([name of dscal after mangling])
echo "$fname"

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([clapack.h vecLib/clapack.h],AC_DEFINE([HAVE_FULL_CLAPACK],[1],[Defined if we have a full clapack.]),AC_DEFINE([HAVE_FULL_CLAPACK],[],[Defined if we have full clapack.]))
AC_CHECK_HEADERS([atlas/clapack.h])
AC_CHECK_HEADERS([string.h float.h stdlib.h stdio.h sys/time.h sys/resource.h limits.h math.h time.h])

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

# Now add remaining libs
# LIBS="$LAPACK_LIBS $LIBS $BLAS_LIBS $FLIBS"

# Output
AC_CONFIG_FILES([tsnnls-libs])    # Write linker info to file. 
AC_CONFIG_FILES([Makefile
		tsnnls/Makefile
		tsnnls/taucs_basic/Makefile])
AC_OUTPUT