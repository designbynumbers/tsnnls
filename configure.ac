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

# Checks for libraries.
#AC_CHECK_LIB(m,sqrt)
#AC_CHECK_LIB(argtable2,arg_parse,[have_argtable2=true],[have_argtable2=false])

# Checks for Darwin.
AC_MSG_CHECKING(whether we should have Accelerate framework)
AC_CANONICAL_HOST
# Since foo is currently untested on any os other
# than solaris, so check the os and quit if not solaris.
UNSUPPORTED_OS="Will try to set up tsnnls with local BLAS/LAPACK"
case "${host_os}" in
  darwin*)
    echo yes
    echo
    echo ===========================================================
    echo Setting up build environment for ${host_cpu}${host_os}
    echo ===========================================================
    AM_CONDITIONAL([DARWIN],[test 1 == 1])
    [HAVE_DARWIN=true]
    ;;
  *)
    echo no
    AM_CONDITIONAL([DARWIN],[test 0 == 1])
    AC_MSG_WARN($UNSUPPORTED_OS)
  esac

#AM_CONDITIONAL([DARWIN],[test "$(uname)" = Darwin],AM_WARN["Found Darwin, will link with Accelerate framework"],AM_WARN["Will attempt to link with system BLAS/LAPACK"])

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Output
AC_CONFIG_FILES([Makefile
		tsnnls/Makefile
		tsnnls/taucs_basic/Makefile])
AC_OUTPUT