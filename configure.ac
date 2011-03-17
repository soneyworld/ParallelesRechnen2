#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_CONFIG_MACRO_DIR([m4])
AC_PREREQ([2.65])
AC_INIT([ParallelPackage2], [0.1], [t.lorentzen@tu-bs.de])
AM_INIT_AUTOMAKE(parallelsuite2, 0.1)
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([config.h])
abi_save_CFLAGS="$CFLAGS"
abi_save_CXXFLAGS="$CXXFLAGS"

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AC_PROG_LN_S

CFLAGS="$abi_save_CFLAGS"
CXXFLAGS="$abi_save_CXXFLAGS"

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday])
AC_ARG_ENABLE(debug,
	[AC_HELP_STRING([--enable-debug], [enable extra debugging output [default=no]]) ],
	[
		CFLAGS="$CFLAGS -ggdb"
		CXXFLAGS="$CXXFLAGS -ggdb"
		AC_DEFINE([DO_DEBUG_OUTPUT], [], [Define to enable extra debugging output.])
	], [
		CFLAGS="$CFLAGS -O3 -ffast-math -mfpmath=sse+387"
		CXXFLAGS="$CXXFLAGS -O3 -ffast-math -mfpmath=sse+387"
	])
AC_ARG_ENABLE(openmp,
 	[AC_HELP_STRING([--disable-openmp], [enable parallel execution with posix threads and openmp [default=yes]]) ],
        [],[
                CFLAGS="$CFLAGS -fopenmp"
                CXXFLAGS="$CXXFLAGS -fopenmp"
                AC_DEFINE([OPENMP], [], [Define to enable openmp support.])
        ])


AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
