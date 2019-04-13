#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([heatFD], [0.2], [stpotter@utexas.edu])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# Check for pdflatex
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
if test -z "$PDFLATEX"; then
	AC_MSG_WARN([Unable to build PDF version of model document])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")

# Checks for libraries.
AX_PATH_GRVY([0.32],[],AC_MSG_ERROR([Unable to locate GRVY or meet minimum version requirement]))
AX_PATH_MASA([0.50],[],AC_MSG_ERROR([Unable to locate MASA or meet minimum version requirement]))
AX_PATH_HDF5([1.8.0],[yes],AC_MSG_ERROR([Unable to locate HD5F or meet minimum version requirement]))
AX_PATH_PETSC([no])

# check for code coverage
AX_CODE_COVERAGE

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])

AC_OUTPUT(Makefile src/Makefile test/Makefile docs/Makefile)
