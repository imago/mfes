dnl Process this file with autoconf to produce a configure script.



AC_PREREQ(2.59)
AC_INIT(mfes, 0.1)

AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()


AC_ARG_WITH([netgen-src],
	[AC_HELP_STRING([--with-netgen-src=dir][Point to netgen source dir])],
    [AC_SUBST([NETGEN_SOURCE], [$withval])],
	AC_MSG_ERROR([You must call configure with the --with-netgen-src option. This tells configure where to find the netgen 5.1 source files e.g. --with-netgen-src=/opt/netgen/sources/netgen-5.1/])
	)
        

MPICXX="mpicxx"
AC_PROG_CXX($MPICXX)
# AC_PROG_CC($MPICC)
# AC_PROG_CXX

AC_PROG_LIBTOOL

AC_CONFIG_FILES(Makefile src/Makefile)


AC_OUTPUT
