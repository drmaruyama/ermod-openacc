#                                               -*- Autoconf -*-

# autoconf <= 2.64 has a bug in AC_OPENMP
AC_PREREQ([2.65])
AC_INIT(erdst-openacc, 1.0.0, Matubayasi Laboratory)
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE

AC_ARG_ENABLE([debug-symbol],
              [AS_HELP_STRING([--enable-debug-symbol],
                        [Enable debug symbols (-g)])],
              [enable_debug=yes], [enable_debug=no])

AS_IF([test x$enable_debug = xno],
      [ac_cv_prog_cc_g=no
       ac_cv_prog_fc_g=no])

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC([nvc gcc cc])
AC_PROG_CC_C_O
AM_PROG_CC_C_O
AC_PROG_FC([nvfortran])
AC_PROG_LN_S

# modified MPI detector (searches mpiifort first)
AC_LANG(Fortran)
ACX_MPI_MY
AC_LANG(Fortran)
ACX_MPI_MY([use_mpi=yes], [use_mpi=no])

AS_IF([test "x$FC" = "x"],
      [AC_MSG_ERROR("Failed to find Fortran Compiler (Fortran 90). Set environment variable FC to work around")]
      [])

AC_ARG_ENABLE([vmdplugin],
	AS_HELP_STRING([--disable-vmdplugin],
		[disable VMD plugin]))

################################

AC_ARG_ENABLE([double],
    [AS_HELP_STRING([--enable-double],
        [Enable conversion of single precision real to double precision (real(8))])],
    [enable_double=yes],
    [enable_double=no])

# Precision (single vs. double) is selected entirely by the -DDP
# preprocessor definition, consumed by precision.F90's "wp" kind
# parameter. Unlike before, no compiler flag that redefines what a
# bare "real" means (e.g. nvfortran's -Mr8) is needed anymore, since
# every real/complex declaration in the source spells out "real(wp)"
# explicitly.
AS_IF([test "x$enable_double" = "xyes"],
      [FCFLAGS="$FCFLAGS -DDP"])

################################

AC_LANG_PUSH(Fortran)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [dnl
      implicit none
      integer(8) :: a = 1
      real(8) :: b = 1.0d0
      
      a = transfer(b, a)
      print *, a
])], [my_check_have_transfer=yes], [my_check_have_transfer=no])

echo -n $my_check_have_transfer
AC_LANG_POP(Fortran)

AM_CONDITIONAL([USE_TRANSFER], [test "x$my_check_have_transfer" = "xyes"])

################################

AC_LANG_PUSH(Fortran)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [dnl
      implicit none

      open(10,file='/tmp/configure.ermod.flush.test')
      write(10,*) 342
      flush(10)
      close(10)
])], [my_check_have_flush=yes], [my_check_have_flush=no])

echo -n $my_check_have_flush
AC_LANG_POP(Fortran)

AM_CONDITIONAL([USE_FLUSH], [test "x$my_check_have_flush" = "xyes"])

################################

# set to f90 as default
AC_FC_SRCEXT(f90)

AC_ARG_ENABLE([mpi],
	AS_HELP_STRING([--disable-mpi],
		[disable MPI even if exists]),
	AS_IF([test x$enable_mpi = xno],
		     [use_mpi=no]),
	[])


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# FIXME: check cpp flags

AC_ARG_ENABLE([defaults],
	AS_HELP_STRING([--disable-defaults],
		[(deprcated)]))

AS_IF([test x$enable_defaults != xno],
CFLAGS="$CFLAGS"
FCFLAGS="$FCFLAGS")

AM_CONDITIONAL(USE_MPI, [test "x$use_mpi" = "xyes"])

AS_IF([test "x$use_mpi" = "xyes"],
	    [mpi_result="enabled"
	     FC_CANONICAL="$FC"
	     FC="$MPIFC"
             FCFLAGS="-DMPI $FCFLAGS"],
	    [mpi_result="disabled"
	     FC_CANONICAL="$FC"])
AC_SUBST([FC_CANONICAL])

# optimization flags
FCFLAGS="$FCFLAGS -fast -Minfo=accel -acc -cudalib=cufft"
AC_DEFUN([MY_FC_OPT_ACCEPT_NATIVE],[
AC_CACHE_CHECK([fortran compiler accepts -march=native],
my_fc_opt_accept_native_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
flag="-march=native"
FCFLAGS="$save_FCFLAGS $flag"
AC_COMPILE_IFELSE([      program main
      end program], [my_fc_opt_accept_native=yes], [my_fc_opt_accept_native=no])
echo -n $my_fc_opt_accept_native
FCFLAGS=$save_FCFLAGS
AC_LANG_POP(Fortran)
])])

# AM_CONDITIONAL([MPI], [test x$use_mpi = xyes])
AM_CONDITIONAL([USE_VMDPLUGIN], [test x$enable_vmdplugin = xyes || test x$enable_vmdplugin = x])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

AC_MSG_RESULT([dnl
MPI: $mpi_result
FC(for ermod)=$FC
FCFLAGS=$FCFLAGS
LDFLAGS=$LDFLAGS
LIBS=$LIBS

Configuration complete; type \"make\" to start compiling])
