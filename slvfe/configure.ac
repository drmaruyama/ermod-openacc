#                                               -*- Autoconf -*-

# autoconf <= 2.64 has a bug in AC_OPENMP
AC_PREREQ([2.65])
#AC_INIT(ermod, 0.3.%rev%, Matubayasi Laboratory)
AC_INIT(ermod, 0.6.0, Matubayasi Laboratory)
AC_CONFIG_MACRO_DIR(../m4)
AC_CONFIG_SRCDIR([slvfe.F90])
AM_INIT_AUTOMAKE

AC_ARG_ENABLE([debug-symbol],
	AS_HELP_STRING([--disable-debug-symbol],
		[Disable debug symbol]))

AS_IF([test x$enable_debug_symbol = xno],
      [ac_cv_prog_cc_g=no
       ac_cv_prog_fc_g=no])

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC([icc fcc gcc cc])
AC_PROG_CC_C_O
AM_PROG_CC_C_O
AC_PROG_FC([ifort xlf90 xlf f95 f90 frt gfortran g95])

# Check for buggy variant of ifort.
AS_IF([test x$FC = xifort],
      [version=`$FC --version`
       case "x$version" in # (
       *12.1.*):
         AC_MSG_ERROR([This program does not run correctly with Intel Compiler 12.1.x])
       ;;
       esac
       ])

# modified MPI detector (searches mpiifort first)
AC_LANG(Fortran)
ACX_MPI_MY
AC_LANG(Fortran)
ACX_MPI_MY([use_mpi=yes], [use_mpi=no])

: ${enable_openmp=no}
AC_ARG_ENABLE([openmp],
              AS_HELP_STRING([--enable-openmp],
                             [use OpenMP (default = no)]))
AC_OPENMP

AS_IF([test "x$FC" = "x"],
      [AC_MSG_ERROR("Failed to find Fortran Compiler (Fortran 90). Set environment variable FC to work around")]
      [])

AC_ARG_ENABLE([vmdplugin],
	AS_HELP_STRING([--disable-vmdplugin],
		[disable VMD plugin]))

AC_DEFUN([MY_CHECK_BLAS], [
  AC_MSG_CHECKING([whether $1 is available])
  save_LIBS=$LIBS
  LIBS="$LIBS $2"
  AC_LINK_IFELSE([AC_LANG_CALL([], [sgemm])],
                 [blas_ok=yes], [blas_ok=no])
  LIBS=$save_LIBS
  AC_MSG_RESULT([$blas_ok])
])

AC_ARG_WITH([mkl],
  AS_HELP_STRING([--with-mkl[[=(MKL version)]]],
    [Automatically set environmental variables for MKL; for a specific version, use this option like --with-mkl=10.2.3]))

AS_IF([test "x$with_mkl" != "x"],
  [case x`uname -m` in # (
  xx86_64) :
  mkl_plat=em64t
  mkl_platnew=intel64
  ;; # (
  xi386) :
  mkl_plat=32
  mkl_platnew=ia32
  ;; # (
  *) :
  ;;
  esac

  AS_IF([test "x$with_mkl" = "xyes"],
    AS_IF([test "x$MKLROOT" = "x"],
      [AC_MSG_FAILURE([could not get MKL version; set MKLROOT or use --with-mkl=(version number)])],
      mklroot="$MKLROOT"),
    mklroot="/opt/intel/mkl/$with_mkl")
  # Check for recent Intel Composer (-mkl=parallel works)
  with_blas="-mkl=parallel"
  MY_CHECK_BLAS([Intel composer default], [$with_blas])
  AS_IF([test "x$blas_ok" = "xyes"],
        [FCFLAGS="-mkl=parallel $FCFLAGS"],
        [
         # Here comes super-duper crazily dirty MKL versioning checker
         # check 10.0 - 10.2
         with_blas="-L$mklroot/lib/$mkl_plat -lmkl_solver_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -lguide -lpthread"
         MY_CHECK_BLAS([MKL 10.0-10.2], 
                       [$with_blas])
         AS_IF([test "x$blas_ok" = "xyes"],
               [],
               # check for 10.3
               [with_blas="-L$mklroot/lib/$mkl_platnew -lmkl_rt $OPENMP_FCFLAGS -lpthread"
                MY_CHECK_BLAS([MKL 10.3+],
                              [$with_blas])
                AS_IF([test "x$blas_ok" = "xyes"],
                      [],
                      [AC_MSG_WARN([Failed to find appropriate MKL. If you are using old mkl, please update])])])
         with_lapack="-L$mklroot/lib/$mkl_plat -lmkl_lapack"
         FCFLAGS="-I$mklroot/include $FCFLAGS"
        ])
])

ACX_BLAS
ACX_LAPACK([lapack_exist=yes],
           [lapack_exist=no])

################################

AC_DEFUN([MY_FC_REAL_DOUBLE],[
AC_CACHE_CHECK([how to read real as real(8) as default],
my_fc_real_double_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
# for default, gfortran, ifort, xlf, Fujitsu f90, Fujitsu f90 (old), Portland, Cray
for flag in "" "-fdefault-real-8 -fdefault-double-8" "-fdefault-real-8" "-autodouble" "-qrealsize=8" "-CcdRR8" "-Ad" "-Mr8" "-s default64"; do
FCFLAGS="$save_FCFLAGS $flag"
# This is very tricky indeed...
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([], [dnl
      implicit none
      real :: a
      integer, parameter :: len = kind(a)
      integer, parameter :: x(len) = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
      print *, "Success"
])], [my_check_realasdouble=yes], [my_check_realasdouble=no])
  if test $my_check_realasdouble = yes; then
    break
  fi
done
if test $my_check_realasdouble = no; then
  AC_MSG_ERROR([unable to find compiler flag to read real as real(8)])
fi
FCFLAGS="$save_FCFLAGS $flag"
if test "x$flag" = x; then
  echo -n "(none)"
else
  echo -n $flag
fi
AC_LANG_POP(Fortran)
])])

MY_FC_REAL_DOUBLE

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

AC_ARG_ENABLE([perf],
	AS_HELP_STRING([--enable-perf],
		[(internal use only)]),
	AS_IF([test x$enable_perf = xyes],
		     [FCFLAGS="$FCFLAGS -fno-omit-frame-pointer -DPERF"
                      CFLAGS="$CFLAGS -fno-omit-frame-pointer"
                      LDFLAGS="$LDFLAGS -lprofiler"]))

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

AC_ARG_ENABLE([opt],
	AS_HELP_STRING([--disable-opt],
		[disable automatic edit to optimization flags]))

AS_IF([test x$enable_opt = xyes || test x$enable_opt = x],
  [FCFLAGS="$FCFLAGS"
   AS_IF([test x$ac_compiler_gnu = xyes],
        [optflag="-O3 -ffast-math -fno-finite-math-only"
         MY_FC_OPT_ACCEPT_NATIVE
         AS_IF([test x$my_fc_opt_accept_native = xyes],
                optflag="$optflag -march=native")],
        [case x$FC in #(
           xmpiifort|xifort)
             optflag="-O3 -no-prec-div -xHost"
             ;; #(
           *)
             optflag="-O"
             ;;
         esac])
   FCFLAGS="$FCFLAGS $optflag"])

# Check OpenMP
AS_IF([test "$enable_openmp" != "no"],
  [FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
   # Check OpenMP version
   AC_MSG_CHECKING([whether OpenMP version >= 3.0])
   AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM([], [dnl
!\$omp task
  print *, "test"
!\$omp end task
])], 
     [AC_MSG_RESULT([yes])], 
     [AC_MSG_RESULT([no])
      AC_MSG_ERROR([OpenMP must be >= 3.0])])])

# AM_CONDITIONAL([MPI], [test x$use_mpi = xyes])
AS_IF([test "x$lapack_exist" != xyes],
	    [AC_MSG_ERROR(
		[dnl
LAPACK library was not found.
Try specifying --with-lapack=<lib> options or set FCFLAGS/FLIBS environment variables to fix the problem.])])

lapack_result="$LAPACK_LIBS $BLAS_LIBS"

AC_CONFIG_FILES([Makefile])
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

AC_MSG_RESULT([dnl
LAPACK: $lapack_result
MPI: $mpi_result
OpenMP: $enable_openmp
FC(for ermod)=$FC
FC_CANONICAL(for slvfe)=$FC_CANONICAL
FCFLAGS=$FCFLAGS

Configuration complete; type \"make\" to start compiling])

