ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = erdst
nobase_pkgdata_SCRIPTS = \
	../tools/AMBER/gen_structure ../tools/AMBER/gen_input \
	../tools/GROMACS/gen_structure ../tools/GROMACS/gen_input \
	../tools/NAMD/gen_structure ../tools/NAMD/gen_input \
	../tools/LAMMPS/gen_structure ../tools/LAMMPS/gen_input \
	../tools/CHARMM/gen_structure ../tools/CHARMM/gen_input \
	../tools/GENESIS/gen_structure ../tools/GENESIS/gen_input \
	../tools/GENESIS/gen_structure_AMBER \
	../tools/GENESIS/gen_structure_CHARMM \
	../tools/GENESIS/gen_structure_GROMACS \
	../tools/refs_converge 

#	tools/GENESIS/gen_structure_GROMACS # <- In GENESIS, GROMACS is unsupported officially so not included in SCRIPTS

dist_pkglibexec_SCRIPTS = $(top_srcdir)/../vmdplugins/libexec/*.so

ersources = enganal.F90 engmain.F90 realcal.F90 spline.F90 fft_iface.F90 mpiproc.F90 engproc.F90 recpcal.F90 insertion.F90 setconf.F90 trajectory.F90 utility.F90 bestfit.F90

EXTRA_DIST = ../parameters_er.example ../parameters_fe.example
EXTRA_DIST += ../extensions/*
EXTRA_DIST += ../trajectory-formats/*
EXTRA_DIST += ../tools/AMBER/*
EXTRA_DIST += ../tools/GROMACS/*
EXTRA_DIST += ../tools/NAMD/*
EXTRA_DIST += ../tools/LAMMPS/*
EXTRA_DIST += ../tools/CHARMM/*
EXTRA_DIST += ../tools/GENESIS/*
EXTRA_DIST += ../tools/mol_dissection
EXTRA_DIST += ../tools/refs_converge ../tools/refs_converge_py2 ../tools/refs_converge_py3
EXTRA_DIST += ../vmdplugins/*.C ../vmdplugins/*.c ../vmdplugins/*.h ../vmdplugins/*.cxx ../vmdplugins/*.hxx ../vmdplugins/Makefile ../vmdplugins/LICENSE
CLEANFILES = *.mod

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = $(FFTWLIBS)
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

erdst_SOURCES = $(ersources)
erdst_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
erdst_FCFLAGS = $(VMDPLUGINFLAGS) $(FFTFLAGS)
erdst_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm

if USE_VMDPLUGIN
erdst_SOURCES += vmdfio.c
endif

if USE_TRANSFER
erdst_FCFLAGS += -DHAVE_TRANSFER
else
erdst_SOURCES += hash_real.c
endif

if USE_FLUSH
erdst_FCFLAGS += -DHAVE_FLUSH
endif

extra_flags = 

# module dependencies
$(top_builddir)/erdst-enganal.$(OBJEXT): $(top_srcdir)/enganal.F90 \
    $(top_builddir)/erdst-engproc.$(OBJEXT) \
    $(top_builddir)/erdst-setconf.$(OBJEXT) \
    $(top_builddir)/erdst-insertion.$(OBJEXT) \
    $(top_builddir)/erdst-trajectory.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-engproc.$(OBJEXT): $(top_srcdir)/engproc.F90 \
    $(top_builddir)/erdst-insertion.$(OBJEXT) \
    $(top_builddir)/erdst-realcal.$(OBJEXT) \
    $(top_builddir)/erdst-recpcal.$(OBJEXT) \
    $(top_builddir)/erdst-bestfit.$(OBJEXT) \
    $(top_builddir)/erdst-utility.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-realcal.$(OBJEXT): $(top_srcdir)/realcal.F90 \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-recpcal.$(OBJEXT): $(top_srcdir)/recpcal.F90 \
    $(top_builddir)/erdst-spline.$(OBJEXT) \
    $(top_builddir)/erdst-fft_iface.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-insertion.$(OBJEXT): $(top_srcdir)/insertion.F90 \
    $(top_builddir)/erdst-setconf.$(OBJEXT) \
    $(top_builddir)/erdst-bestfit.$(OBJEXT) \
    $(top_builddir)/erdst-trajectory.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-setconf.$(OBJEXT): $(top_srcdir)/setconf.F90 \
    $(top_builddir)/erdst-utility.$(OBJEXT) \
    $(top_builddir)/erdst-trajectory.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-bestfit.$(OBJEXT): $(top_srcdir)/bestfit.F90 \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT)
$(top_builddir)/erdst-trajectory.$(OBJEXT): $(top_srcdir)/trajectory.F90 \
    $(top_builddir)/erdst-utility.$(OBJEXT)
$(top_builddir)/erdst-utility.$(OBJEXT): $(top_srcdir)/utility.F90 \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-mpiproc.$(OBJEXT): $(top_srcdir)/mpiproc.F90 \
    $(top_builddir)/erdst-engmain.$(OBJEXT)

