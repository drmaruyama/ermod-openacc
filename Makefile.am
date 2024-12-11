ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod slvfe
nobase_pkgdata_SCRIPTS = \
	tools/AMBER/gen_structure tools/AMBER/gen_input \
	tools/GROMACS/gen_structure tools/GROMACS/gen_input \
	tools/NAMD/gen_structure tools/NAMD/gen_input \
	tools/LAMMPS/gen_structure tools/LAMMPS/gen_input \
	tools/CHARMM/gen_structure tools/CHARMM/gen_input \
	tools/GENESIS/gen_structure tools/GENESIS/gen_input \
	tools/GENESIS/gen_structure_AMBER tools/GENESIS/gen_structure_CHARMM \
	tools/GENESIS/gen_structure_GROMACS \
	tools/refs_converge 

#	tools/GENESIS/gen_structure_GROMACS # <- In GENESIS, GROMACS is unsupported officially so not included in SCRIPTS

dist_pkglibexec_SCRIPTS = $(top_srcdir)/vmdplugins/libexec/*.so

ersources = enganal.F90 engmain.F90 realcal.F90 spline.F90 fft_iface.F90 mpiproc.F90 engproc.F90 recpcal.F90 insertion.F90 setconf.F90 trajectory.F90 utility.F90 bestfit.F90

EXTRA_DIST = parameters_er.example parameters_fe.example
EXTRA_DIST += extensions/*
EXTRA_DIST += trajectory-formats/*
EXTRA_DIST += tools/AMBER/*
EXTRA_DIST += tools/GROMACS/*
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/LAMMPS/*
EXTRA_DIST += tools/CHARMM/*
EXTRA_DIST += tools/GENESIS/*
EXTRA_DIST += tools/mol_dissection
EXTRA_DIST += tools/refs_converge tools/refs_converge_py2 tools/refs_converge_py3
EXTRA_DIST += vmdplugins/*.C vmdplugins/*.c vmdplugins/*.h vmdplugins/*.cxx vmdplugins/*.hxx vmdplugins/Makefile vmdplugins/LICENSE
CLEANFILES = *.mod

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = $(FFTWLIBS)
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

ermod_SOURCES = $(ersources)
ermod_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
ermod_FCFLAGS = $(VMDPLUGINFLAGS) $(FFTFLAGS)
ermod_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm 

if USE_VMDPLUGIN
ermod_SOURCES += vmdfio.c
endif

if USE_TRANSFER
ermod_FCFLAGS += -DHAVE_TRANSFER
else
ermod_SOURCES += hash_real.c
endif

if USE_FLUSH
ermod_FCFLAGS += -DHAVE_FLUSH
endif

extra_flags = 

# module dependencies
$(top_builddir)/ermod-enganal.$(OBJEXT): $(top_srcdir)/enganal.F90 \
    $(top_builddir)/ermod-engproc.$(OBJEXT) \
    $(top_builddir)/ermod-setconf.$(OBJEXT) \
    $(top_builddir)/ermod-insertion.$(OBJEXT) \
    $(top_builddir)/ermod-trajectory.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-engproc.$(OBJEXT): $(top_srcdir)/engproc.F90 \
    $(top_builddir)/ermod-insertion.$(OBJEXT) \
    $(top_builddir)/ermod-realcal.$(OBJEXT) \
    $(top_builddir)/ermod-recpcal.$(OBJEXT) \
    $(top_builddir)/ermod-bestfit.$(OBJEXT) \
    $(top_builddir)/ermod-utility.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-realcal.$(OBJEXT): $(top_srcdir)/realcal.F90 \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-recpcal.$(OBJEXT): $(top_srcdir)/recpcal.F90 \
    $(top_builddir)/ermod-spline.$(OBJEXT) \
    $(top_builddir)/ermod-fft_iface.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-insertion.$(OBJEXT): $(top_srcdir)/insertion.F90 \
    $(top_builddir)/ermod-setconf.$(OBJEXT) \
    $(top_builddir)/ermod-bestfit.$(OBJEXT) \
    $(top_builddir)/ermod-trajectory.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-setconf.$(OBJEXT): $(top_srcdir)/setconf.F90 \
    $(top_builddir)/ermod-utility.$(OBJEXT) \
    $(top_builddir)/ermod-trajectory.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-bestfit.$(OBJEXT): $(top_srcdir)/bestfit.F90 \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT)
$(top_builddir)/ermod-trajectory.$(OBJEXT): $(top_srcdir)/trajectory.F90 \
    $(top_builddir)/ermod-utility.$(OBJEXT)
$(top_builddir)/ermod-utility.$(OBJEXT): $(top_srcdir)/utility.F90 \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-mpiproc.$(OBJEXT): $(top_srcdir)/mpiproc.F90 \
    $(top_builddir)/ermod-engmain.$(OBJEXT)

# slvfe part
slvfe_SOURCES = slvfe.F90 sfemain.F90 sfecorrect.F90
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
slvfe_LINK = $(FC_CANONICAL) $(FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
slvfe.$(OBJEXT): slvfe.F90 sfemain.$(OBJEXT) sfecorrect.$(OBJEXT)
sfemain.$(OBJEXT): sfemain.F90
sfecorrect.$(OBJEXT): sfecorrect.F90 sfemain.$(OBJEXT)

# these are necessary to compile slvfe's sources without mpi
$(top_builddir)/slvfe.$(OBJEXT): $(top_srcdir)/slvfe.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfemain.$(OBJEXT): $(top_srcdir)/sfemain.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfecorrect.$(OBJEXT): $(top_srcdir)/sfecorrect.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
