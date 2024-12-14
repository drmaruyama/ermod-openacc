ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = slvfe

CLEANFILES = *.mod

slvfe_SOURCES = slvfe.F90 sfemain.F90 sfecorrect.F90
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
slvfe_LINK = $(FC) $(FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
slvfe.$(OBJEXT): slvfe.F90 sfemain.$(OBJEXT) sfecorrect.$(OBJEXT)
sfemain.$(OBJEXT): sfemain.F90
sfecorrect.$(OBJEXT): sfecorrect.F90 sfemain.$(OBJEXT)

# these are necessary to compile slvfe's sources without mpi
$(top_builddir)/slvfe.$(OBJEXT): $(top_srcdir)/slvfe.F90
	$(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfemain.$(OBJEXT): $(top_srcdir)/sfemain.F90
	$(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfecorrect.$(OBJEXT): $(top_srcdir)/sfecorrect.F90
	$(FC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
