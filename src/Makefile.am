## process this file with automake to produce Makefile.in
bin_PROGRAMS = trisurf tsmeasure tspoststat
trisurf_SOURCES = main.c
GITVERSION:=$(shell git --no-pager describe --tags --always --dirty)
AM_CFLAGS = -Wall -Wextra -DTS_VERSION=\"$(GITVERSION)\" -Wno-error=date-time # -Werror
AM_CPPFLAGS = ${libxml2_CFLAGS} -fgnu89-inline
trisurf_LDADD = libtrisurf.la ${libcurl_LIBS} ${libxml2_LIBS}

tsmeasure_SOURCES = tsmeasure.c
tsmeasure_LDADD = libtrisurf.la ${libcurl_LIBS} ${libxml2_LIBS}

tspoststat_SOURCES = tspoststat.c
tspoststat_LDADD = libtrisurf.la ${libcurl_LIBS} ${libxml2_LIBS}

#------------- LIBS ----------

lib_LTLIBRARIES= libtrisurf.la
libtrisurf_la_SOURCES= general.c vertex.c bond.c triangle.c cell.c vesicle.c initial_distribution.c io.c dumpstate.c frame.c energy.c timestep.c vertexmove.c bondflip.c poly.c stats.c sh.c shcomplex.c constvol.c snapshot.c restore.c cluster.c
#libtrisurf_ey_la_CPPFLAGS = ${libxml2_CFLAGS}
libtrisurf_la_LIBADD=${libxml2_LIBS}
pkginclude_HEADERS=general.h vertex.h bond.h triangle.h cell.h vesicle.h initial_distribution.h io.h dumpstate.h frame.h energy.h timestep.h vertexmove.h bondflip.h poly.h stats.h sh.h shcomplex.h constvol.h snapshot.h restore.h cluster.h
