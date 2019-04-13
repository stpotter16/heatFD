AUTOMAKE_OPTIONS = foreign

if PETSC_ENABLED
	include $(PETSC_DIR)/lib/petsc/conf/variables
	AM_CXXFLAGS += -DINCLUDE_PETSC $(PETSC_CC_INCLUDES)
endif

bin_PROGRAMS = heatFD

heatFD_SOURCES = main.c utilities.c utilities.h

AM_CFLAGS = $(GRVY_CFLAGS)
AM_CFLAGS += $(MASA_CXXFLAGS)
AM_CFLAGS += $(HDF5_CFLAGS)
LIBS = $(GRVY_LIBS)
LIBS += $(MASA_LIBS)
LIBS += $(HDF5_LIBS)
LIBS += -lm

if PETSC_ENABLED
	heatFD_LDADD += ${PETSC_LIB}
endif
