BINDIR = ../bin
OBJDIR = ../obj
vpath %.o $(OBJDIR)
VPATH = $(BINDIR)

HEALPIXPATH =  /usr/local/Healpix_Partpix/src/cxx/osx
HDFDIR     =  /usr/local/hdf5-1.8.9/hdf5
HDF        =  /usr/local/hdf5-1.8.9/hdf5/lib/libhdf5.a /usr/local/hdf5-1.8.9/hdf5/lib/libhdf5_hl.a /usr/local/hdf5-1.8.9/hdf5/lib/libhdf5_cpp.a /usr/local/hdf5-1.8.9/hdf5/lib/libhdf5_hl_cpp.a
OPENCLDIR   = /System/Library/Frameworks/OpenCL.framework/Versions/A

CC         = g++
CFLAGS     = -Wall -MD -O3 -I$(HEALPIXPATH)/include/ -I/opt/local/include/ -I/usr/local -I$(OPENCLDIR)/Headers/ -I$(HDFDIR)/include
LFLAGS     = -L$(HEALPIXPATH)/lib -L/opt/local/lib -lcfitsio -lcxxsupport -lhealpix_cxx -lpsht -lc_utils -L$(HDFDIR)/lib -fopenmp -L$(OPENCLDIR)/Libraries/ -framework OpenCL -lz


DEP = ./Dependency_files

BINARIES = correlate make_maps

all:    $(DEP) \
        $(BINARIES) \
        TIDY

$(BINARIES):
	$(CC) $(CFLAGS) $^  -o $(BINDIR)/$@ $(LFLAGS) -lm

# Disable making of c files
%.o : %.c

# Disable .cc implicit rule
%.o : %.cc

# Produce object files from .C files only.
%.o : %.C
	$(CC) $(MPFLAG) $(CFLAGS) -c $< -o $@ 
	@  mv *.d $(DEP)

make_maps: make_maps.o $(HDF)
correlate: correlate.o $(HDF)

$(DEP):
	mkdir $@
	touch $@/default_dependency.d

.PHONY: TIDY
TIDY:
	-@      mv *.o $(OBJDIR)
#       -@      mv *.d $(DEP)

-include $(DEP)/*.d *.d
