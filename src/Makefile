HEALPIXPATH =  ../Healpix_Partpix/src/cxx/osx
HDFDIR     =  /usr/local/hdf5-1.8.11/hdf5
HDF        =  /usr/local/hdf5-1.8.11/hdf5/lib/libhdf5.a /usr/local/hdf5-1.8.11/hdf5/lib/libhdf5_hl.a /usr/local/hdf5-1.8.11/hdf5/lib/libhdf5_cpp.a /usr/local/hdf5-1.8.11/hdf5/lib/libhdf5_hl_cpp.a
OPENCLDIR   = /System/Library/Frameworks/OpenCL.framework/Versions/A
#PYTHONDIR   = /Library/Frameworks/Python.framework/Versions/2.7/include/python2.7

CC         = g++
SWIG       = /usr/local/bin/swig
CFLAGS     = -Wall -MD -O3 -fopenmp -I$(HEALPIXPATH)/include/ -I/opt/local/include/ -I/usr/local -I$(OPENCLDIR)/Headers/ -I$(HDFDIR)/include 
LFLAGS     = -lgomp -L/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/x86_64 -L$(HEALPIXPATH)/lib -L/opt/local/lib -lcfitsio -lcxxsupport -lhealpix_cxx -lpsht -lc_utils -L$(HDFDIR)/lib -L$(OPENCLDIR)/Libraries/ -framework OpenCL -lz 

all:    MK_LIBS

MK_LIBS:
	$(SWIG) -c++ -python make_maps.i 
	$(CC) -fpic -c make_maps.C make_maps_wrap.cxx $(CFLAGS) -I/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/ -I/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/
	ld -bundle -flat_namespace -undefined suppress -lm -lpython2.7 -L/Library/Frameworks/Python.framework/Versions/2.7/lib $(LFLAGS) $(HDF) -o _make_maps.so make_maps.o make_maps_wrap.o
	$(CC) make_maps.o $(HDF) -o ../bin/make_maps $(LFLAGS) -lm
	$(SWIG) -c++ -python correlate.i 
	$(CC) -fpic -c correlate.C correlate_wrap.cxx $(CFLAGS) -I/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/ -I/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/
	ld -bundle -flat_namespace -undefined suppress -lm -lpython2.7 -L/Library/Frameworks/Python.framework/Versions/2.7/lib $(LFLAGS) $(HDF) -o _correlate.so correlate.o correlate_wrap.o
	$(CC) correlate.o $(HDF) -o ../bin/correlate $(LFLAGS) -lm
	$(SWIG) -c++ -python make_maps3pt.i 
	$(CC) -fpic -c make_maps3pt.C make_maps3pt_wrap.cxx $(CFLAGS) -I/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/ -I/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/
	ld -bundle -flat_namespace -undefined suppress -lm -lpython2.7 -L/Library/Frameworks/Python.framework/Versions/2.7/lib $(LFLAGS) $(HDF) -o _make_maps3pt.so make_maps3pt.o make_maps3pt_wrap.o
	$(CC) make_maps3pt.o $(HDF) -o ../bin/make_maps3pt $(LFLAGS) -lm
	$(SWIG) -c++ -python correlate3pt.i 
	$(CC) -fpic -c correlate3pt.C correlate3pt_wrap.cxx $(CFLAGS) -I/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/ -I/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/
	ld -bundle -flat_namespace -undefined suppress -lm -lpython2.7 -L/Library/Frameworks/Python.framework/Versions/2.7/lib $(LFLAGS) $(HDF) -o _correlate3pt.so correlate3pt.o correlate3pt_wrap.o
	$(CC) correlate3pt.o $(HDF) -o ../bin/correlate3pt $(LFLAGS) -lm

