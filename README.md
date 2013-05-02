
Description of the pxcorr inputs, outputs, and functions.

Anne Bauer
May 2013

This package includes code to calculate information about a galaxy catalog, make Partpix maps from the catalog, and calculate correlations between the maps.

The idea is to call scripts/wrapper.py.  Right now there are some configuration details in the top part of wrapper.py such as the desired angular binning and the name of the input catalog (there are 2 examples in the example_input directory).  The wrapper calls some python scripts to figure out some metadata and then calls C++ code to make the maps and then to calculate the correlations.  Then a python script packages up the results into an hdf5 file.  Also, the results are also written out (currently in JSON format) in ascii.

The C++ functions can be called as executables, but that is no longer the way they are meant to be called.  For testing, however, you can call them like

make_maps metadata.txt catalog1 None counts mags

to make a map from catalog1 with no given mask (meaning that the code will make some approximate mask itself).  actually this will make 2 maps: one for delta_counts and one for delta_mags, as indicated in the last two arguments.  they'll be called dc_map.h5 and dm_map.h5.  metadata.txt and catalog0 and catalog1 are examples given in the example_input directory, for testing.

correlate map1.h5 map2.h5

will correlate the two given maps.  The output will be in JSON format in ascii files, one file for the correlation and another for the covariance.  I plan to change this in the future to output the jackknife results explicitly rather than the pre-computed covariance matrix, since as it is now it is impossible to collect a lot of different cross-correlations and compile a complete covariance.

There is a makefile in the source directory that compiles make_maps and correlate into executables and also into shared libraries that are then wrapped with python scripts so that they can be called as modules from wrapper.py.  This process is done using swig.  You will probably have to change some paths in the top of the Makefile, e.g. to your python directory.

The hdf5 output file (pxcorr_out.h5) that is the product of wrapper.py can be given directly as input to Martin Eriksen's modeling & likelihood code.

