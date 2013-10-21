
Description of the pxcorr inputs, outputs, and functions.

Anne Bauer
October 2013

This package includes code to calculate information about a galaxy catalog, make Partpix maps from the catalog, and calculate correlations between the maps.

The idea is to call scripts/wrapper.py.  You need to give it the name of a config file;  an example is copied below.

The wrapper figures out some metadata and then calls C++ code to make the maps and then to calculate the correlations.  Then the wrapper packages up the results into an hdf5 file.  Also, the results are also written out (currently in JSON format) in ascii.


There is a makefile in the source directory that compiles make_maps and correlate into executables and also into shared libraries that are then wrapped with python scripts so that they can be called as modules from wrapper.py.  This process is done using swig.  The src/Makefile calls the src/make/Makefile that is appropriate for your architecture.  So far there are Makefiles available that work on my mac (Makefile.Darwin) and the linux setup at PIC (Makefile.Linux).

You will probably have to change some paths in the top of the Makefile, e.g. to your python directory.

The hdf5 output file (pxcorr_out.h5) that is the product of wrapper.py can be given directly as input to Martin Eriksen's modeling & likelihood code.

Below is an example config file for doing auto and cross-correlations of two catalogs.  

sample_catalog has format
ra dec z_phot z_spec mag
while the ".fits" file is a Healpix map.  My code will understand a healpix map, but if you do correlations with a healpix map you can't give the resulting hdf5 file to martin erikson's modeling code because the metadata about the sample isn't there (which is appropriate, since that info is appropriate to a galaxy catalog rather than the stuff that will be in the maps, like systematics).

There needs to be one mask per catalog.  If the catalog is a Healpix mask then the code takes any value that doesn't equal 0 or Healpix_Undef to be a valid pixel.  

There also should be one n_of_z file, with values 
z_spec z_phot weight
per catalog.  If the catalog is a Healpix file then it doesn't use this info.  If the file is not given, as below, then it takes the z_spec and z_phot values from the input catalog.

There must be one mag cut and z bin mean & width per catalog.  If you want to do several redshift bins from one catalog, you have to put the catalog filename each time.

**** example config file ****

# catalogs
catalog_filenames: ["sample_catalog", "../../mymasks/sva1_spte_holymolys_maglims_i_cut_23.4_r9.fits"]
mask_filenames: ["../../mymasks/sva1_spte_holymolys_maglims_i_cut_23.4_r9.dat", None]
nofz_filenames: [None, None]
mag_cuts: [ 99, None ]
z_means: [ 0.6, None ]
z_widths: [ 0.2, None ]

# angular bins
ang_mean: [ 0.5, 1.0, 1.5 ]
ang_width: [ 0.5, 0.5, 0.5 ]

# config
use_counts: True
use_mags: False
only_auto: False
only_makemaps: False
pop: 'faint'
