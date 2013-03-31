#!/usr/bin/env python
# encoding: UTF8

import sys
import json
import os
import pdb
import numpy as np
import tables
import ctypes

from make_meta import make_metadata
from parse_sample import parse_data
from packobs import add_corr, add_cov

sys.path.append("../src/")
#import makemaps

#make_maps.hello_world

catalog_filename = "/Users/bauer/correlations/pic_test/catalogs/radeczpzsr_wmu_small"
mag_cuts = [ 21.2, 23.2, 24.1, 24.1 ]
use_counts = False
use_mags = False
pop = 'faint'
z_mean = [0.3, 0.6, 0.9, 1.2]
z_width = [0.2, 0.2, 0.2, 0.2]
ang_mean = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
ang_width = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]


# make some metadata to tell us what to do
hdf5file = tables.openFile('pxcorr_out.hdf5', 'w')
ang_info = make_metadata(hdf5file, z_mean, z_width, ang_mean, ang_width, pop)
print ang_info

# parse the catalog: construct N(z) data, calculate slope, etc.
subcat_filenames = parse_data( catalog_filename, mag_cuts, hdf5file, sparse=True )
print subcat_filenames

exit(1)

ang_means_carray = []
ang_widths_carray = []

for filename in subcat_filenames:
    pass
    # makemaps.make_maps( filename, 'None', ang_means_carray, ang_widths_carray, len(ang_info[0]), use_counts, use_mags );

    # correlate
    
# packobs
add_corr(hdf5file)
add_cov(hdf5file)

hdf5file.close()

# the end.
