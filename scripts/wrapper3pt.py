#!/usr/bin/env python
# encoding: UTF8

import sys
import json
import os
import pdb
import numpy as np
import tables
import ctypes

from make_meta import make_metadata3pt
from parse_sample import parse_data3pt
from packobs import add_corr, add_cov

sys.path.append("/Users/bauer/software/pxcorr/src/")
import make_maps3pt
import correlate3pt

catalog_filename = "/Users/bauer/surveys/SDSS/megaz/pol/dr7/radecz1z1i"
mask_filename = "/Users/bauer/surveys/SDSS/megaz/dr7/megaz_dr7.1024.mask.dat"

#catalog_filename = "../example_input/catalog0"
#mag_cuts = [ 20.3, 22.4, 23.8, 24.2 ]
mag_cuts = [ 19.8 ] #, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0 ]
#mag_cuts = [ 23, 24 ]
use_counts = True
use_mags = False
pop = 'faint'
z_mean = [ 0.525 ] #, 0.425, 0.575, 0.725, 0.875, 1.025, 1.175, 1.325 ] 
z_width = [ 0.05 ] #, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15 ]

angles = [ 10, 30, 50, 70, 90, 110, 130, 150, 170 ]
r_12s = [ 3.0, 2.0, 1.0, 0.5 ]
r_23s = [ 1.5, 1.0, 0.5, 0.25 ]
only_auto = True

# make some metadata to tell us what to do
hdf5file = tables.openFile('pxcorr_out.h5', 'w')
md_info = make_metadata3pt( hdf5file, z_mean, z_width, angles, r_12s, r_23s, pop)
print md_info

# parse the catalog: construct N(z) data, calculate slope, etc.
subcat_filenames = parse_data3pt( catalog_filename, mag_cuts, hdf5file, sparse=True )
# print subcat_filenames


for fi, filename in enumerate(subcat_filenames):
    suffix = "_" + str(fi)
    make_maps3pt.make_maps3pt( filename, mask_filename, angles, r_12s, r_23s, suffix );
    print "Finished making map from %s" %filename
    

nf = len( subcat_filenames )

# only do one map for now
for i in range(nf):
    name1 = "dc_map3pt_" + str(i) + ".h5"
    print "Correlating %s" %(name1)
    suffix = "." + str(i) + "c"
    correlate3pt.correlate3pt( name1, name1, name1, suffix )
        
# skip packobs for now.
exit(0)



# packobs

print "Adding correlations"

# add correlations to the file
hdf5file.createGroup('/', 'corr')

corr = np.ones((len(ang_mean), len(z_mean), len(z_mean)))
for i in range(nf):
    for j in range(i,nf):
        file_path = "corr." + str(i) + "c." + str(j) + "c"
        datafile = open(file_path)
        data = json.load(datafile)
        datafile.close()
        corr[:,i,j] = data[0]
        if i != j:
            corr[:,j,i] = data[0]
corrobj = hdf5file.createArray('/corr', 'corr1', corr)
corrobj.setAttr('ftype0', json.dumps('counts'))
corrobj.setAttr('pop0', json.dumps('faint'))
corrobj.setAttr('ftype1', json.dumps('counts'))
corrobj.setAttr('pop1', json.dumps('faint'))

# add covariance to the file
cov = np.ones((len(ang_mean), len(ang_mean), len(z_mean), len(z_mean), len(z_mean), len(z_mean)))

# TBD: read in the jackknife results and make the big covariance matrix.
# need to change the correlate output...

print "Adding covariances"

hdf5file.createGroup('/', 'cov')

for i in range(nf):
    for j in range(i,nf):
        
        file_path = "cov." + str(i) + "c." + str(j) + "c"
        datafile = open(file_path)
        data = json.load(datafile)
        datafile.close()

        cov[:,:,i,j,i,j] = data
        cov[:,:,j,i,i,j] = data
        cov[:,:,i,j,j,i] = data
        cov[:,:,j,i,j,i] = data

corrobj = hdf5file.createArray('/cov', 'cov1', cov)
corrobj.setAttr('ftype0', json.dumps('counts'))
corrobj.setAttr('pop0', json.dumps('faint'))
corrobj.setAttr('ftype1', json.dumps('counts'))
corrobj.setAttr('pop1', json.dumps('faint'))
corrobj.setAttr('ftype2', json.dumps('counts'))
corrobj.setAttr('pop2', json.dumps('faint'))
corrobj.setAttr('ftype3', json.dumps('counts'))
corrobj.setAttr('pop3', json.dumps('faint'))

hdf5file.close()

# the end.
