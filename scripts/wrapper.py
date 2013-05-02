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

sys.path.append("/Users/bauer/software/pxcorr/src/")
import make_maps
import correlate

#catalog_filename = "/Users/bauer/mice/des_mice/mice-des-v0.4-r1.4/v0.4_HOD_bcnz_with_id/radeczpzsr_wmu"
catalog_filename = "../example_input/catalog0"
#mag_cuts = [ 20.3, 22.4, 23.8, 24.2 ]
mag_cuts = [ 24.0 ]#, 24.0, 24.0, 24.0 ]
#mag_cuts = [ 22.7, 22.7, 22.7 ]
use_counts = True
use_mags = False
pop = 'faint'
z_mean = [0.7] #, 0.6, 0.9, 1.2]
z_width = [0.2] #, 0.2, 0.2, 0.2]
ang_mean = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8] #, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
ang_width = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1] #, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
#z_mean = [0.205, 0.605, 1.005]
#z_width = [0.4, 0.4, 0.4]
#ang_mean = [0.165, 0.225, 0.285] #[0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.095, 0.115, 0.165, 0.225, 0.285]
#ang_width = [0.06, 0.06, 0.06] #[0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0.01, 0.02, 0.02, 0.04, 0.04, 0.06, 0.06, 0.06]

# make some metadata to tell us what to do
hdf5file = tables.openFile('pxcorr_out.h5', 'w')
ang_info = make_metadata( hdf5file, z_mean, z_width, ang_mean, ang_width, pop)
print ang_info

# parse the catalog: construct N(z) data, calculate slope, etc.
subcat_filenames = parse_data( catalog_filename, mag_cuts, hdf5file, sparse=True )
# print subcat_filenames


for fi, filename in enumerate(subcat_filenames):
    suffix = "_" + str(fi)
    make_maps.make_maps( filename, 'None', ang_mean, ang_width, use_counts, use_mags, suffix );
    print "Finished making map from %s" %filename

nf = len( subcat_filenames )

for i in range(nf):
    for j in range(i, nf):
        name1 = "dc_map_" + str(i) + ".h5"
        name2 = "dc_map_" + str(j) + ".h5"
        print "Correlating %s and %s" %(name1, name2)
        suffix = "." + str(i) + "c." + str(j) + "c"
        correlate.correlate( name1, name2, suffix )
        
        if use_mags:
            name2 = "dm_map_" + str(j) + ".h5"
            suffix = "." + str(i) + "c." + str(j) + "m"
            print "Correlating %s and %s" %(name1, name2)
            correlate.correlate( name1, name2, suffix )
            name1 = "dm_map_" + str(i) + ".h5"
            suffix = "." + str(i) + "m." + str(j) + "m"
            print "Correlating %s and %s" %(name1, name2)
            correlate.correlate( name1, name2, suffix )
            name2 = "dc_map_" + str(j) + ".h5"
            suffix = "." + str(i) + "m." + str(j) + "c"
            print "Correlating %s and %s" %(name1, name2)
            correlate.correlate( name1, name2, suffix )

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
corrobj.setAttr('ftype0', json.dumps(['counts']))
corrobj.setAttr('ftype1', json.dumps(['counts']))

if( use_mags ):
    for i in range(nf):
        for j in range(i,nf):
            file_path = "corr." + str(i) + "c." + str(j) + "m"
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()
            corr[:,i,j] = data[0]
            if i != j:
                file_path = "corr." + str(i) + "m." + str(j) + "c"
                datafile = open(file_path)
                data = json.load(datafile)
                datafile.close()
                corr[:,j,i] = data[0]
    corrobj = hdf5file.createArray('/corr', 'corr2', corr)
    corrobj.setAttr('ftype0', json.dumps(['counts']))
    corrobj.setAttr('ftype1', json.dumps(['mags']))
    
    for i in range(nf):
        for j in range(i,nf):
            file_path = "corr." + str(i) + "m." + str(j) + "m"
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()
            corr[:,i,j] = data[0]
            if i != j:
                corr[:,j,i] = data[0]
    corrobj = hdf5file.createArray('/corr', 'corr3', corr)
    corrobj.setAttr('ftype0', json.dumps(['mags']))
    corrobj.setAttr('ftype1', json.dumps(['mags']))

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
corrobj.setAttr('ftype0', json.dumps(['counts']))
corrobj.setAttr('ftype1', json.dumps(['counts']))
corrobj.setAttr('ftype2', json.dumps(['counts']))
corrobj.setAttr('ftype3', json.dumps(['counts']))

if( use_mags ):
    for i in range(nf):
        for j in range(i,nf):
        
            file_path = "cov." + str(i) + "c." + str(j) + "m"
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()

            cov[:,:,i,j,i,j] = data

    corrobj = hdf5file.createArray('/cov', 'cov2', cov)
    corrobj.setAttr('ftype0', json.dumps(['counts']))
    corrobj.setAttr('ftype1', json.dumps(['mags']))
    corrobj.setAttr('ftype2', json.dumps(['counts']))
    corrobj.setAttr('ftype3', json.dumps(['mags']))
    
    for i in range(nf):
        for j in range(i,nf):
        
            file_path = "cov." + str(i) + "m." + str(j) + "c"
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()

            cov[:,:,i,j,i,j] = data

    corrobj = hdf5file.createArray('/cov', 'cov3', cov)
    corrobj.setAttr('ftype0', json.dumps(['mags']))
    corrobj.setAttr('ftype1', json.dumps(['counts']))
    corrobj.setAttr('ftype2', json.dumps(['mags']))
    corrobj.setAttr('ftype3', json.dumps(['counts']))
    
    for i in range(nf):
        for j in range(i,nf):

            file_path = "cov." + str(i) + "m." + str(j) + "m"
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()

            cov[:,:,i,j,i,j] = data
            cov[:,:,j,i,i,j] = data
            cov[:,:,i,j,j,i] = data
            cov[:,:,j,i,j,i] = data

    corrobj = hdf5file.createArray('/cov', 'cov4', cov)
    corrobj.setAttr('ftype0', json.dumps(['mags']))
    corrobj.setAttr('ftype1', json.dumps(['mags']))
    corrobj.setAttr('ftype2', json.dumps(['mags']))
    corrobj.setAttr('ftype3', json.dumps(['mags']))


hdf5file.close()

# the end.
