#!/usr/bin/env python
# encoding: UTF8

import sys
import json
import os
import pdb
import numpy as np
import tables
import ctypes
import yaml

from make_meta import make_metadata
from parse_sample import parse_data, add_nofz
from packobs import add_corr, add_cov

sys.path.append("/Users/bauer/software/pxcorr/src/")
import make_maps
import correlate

if len(sys.argv) != 2:
    print "Usage: wrapper.py config_file"
    exit(0)

config_file = sys.argv[1]

config = yaml.load(open(config_file, 'r'))
catalog_filename = config["catalog_filename"]
mask_filename = config["mask_filename"]
nofz_filename = config["nofz_filename"]
mag_cuts = config["mag_cuts"]
use_counts = config["use_counts"]
use_mags = config["use_mags"]
pop = config["pop"]
z_mean = config["z_mean"]
z_width = config["z_width"]
ang_mean = config["ang_mean"]
ang_width = config["ang_width"]
only_auto = config["only_auto"]
only_makemaps = config["only_makemaps"]

# catalog_filename = "/Users/bauer/correlations/poldr7/cross_speclrgs/speclrg_notdr9sz_radeczzmag" #"/Users/bauer/surveys/SDSS/megaz/pol/dr7/radecz1z1i"
# mask_filename = "/Users/bauer/correlations/poldr7/cross_speclrgs/megaz_dr7.1024.mask.dat"
# nofz_filename = None
# 
# 
# mag_cuts = [ 99. ]
# use_counts = True
# use_mags = False
# pop = 'faint'
# z_mean = [ 0.25 ]
# z_width = [ 0.15 ]
# ang_mean = [ 0.02, 0.04, 0.06, 0.08, 0.11, 0.15, 0.19, 0.23, 0.27, 0.31, 0.36, 0.42, 0.5, 0.6, 0.7, 0.8, 0.9 ]
# ang_width = [ 0.02, 0.02, 0.02, 0.02, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.06, 0.06, 0.1, 0.1, 0.1, 0.1, 0.1 ]
# only_auto = False
# only_makemaps = True

# make some metadata to tell us what to do
hdf5file = tables.openFile('pxcorr_out.h5', 'w')
ang_info = make_metadata( hdf5file, z_mean, z_width, ang_mean, ang_width, pop)
print ang_info

# do we want to figure out the N(z) from the catalog or a separate training set file?
nofz_from_data = False
if nofz_filename == catalog_filename or nofz_filename is None or nofz_filename == 'None':
    nofz_from_data = True

# parse the catalog: construct N(z) data, calculate slope, etc.
subcat_filenames = parse_data( catalog_filename, mag_cuts, hdf5file, add_nofz=False, sparse=True )
# print subcat_filenames

# if the N(z) is from a separate file, put that in the hdf5 file now.
# note that this info must correspond to (have the same cuts as) the data catalog!
if not nofz_from_data:
    add_nofz( nofz_filename, hdf5file, pop )

for fi, filename in enumerate(subcat_filenames):
    suffix = "_" + str(fi)
    make_maps.make_maps( filename, mask_filename, ang_mean, ang_width, use_counts, use_mags, suffix );
    print "Finished making map from %s" %filename

if only_makemaps:
    hdf5file.close()
    exit(0)

nf = len( subcat_filenames )

for i in range(nf):
    for j in range(i, nf):
        name1 = "dc_map_" + str(i) + ".h5"
        name2 = "dc_map_" + str(j) + ".h5"
        if( only_auto and i != j ):
            continue
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
corrobj.setAttr('ftype0', json.dumps('counts'))
corrobj.setAttr('pop0', json.dumps('faint'))
corrobj.setAttr('ftype1', json.dumps('counts'))
corrobj.setAttr('pop1', json.dumps('faint'))

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
    corrobj.setAttr('ftype0', json.dumps('counts'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    
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
    corrobj.setAttr('ftype0', json.dumps('mags'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
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

if( use_mags ):
    for i in range(nf):
        for j in range(i,nf):
        
            file_path = "cov." + str(i) + "c." + str(j) + "m"
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()

            cov[:,:,i,j,i,j] = data

    corrobj = hdf5file.createArray('/cov', 'cov2', cov)
    corrobj.setAttr('ftype0', json.dumps('counts'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    corrobj.setAttr('ftype2', json.dumps('counts'))
    corrobj.setAttr('pop2', json.dumps('faint'))
    corrobj.setAttr('ftype3', json.dumps('mags'))
    corrobj.setAttr('pop3', json.dumps('faint'))
    
    for i in range(nf):
        for j in range(i,nf):
        
            file_path = "cov." + str(i) + "m." + str(j) + "c"
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()

            cov[:,:,i,j,i,j] = data

    corrobj = hdf5file.createArray('/cov', 'cov3', cov)
    corrobj.setAttr('ftype0', json.dumps('mags'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('counts'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    corrobj.setAttr('ftype2', json.dumps('mags'))
    corrobj.setAttr('pop2', json.dumps('faint'))
    corrobj.setAttr('ftype3', json.dumps('counts'))
    corrobj.setAttr('pop3', json.dumps('faint'))
    
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
    corrobj.setAttr('ftype0', json.dumps('mags'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    corrobj.setAttr('ftype2', json.dumps('mags'))
    corrobj.setAttr('pop2', json.dumps('faint'))
    corrobj.setAttr('ftype3', json.dumps('mags'))
    corrobj.setAttr('pop3', json.dumps('faint'))


hdf5file.close()

# the end.
