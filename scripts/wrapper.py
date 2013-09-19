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


# make some metadata to tell us what to do
hdf5file = tables.openFile('pxcorr_out.h5', 'w')
ang_info = make_metadata( hdf5file, z_mean, z_width, ang_mean, ang_width, pop)
print ang_info

# do we want to figure out the N(z) from the catalog or a separate training set file?
nofz_from_data = False
if nofz_filename == catalog_filename or nofz_filename is None or nofz_filename == 'None':
    nofz_from_data = True

# parse the catalog: construct N(z) data, calculate slope, etc.
subcat_filenames = parse_data( catalog_filename, mag_cuts, hdf5file, nofz_from_data, sparse=True )
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

# FOR JORGE:
# there are 3 loops here: i & j are over redshift bins, and r is over radial bins.
# now i'm just looping serially through them all, but i would like you to send 
# each ijr combination as a separate job to a node.  
# basically, each call to correlate.correlate could be replaced by a something 
# that calls correlate.correlate on a node.  the other stuff about defining 
# "name1" and "name2" and "suffix" are just to make the input parameters for correlate.  
# then you'd have to wait for all of the jobs to be done before going on to the 
# "Adding correlations" part, which collects all of the little text files written by 
# the correlate jobs (called "file_path" below;  i guess you'll have to do something 
# about the file locations too.  correlate.correlate currently puts those output files 
# in the directory from which it is called.).

for i in range(nf):
    for j in range(i, nf):
        name1 = "dc_map_" + str(i) + ".h5"
        name2 = "dc_map_" + str(j) + ".h5"
        if( only_auto and i != j ):
            continue
        print "Correlating %s and %s" %(name1, name2)
        
        for r in range(len(config["ang_mean"])):
            suffix = "." + str(i) + "c." + str(j) + "c." + str(r)
            correlate.correlate( name1, name2, suffix, r )
        
        if use_mags:
            name2 = "dm_map_" + str(j) + ".h5"
            print "Correlating %s and %s" %(name1, name2)
            for r in range(len(config["ang_mean"])):
                suffix = "." + str(i) + "c." + str(j) + "m." + str(r)
                correlate.correlate( name1, name2, suffix, r )
            name1 = "dm_map_" + str(i) + ".h5"
            print "Correlating %s and %s" %(name1, name2)
            for r in range(len(config["ang_mean"])):
                suffix = "." + str(i) + "m." + str(j) + "m." + str(r)
                correlate.correlate( name1, name2, suffix, r )
            name2 = "dc_map_" + str(j) + ".h5"
            print "Correlating %s and %s" %(name1, name2)
            for r in range(len(config["ang_mean"])):
                suffix = "." + str(i) + "m." + str(j) + "c." + str(r)
                correlate.correlate( name1, name2, suffix, r )

# packobs

print "Adding correlations"

# add correlations to the file
hdf5file.createGroup('/', 'corr')

corr = np.ones((len(ang_mean), len(z_mean), len(z_mean)))
for i in range(nf):
    for j in range(i,nf):
        for r in range(len(config["ang_mean"])):
            file_path = "corr." + str(i) + "c." + str(j) + "c." + str(r)
            datafile = open(file_path)
            data = json.load(datafile)
            datafile.close()
            corr[r,i,j] = data[0] # correlation value at r.
            if i != j:
                corr[r,j,i] = data[0]
corrobj = hdf5file.createArray('/corr', 'corr1', corr)
corrobj.setAttr('ftype0', json.dumps('counts'))
corrobj.setAttr('pop0', json.dumps('faint'))
corrobj.setAttr('ftype1', json.dumps('counts'))
corrobj.setAttr('pop1', json.dumps('faint'))

if( use_mags ):
    for i in range(nf):
        for j in range(i,nf):
            for r in range(len(config["ang_mean"])):
                file_path = "corr." + str(i) + "c." + str(j) + "m." + str(r)
                datafile = open(file_path)
                data = json.load(datafile)
                datafile.close()
                corr[r,i,j] = data[0]
                if i != j:
                    file_path = "corr." + str(i) + "m." + str(j) + "c." + str(r)
                    datafile = open(file_path)
                    data = json.load(datafile)
                    datafile.close()
                    corr[r,j,i] = data[0]
    corrobj = hdf5file.createArray('/corr', 'corr2', corr)
    corrobj.setAttr('ftype0', json.dumps('counts'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    
    for i in range(nf):
        for j in range(i,nf):
            for r in range(len(config["ang_mean"])):
                file_path = "corr." + str(i) + "m." + str(j) + "m." + str(r)
                datafile = open(file_path)
                data = json.load(datafile)
                datafile.close()
                corr[r,i,j] = data[0]
                if i != j:
                    corr[r,j,i] = data[0]
    corrobj = hdf5file.createArray('/corr', 'corr3', corr)
    corrobj.setAttr('ftype0', json.dumps('mags'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))

# add covariance to the file
cov = np.zeros((len(ang_mean), len(ang_mean), len(z_mean), len(z_mean), len(z_mean), len(z_mean)))


print "Calculating covarinances"
njk = 0
for i in range(nf):
    for j in range(i,nf):
        for r1 in range(len(ang_mean)):
            file_path = "jks." + str(i) + "c." + str(j) + "c." + str(r1)
            datafile = open(file_path)
            jk1 = json.load(datafile) # list of jackknife values
            datafile.close()
            jk1_mean = np.mean(jk1)
            # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
            njk = len(jk1)
            for r2 in range(len(ang_mean)):
                file_path = "jks." + str(i) + "c." + str(j) + "c." + str(r2)
                datafile = open(file_path)
                jk2 = json.load(datafile) # list of jackknife values
                datafile.close()
                jk2_mean = np.mean(jk2)
                # print "r2 %d: read in %d jk values, mean %e (first %e)" %(r2, len(jk2), jk2_mean, jk2[0])
                if njk != len(jk2):
                    print "Error, %d and %d jks are sizes %d != %d" %(r1,r2,len(jk1),len(jk2))
                for k in range(njk):
                    cov[r1,r2,i,j,i,j] = cov[r1,r2,i,j,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                    if i != j:
                        cov[r1,r2,j,i,i,j] = cov[r1,r2,j,i,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                        cov[r1,r2,i,j,j,i] = cov[r1,r2,i,j,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                        cov[r1,r2,j,i,j,i] = cov[r1,r2,j,i,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
cov = (float(njk-1.0))/float(njk) * cov

# TBD: include area fraction to account for jackknife area being a bit smaller than the total?
# Would need to get that info from the correlation code...  
# Should be a small correction, and leaving it out makes the errors a bit bigger (conservative).

print "Adding covariances"

hdf5file.createGroup('/', 'cov')

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

    cov = np.zeros((len(ang_mean), len(ang_mean), len(z_mean), len(z_mean), len(z_mean), len(z_mean)))

    njk = 0
    for i in range(nf):
        for j in range(i,nf):
            for r1 in range(len(ang_mean)):
                file_path = "jks." + str(i) + "c." + str(j) + "m." + str(r1)
                datafile = open(file_path)
                jk1 = json.load(datafile) # list of jackknife values
                datafile.close()
                jk1_mean = np.mean(jk1)
                # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
                njk = len(jk1)
                for r2 in range(len(ang_mean)):
                    file_path = "jks." + str(i) + "c." + str(j) + "m." + str(r2)
                    datafile = open(file_path)
                    jk2 = json.load(datafile) # list of jackknife values
                    datafile.close()
                    jk2_mean = np.mean(jk2)
                    # print "r2 %d: read in %d jk values, mean %e (first %e)" %(r2, len(jk2), jk2_mean, jk2[0])
                    if njk != len(jk2):
                        print "Error, %d and %d jks are sizes %d != %d" %(r1,r2,len(jk1),len(jk2))
                    for k in range(njk):
                        cov[r1,r2,i,j,i,j] = cov[r1,r2,i,j,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                        if i != j:
                            cov[r1,r2,j,i,i,j] = cov[r1,r2,j,i,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                            cov[r1,r2,i,j,j,i] = cov[r1,r2,i,j,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                            cov[r1,r2,j,i,j,i] = cov[r1,r2,j,i,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
    cov = (float(njk-1.0))/float(njk) * cov
    
    corrobj = hdf5file.createArray('/cov', 'cov2', cov)
    corrobj.setAttr('ftype0', json.dumps('counts'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    corrobj.setAttr('ftype2', json.dumps('counts'))
    corrobj.setAttr('pop2', json.dumps('faint'))
    corrobj.setAttr('ftype3', json.dumps('mags'))
    corrobj.setAttr('pop3', json.dumps('faint'))
    
    
    cov = np.zeros((len(ang_mean), len(ang_mean), len(z_mean), len(z_mean), len(z_mean), len(z_mean)))

    njk = 0
    for i in range(nf):
        for j in range(i,nf):
            for r1 in range(len(ang_mean)):
                file_path = "jks." + str(i) + "m." + str(j) + "c." + str(r1)
                datafile = open(file_path)
                jk1 = json.load(datafile) # list of jackknife values
                datafile.close()
                jk1_mean = np.mean(jk1)
                # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
                njk = len(jk1)
                for r2 in range(len(ang_mean)):
                    file_path = "jks." + str(i) + "m." + str(j) + "c." + str(r2)
                    datafile = open(file_path)
                    jk2 = json.load(datafile) # list of jackknife values
                    datafile.close()
                    jk2_mean = np.mean(jk2)
                    # print "r2 %d: read in %d jk values, mean %e (first %e)" %(r2, len(jk2), jk2_mean, jk2[0])
                    if njk != len(jk2):
                        print "Error, %d and %d jks are sizes %d != %d" %(r1,r2,len(jk1),len(jk2))
                    for k in range(njk):
                        cov[r1,r2,i,j,i,j] = cov[r1,r2,i,j,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                        if i != j:
                            cov[r1,r2,j,i,i,j] = cov[r1,r2,j,i,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                            cov[r1,r2,i,j,j,i] = cov[r1,r2,i,j,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                            cov[r1,r2,j,i,j,i] = cov[r1,r2,j,i,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
    cov = (float(njk-1.0))/float(njk) * cov

    corrobj = hdf5file.createArray('/cov', 'cov3', cov)
    corrobj.setAttr('ftype0', json.dumps('mags'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('counts'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    corrobj.setAttr('ftype2', json.dumps('mags'))
    corrobj.setAttr('pop2', json.dumps('faint'))
    corrobj.setAttr('ftype3', json.dumps('counts'))
    corrobj.setAttr('pop3', json.dumps('faint'))
    
    
    cov = np.zeros((len(ang_mean), len(ang_mean), len(z_mean), len(z_mean), len(z_mean), len(z_mean)))

    njk = 0
    for i in range(nf):
        for j in range(i,nf):
            for r1 in range(len(ang_mean)):
                file_path = "jks." + str(i) + "m." + str(j) + "m." + str(r1)
                datafile = open(file_path)
                jk1 = json.load(datafile) # list of jackknife values
                datafile.close()
                jk1_mean = np.mean(jk1)
                # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
                njk = len(jk1)
                for r2 in range(len(ang_mean)):
                    file_path = "jks." + str(i) + "m." + str(j) + "m." + str(r2)
                    datafile = open(file_path)
                    jk2 = json.load(datafile) # list of jackknife values
                    datafile.close()
                    jk2_mean = np.mean(jk2)
                    # print "r2 %d: read in %d jk values, mean %e (first %e)" %(r2, len(jk2), jk2_mean, jk2[0])
                    if njk != len(jk2):
                        print "Error, %d and %d jks are sizes %d != %d" %(r1,r2,len(jk1),len(jk2))
                    for k in range(njk):
                        cov[r1,r2,i,j,i,j] = cov[r1,r2,i,j,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                        if i != j:
                            cov[r1,r2,j,i,i,j] = cov[r1,r2,j,i,i,j] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                            cov[r1,r2,i,j,j,i] = cov[r1,r2,i,j,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
                            cov[r1,r2,j,i,j,i] = cov[r1,r2,j,i,j,i] + (jk1[k]-jk1_mean)*(jk2[k]-jk2_mean)
    cov = (float(njk-1.0))/float(njk) * cov
    
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
