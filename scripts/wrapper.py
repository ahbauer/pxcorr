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
from parse_sample import nofz_to_hdf5, noise_to_hdf5, slopes_to_hdf5, construct_inputs
from packobs import add_corr, add_cov

sys.path.append("/Users/bauer/software/pxcorr/src/")
import correlate

if len(sys.argv) != 2:
    print "Usage: wrapper.py config_file"
    exit(0)

config_file = sys.argv[1]

config = yaml.load(open(config_file, 'r'))
catalog_filenames = config["catalog_filenames"]
mask_filenames = config["mask_filenames"]
nofz_filenames = config["nofz_filenames"]
mag_cuts = config["mag_cuts"]
z_means = config["z_means"]
z_widths = config["z_widths"]

ang_means = config["ang_mean"]
ang_widths = config["ang_width"]

use_counts = config["use_counts"]
use_mags = config["use_mags"]
only_auto = config["only_auto"]
only_makemaps = config["only_makemaps"]
only_correlate = config["only_correlate"]
pop = config["pop"]



nf = len(catalog_filenames)
nr = len(ang_means)

nobjs = np.ones(nf)
slopes_filenames = ['']*nf # this is ok because strings are immutable...

if not only_correlate:
    
    # do the following for each map
    # note: this can be parallelized and done separately for each map
    # the maps need to be copied back from the nodes, though.
    
    for c in range(nf):
        # the suffix is the identifier we add to the filenames saved inside construct_inputs to know that it's from this catalog, etc.
        # it will be something associated with the job id or whatever.
        suffix = "_" + str(c)
        nobjs[c], slopes_filenames[c], nofz_filenames[c] = construct_inputs( suffix, catalog_filenames[c], nofz_filenames[c], mask_filenames[c], ang_means, ang_widths, z_means[c], z_widths[c], mag_cuts[c], use_counts, use_mags )    
    
    if only_makemaps:
        exit(0)



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

corr = np.ones((nr, nf, nf))
corrcm = None
corrmc = None
corrmm = None
if use_mags:
    corrcm = np.ones((nr, nf, nf))
    corrmc = np.ones((nr, nf, nf))
    corrmm = np.ones((nr, nf, nf))

jks = []
jkscm = []
jksmc = []
jksmm = []
for r in range(nr):
    jks.append([])
    jkscm.append([])
    jksmc.append([])
    jksmm.append([])
    for i in range(nf):
        jks[r].append([])
        jkscm[r].append([])
        jksmc[r].append([])
        jksmm[r].append([])
        for j in range(nf):
            jks[r][i].append([None])
            jkscm[r][i].append([None])
            jksmc[r][i].append([None])
            jksmm[r][i].append([None])

for i in range(nf):
    for j in range(i, nf):
        name1 = "dc_map_" + str(i) + ".h5"
        name2 = "dc_map_" + str(j) + ".h5"
        if( only_auto and i != j ):
            continue
        print "Correlating %s and %s" %(name1, name2)
        
        for r in range(nr):
            suffix = "." + str(i) + "c." + str(j) + "c." + str(r)
            print "Calling for r=%d" %r
            result = correlate.correlate( name1, name2, suffix, r )
            corr[r,i,j] = result[0]
            jks[r][i][j]=result[1:]
            if i != j:
                corr[r,j,i] = corr[r,i,j]
                jks[r][j][i] = jks[r][i][j]
        
        if use_mags:
            name2 = "dm_map_" + str(j) + ".h5"
            print "Correlating %s and %s" %(name1, name2)
            for r in range(nr):
                suffix = "." + str(i) + "c." + str(j) + "m." + str(r)
                result = correlate.correlate( name1, name2, suffix, r )
                corrcm[r,i,j] = result[0]
                jkscm[r][i][j] = result[1:]
                if i != j:
                    corrcm[r,j,i] = corrcm[r,i,j]
                    jkscm[r][j][i] = jkscm[r][i][j]
            name1 = "dm_map_" + str(i) + ".h5"
            print "Correlating %s and %s" %(name1, name2)
            for r in range(nr):
                suffix = "." + str(i) + "m." + str(j) + "m." + str(r)
                result = correlate.correlate( name1, name2, suffix, r )
                corrmm[r,i,j] = result[0]
                jksmm[r][i][j] = result[1:]
                if i != j:
                    corrmm[r,j,i] = corrmm[r,i,j]
                    jksmm[r][j][i] = jksmm[r][i][j]
            name2 = "dc_map_" + str(j) + ".h5"
            print "Correlating %s and %s" %(name1, name2)
            for r in range(nr):
                suffix = "." + str(i) + "m." + str(j) + "c." + str(r)
                result = correlate.correlate( name1, name2, suffix, r )
                corrmc[r,i,j] = result[0]
                jksmc[r][i][j] = result[1:]
                if i != j:
                    corrmc[r,j,i] = corrmc[r,i,j]
                    jksmc[r][j][i] = jksmc[r][i][j]

# packobs



print "Starting hdf5 file stuff"
hdf5file = tables.openFile('pxcorr_out.h5', 'w')
make_metadata( hdf5file, z_means, z_widths, ang_means, ang_widths, pop)
noise_to_hdf5( hdf5file, pop, nobjs )
for index in range(nf):
    slopes_to_hdf5( hdf5file, slopes_filenames[index], index, pop )
    nofz_to_hdf5( hdf5file, nofz_filenames[index], index, pop )



print "Adding correlations"

# add correlations to the file
hdf5file.createGroup('/', 'corr')

corrobj = hdf5file.createArray('/corr', 'corr1', corr)
corrobj.setAttr('ftype0', json.dumps('counts'))
corrobj.setAttr('pop0', json.dumps('faint'))
corrobj.setAttr('ftype1', json.dumps('counts'))
corrobj.setAttr('pop1', json.dumps('faint'))

if( use_mags ):
    
    corrobj = hdf5file.createArray('/corr', 'corr2', corrcm)
    corrobj.setAttr('ftype0', json.dumps('counts'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    
    corrobj = hdf5file.createArray('/corr', 'corr3', corrmc)
    corrobj.setAttr('ftype0', json.dumps('mags'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('counts'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    
    corrobj = hdf5file.createArray('/corr', 'corr4', corrmm)
    corrobj.setAttr('ftype0', json.dumps('mags'))
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('mags'))
    corrobj.setAttr('pop1', json.dumps('faint'))

# add covariance to the file

cov = np.zeros((nr, nr, len(z_means), len(z_means), len(z_means), len(z_means)))

print "Calculating covarinances"
njk = 0
for i in range(nf):
    for j in range(nf):
        for k in range(nf):
            for l in range(nf):
                for r1 in range(nr):
                    jk1 = jks[r1][i][j]
                    jk1_mean = np.mean(jk1)
                    # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
                    njk = len(jk1)
                    for r2 in range(nr):
                        jk2 = jks[r2][k][l]
                        jk2_mean = np.mean(jk2)
                        # print "r2 %d: read in %d jk values, mean %e (first %e)" %(r2, len(jk2), jk2_mean, jk2[0])
                        if njk != len(jk2):
                            print "Error, %d and %d jks are sizes %d != %d" %(r1,r2,len(jk1),len(jk2))
                        for jk in range(njk):
                            cov[r1,r2,i,j,k,l] = cov[r1,r2,i,j,k,l] + (jk1[jk]-jk1_mean)*(jk2[jk]-jk2_mean)

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

print corr[:,0,0]
print cov[:,:,0,0,0,0]

if( use_mags ):

    cov = np.zeros((nr, nr, len(z_means), len(z_means), len(z_means), len(z_means)))

    njk = 0
    for i in range(nf):
        for j in range(i,nf):
            for r1 in range(nr):
                jk1 = jkscm[r1][i][j]
                jk1_mean = np.mean(jk1)
                # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
                njk = len(jk1)
                for r2 in range(nr):
                    jk2 = jkscm[r2][i][j]
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
    
    
    cov = np.zeros((nr, nr, len(z_means), len(z_means), len(z_means), len(z_means)))

    njk = 0
    for i in range(nf):
        for j in range(i,nf):
            for r1 in range(nr):
                jk1 = jksmc[r1][i][j]
                jk1_mean = np.mean(jk1)
                # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
                njk = len(jk1)
                for r2 in range(nr):
                    jk2 = jksmc[r2][i][j]
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
    
    
    cov = np.zeros((nr, nr, len(z_means), len(z_means), len(z_means), len(z_means)))

    njk = 0
    for i in range(nf):
        for j in range(i,nf):
            for r1 in range(nr):
                jk1 = jksmm[r1][i][j]
                jk1_mean = np.mean(jk1)
                # print "r1 %d: read in %d jk values, mean %e (first %e)" %(r1, len(jk1), jk1_mean, jk1[0])
                njk = len(jk1)
                for r2 in range(nr):
                    jk2 = jksmm[r2][i][j]
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
