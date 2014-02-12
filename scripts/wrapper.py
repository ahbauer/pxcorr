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

import partpix_map
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

# use_counts = config["use_counts"]
# use_mags = config["use_mags"]
only_auto = config["only_auto"]
only_makemaps = config["only_makemaps"]
only_correlate = config["only_correlate"]
pops = config["pops"]
ftypes = config["ftypes"]


nf = len(catalog_filenames)
nr = len(ang_means)

nobjs = np.ones(nf)
slopes_filenames = ['']*nf # this is ok because strings are immutable...

if not only_correlate:
    
    # do the following for each map
    # note: this can be parallelized and done separately for each map
    # the maps need to be copied back from the nodes, though.
    
    for c in range(nf):
        for bin in range(len(z_means[c])):
            # the suffix is the identifier we add to the filenames saved inside construct_inputs to know that it's from this catalog, etc.
            # it will be something associated with the job id or whatever.
            suffix = str(c) + "z" + str(bin)
            nobjs[c], slopes_filenames[c], nofz_filenames[c] = construct_inputs( suffix, catalog_filenames[c], ftypes[c], ang_means, ang_widths, nofz_filenames[c], mask_filenames[c][bin], z_means[c][bin], z_widths[c][bin], mag_cuts[c][bin] )
    
    if only_makemaps:
        exit(0)



# initialize the lists that will hold all the correlation function results
# indices: file1 ftype_for_file1 file2 ftype_for_file2 radial_bin z_bin_file1 z_bin_file2
corrs_all = []
jks_all = []
for i in range(nf):
    corrs_all.append([])
    jks_all.append([])
    nzi = len(z_means[i])
    for fi in range(len(ftypes[i])):
        corrs_all[i].append([])
        jks_all[i].append([])
        for j in range(nf):
            corrs_all[i][fi].append([])
            jks_all[i][fi].append([])
            nzj = len(z_means[j])
            for fj in range(len(ftypes[j])):
                corrs_all[i][fi][j].append(np.ones((nr, nzi, nzj)))
                jks_all[i][fi][j].append([])
                for r in range(nr):
                    jks_all[i][fi][j][fj].append([])
                    for zi in range(nzi):
                        jks_all[i][fi][j][fj][r].append([])
                        for zj in range(nzj):
                            # jks_all[i][fi][j][fj][r][zi].append([None])
                            jk_tmp = partpix_map.Partpix_Map()
                            jks_all[i][fi][j][fj][r][zi].append(jk_tmp)


for i in range(nf):
    nzi = len(z_means[i])
    for fi in range(len(ftypes[i])):
        fi_char = ""
        if ftypes[i][fi] == 'counts':
            fi_char = 'c'
        elif ftypes[i][fi] == 'mags':
            fi_char = 'm'
        else:
            err_msg = "ftype isn't understood: " + ftypes[i][fi]
            raise Exception, err_msg
        
        for j in range(i, nf):
            nzj = len(z_means[j])
            for fj in range(len(ftypes[j])):
                fj_char = ""
                if ftypes[j][fj] == 'counts':
                    fj_char = 'c'
                elif ftypes[j][fj] == 'mags':
                    fj_char = 'm'
                else:
                    err_msg = "ftype isn't understood: " + ftypes[j][fj]
                    raise Exception, err_msg
                    
                for zi in range(nzi):
                    for zj in range(nzj):
                        
                        name1 = "d" + fi_char + "_map" + str(i) + "z" + str(zi) + ".h5"
                        name2 = "d" + fj_char + "_map" + str(j) + "z" + str(zj) + ".h5"
                        if( only_auto and not (i == j and fi == fj and zi == zj) ):
                            continue
                        print "Correlating %s and %s" %(name1, name2)
        
                        for r in range(nr):
                            suffix = "." + str(i) + "z" + str(zi) + fi_char + "." + str(j) + "z" + str(zj) + fj_char + "." + str(r)
                            print "Calling for r=%d" %r
                            result = correlate.correlate( name1, name2, suffix, r )
                            corrs_all[i][fi][j][fj][r,zi,zj] = result[0]
                            # jks_all[i][fi][j][fj][r][zi][zj] = result[1:]
                            jks_all[i][fi][j][fj][r][zi][zj].read_from_ascii('jks' + suffix + '.dat')


# packobs



print "Starting hdf5 file stuff"
hdf5file = tables.openFile('pxcorr_out.h5', 'w')
make_metadata( hdf5file, z_means, z_widths, ang_means, ang_widths, mag_cuts, pops, ftypes)
for index in range(nf):
    slopes_to_hdf5( hdf5file, slopes_filenames[index], pops[index] )
    nofz_to_hdf5( hdf5file, nofz_filenames[index], pops[index] )



print "Adding correlations and covariances"

# add correlations to the file
hdf5file.createGroup('/', 'corr')
hdf5file.createGroup('/', 'cov')

corr_count = 1
for i in range(nf):
    for fi in range(len(ftypes[i])):
        for j in range(i,nf):
            for fj in range(len(ftypes[j])):
                
                corrobj = hdf5file.createArray('/corr', 'corr' + str(corr_count), corrs_all[i][fi][j][fj])
                corrobj.setAttr('ftype0', json.dumps(ftypes[i][fi]))
                corrobj.setAttr('ftype1', json.dumps(ftypes[j][fj]))
                corrobj.setAttr('pop0', json.dumps(pops[i]))
                corrobj.setAttr('pop1', json.dumps(pops[j]))

                # add covariance to the file

                print "Calculating covariances"
                nzi = len(z_means[i])
                nzj = len(z_means[j])
                cov = np.zeros((nr, nr, nzi, nzj, nzi, nzj))
                njk = 0
                for zi in range(nzi):
                    for zj in range(nzj):
                        for zk in range(nzi):
                            for zl in range(nzj):
                                for r1 in range(nr):
                                    jk1 = jks_all[i][fi][j][fj][r1][zi][zj]
                                    jk1_mean = jk1.mean() # np.mean(jk1)
                                    for r2 in range(nr):
                                        jk2 = jks_all[i][fi][j][fj][r2][zk][zl]
                                        jk2_mean = jk2.mean() # np.mean(jk2)
                                        try:
                                            jk1_intersection, jk2_intersection = jk1.intersection(jk2)
                                            if jk1_intersection.npartpix < 0.5*jk1.npartpix:
                                                print "Warning, small jackknife overlap area!! %d vs %d pixels" %(jk1_intersection.npartpix, jk1.npartpix)
                                            # if njk != len(jk2):
                                            #     print "Error, %d and %d jks are sizes %d != %d" %(r1,r2,len(jk1),len(jk2))
                                            #     raise AssertionError, "Jackknife region number mismatch"
                                            njk = jk1_intersection.npartpix
                                            for jk in range(njk):
                                                # cov[r1,r2,zi,zj,zk,zl] = cov[r1,r2,zi,zj,zk,zl] + (jk1[jk]-jk1_mean)*(jk2[jk]-jk2_mean)
                                                cov[r1,r2,zi,zj,zk,zl] = cov[r1,r2,zi,zj,zk,zl] + (jk1_intersection.partmap[jk]-jk1_mean)*(jk2_intersection.partmap[jk]-jk2_mean)
                                            cov[r1,r2,zi,zj,zk,zl] = (float(njk-1.0))/float(njk) * cov[r1,r2,zi,zj,zk,zl]
                                            print "Successfully set covariance = %e for zi, zj, zk, zl, r1, r2: %d %d %d %d %d %d" %(cov[r1,r2,zi,zj,zk,zl], zi, zj, zk, zl, r1, r2)
                                        except AssertionError:
                                            print "Setting covariance to 0 for zi, zj, zk, zl, r1, r2: %d %d %d %d %d %d" %(zi, zj, zk, zl, r1, r2)
                                        
                covobj = hdf5file.createArray('/cov', 'cov' + str(corr_count), cov)
                covobj.setAttr('ftype0', json.dumps(ftypes[i][fi]))
                covobj.setAttr('ftype1', json.dumps(ftypes[j][fj]))
                covobj.setAttr('ftype2', json.dumps(ftypes[i][fi]))
                covobj.setAttr('ftype3', json.dumps(ftypes[j][fj]))
                covobj.setAttr('pop0', json.dumps(pops[i]))
                covobj.setAttr('pop1', json.dumps(pops[j]))
                covobj.setAttr('pop2', json.dumps(pops[i]))
                covobj.setAttr('pop3', json.dumps(pops[j]))
                
                if( ftypes[i] == 'counts' and ftypes[j] == 'counts' ):
                    noise_to_hdf5( hdf5file, [pops[i],pops[j]], ['counts','counts'], corr_count, nobjs ) # only for counts autocorr for now...
                
                corr_count += 1
                

# TBD: include area fraction to account for jackknife area being a bit smaller than the total?
# Would need to get that info from the correlation code...  
# Should be a small correction, and leaving it out makes the errors a bit bigger (conservative).


# print corr[:,0,0]
# print cov[:,:,0,0,0,0]

hdf5file.close()

# the end.
