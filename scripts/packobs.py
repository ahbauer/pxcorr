#!/usr/bin/env python
# encoding: UTF8

import json
import os
import pdb
import numpy as np
import tables

des_dir = '/Users/bauer/correlations/pic_test/catalogs'
toadd = [\
(0,0,'results00'),
(0,1,'results01'),
(0,2,'results02'),
(0,3,'results03'),
(1,1,'results11'),
(1,2,'results12'),
(1,3,'results13'),
(2,2,'results22'),
(2,3,'results23'),
(3,3,'results33')]

nbins = 4
nangles = 18

# no longer need:
metafiles = ['metadata_faint2.txt', 'metadata_faint6.txt', 'metadata_faint10.txt', 'metadata_faint14.txt']
z_mean = [0.24, 0.45, 0.70, 1.0]
z_width = [0.1, 0.15, 0.2, 0.22]


def define_meta(f):
    f.createGroup('/', 'meta')
    meta = f.createArray('/meta', 'meta', np.ones(0))

    metafile = open(os.path.join(des_dir, metafiles[0]))
    metadata = json.load(metafile)
    metafile.close()
    u_mean = metadata["ang_mean"]
    assert len(u_mean) == nangles
    u_width = metadata["ang_width"]

    #print u_mean
    #print u_width

    meta.setAttr("ang_mean", json.dumps(list(u_mean)))
    meta.setAttr('ang_width', json.dumps(list(u_width)))
    meta.setAttr('fourier', json.dumps(False))
    meta.setAttr('pop', json.dumps(['faint']))

    f.createGroup('/meta', 'faint')
    metafaint = f.createArray('/meta/faint', 'meta', np.ones(0))
    metafaint.setAttr('z_mean', json.dumps(z_mean))
    metafaint.setAttr('z_width', json.dumps(z_width))

def add_corr(f):
    corr = np.ones((nangles, nbins, nbins))

    for i,j, file_name in toadd:
        file_path = os.path.join(des_dir, file_name, "correlation")
        datafile = open(file_path)
        data = json.load(datafile)
        datafile.close()

        #print data[0]
        corr[:,i,j] = data[0]

    f.createGroup('/', 'corr')
    corrobj = f.createArray('/corr', 'corr1', corr)
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype0', json.dumps('counts'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('counts'))

def add_cov(f):
    cov = np.ones((nangles, nangles, nbins, nbins))

    for i,j, file_name in toadd:
        file_path = os.path.join(des_dir, file_name, "covariance")
        datafile = open(file_path)
        data = json.load(datafile)
        datafile.close()

        cov[:,:,i,j] = data

    f.createGroup('/', 'cov')
    corrobj = f.createArray('/cov', 'cov1', cov)
    corrobj.setAttr('pop0', json.dumps('faint'))
    corrobj.setAttr('ftype0', json.dumps('counts'))
    corrobj.setAttr('pop1', json.dumps('faint'))
    corrobj.setAttr('ftype1', json.dumps('counts'))
    corrobj.setAttr('pop2', json.dumps('faint'))
    corrobj.setAttr('ftype2', json.dumps('counts'))
    corrobj.setAttr('pop3', json.dumps('faint'))
    corrobj.setAttr('ftype3', json.dumps('counts'))

def main():
    f = tables.openFile('mycorr.hdf5', 'w')
    define_meta(f)
    add_corr(f)
    add_cov(f)
    f.close()

if __name__ == '__main__':
    main()
