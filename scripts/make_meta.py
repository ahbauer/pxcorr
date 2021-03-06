#!/usr/bin/env python
# encoding: UTF8

import json
import os
import pdb
import numpy as np
import tables

# write an hdf5 file with metadata, for input to other code.

def make_metadata(f, z_mean, z_width, ang_mean, ang_width, mag_cuts, pops, ftypes):
    # f = tables.openFile('metadata.hdf5', 'w')
    
    f.createGroup('/', 'meta')
    meta = f.createArray('/meta', 'meta', np.ones(0))

    meta.setAttr("ang_mean", json.dumps(ang_mean))
    meta.setAttr('ang_width', json.dumps(ang_width))
    meta.setAttr('fourier', json.dumps(False))
    meta.setAttr('pop', json.dumps(pops))
    # meta.setAttr('ftype', json.dumps(ftypes))

    for i in range(len(pops)):
        f.createGroup('/meta', pops[i])
        groupname = '/meta/' + pops[i]
        metapop = f.createArray(groupname, 'meta', np.ones(0))
        metapop.setAttr('z_mean', json.dumps(z_mean[i]))
        metapop.setAttr('z_width', json.dumps(z_width[i]))
        metapop.setAttr('mag_limit', json.dumps(mag_cuts[i]))
        metapop.setAttr('ftype', json.dumps(ftypes[i]))

def make_metadata_task(f, catalogs, ang_mean, ang_width):
    
    f.createGroup('/', 'meta')
    meta = f.createArray('/meta', 'meta', np.ones(1))

    meta.setAttr("ang_mean", json.dumps(ang_mean))
    meta.setAttr('ang_width', json.dumps(ang_width))
    meta.setAttr('fourier', json.dumps(False))
    
    for catalog in catalogs:
        f.createGroup('/meta', catalog['pop'])
        groupname = '/meta/' + catalog['pop']
        metapop = f.createArray(groupname, 'meta', np.ones(1))
        metapop.setAttr('ftype', json.dumps(catalog['ftype']))
        z_means = []
        z_widths = []
        mag_cuts = []
        ftypes = catalog['ftype']
        for cut in catalog['cuts']:
            z_means.append(cut.get('z_mean',0.01))
            z_widths.append(cut.get('z_width',0.01))
            mag_cuts.append(cut.get('mag',999))
        metapop.setAttr('z_mean', json.dumps(z_means))
        metapop.setAttr('z_width', json.dumps(z_widths))
        metapop.setAttr('mag_limit', json.dumps(mag_cuts))
        metapop.setAttr('ftype', json.dumps(ftypes))
    
def make_metadata3pt(f, z_mean, z_width, angles, r_12s, r_23s, pop):
    # f = tables.openFile('metadata.hdf5', 'w')

    f.createGroup('/', 'meta')
    meta = f.createArray('/meta', 'meta', np.ones(0))

    meta.setAttr("angles", json.dumps(angles))
    meta.setAttr('r_12s', json.dumps(r_12s))
    meta.setAttr('r_23s', json.dumps(r_23s))
    meta.setAttr('fourier', json.dumps(False))
    meta.setAttr('pop', json.dumps([pop]))

    f.createGroup('/meta', pop)
    groupname = '/meta/' + pop
    metapop = f.createArray(groupname, 'meta', np.ones(0))
    metapop.setAttr('z_mean', json.dumps(z_mean))
    metapop.setAttr('z_width', json.dumps(z_width))

    # f.close()

    info = []
    info.append(angles)
    info.append(r_12s)
    info.append(r_23s)
    
    return info


if __name__ == '__main__':
    make_metadata()