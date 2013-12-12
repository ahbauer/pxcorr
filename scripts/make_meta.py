#!/usr/bin/env python
# encoding: UTF8

import json
import os
import pdb
import numpy as np
import tables

# write an hdf5 file with metadata, for input to other code.

def make_metadata(f, z_mean, z_width, ang_mean, ang_width, pop):
    # f = tables.openFile('metadata.hdf5', 'w')
    
    f.createGroup('/', 'meta')
    meta = f.createArray('/meta', 'meta', np.ones(0))

    meta.setAttr("ang_mean", json.dumps(ang_mean))
    meta.setAttr('ang_width', json.dumps(ang_width))
    meta.setAttr('fourier', json.dumps(False))
    meta.setAttr('pop', json.dumps([pop]))

    f.createGroup('/meta', pop)
    groupname = '/meta/' + pop
    metapop = f.createArray(groupname, 'meta', np.ones(0))
    metapop.setAttr('z_mean', json.dumps(z_mean))
    metapop.setAttr('z_width', json.dumps(z_width))
    
    # f.close()
    
    # ang_info = np.zeros((2,len(ang_mean)))
    # ang_info[0:] = ang_mean
    # ang_info[1:] = ang_width
    # return ang_info

    
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