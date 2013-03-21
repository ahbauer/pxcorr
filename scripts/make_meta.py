#!/usr/bin/env python
# encoding: UTF8

import json
import os
import pdb
import numpy as np
import tables

# write an hdf5 file with metadata, for input to other code.
pop = 'faint'

z_mean = [0.24, 0.45, 0.70, 1.0]
z_width = [0.1, 0.15, 0.2, 0.22]

ang_mean = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
ang_width = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]


def make_metadata(f):
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
    
    ang_info = np.zeros((2,len(ang_mean)))
    ang_info[0:] = ang_mean
    ang_info[1:] = ang_width
    return ang_info


if __name__ == '__main__':
    make_metadata()