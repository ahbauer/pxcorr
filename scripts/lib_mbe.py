#!/usr/bin/env python
# encoding: UTF8
# Library to ease setting and reading json
# formated attrictures in HDF5 files.


# this file is copied from martin eriksen, so that we don't have to do
# from sci.lib.libh5attr import h5getattr, h5setattr

import pdb
import json

import numpy as np
try:
    import tables
except ImportError:
    pass

# - All attributes are assumed to be JSON object.
# - Lists in JSON are converted to numpy arrays.

def h5getattr(obj, attr):
    """Read in attribute and convert JSON."""

    obj_attrs = getattr(obj, 'attrs')

    assert attr in obj_attrs, 'Missing attribute.'
    try:
        injson = getattr(obj_attrs, attr)
    except TypeError:
        pdb.set_trace()
        #raise TypeError, 'Error loading attribute: '+ attr

    try:
        val = json.loads(injson)
    except TypeError:
        pdb.set_trace()

    if not isinstance(val, list):
        return val

    try:
        float(val[0])
        return np.array(val)
    except ValueError:
        return val

def h5setattr(obj, attr, val):
    """Set attribute in JSON format."""

    if isinstance(val, np.ndarray):
        kind = val.dtype.kind
        if kind == 'i':
            val = map(int, val)
        elif kind == 'f':
            val = map(float, val)
        else:
            raise ValueError('Not implemented')

    val = json.dumps(val)
 
    setattr(getattr(obj, 'attrs'), attr, val)  
