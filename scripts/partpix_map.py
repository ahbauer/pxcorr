#!/usr/bin/env python
# encoding: UTF8

"""

Anne Bauer
01/2014

This is THE BEGINNING of a partpix python class.

It only has A FEW functions and is not equivalent to the C++ version, but is a bit 
useful for dealing with partpix maps and I will add to it as necessary.

"""

import sys
import json
import os
import numpy as np
import tables

class Partpix_Map:
    def __init__(self, o=None, s=0):
        self.order = o
        self.scheme = s
        self.partmap = None
        self.pixel_mapping_arraytohigh = None
        self.npartpix = 0
        
    def Npartpix(self):
        return self.npartpix
        
    def mean(self):
        if self.partmap is None:
            err_msg = "Trying to take the mean of an undefined Partpix_Map"
            raise Exception, err_msg
        return np.mean(self.partmap)
        
    def intersection( self, map2 ):
        
        if self.order != map2.order:
            err_msg = "Trying to intersect two maps with different orders: {0} and {1}".format(self.order, map2.order)
            print err_msg
            raise AssertionError, err_msg
        if self.scheme != map2.scheme:
            err_msg = "Trying to intersect two maps with different schemes: {0} and {1}".format(self.scheme, map2.scheme)
            print err_msg
            raise AssertionError, err_msg
        
        i = 0
        j = 0
        partmap1 = []
        partmap2 = []
        pix = []
        while i<self.npartpix and j<map2.npartpix:
            while i<self.npartpix and self.pixel_mapping_arraytohigh[i] < map2.pixel_mapping_arraytohigh[j]:
                i += 1
            while j<map2.npartpix and self.pixel_mapping_arraytohigh[i] > map2.pixel_mapping_arraytohigh[j]:
                j += 1
            while i<self.npartpix and j<map2.npartpix and self.pixel_mapping_arraytohigh[i] == map2.pixel_mapping_arraytohigh[j]:
                partmap1.append(self.partmap[i])
                partmap2.append(map2.partmap[j])
                pix.append(self.pixel_mapping_arraytohigh[i])
                i += 1
                j += 1
        outmap1 = Partpix_Map( self.order, self.scheme )
        outmap2 = Partpix_Map( self.order, self.scheme )
        outmap1.pixel_mapping_arraytohigh = np.array(pix)
        outmap2.pixel_mapping_arraytohigh = np.array(pix)
        outmap1.partmap = np.array(partmap1)
        outmap2.partmap = np.array(partmap2)
        outmap1.npartpix = len(partmap1)
        outmap2.npartpix = len(partmap2)
        
        return outmap1, outmap2
        
        
    def read_from_ascii( self, filename ):
        file = open(filename, 'r')
        filelines = file.readlines()
        (o, s) = filelines.pop(0).split()
        self.order = int(o)
        if s == 'RING' or s == 0:
            self.scheme = 0
        elif s == 'NEST' or s == 1:
            self.scheme = 1
        else:
            err_msg = "Problem understanding scheme {0}".format(s)
            raise Exception, err_msg
        self.npartpix = len(filelines)
        self.pixel_mapping_arraytohigh = np.zeros(self.npartpix, dtype='int64')
        self.partmap = np.zeros(self.npartpix, dtype='float')
        c = 0
        for line in filelines:
            (p, v) = line.split()
            self.pixel_mapping_arraytohigh[c] = int(p)
            self.partmap[c] = float(v)
            c += 1
        
        
    def read_from_hdf5table( self, table ):
        
        assert isinstance(table, tables.table.Table), "This is not an hdf5 table!"
        colset = set(table.cols._v_colnames)
        assert set(['pixel', 'value']).issubset(colset), msg_cols
        
        ps = table.read(field='pixel')
        vs = table.read(field='value')
        
        if len(ps) != len(vs):
            print "Error, hdf5 partpix map table has a different number of pixels from values (%d vs. %d)" %(len(ps), len(vs))
            
        obj_attrs = getattr(table, 'attrs')
        self.order = getattr(obj_attrs, 'order')
        self.scheme = getattr(obj_attrs, 'scheme')
        self.npartpix = len(ps)
        
        self.pixel_mapping_arraytohigh = np.zeros(self.npartpix, dtype='int64')
        self.partmap = np.zeros(self.npartpix, dtype='float')
        for i in range(self.npartpix):
            self.pixel_mapping_arraytohigh[i] = ps[i]
            self.partmap[i] = vs[i]
    
    
        