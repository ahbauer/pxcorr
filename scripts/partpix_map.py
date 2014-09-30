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
from scipy import weave
from scipy.weave import converters

class partpix_entry(tables.IsDescription):
    pixel = tables.Int32Col(dflt=0, pos=0)
    value = tables.Float32Col(dflt=0.0, pos=1)


class Partpix_Map:
    def __init__(self, o=None, s=0):
        self.order = o
        self.scheme = s
        self.partmap = None
        self.pixel_mapping_arraytohigh = None
        self.npartpix = 0
        self.npix = 0
    
    def update(self):
        self.npartpix = len(self.partmap)
        nside = int(2.0**self.order)
        self.npix = int(12*nside*nside)
        indices = np.argsort(self.pixel_mapping_arraytohigh)
        self.pixel_mapping_arraytohigh = self.pixel_mapping_arraytohigh[indices]
        self.partmap = self.partmap[indices]
        
    def Npartpix(self):
        return self.npartpix
    
    def Npix(self):
        return self.npix
        
    def Nside(self):
        return int(2.0**self.order)
    
    def __getitem__(self, i):
        indices = np.searchsorted(self.pixel_mapping_arraytohigh,i)
        assert self.pixel_mapping_arraytohigh[indices] == i, "partpix_map:__getitem__: all indices not found in the map"
        return self.partmap[indices]
        
    def mean(self):
        if self.partmap is None:
            err_msg = "Trying to take the mean of an undefined Partpix_Map"
            raise Exception, err_msg
        return np.mean(self.partmap)
        
    def sum(self):
        if self.partmap is None:
            err_msg = "Trying to take the sum of an undefined Partpix_Map"
            raise Exception, err_msg
        return np.sum(self.partmap)
        
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
        while i<self.npartpix-1 and j<map2.npartpix-1:
            while (i<self.npartpix-1) and (self.pixel_mapping_arraytohigh[i] < map2.pixel_mapping_arraytohigh[j]):
                i += 1
            while (j<map2.npartpix-1) and (self.pixel_mapping_arraytohigh[i] > map2.pixel_mapping_arraytohigh[j]):
                j += 1
            while (i<self.npartpix) and (j<map2.npartpix) and (self.pixel_mapping_arraytohigh[i] == map2.pixel_mapping_arraytohigh[j]):
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
        
        outmap1.update()
        outmap2.update()
        
        return outmap1, outmap2
        
    def bool( self, cutoff=0.5 ):
        if self.npartpix == 0:
            return
        for i in range(self.npartpix):
            if self.partmap[i] >= cutoff:
                self.partmap[i] = 1.0
            else:
                self.partmap[i] = 0.0
    
    def read_from_ascii( self, filename ):
        file = open(filename, 'r')
        filelines = file.readlines()
        (s, o) = filelines.pop(0).split()
        self.order = int(o)
        nside = int(2.0**self.order)
        self.npix = int(12*nside*nside)
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
        self.update()
        
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
        self.update()
    
    def write_to_pytable( self, table ):
        
        assert isinstance(table, tables.table.Table), "This is not an hdf5 table!"
        table.setAttr('order', self.order)
        table.setAttr('scheme', self.scheme)
        row = table.row
        for i in range(self.npartpix):
            row['pixel']  = self.pixel_mapping_arraytohigh[i]
            row['value'] = self.partmap[i]
            row.append()
        table.flush()
    
    def to_healpix(self):
        healpix_map = np.zeros(self.npix)
        healpix_map[self.pixel_mapping_arraytohigh] = self.partmap
        return healpix_map
    
    def from_healpix(self, hmap ):
        nside = np.sqrt(len(hmap)/12.)
        nside = int(np.floor(nside))
        self.order = np.log2(nside)
        self.pixel_mapping_arraytohigh = np.arange(len(hmap))[(hmap != 0.) & (hmap > -1.e10)] # [i for i in range(len(hmap)) if (hmap[i] != 0. and hmap[i] > -1.e10)]
        self.partmap = hmap[(hmap != 0.) & (hmap > -1.e10)]
        self.update()
        
def read_from_hdf5file( map_filename, read_map=True, read_mask=True ):
    
    mapfile = tables.openFile(map_filename)
    map1 = None
    mask1 = None
    if read_mask:
        mask1 = Partpix_Map()
        mask_pixels = mapfile.getNode('/data/', 'mask_pixels')
        mask1.pixel_mapping_arraytohigh = mask_pixels.read() # a numpy array
        # this is sloppy;  the order & scheme should be metadata, not the first values in the array.
        mask1.order = mask1.pixel_mapping_arraytohigh[0]
        mask1.pixel_mapping_arraytohigh = np.delete(mask1.pixel_mapping_arraytohigh, 0) 
        mask_data = mapfile.getNode('/data/', 'mask_values')
        mask1.partmap = mask_data.read()
        mask1.scheme = mask1.partmap[0]
        mask1.partmap = np.delete(mask1.partmap, 0)
        mask1.update()
    if read_map:
        map1 = Partpix_Map()
        map_pixels = mapfile.getNode('/data/', 'pixels')
        map1.pixel_mapping_arraytohigh = map_pixels.read() # a numpy array
        map1.order = map1.pixel_mapping_arraytohigh[0]
        map1.pixel_mapping_arraytohigh = np.delete(map1.pixel_mapping_arraytohigh, 0)
        map_data = mapfile.getNode('/data/', 'values')
        map1.partmap = map_data.read()
        map1.scheme = map1.partmap[0]
        map1.partmap = np.delete(map1.partmap, 0)
        map1.update()
    mapfile.close()
    return map1, mask1


code_degrade_nests="""
using namespace std;

short ctab[0x100];
short utab[0x100];

for (int m=0; m<0x100; ++m){
  ctab[m] = short(
       (m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
    | ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4));
  utab[m] = short(
       (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
    | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7));
}
  

int num_pixels = n_p;
for( unsigned int ppp=0; ppp<num_pixels; ppp++ ){
    
    int pix = pixs[ppp];
    int ix, iy, face_num;
    int npface_;
    
    int nside_in  = 1<<order_in;
    int nside_out  = 1<<order_out;
    npface_ = nside_in*nside_in;
    
  face_num = pix>>(2*order_in);
  pix &= (npface_-1);
  int raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  ix = ctab[raw&0xff] | (ctab[raw>>8]<<4);
  pix >>= 1;
  raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  iy = ctab[raw&0xff] | (ctab[raw>>8]<<4);


  // degrade
  float fact=nside_in/nside_out;
  ix = int(ix/fact);
  iy = int(iy/fact);
  
  
  pixs[ppp] = (face_num<<(2*order_out)) +
      (utab[ix&0xff] | (utab[ix>>8]<<16)
    | (utab[iy&0xff]<<1) | (utab[iy>>8]<<17));

}

"""

code_degrade_rings = """
using namespace std;

const int jrll[12] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
const int jpll[12] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

int num_pixels = n_p;
for( unsigned int ppp=0; ppp<num_pixels; ppp++ ){
    
    int pix = pixs[ppp];
    int ix, iy, face_num;
    int nside_in;
    {
    int iring, iphi, kshift, nr;
    int nside_;
    int npface_, ncap_, npix_;
    double fact1_, fact2_;
        
    nside_  = 1<<order_in;
    npface_ = nside_*nside_;
    ncap_   = (npface_-nside_)<<1;
    npix_   = 12*npface_;
    fact2_  = 4./npix_;
    fact1_  = (nside_<<1)*fact2_;;
    int nl2 = 2*nside_;
    
    nside_in = nside_;
    
if (pix<ncap_) // North Polar cap
  {
  iring = int(0.5*(1+int(sqrt(1+2*pix+0.5)))); //counted from North pole
  iphi  = (pix+1) - 2*iring*(iring-1);
  kshift = 0;
  nr = iring;
  face_num=0;
  int tmp = iphi-1;
  if (tmp>=(2*iring))
    {
    face_num=2;
    tmp-=2*iring;
    }
  if (tmp>=iring) ++face_num;
  }
else if (pix<(npix_-ncap_)) // Equatorial region
  {
  int ip = pix - ncap_;
  if (order_in>=0)
    {
    iring = (ip>>(order_in+2)) + nside_; // counted from North pole
    iphi  = (ip&(4*nside_-1)) + 1;
    }
  else
    {
    iring = (ip/(4*nside_)) + nside_; // counted from North pole
    iphi  = (ip%(4*nside_)) + 1;
    }
  kshift = (iring+nside_)&1;
  nr = nside_;
  unsigned int ire = iring-nside_+1;
  unsigned int irm = nl2+2-ire;
  int ifm, ifp;
  if (order_in>=0)
    {
    ifm = (iphi - ire/2 + nside_ -1) >> order_in;
    ifp = (iphi - irm/2 + nside_ -1) >> order_in;
    }
  else
    {
    ifm = (iphi - ire/2 + nside_ -1) / nside_;
    ifp = (iphi - irm/2 + nside_ -1) / nside_;
    }
  if (ifp == ifm) // faces 4 to 7
    face_num = (ifp==4) ? 4 : ifp+4;
  else if (ifp<ifm) // (half-)faces 0 to 3
    face_num = ifp;
  else // (half-)faces 8 to 11
    face_num = ifm + 8;
  }
else // South Polar cap
  {
  int ip = npix_ - pix;
  iring = int(0.5*(1+int(sqrt(2*ip-1+0.5)))); //counted from South pole
  iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
  kshift = 0;
  nr = iring;
  iring = 2*nl2-iring;
  face_num=8;
  int tmp = iphi-1;
  if (tmp>=(2*nr))
    {
    face_num=10;
    tmp-=2*nr;
    }
  if (tmp>=nr) ++face_num;
  }

    int irt = iring - (jrll[face_num]*nside_) + 1;
    int ipt = 2*iphi- jpll[face_num]*nr - kshift -1;
    if (ipt>=nl2) ipt-=8*nside_;

    ix =  (ipt-irt) >>1;
    iy =(-(ipt+irt))>>1;

    }
    // end ring_to_xyf, start xyf_to_ring
    {
    int nside_;
    int npface_, ncap_, npix_;
    nside_  = 1<<order_out;
    
    npface_ = nside_*nside_;
    ncap_   = (npface_-nside_)<<1;
    npix_   = 12*npface_;
    int nl4 = 4*nside_;
    
    // degrade!!
    float fact=nside_in/nside_;
    ix = int(ix/fact);
    iy = int(iy/fact);
    
    int jr = (jrll[face_num]*nside_) - ix - iy  - 1;
    
    int nr, kshift, n_before;
    if (jr<nside_)
      {
      nr = jr;
      n_before = 2*nr*(nr-1);
      kshift = 0;
      }
    else if (jr > 3*nside_)
      {
      nr = nl4-jr;
      n_before = npix_ - 2*(nr+1)*nr;
      kshift = 0;
      }
    else
      {
      nr = nside_;
      n_before = ncap_ + (jr-nside_)*nl4;
      kshift = (jr-nside_)&1;
      }

    int jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
    if (jp>nl4)
      jp-=nl4;
    else
      if (jp<1) jp+=nl4;

    pixs[ppp] = n_before + jp - 1;
    }
    }

"""

def degrade( map_in, order_out ):
    
    pixs = list(map_in.pixel_mapping_arraytohigh.copy())
    order_in = int(map_in.order)
    order_out = int(order_out)
    n_p = len(pixs)
    new_pixs = list(np.zeros(n_p))
    args = ['pixs', 'n_p', 'order_in', 'order_out'],
    # print "before: %d %d %d %d" %(len(pixs), pixs[0], pixs[1], pixs[-1])
    if map_in.scheme == 0:
        weave.inline(code_degrade_rings, *args, type_converters=converters.blitz)
    else:
        weave.inline(code_degrade_nests, *args, type_converters=converters.blitz)
    # print "after: %d %d %d %d" %(len(pixs), pixs[0], pixs[1], pixs[-1])
    lrpixs = pixs
    
    map_out = Partpix_Map()
    unique_pixels, unique_indices = np.unique(lrpixs, return_inverse=True)
    nu = len(unique_pixels)
    map_out.pixel_mapping_arraytohigh = unique_pixels
    map_out.partmap = np.zeros(nu)
    for p in range(len(lrpixs)):
        map_out.partmap[unique_indices[p]] += map_in.partmap[p]
        
    map_out.order = order_out
    map_out.update()
    fact = map_in.Nside()/map_out.Nside()
    map_out.partmap /= (fact*fact)
    
    # hpmap = map_out.to_healpix()
    # healpy.fitsfunc.write_map("temp.fits",hpmap)
     
    return map_out

# def upgrade( map_in, order_out ):
#     
#     assert map_in.scheme == 0, "partpix_mapupgrade error: only implemented for RING"
#     
#     weave.inline(code_upgrade_rings, *args, type_converters=converters.blitz)
#     
#     fact = int(2.0**order_out)/map_in.Nside();
#     for m0 in range(map_in.Npartpix()):
#         m = map_in.pixel_mapping_arraytohigh[m0]
#         x,y,f = ring_to_xyf(m)
#         for j in range(fact*y,fact*(y+1),1):
#             for i in range(fact*x,fact*(x+1),1):
#                 p = xyf_to_ring(i,j,f)
#                 map_out[p] = map_in.partmap[m0]
      
      