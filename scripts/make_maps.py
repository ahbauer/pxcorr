# to replace the C++ make_maps
# since this whole swig thing kind of sucks and is likely unnecessary

import os
import sys
import math
from array import array
import numpy as np
import tables
import healpy
from partpix_map import Partpix_Map
from partpix_map import write_to_hdf5file
from partpix_map import degrade

def choose_order(ang_widths):
    min_footprint_order = 3
    max_order = 18
    pixel_multiple = 2.5
    
    min_width = np.min(ang_widths)
    orders = np.arange(min_footprint_order,max_order+1,1)
    pixel_sizes = np.sqrt(41253./(12.0*(2.0**orders)**2.0))
    good_orders = orders[pixel_multiple*pixel_sizes < min_width]
    if len(good_orders) == 0:
        msg = 'make_maps.py:choose_order: Trying to use an order > max_order={0}'.format(max_order)
        raise Exception(msg)
    order = good_orders[0]
    print 'Using order {0} for minimum angular bin width {1}'.format(order, min_width)
    return order


def construct_footprint(input_map, min_footprint_order=3):
    # degrade the mask to each res and count non-zero pixels
    footprint_order = min_footprint_order;
    footprint_area = 50000.
    footprint_map = None
    while footprint_order < order:
        footprint_map = input_partpix.degrade(input_partpix, footprint_order)
        footprint_area_new = np.sum(footprint_map.partmap>0.)*(41253./footprint_map.Npix())
        if footprint_area_new <= 0.8*footprint_area:
            footprint_area = footprint_area_new
            print 'Footprint order {0} uses area {1} deg2... trying next order'.format(footprint_order, footprint_area)
            footprint_order += 1
        else:
            break
    # the previous order was good
    footprint_order -= 1;
    footprint_map = input_partpix.degrade(input_partpix, footprint_order)
    footprint_map.partmap[footprint_map.partmap!=0.] = 1
    print 'Using footprint order {0} with {1} deg2'.format(footprint_order, footprint_area)
    return footprint_map


def check_filesize(mask):
    max_mapsize = 0.8 # GB
    n_maskpix = mask.sum()
    n_highrespix = np.sum(mask.partmap[mask.partmap>0.5])
    expected_mapsize = ( (8 + 4)*n_highrespix + (8 + 1)*n_maskpix )/1000000.
    print 'Your map will be about {0} MB'.format(expected_mapsize)
    if expected_mapsize > max_mapsize*1000:
        msg = 'Your map will be too big: {0} MB.  Please run with a larger angular bin size or a smaller area.  Sorry... blame PIC!'.format(expected_mapsize)
        raise Exception(msg)
    return


def make_maps_2pt(subcat_filename, mask_filename, ang_means, ang_widths, use_counts, use_mags, suffix, pop, zbin):
    
    if (not use_counts) and (not use_mags):
        msg = 'Not using either counts or mags... doing nothing!'
        raise Exception(msg)
    
    # from the angular widths, figure out the resolution
    order = choose_order(ang_widths)
    nside = 2**order
    
    # is the input a fits file?
    if subcat_filename.endswith(".fits") or subcat_filename.endswith(".dat"):
        
        if use_mags:
            msg = 'Error: The input is a map, but you want to use magnitudes too?  What does that mean?'
            raise Exception(msg)
        
        input_partpix = Partpix_Map()
        if subcat_filename.endswith(".fits"):        
            # read it in
            input_map = healpy.fitsfunc.read_map(subcat_filename)
            # convert to partpix
            input_partpix.from_healpix(input_map)
        else:
            input_partpix.read_from_ascii(subcat_filename)
        
        # upgrade the input map if necessary to the min res
        if input_partpix.order < order:
            msg = 'make_maps.py: Need to implement upgrade in partpix_map.py!'
            raise Exception(msg)
        
        input_nside = healpy.pixelfunc.get_nside(input_map)
        input_order = np.log2(input_nside)
        
        # make a mask and footprint if there isn't a mask
        dcmask = Partpix_Map()
        if( mask_filename.lower() in ['none', 'null'] ):
            
            # Construct the mask from the input map
            dcmask.copy(input_partpix)
            dcmask.partmap[dcmask.partmap <-1.e30] = 0
            dcmask.partmap[np.isnan(dcmask.partmap)] = 0
            dcmask.partmap[dcmask.partmap != 0.] = 1
            print 'Finished filling in the mask from the Healpix input itself'
        
        # if there is an input mask, read it in
        else:
            dcmask.read_from_ascii(mask_filename)
            
            # upgrade the mask if necessary
            if dcmask.order < order:
                msg = 'make_maps.py: Need to implement upgrade in partpix_map.py!'
                raise Exception(msg)
            
        # check the file size
        check_filesize(dcmask)
        
        # write it out
        write_to_hdf5file(outfilename, input_partpix, dcmask, ang_means, ang_widths, use_counts, use_mags, suffix, 'counts', pop, zbin, 0)
        
    # else, it is an ascii galaxy catalog.
    else:
        # make a list of pointings, mags, and weights.
        input_all = np.fromfile(subcat_filename)
        phis = math.pi/180. * input_all[0::4]  # ra = mypointing.phi*180./3.1415926;
        thetas = math.pi/180. * (90.0-input_all[1::4])  # dec = 90. - mypointing.theta*180./3.1415926;
        mags = input_all[2::4]
        weights = input_all[3::4]
        print 'Read in {0} objects'.format(len(phis))
        
        # read in the mask.  if there isn't one, make it from the catalog but warn.  (same resolution as the map.)
        mask = Partpix_Map()
        if( mask_filename.lower() in ['none', 'null'] ):
            print 'make_maps.py WARNING: constructing the mask from the catalog;  this is usually a bad idea.'
            mask.order = order
            mask_pix = np.zeros(len(phis),dtype='int')
            for i, pix in zip(range(len(thetas)), healpy.pixelfunc.ang2pix(nside,thetas,phis)):
                mask_pix[i] = pix
            mask_pix = np.unique(mask_pix)
            mask.pixel_mapping_arraytohigh = mask_pix
            mask.partmap = np.ones(len(mask_pix),dtype='int')
            mask.update()
        else:
            mask.read_from_ascii(mask_filename)
            # upgrade the mask if necessary
            if mask.order < order:
                msg = 'make_maps.py: Need to implement upgrade in partpix_map.py!'
                raise Exception(msg)
            # downgrade the mask if necessary!
            elif mask.order > order:
                mask = degrade( mask, order )
        
        # make an identical mag mask if desired.  (same resolution as the map.)
        mag_mask = Partpix_Map()
        if use_mags:
            mag_mask.copy(mask)
        
        # check file size
        check_filesize(mask)
        
        # make the counts map
        # also make mag maps if desired
        dcmap = Partpix_Map()
        dcmap.copy(mask)
        dcmap.partmap = np.zeros(len(dcmap.partmap))
        dmmap = Partpix_Map()
        dmmask = Partpix_Map()
        if use_mags:
            dmmap.copy(mask)
            dmmap.partmap = np.zeros(len(dmmap.partmap))
            dmmask.copy(mask)
        
        pixs = healpy.pixelfunc.ang2pix(nside,thetas,phis)
        # sort by pixel order, to match the mask ordering
        indices = np.argsort(pixs)
        pixs = pixs[indices]
        mags = mags[indices]
        weights = weights[indices]
        # this is very c-ish!  but my pythonic go-to was way too slow...
        i=0
        j=0
        nobjs = len(pixs)
        npx = dcmap.Npartpix()
        while(1):
            if pixs[i] == dcmap.pixel_mapping_arraytohigh[j]:
                dcmap.partmap[j] += weights[i]
                if use_mags:
                    dmmap.partmap[j] += mags[i]*weights[i]
                i += 1
            elif pixs[i] < dcmap.pixel_mapping_arraytohigh[j]:
                i += 1
            else:
                j += 1
            if (i == nobjs) or (j == npx):
                break
                
        # for pix, mag, weight in zip(pixs, mags, weights):
        #     if pix in mask.pixel_mapping_arraytohigh:
        #         dcmap[pix] += weight
        #         if use_mags:
        #             dmmap[pix] += mag*weight
        
        # make the mag -> deta_mag map
        # only valid for pixels with counts in them
        if use_mags:
            print 'Normalizing the magnitude map by the mean {0}'.format(np.mean(mags))
            dmmap.partmap[dcmap.partmap>0] = dmmap.partmap[dcmap.partmap>0]/dcmap.partmap[dcmap.partmap>0] - np.mean(mags)
            dmmask.partmap[dcmap.partmap==0] = 0
        
        # convert the counts map into a delta map
        nobjs = np.sum(dcmap.partmap)
        if nobjs == 0:
            msg = 'make_maps: number of objects in the dc_map = 0!?'
            raise Exception(msg)
        dcmap.partmap = dcmap.partmap/np.mean(dcmap.partmap) - 1.0
        
        # write out the maps
        print 'Writing out the map(s)'
        write_to_hdf5file('dc_map' + suffix + '.h5', dcmap, mask, ang_means, ang_widths, 'counts', pop, zbin, nobjs)
        
        if use_mags:
            write_to_hdf5file('dm_map' + suffix + '.h5', dmmap, dmmask, ang_means, ang_widths, 'mags', pop, zbin, nobjs)
        
    return

