import os
import sys
import ast
import itertools
import math
from array import array
import numpy as np
from scipy import weave
from scipy.weave import converters
import tables
import healpy
import partpix_map

def read_map( mapname ): # written in make_maps.py
    
    mapfile = tables.open_file(mapname)
    
    # read in the map data
    mask_pixels = None
    try:
        mask_pixels = mapfile.root.data.mask_pixels.read()
    except tables.exceptions.NoSuchNodeError:
        msg = 'Error, data/mask_pixels not in the map file {0}'.format(mapname)
        raise Exception(msg)
    mask_values = None
    try:
        mask_values = mapfile.root.data.mask_values.read()
    except tables.exceptions.NoSuchNodeError:
        msg = 'Error, data/mask_values not in the map file {0}'.format(mapname)
        raise Exception(msg)
    pixels = None
    try:
        pixels = mapfile.root.data.pixels.read()
    except tables.exceptions.NoSuchNodeError:
        msg = 'Error, data/pixels not in the map file {0}'.format(mapname)
        raise Exception(msg)
    values = None
    try:
        values = mapfile.root.data.values.read()
    except tables.exceptions.NoSuchNodeError:
        msg = 'Error, data/values not in the map file {0}'.format(mapname)
        raise Exception(msg)
    
    # read in the metadata
    attrs = None
    try:
        attrs = mapfile.root.meta.meta.attrs
    except AttributeError:
        msg = 'Error, meta/meta does not have attributes in the map file {0}'.format(mapname)
        raise Exception(msg)
    
    ftype = attrs['ftype']
    pop = attrs['pop']
    nobj = ast.literal_eval(attrs['nobj']) # parse to numbers
    zbin = ast.literal_eval(attrs['zbin'])
    u_mean = ast.literal_eval(attrs['u_mean'])
    u_width = ast.literal_eval(attrs['u_width'])
    order = ast.literal_eval(attrs['order'])
    scheme = ast.literal_eval(attrs['scheme'])
    
    mapfile.close()
    
    dmap = partpix_map.Partpix_Map(order, scheme)
    dmap.pixel_mapping_arraytohigh = pixels
    dmap.partmap = values
    dmap.update()
    
    dmask = partpix_map.Partpix_Map(order, scheme)
    dmask.pixel_mapping_arraytohigh = mask_pixels
    dmask.partmap = mask_values
    dmask.update()
    
    # only keep the mask pixels that are good
    mask_indices = (dmask.partmap>0.5)
    dmask.pixel_mapping_arraytohigh = dmask.pixel_mapping_arraytohigh[mask_indices]
    dmask.partmap = dmask.partmap[mask_indices]
    dmask.update()
    
    # only return the part of the map that lies within the mask
    dmap2, dmask2 = dmap.intersection(dmask)
    
    return dmap2, dmask, ftype, pop, nobj, zbin, u_mean, u_width



def correlate( mapname1, mapname2, r, jk_order, degrade_maps, outfilename, exclude_0dist = True ):
    
    pixel_multiple = 2.5
    jk_threshold_frac = 0.5
    
    # read in the two maps into partpix maps
    (map1, mask1, ftype1, pop1, nobj1, zbin1, u_mean1, u_width1) = read_map(mapname1)
    (map2, mask2, ftype2, pop2, nobj2, zbin2, u_mean2, u_width2) = read_map(mapname2)
    print 'Read in the maps, orders {0} & {1}'.format(map1.order, map2.order)
    print 'Using angular bin {0}: {1} degrees'.format(r, u_mean1[r])
    
    # map orders must be the same!
    assert map1.order == map2.order, 'correlate: maps must have the same order!  {0} != {1}'.format(map1.order, map2.order)
    
    # make the jackknife map at the given resolution
    # using the jk_threshold if the jackknife map is lower res than the mask
    jkmap = partpix_map.Partpix_Map(jk_order, 0)
    if jk_order > map1.order:  # weird
        msg = 'Error, your jackknife resolution is higher than your data resolution!'
        raise Exception(msg)
    elif jk_order == map1.order:
        jk_map, jk_map2 = mask1.intersection(mask2)
    else: # jackknife map lower resolution than the data
        mask1d = partpix_map.degrade(mask1, jk_order)
        mask2d = partpix_map.degrade(mask2, jk_order)
        mask1i, mask2i = mask1d.intersection(mask2d)
        mask1i.partmap[mask1i.partmap<jk_threshold_frac] = 0
        mask1i.partmap[mask2i.partmap<jk_threshold_frac] = 0
        jkmap.copy(mask1i)
    jk_nside = jkmap.Nside()
    # condense the jk map to only the good pixels
    indices = jkmap.partmap>=jk_threshold_frac
    jkmap.partmap = jkmap.partmap[indices]
    jkmap.pixel_mapping_arraytohigh = jkmap.pixel_mapping_arraytohigh[indices]
    jkmap.update()
    njk = int(jkmap.Npartpix())
    print 'Constructed a jackknife map of order {0}, {1} good pixels'.format(jk_order, njk)
    
    # make jk_index_array to save time in the loops
    # takes some space, but the jkmap should be pretty low resolution
    jkindex_array = njk*np.ones(jkmap.Npix(),dtype='int')
    for i in range(jkmap.Npartpix()):
        pix = jkmap.pixel_mapping_arraytohigh[i]
        jkindex_array[pix] = i
        
    # figure out the mask area (unused)
    area_mask1 = 41252.962*np.sum(mask1.partmap>0.5)/mask1.Npix();
    area_mask2 = 41252.962*np.sum(mask2.partmap>0.5)/mask2.Npix();
    area_jk = 41252.962*np.sum(jkmap.partmap>=jk_threshold_frac)/jkmap.Npix();
    # area_fraction = area_mask/area_jk;
    print 'Mask areas: {0:4.1f} & {1:4.1f} deg2, Jackknife area {2:4.1f} deg2'.format(area_mask1, area_mask2, area_jk)
    
    # degrade the maps if we can, for the resolution of the angular bin
    if degrade_maps:
        order = jk_order
        while(1):
            pixel_size = np.sqrt(41253./(12.0*((2.0**order)**2.0)))
            if pixel_multiple*pixel_size < u_width1[r]:
                break
            order += 1
        if order < map1.order:
            print 'Degrading the maps from order {0} to {1}'.format(map1.order, order)
            map1 = partpix_map.degrade(map1, order)
            mask1 = partpix_map.degrade(mask1, order)
            map2 = partpix_map.degrade(map2, order)
            mask2 = partpix_map.degrade(mask2, order)
    
    sumij = np.zeros(njk+1)
    nij = np.zeros(njk+1)
    sumij_injk = np.zeros( (njk+1, njk+1) )
    nij_injk = np.zeros( (njk+1, njk+1) )
    
    # for pixel in map1
    nside1 = map1.Nside()
    nside2 = map2.Nside()
    print 'Beginning correlation'
    vxs, vys, vzs = healpy.pix2vec(nside1, map1.pixel_mapping_arraytohigh)
    jki1s = jkindex_array[healpy.vec2pix(jk_nside,vxs, vys, vzs)]
    u_high = (u_mean1[r]+0.5*u_width1[r])*math.pi/180.
    u_low = (u_mean1[r]-0.5*u_width1[r])*math.pi/180.
    
    for val1, p1, vx, vy, vz, jkpix_index1 in zip(map1.partmap,map1.pixel_mapping_arraytohigh,vxs, vys, vzs,jki1s):
        
        # get annulus of pixels around p
        p2s_low = healpy.query_disc(nside1, (vx, vy, vz), u_low, inclusive=False)
        p2s_high = healpy.query_disc(nside1, (vx, vy, vz), u_high, inclusive=False)
        annulus = np.setdiff1d(p2s_high, p2s_low, assume_unique=True)
        
        # take only the pixels in map 2
        nobjs = len(annulus)
        npx = mask2.Npartpix()
        out_array = np.zeros(nobjs,dtype='int')
        
        code_intersect = """
        int i=0;
        int j=0;
        int k=0;
        while(1){
            if( annulus[i] == pmapping[j] ){
                out_array[k] = annulus[i];
                ++k;
                ++j;
                ++i;
            }
            else if( annulus[i] < pmapping[j] ){
                ++i;
            }
            else{
                ++j;
            }
            if( (i == nobjs) || (j == npx) ){
                break;
            }
        }
        return_val = k;
        """
        pmapping = mask2.pixel_mapping_arraytohigh
        k = weave.inline(code_intersect, ['annulus','pmapping','out_array','nobjs','npx'])
        annulus = out_array[:k]
        # intersect1d is slow:
        # annulus = np.intersect1d(annulus, mask2.pixel_mapping_arraytohigh, assume_unique=True)
        
        # for pixel j in annulus
        theta2s, phi2s = healpy.pix2ang(nside2, annulus)
        jki2s = jkindex_array[healpy.ang2pix(jk_nside,theta2s,phi2s)]
        
        # do not include the same pixel in the two maps?
        if exclude_0dist:
            indices = (annulus != p1)
            annulus = annulus[indices]
            jki2s = jki2s[indices]
        
        map2_vals = map2[annulus]
        npix2 = len(annulus)
        code_add = """
        for(int i=0; i<npix2; ++i){
            int pos = int(int(jkpix_index1)*(njk+1) + jki2s[i]);
            sumij_injk[pos] += float(val1)*float(map2_vals[i]);
            nij_injk[pos] += 1.0;
        }
        """
        weave.inline(code_add, ['sumij_injk','nij_injk','njk','npix2','jkpix_index1','jki2s','map2_vals','val1'])
        
        # for p2 
            
    # for p1
    
    # rearrange jackknife info to exclude the pixel
    print 'Compiling the results'
    
    code_jk = """
    for(int i = 0; i < njk+1; i++){
        for(int j = 0; j < njk+1; j++){
            int pos = i*(njk+1) + j;
            for(int k = 0; k < njk; k++){
                if( (k != i) && (k != j) ){
                    sumij[k] += sumij_injk[pos];
                    nij[k] += nij_injk[pos];
                }
            }
            sumij[njk] += sumij_injk[pos];
            nij[njk] += nij_injk[pos];
        }
    }
    """
    weave.inline(code_jk, ['sumij','nij','sumij_injk','nij_injk','njk'])
    
    if not np.all(nij>0):
        empties = np.sum(nij==0)
        print 'correlate: {0} empty JK pixels!'.format(empties)
    indices = nij>0
    corr_wo1 = np.zeros(njk+1)
    corr_wo1[indices] = sumij[indices]/nij[indices]
    correlation = corr_wo1[-1]
    jkmap.partmap = corr_wo1[:njk]
    
    # write out the hdf5 file
    print 'Writing output'
    outfile = tables.openFile(outfilename, 'w')
    outfile.createGroup('/', 'meta', 'Metadata')
    meta_array = outfile.createArray('/meta', 'meta', 'Dummy')
    meta_array.setAttr('fourier', 'false')
    meta_array.setAttr('ftype0', ftype1)
    meta_array.setAttr('pop0', pop1)
    meta_array.setAttr('zbin0', str(zbin1))
    meta_array.setAttr('nobj0', str(nobj1))
    meta_array.setAttr('ftype1', ftype2)
    meta_array.setAttr('pop1', pop2)
    meta_array.setAttr('zbin1', str(zbin2))
    meta_array.setAttr('nobj1', str(nobj2))
    
    corr_group = outfile.createGroup('/', 'corr', 'corr')
    corr0 = outfile.createArray('/corr', 'corr0', [correlation], 'correlation')
    
    jk_group = outfile.createGroup('/', 'JK', 'JK')
    jk_table = outfile.create_table(jk_group, 'JK0', partpix_map.partpix_entry, 'JK0')
    jk_table.setAttr('order', str(jk_order))
    jk_table.setAttr('scheme', '0')
    jk_row = jk_table.row
    for k in range(jkmap.Npartpix()):
        jk_row['pixel'] = jkmap.pixel_mapping_arraytohigh[k]
        jk_row['value'] = jkmap.partmap[k]
        jk_row.append()
    jk_table.flush()
    
    outfile.close()
    

def main():
    
    if (len(sys.argv) != 3) and (len(sys.argv) != 5):
        print 'Usage: correlate.py map1 map2 (radial_bin=0 jk_order=7)'
        print '       where the maps are hdf5 files produced by make_maps'
        exit(1)
    
    r = 0
    jk_order = 7
    if len(sys.argv) == 5:
        r = int(sys.argv[3])
        jk_order = int(sys.argv[4])
    
    print 'Calling correlate'
    correlate( sys.argv[1], sys.argv[2], r, jk_order, True, 'correlation.h5' )
    

if __name__ == '__main__':
    main()

    