#!/usr/bin/env python
# encoding: UTF8

import sys
import numpy as np
import partpix_map
import healpy

def main():
    
    if len(sys.argv) < 3:
        print "Usage: svd_maps.py map1 map2 ..."
        exit(1)
    input_files = sys.argv[1:]
    print "Orthogonalizing maps {0}".format(input_files)
    
    maps = []
    for infile in input_files:
        if infile[-4:] == '.dat':
            pmap = partpix_map.Partpix_Map()
            pmap.read_from_ascii(infile)
            maps.append(pmap)
            # hmap = pmap.to_healpix()
            # maps.append(hmap)
            print "Read in {0}".format(infile)
        elif infile[-4:] == 'fits':
            hmap = healpy.fitsfunc.read_map(infile)
            pmap = partpix_map.Partpix_Map()
            pmap.from_healpix(hmap)
            maps.append(pmap)
            # maps.append(hmap)
            print "Read in {0}".format(infile)
        
    print "Read in the maps"
    
    # # make them all the resolution of the highest one
    # # also, at this point convert to partpix.
    # max_nside = max([healpy.pixelfunc.get_nside(m) for m in maps])
    # for m in range(len(maps)):
    #     if healpy.pixelfunc.get_nside(maps[m]) < max_nside:
    #         # maps[m] = upgrade(maps[m], max_order)
    #         #print "Need to upgrade {0} to order {1}!".format(input_files[m],max_order)
    #         #exit(1)
    #         maps[m] = healpy.pixelfunc.ud_grade(maps[m], max_nside)
    #     pmap = partpix_map.Partpix_Map()
    #     pmap.from_healpix(maps[m])
    #     maps[m] = pmap
        
    
    # this should intersection them all?
    for p1, pmap1 in enumerate(maps):
        map1 = pmap1
        for pmap2 in maps:
            if pmap1 == pmap2:
                continue
            map1, map2 = map1.intersection(pmap2)
        maps[p1] = map1
    
    print "Intersected the maps"
    
    map_matrix = np.zeros((maps[0].npartpix, len(maps)))
    for p, pmap in enumerate(maps):
        map_matrix[:,p] = pmap.partmap
    
    print "About to decompose"
    
    U,s,V = np.linalg.svd(map_matrix, full_matrices=False)
    
    print "Decomposed!"
    
    for p in range(len(s)):
        omap = partpix_map.Partpix_Map()
        omap.pixel_mapping_arraytohigh = maps[0].pixel_mapping_arraytohigh
        omap.partmap = U[:,p]
        omap.order = maps[0].order
        omap.scheme = maps[0].scheme
        omap.update()
        hmap = omap.to_healpix()
        ofname = 'svd_map_' + str(s[p]) + '.fits'
        healpy.fitsfunc.write_map(ofname, hmap)
    
if __name__ == '__main__':
    main()
