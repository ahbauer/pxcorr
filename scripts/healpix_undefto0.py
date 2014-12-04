#!/usr/bin/env python
# encoding: UTF8

import sys
import healpy
import partpix_map

if len(sys.argv) != 3:
    print 'Usage: healpix_undefto0.py input_mapname output_mapname'
    exit(1)

infilename = sys.argv[1]
outfilename = sys.argv[2]

hmap = healpy.fitsfunc.read_map(infilename)
pmap = partpix_map.Partpix_Map()
pmap.from_healpix(hmap)
hmap = pmap.to_healpix()

healpy.fitsfunc.write_map(outfilename, hmap)

