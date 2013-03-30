#!/usr/bin/env python
# encoding: UTF8

import sys
import os
import math
import numpy
import tables
import json
from scipy.stats import gaussian_kde

class photoz_entry(tables.IsDescription):
    z_s = tables.Float32Col(dflt=0.0, pos=0) 
    z_p = tables.Float32Col(dflt=0.0, pos=1)

class slopes_entry(tables.IsDescription):
    z = tables.Float32Col(dflt=0.0, pos=0) 
    counts = tables.Float32Col(dflt=0.0, pos=1)
    mag = tables.Float32Col(dflt=0.0, pos=1)
    size = tables.Float32Col(dflt=0.0, pos=1)

def parse_data( filename, mag_cuts, f ):

    # read the needed metadata
    meta = f.getNode('/meta', 'meta')
    injson = getattr(getattr(meta, 'attrs'), 'pop')
    pop = (json.loads(injson))[0]
    #print "pop: %s" %pop
    injson = getattr(getattr(meta, 'attrs'), 'ang_mean')
    ang_means = json.loads(injson)
    #print "ang_means length %d:" %len(ang_means)
    #print ang_means
    metapopname = '/meta/' + pop
    metapop = f.getNode(metapopname, 'meta')
    injson = getattr(getattr(metapop, 'attrs'), 'z_mean')
    z_means = json.loads(injson)
    injson = getattr(getattr(metapop, 'attrs'), 'z_width')
    z_widths = json.loads(injson)
    #print "z_means length %d:" %len(z_means)
    #print z_means
    
    # arrays for output later
    z_phot = []
    z_spec = []

    # what kind of binning do we want for the slopes?
    mag_nzs = 100
    mag_maxz = 2.0

    # make a matrix with sides z_phot and mag.
    # the noise will be determined from the counts vs z 
    nbins_z = len(z_means)
    min_z = 0.0
    max_z = z_means[nbins_z-1] + z_widths[nbins_z-1]
    data = numpy.zeros(nbins_z)

    # open the output files
    # these at some point will turn into FIFOs...?
    filehandles = []
    outcat_filenames = []
    for zbin in range(nbins_z):
        outfilename = 'sample_catalog' + str(zbin)
        outcat_filenames.append(outfilename)
        f1=open(outfilename, 'w')
        filehandles.append(f1)

    mags = []
    if( mag_maxz < max_z ):
        mag_maxz = max_z
    for z in range(mag_nzs):
        mags.append([])
    
    # read the input catalog!
    # this SHOULD be an sql query to the DB...
    # but, for now it is just an input file.    catalog = open(filename, "r")
    catalog = open(filename, "r")
    for line in catalog:
        splitted_line = line.split()
        
        ra = float(splitted_line[0])
        dec = float(splitted_line[1])

        # for the N(z) hdf5 file to be given to the modelling code
        specz = float(splitted_line[3])
        photz = float(splitted_line[2])
        if( specz > 0.0 and specz < 10.0 ):
            z_spec.append(specz)
            z_phot.append(photz)
        
        mag = float(splitted_line[4])
        
        # keep the info for the slopes even if outside the bin ranges
        zbin_for_mag = int(math.floor(photz/(mag_maxz/mag_nzs)))
        if( zbin_for_mag > mag_nzs ):
            continue
        mags[zbin_for_mag].append(mag)
        
        # only keep objs within the z binning range
        if( photz < min_z or photz > max_z ):
            continue
        
        bin_z = -1
        for zbin in range(nbins_z):
            if( photz >= z_means[zbin]-z_widths[zbin] and photz < z_means[zbin]+z_widths[zbin] ):
                bin_z = zbin
                break
        if bin_z < 0:
            continue
        data[bin_z] += 1
    
        # if brighter than the mag cut, we want to correlate these!
        if( mag < mag_cuts[bin_z] ):
            # print to output "stream"
            outstring = "%f %f %f\n" %(ra, dec, mag)
            filehandles[zbin].write(outstring)

    
    # print >> sys.stderr, "read in catalog.  %d spectroscopic zs." %(len(z_spec))
    for zbin in range(nbins_z):
        filehandles[zbin].close()

    # loop over z bins to calculate one slope per each
    slope_array = numpy.zeros(mag_nzs)
    for zbin in range(mag_nzs):

        # don't bother trying if there are very few objs
        if( len(mags[zbin]) < 100 ):
            slope_array[zbin] = 0.
            continue

        magz = (zbin+0.5)*mag_maxz/mag_nzs
        if magz >= max_z:
            zbin2 = nbins_z-1
        else:
            zbin2 = -1
            for z in range(nbins_z):
                if( magz >= z_means[z]-z_widths[z] and magz < z_means[z]+z_widths[z] ):
                    zbin2 = z
                    break
            if zbin2 < 0:
                zbin2 = 0
        mag_cut = mag_cuts[zbin2]

        # kde!
        # NOTE: i'm setting alpha=0 for bins with no data!
        kde = gaussian_kde(mags[zbin])
        y_kde1 = kde.evaluate(mag_cut)
        y_kde2 = kde.evaluate(mag_cut+0.01)
        if( y_kde1 == 0 or y_kde2 == 0 ):
            slope_kde = 0.4
        else:
            slope_kde = (numpy.log10(y_kde2)-numpy.log10(y_kde1))/0.01
        # print >> sys.stderr, "kde slope %f" %slope_kde

        # for the actual result, let's use the fit slope, and use the other two to estimate an error.
        slope_array[zbin] = 2.5*slope_kde - 1
        print >> sys.stderr, "z bin %d at %f s = %f alpha = %f" %(zbin, (zbin+0.5)*mag_maxz/mag_nzs, slope_kde, slope_array[zbin])
    
    
    
    # now, write outputs...
    # N(z) and slopes in an hdf5 file.
    # f = tables.openFile('sample_info.hdf5', 'w')

    # photoz: a table
    f.createGroup('/', 'photoz')
    f.createGroup('/photoz', 'catalog')
    photoz_table = f.createTable("/photoz/catalog", pop, photoz_entry)
    photoz_table.setAttr('pop', json.dumps(pop))
    row = photoz_table.row
    for i in range(len(z_phot)):
        # First, assign the values to the Particle record
        row['z_p']  = z_phot[i]
        row['z_s'] = z_spec[i]
        row.append()
    photoz_table.flush()

    # noise: an array of one value per redshift bin.
    f.createGroup('/', 'noise')
    noise_array = numpy.zeros(nbins_z)
    for z in range(nbins_z):
        if data[z] > 0:
            noise_array[z] = 1.0/data[z]

    noise = f.createArray('/noise', 'noise1', numpy.diag(noise_array))
    noise.setAttr('ftype0', json.dumps('counts'))
    noise.setAttr('pop0', json.dumps(pop))
    noise.setAttr('ftype1', json.dumps('counts'))
    noise.setAttr('pop1', json.dumps(pop))

    # slopes: a table, for lots of redshift values.
    f.createGroup('/', 'slopes')
    slopes_table = f.createTable('/slopes', 'slope1', slopes_entry)
    slopes_table.setAttr('ftype', json.dumps('counts'))
    slopes_table.setAttr('pop', json.dumps(pop))
    row = slopes_table.row
    for i in range(mag_nzs):
        row['z'] = (i+0.5)*mag_maxz/mag_nzs
        row['counts'] = slope_array[i]
        row['mag'] = 0. # TODO!
        row['size'] = 0. # TODO?
        row.append()
    slopes_table.flush()
    
    return outcat_filenames
    


def main():

    # Right now this is a mess, because i've tried a number of ways of doing this and none seems very good...
    # But let's leave this as it is for now, and improve it in the future.
    # I think that the best approach might be to fit an appropriate functional form to the curve.  (better than a 3rd order polynomial?)

    if len(sys.argv) != 4:
        print "Usage: parse_sample.py catalog_filename mag_cut metadata_filename"
        print "       where the galaxy catalog is like: ra dec z_phot z_spec magnitude"
        print "       z_spec<0 or z_spec>10 will be treated as a null value"
        print "       metadata file is hdf5"

    filename = sys.argv[1]
    mag_cut = float(sys.argv[2])
    metadata_filename = sys.argv[3]

    parse_data( filename, mag_cut, metadata_filename )


if __name__ == '__main__':
    main()
