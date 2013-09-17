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
    z_s = tables.Float64Col(dflt=0.0, pos=0) 
    z_p = tables.Float64Col(dflt=0.0, pos=1)
    weight = tables.Float64Col(dflt=0.0, pos=2)

class slopes_entry(tables.IsDescription):
    z = tables.Float64Col(dflt=0.0, pos=0) 
    counts = tables.Float64Col(dflt=0.0, pos=1)
    mag = tables.Float64Col(dflt=0.0, pos=1)
    size = tables.Float64Col(dflt=0.0, pos=1)

def measure_slopes( mags, f, mag_maxz, mag_nzs, mag_cuts, z_edges, pop, use_mags ):
    
    # loop over z bins to calculate one slope per each
    slope_array = numpy.zeros(mag_nzs)
    if not use_mags:
        print >> sys.stderr, "setting all alphas to zero!"
        for zbin in range(mag_nzs):
            slope_array[zbin] = 0.
        
    else:
        for zbin_fine in range(mag_nzs):

            # don't bother trying if there are very few objs
            if( len(mags[zbin_fine]) < 100 ):
                slope_array[zbin_fine] = 0.
                continue

            magz = (zbin_fine+0.5)*mag_maxz/mag_nzs
            if magz >= z_edges[len(z_edges)-1]:
                zbin2 = len(z_edges)-2 # zbin2 is a bin in the correlation z binning
            else:
                inds = numpy.digitize([magz], z_edges)
                zbin2 = inds[0]-1
                if zbin2 < 0:
                    zbin2 = 0
            mag_cut = mag_cuts[zbin2]

            # kde!
            # NOTE: i'm setting alpha=0 for bins with no data!
            kde = gaussian_kde(mags[zbin_fine])
            y_kde1 = kde.evaluate(mag_cut)
            y_kde2 = kde.evaluate(mag_cut+0.01)
            if( y_kde1 == 0 or y_kde2 == 0 ):
                slope_kde = 0.4
            else:
                slope_kde = (numpy.log10(y_kde2)-numpy.log10(y_kde1))/0.01
            # print >> sys.stderr, "kde slope %f" %slope_kde

            # for the actual result, let's use the fit slope, and use the other two to estimate an error.
            slope_array[zbin_fine] = 2.5*slope_kde - 1
            print >> sys.stderr, "z bin %d at %f s = %f alpha = %f" %(zbin_fine, (zbin_fine+0.5)*mag_maxz/mag_nzs, slope_kde, slope_array[zbin_fine])

    # save slopes: a table, for lots of redshift values.
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
    
def parse_data3pt( filename, mag_cuts, f, sparse=True ):

    # read the needed metadata
    meta = f.getNode('/meta', 'meta')
    injson = getattr(getattr(meta, 'attrs'), 'pop')
    pop = (json.loads(injson))[0]
#    injson = getattr(getattr(meta, 'attrs'), 'angles')
#    angles = json.loads(injson)
    metapopname = '/meta/' + pop
    metapop = f.getNode(metapopname, 'meta')
    injson = getattr(getattr(metapop, 'attrs'), 'z_mean')
    z_means = json.loads(injson)
    injson = getattr(getattr(metapop, 'attrs'), 'z_width')
    z_widths = json.loads(injson)

    z_edges = []
    z_edges.append(z_means[0]-z_widths[0]/2.0)
    for i in range(len(z_means)):
        z_edges.append(z_means[i]+z_widths[i]/2.0)

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
    max_z = z_means[nbins_z-1] + z_widths[nbins_z-1]/2.0
    data = numpy.zeros(nbins_z)

    # open the output files
    # these at some point will turn into FIFOs...?
    # samples = []
    # for zbin in range(nbins_z):
    #     samples.append(dict())
    #     samples[zbin]['ra'] = []
    #     samples[zbin]['dec'] = []
    #     samples[zbin]['mag'] = []

    filehandles = []
    outcat_filenames = []
    for zbin in range(nbins_z):
        outfilename = 'sample_catalog' + str(zbin)
        outcat_filenames.append(outfilename)
        f1=open(outfilename, 'wb')
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
    nsample = numpy.zeros(nbins_z)
    use_mags = True
    for line in catalog:
        splitted_line = line.split()

        # for the N(z) hdf5 file to be given to the modelling code
        if( len(splitted_line) > 4 ):
            photz = float(splitted_line[2])
            specz = float(splitted_line[3])
            mag = float(splitted_line[4])

        elif( len(splitted_line) == 4 ):
            photz = float(splitted_line[2])
            specz = float(splitted_line[3])
            mag = 0.
            use_mags = False

        elif( len(splitted_line) == 3 ):
            photz = float(splitted_line[2])
            specz = photz
            mag = 0.
            use_mags = False

        # what correlation redshift bin are we in?
        if( photz <= z_edges[0] or photz >= z_edges[-1] ):
            continue
        inds = numpy.digitize([photz], z_edges)
        bin_z = inds[0]-1

        # keep the info for the slopes even if outside the bin ranges
        # but get rid of objects well below the mag limit since those will 
        # probably be bad.
        zbin_for_mag = int(math.floor(photz/(mag_maxz/mag_nzs)))
        if( zbin_for_mag > mag_nzs-1 ):
            continue
        interp_mag_cut = numpy.interp(photz, z_means, mag_cuts)
        if mag < interp_mag_cut+1:
            mags[zbin_for_mag].append(mag)

        # only keep objs within the z binning range
        if( photz < min_z or photz > max_z ):
            continue

        if bin_z < 0:
            continue
        data[bin_z] += 1

        # if brighter than the mag cut, we want to correlate these!
        if( mag < mag_cuts[bin_z] ):
            # print to output "stream"
            # outstring = "%s %s %f\n" %(splitted_line[0], splitted_line[1], mag)
            # filehandles[bin_z].write(outstring)
            if( specz > 0.0 and specz < 10.0 ):
                z_spec.append(specz)
                z_phot.append(photz)

            outarray = numpy.array([float(splitted_line[0]), float(splitted_line[1]), mag])
            outarray.tofile(filehandles[bin_z])
            nsample[bin_z] += 1

            # samples[bin_z]['ra'].append(float(splitted_line[0]))
            # samples[bin_z]['dec'].append(float(splitted_line[1]))
            # samples[bin_z]['mag'].append(mag)


    print >> sys.stderr, "read in catalog.  %d spectroscopic zs." %(len(z_spec))
    for zbin in range(nbins_z):
        print >> sys.stderr, "bin %d: %d objs" %(zbin, nsample[zbin])
        filehandles[zbin].close()


    # save N(z) info
    n_sparse = 20000
    if sparse and len(z_phot)>n_sparse:
        inds = numpy.random.randint(0, len(z_phot), n_sparse)
        z_phot2 = numpy.array(z_phot)
        z_phot2 = z_phot2[inds]
        z_phot = z_phot2
        z_spec2 = numpy.array(z_spec)
        z_spec2 = z_spec2[inds]
        z_spec = z_spec2

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
        row['weight'] = 1.0
        row.append()
    photoz_table.flush()

    # save noise info: an array of one value per redshift bin.
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

    # measure slopes of the number counts, etc.
    measure_slopes( mags, f, mag_maxz, mag_nzs, mag_cuts, z_edges, pop, use_mags )

    return outcat_filenames


def parse_data( filename, mag_cuts, f, add_nofz=True, sparse=True ):

    # read the needed metadata
    meta = f.getNode('/meta', 'meta')
    injson = getattr(getattr(meta, 'attrs'), 'pop')
    pop = (json.loads(injson))[0]
    injson = getattr(getattr(meta, 'attrs'), 'ang_mean')
    ang_means = json.loads(injson)
    metapopname = '/meta/' + pop
    metapop = f.getNode(metapopname, 'meta')
    injson = getattr(getattr(metapop, 'attrs'), 'z_mean')
    z_means = json.loads(injson)
    injson = getattr(getattr(metapop, 'attrs'), 'z_width')
    z_widths = json.loads(injson)
    
    z_edges = []
    z_edges.append(z_means[0]-z_widths[0]/2.0)
    for i in range(len(z_means)):
        z_edges.append(z_means[i]+z_widths[i]/2.0)
    
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
    max_z = z_means[nbins_z-1] + z_widths[nbins_z-1]/2.0
    data = numpy.zeros(nbins_z)

    # open the output files
    # these at some point will turn into FIFOs...?
    # samples = []
    # for zbin in range(nbins_z):
    #     samples.append(dict())
    #     samples[zbin]['ra'] = []
    #     samples[zbin]['dec'] = []
    #     samples[zbin]['mag'] = []

    filehandles = []
    outcat_filenames = []
    for zbin in range(nbins_z):
        outfilename = 'sample_catalog' + str(zbin)
        outcat_filenames.append(outfilename)
        f1=open(outfilename, 'wb')
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
    nsample = numpy.zeros(nbins_z)
    use_mags = True
    for line in catalog:
        splitted_line = line.split()

        # for the N(z) hdf5 file to be given to the modelling code
        if( len(splitted_line) > 4 ):
            photz = float(splitted_line[2])
            specz = float(splitted_line[3])
            mag = float(splitted_line[4])
        
        elif( len(splitted_line) == 4 ):
            photz = float(splitted_line[2])
            specz = float(splitted_line[3])
            mag = 0.
            use_mags = False
            
        elif( len(splitted_line) == 3 ):
            photz = float(splitted_line[2])
            specz = photz
            mag = 0.
            use_mags = False
        
        # what correlation redshift bin are we in?
        if( photz <= z_edges[0] or photz >= z_edges[-1] ):
            continue
        inds = numpy.digitize([photz], z_edges)
        bin_z = inds[0]-1

        # keep the info for the slopes even if outside the bin ranges
        # but get rid of objects well below the mag limit since those will 
        # probably be bad.
        zbin_for_mag = int(math.floor(photz/(mag_maxz/mag_nzs)))
        if( zbin_for_mag > mag_nzs-1 ):
            continue
        interp_mag_cut = numpy.interp(photz, z_means, mag_cuts)
        if mag < interp_mag_cut+1:
            mags[zbin_for_mag].append(mag)
        
        # only keep objs within the z binning range
        if( photz < min_z or photz > max_z ):
            continue
        
        if bin_z < 0:
            continue
        data[bin_z] += 1
    
        # if brighter than the mag cut, we want to correlate these!
        if( mag < mag_cuts[bin_z] ):
            # print to output "stream"
            # outstring = "%s %s %f\n" %(splitted_line[0], splitted_line[1], mag)
            # filehandles[bin_z].write(outstring)
            if( add_nofz and specz > 0.0 and specz < 10.0 ):
                z_spec.append(specz)
                z_phot.append(photz)
            
            outarray = numpy.array([float(splitted_line[0]), float(splitted_line[1]), mag])
            outarray.tofile(filehandles[bin_z])
            nsample[bin_z] += 1

            # samples[bin_z]['ra'].append(float(splitted_line[0]))
            # samples[bin_z]['dec'].append(float(splitted_line[1]))
            # samples[bin_z]['mag'].append(mag)

    
    print >> sys.stderr, "read in catalog."
    if add_nofz:
        print >> sys.stderr, "%d spectroscopic zs." %(len(z_spec))

    for zbin in range(nbins_z):
        print >> sys.stderr, "bin %d: %d objs" %(zbin, nsample[zbin])
        filehandles[zbin].close()

    
    # save N(z) info
    if add_nofz:
        n_sparse = 20000
        if sparse and len(z_phot)>n_sparse:
            inds = numpy.random.randint(0, len(z_phot), n_sparse)
            z_phot2 = numpy.array(z_phot)
            z_phot2 = z_phot2[inds]
            z_phot = z_phot2
            z_spec2 = numpy.array(z_spec)
            z_spec2 = z_spec2[inds]
            z_spec = z_spec2
        
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
            row['weight'] = 1.0
            row.append()
        photoz_table.flush()

    # save noise info: an array of one value per redshift bin.
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
    
    # measure slopes of the number counts, etc.
    measure_slopes( mags, f, mag_maxz, mag_nzs, mag_cuts, z_edges, pop, use_mags )

    return outcat_filenames


def add_nofz( filename, f, pop ):
    
    catalog = open(filename, "r")
    z_phot = []
    z_spec = []
    weight = []
    for line in catalog:
        splitted_line = line.split()
        if( splitted_line[0] == '#' ):
            continue
        z_spec.append(splitted_line[0])
        z_phot.append(splitted_line[1])
        weight.append(splitted_line[2])
    
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
        row['weight'] = weight[i]
        row.append()
    photoz_table.flush()
        
        
def main():

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
