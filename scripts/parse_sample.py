#!/usr/bin/env python
# encoding: UTF8

import sys
import os
import math
import numpy
import tables
import json
from array import array
from scipy.stats import gaussian_kde
sys.path.append("/Users/bauer/software/pxcorr/src/")
import make_maps

class photoz_entry(tables.IsDescription):
    zbin = tables.Int32Col(dflt=0, pos=0)
    z_s = tables.Float64Col(dflt=0.0, pos=1) 
    z_p = tables.Float64Col(dflt=0.0, pos=2)
    weight = tables.Float64Col(dflt=0.0, pos=3)

class slopes_entry(tables.IsDescription):
    zbin = tables.Int32Col(dflt=0, pos=0) 
    z = tables.Float64Col(dflt=0.0, pos=1) 
    counts = tables.Float64Col(dflt=0.0, pos=2)
    mag = tables.Float64Col(dflt=0.0, pos=3)
    size = tables.Float64Col(dflt=0.0, pos=4)



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

            outarray = array('d', [float(splitted_line[0]), float(splitted_line[1]), mag])
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


def parse_data( filename, mag_cut, suffix, z_mean, z_width, add_nofz=True, sparse=True ):
    
    z_edges = []
    z_edges.append(z_mean-z_width/2.0)
    z_edges.append(z_mean+z_width/2.0)
    
    # arrays for output later
    z_phot = []
    z_spec = []
    
    # what kind of binning do we want for the slopes?
    mag_nzs = 100
    mag_maxz = 2.0


    outfilename = 'sample_catalog_' + suffix
    f1=open(outfilename, 'wb')

    mags = []
    if( mag_maxz < z_edges[1] ):
        mag_maxz = z_edges[1]
    for z in range(mag_nzs):
        mags.append([])
    
    # read the input catalog!
    # this SHOULD be an sql query to the DB...
    # but, for now it is just an input file.    catalog = open(filename, "r")
    catalog = open(filename, "r")
    nsample = 0
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
        
        
        # are we in the correlation z bin?
        if( photz <= z_edges[0] or photz >= z_edges[-1] ):
            continue
        
        # keep the info for the slopes even if outside the mag bin ranges
        # but get rid of objects well below the mag limit since those will 
        # probably be bad.
        zbin_for_mag = int(math.floor(photz/(mag_maxz/mag_nzs)))
        if( zbin_for_mag > mag_nzs-1 ):
            continue
        if mag < mag_cut+1:
            mags[zbin_for_mag].append(mag)
        
        
        # if brighter than the mag cut, we want to correlate these!
        if( mag < mag_cut ):
            # print to output "stream"
            if( add_nofz and specz > 0.0 and specz < 10.0 ):
                z_spec.append(specz)
                z_phot.append(photz)
            
            outarray = array('d', [float(splitted_line[0]), float(splitted_line[1]), mag]) # not numpy.array
            outarray.tofile(f1)
            nsample += 1
    
    print >> sys.stderr, "read in catalog."
    if add_nofz:
        print >> sys.stderr, "%d spectroscopic zs." %(len(z_spec))

    print >> sys.stderr, "%d objs" %(nsample)
    f1.close()
    
    # save N(z) info to a file
    nofz_filename = "nofz_" + suffix + ".ssv"
    if add_nofz:
        n_sparse = 10000
        if sparse and len(z_phot)>n_sparse:
            inds = numpy.random.randint(0, len(z_phot), n_sparse)
            z_phot2 = numpy.array(z_phot)
            z_phot2 = z_phot2[inds]
            z_phot = z_phot2
            z_spec2 = numpy.array(z_spec)
            z_spec2 = z_spec2[inds]
            z_spec = z_spec2
        
        # write to output file, so the epilogue can read it in and save it to the hdf5 file.
        outfile = open( nofz_filename, 'w' )
        for i in range(n_sparse):
            outfile.write( "{0:.4f} {0:.4f} 1.0\n".format(z_phot[i], z_spec[i]) )
        outfile.close()
    
    
    # measure slopes of the number counts, etc.
    slopes_filename = measure_slopes( mags, mag_maxz, mag_nzs, mag_cut, suffix, use_mags )
    
    return outfilename, slopes_filename, nofz_filename



def measure_slopes( mags, mag_maxz, mag_nzs, mag_cut, suffix, use_mags ):

    delta_mag = 0.05

    # loop over z bins to calculate one slope per each
    slope_array = numpy.zeros(mag_nzs)
    slope_m_array = numpy.zeros(mag_nzs)
    if not use_mags:
        print >> sys.stderr, "setting all alphas to zero!"

    else:
        for zbin_fine in range(mag_nzs):

            # don't bother trying if there are very few objs
            nm = len(mags[zbin_fine])
            if( nm < 100 ):
                continue
            mag_array = numpy.array(mags[zbin_fine])
            minmag = numpy.min(mag_array)
            # print magarray
            # magz = (zbin_fine+0.5)*mag_maxz/mag_nzs

            # counts! kde!
            # NOTE: i'm setting alpha=0 for bins with no data!
            # kde = gaussian_kde(mag_array)
            # y_kde1 = nm*kde.integrate_box_1d(minmag, mag_cut) #kde.evaluate(mag_cut)
            # y_kde2 = nm*kde.integrate_box_1d(minmag, mag_cut+delta_mag) #kde.evaluate(mag_cut+delta_mag)
            # if( y_kde1 == 0 or y_kde2 == 0 ):
            #     slope_kde = 0.4
            # else:
            #     slope_kde = (numpy.log10(y_kde2)-numpy.log10(y_kde1))/delta_mag
            # slope_array[zbin_fine] = 2.5*slope_kde - 1
            
            
            # mags!
            mags1 = mag_array[(mag_array < mag_cut)]
            mags2 = mag_array[(mag_array < mag_cut+delta_mag)]
            slope_m_array[zbin_fine] = -1.0857*(1.0-(numpy.mean(mags2)-numpy.mean(mags1))/delta_mag)
            
            # no kde?
            n1 = len(mags1)
            n2 = len(mags2)
            slope_counts = 0.4
            if n1> 0 and n2 > 0:
                slope_counts = (numpy.log10(n2)-numpy.log10(n1))/delta_mag
            slope_array[zbin_fine] = 2.5*slope_counts - 1

            # print >> sys.stderr, "z bin %d at %f s = %f alpha_c = %f alpha_m = %f" %(zbin_fine, (zbin_fine+0.5)*mag_maxz/mag_nzs, slope_kde, slope_array[zbin_fine], slope_m_array[zbin_fine])
            print >> sys.stderr, "z bin %d at %f s = %f alpha_c = %f alpha_m = %f" %(zbin_fine, (zbin_fine+0.5)*mag_maxz/mag_nzs, slope_counts, slope_array[zbin_fine], slope_m_array[zbin_fine])

    assert len(slope_array), 'Error parse_sample.py:measure_slopes slope_array has length zero!'
    
    # write to output file, so the epilogue can read it in and save it to the hdf5 file.
    outfilename = "slopes_" + suffix + ".ssv"
    outfile = open( outfilename, 'w' )
    for s in range(len(slope_array)):
        outfile.write( "{0:.4f} {0:.4f} {0:.4f}\n".format((s+0.5)*mag_maxz/mag_nzs, slope_array[s], slope_m_array[s]) )
    outfile.close()
    return outfilename



def slopes_to_hdf5( f, maps_list ):
    
    try:
        g = f.getNode('/', 'slopes', 'Group')
    except tables.exceptions.NoSuchNodeError:
        f.createGroup('/', 'slopes')
    
    slope_info = dict()
    for map1 in maps_list:
        slope_info[map1['pop']] = { 'zbin': map1['cut_index'], 'filename': map1['slope_path'] }
        
    for pop in slope_info.keys():
        try:
            slopes_table = f.createTable('/slopes', pop, slopes_entry)
        except tables.exceptions.NodeError:
            # we've already written the slopes for this population.  return!
            return
        
        slopes_table.setAttr('pop', json.dumps(pop))
        row = slopes_table.row
        
        # read in the slopes from the file
        infile = open( slope_info[pop]['filename'], 'r' )
        for line in infile:
            entries = line.split()
            row['zbin'] = slope_info[pop]['zbin']
            row['z'] = float(entries[0])
            row['counts'] = float(entries[1])
            row['mag'] = float(entries[2])
            row['size'] = 0. # TODO?
            row.append()
        infile.close()
        
        slopes_table.flush()


def nofz_to_hdf5( f, maps_list ):
    
    try:
        g = f.getNode('/', 'photoz', 'Group')
    except tables.exceptions.NoSuchNodeError:
        f.createGroup('/', 'photoz')
    
    try:
        g = f.getNode('/photoz', 'catalog', 'Group')
    except tables.exceptions.NoSuchNodeError:
        f.createGroup('/photoz', 'catalog')
            
    nofz_info = dict()
    for map1 in maps_list:
        nofz_info[map1['pop']] = { 'zbin': map1['cut_index'], 'filename': map1['nofz_path'] }
        
    for pop in nofz_info.keys():
        try:
            photoz_table = f.createTable('/photoz/catalog', pop, photoz_entry)
        except tables.exceptions.NodeError:
            # we've already written the nofz for this population.  return!
            return
        
        photoz_table.setAttr('pop', json.dumps(pop))
        row = photoz_table.row
        
        # read in the slopes from the file
        nofz_filename = nofz_info[pop]['filename']
        if not (nofz_filename is None or nofz_filename == 'None' or nofz_filename == 'none'):
            infile = open( nofz_filename, 'r' )
            for line in infile:
                entries = line.split()
                # First, assign the values to the Particle record
                row['zbin'] = nofz_info[pop]['zbin']
                row['z_p']  = float(entries[0])
                row['z_s'] = float(entries[1])
                row['weight'] = float(entries[2])
                row.append()
            infile.close()
        
        photoz_table.flush()
    
    


def construct_inputs( suffix, catalog_filename, ftypes, ang_means, ang_widths, nofz_filename, mask_filename, z_mean, z_width, zbin, mag_cut, pop ):
    
    # do we want to figure out the N(z) from the catalog or a separate training set file?
    nofz_from_data = False
    if nofz_filename[0] == catalog_filename or nofz_filename[0] is None or nofz_filename[0] == 'None':
        nofz_from_data = True

    subcat_filename = None
    slopes_filename = None
    nobj = 1
    if catalog_filename.endswith(".fits"):
        print "Found an input Healpix map!"
        subcat_filename = catalog_filename
        
        # make sure there exists a (dummy) nofz_file to read in later
        if nofz_from_data:
            # save N(z) info to a file
            nofz_filename = "nofz_" + suffix + ".ssv"
            # write to output file, so the epilogue can read it in and save it to the hdf5 file.
            outfile = open( nofz_filename, 'w' )
            outfile.write( "0.0 0.0 1.0\n" )
            outfile.close()
            nofz_filename = [nofz_filename]
            
        # make sure there exists a (dummy) slopes file to read in later
        slopes_filename = "slopes_" + suffix + ".ssv"
        # write to output file, so the epilogue can read it in and save it to the hdf5 file.
        outfile = open( slopes_filename, 'w' )
        outfile.write( "0.0 0.0 0.0\n" )
        outfile.close()
                
    else:
        # parse the catalog: construct N(z) data, calculate slope, etc.
        # these things are necessary for interpretation with martin erikson's modeling code.
        subcat_filename, slopes_filename, nofz_filename2 = parse_data( catalog_filename, mag_cut, suffix, z_mean, z_width, nofz_from_data, sparse=True )
        if nofz_from_data:
            nofz_filename = [nofz_filename2]
    
    use_counts = False
    if 'counts' in ftypes:
        use_counts = True
    use_mags = False
    if 'mags' in ftypes:
        use_mags = True
    make_maps.make_maps( subcat_filename, mask_filename, ang_means, ang_widths, use_counts, use_mags, suffix, pop, zbin );
    print "Finished making map from %s" %subcat_filename

    return slopes_filename, nofz_filename


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
