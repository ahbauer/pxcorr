#!/usr/bin/env python
# encoding: UTF8

import sys
import os
import math
import numpy
import tables
import json

class photoz_entry(tables.IsDescription):
    z_s = tables.Float32Col(dflt=0.0, pos=0) 
    z_p = tables.Float32Col(dflt=0.0, pos=1)

    
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

    # make a matrix with sides z_phot and mag.
    # the slope of the number counts will be determined from the counts vs mag, 
    # but we also would like to keep vs z for diagnostics.
    nbins_z = len(z_means)
    min_z = 0.0
    max_z = z_means[nbins_z-1] + z_widths[nbins_z-1]
    nbins_m = 100 # must be even
    binwidth_m = 0.1
    min_m = 20.05
    max_m = min_m + nbins_m*binwidth_m
    for mag_cut in mag_cuts:
        if( mag_cut < min_m+binwidth_m or mag_cut > max_m-binwidth_m ):
            print "Error, mag cut of %f is not within the expected range of %f to %f (with buffer %f)" %(mag_cut, min_m, max_m, binwidth_m)
            return;
    data = numpy.zeros((nbins_z,nbins_m))

    # open the output files
    # these at some point will turn into FIFOs...?
    filehandles = []
    outcat_filenames = []
    for zbin in range(nbins_z):
        outfilename = 'sample_catalog' + str(zbin)
        outcat_filenames.append(outfilename)
        f1=open(outfilename, 'w')
        filehandles.append(f1)

    # read the input catalog!
    # this SHOULD be an sql query to the DB...
    # but, for now it is just an input file.
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
        
        if( photz < min_z or photz > max_z or mag < min_m or mag > max_m ):
            continue
        
        bin_m = math.floor((mag-min_m)/binwidth_m)
        bin_z = -1
        for zbin in range(nbins_z):
            if( photz >= z_means[zbin]-z_widths[zbin] and photz < z_means[zbin]+z_widths[zbin] ):
                bin_z = zbin
                break
        if bin_z < 0:
            continue
        data[bin_z][bin_m] += 1
    
        if( mag < mag_cuts[bin_z] ):
            # print to output "stream"
            outstring = "%f %f %f\n" %(ra, dec, mag)
            filehandles[zbin].write(outstring)
    
    # print >> sys.stderr, "read in catalog.  %d spectroscopic zs." %(len(z_spec))
    for zbin in range(nbins_z):
        filehandles[zbin].close()

    # loop over z bins to calculate one slope per each
    slope_array = numpy.zeros(nbins_z)
    slope_error_array = numpy.zeros(nbins_z)
    for zbin in range(nbins_z):
        
        data_array0 = data[zbin].copy()
        if( data_array0.sum() == 0 ):
            print >> sys.stderr, "skipping redshift bin %d (%f)" %(zbin, z_means[zbin])
            continue

        # cheating, to avoid stupid log10(0) errors.
        for i in range(len(data_array0)):
            if data_array0[i] == 0:
                data_array0[i] = 1            
        data_array = numpy.log10(data_array0)

        # what is the bin (fractional) index over which we want the slope?
        bin_mag_cut_float = (mag_cuts[zbin]-min_m)/binwidth_m
        bin_mag_cut = math.floor(bin_mag_cut_float + 0.0001)
        # print >> sys.stderr, "bin of the mag cut = %f -> %d" %(bin_mag_cut_float, bin_mag_cut)

        # check to see that our bins are big enough to yield decent statistics
        while data_array[bin_mag_cut] < 2.5:
            # print >> sys.stderr, "the statistics in our desired bin are too small (N = %d)... rebinning." %data_array0[bin_mag_cut]
            binwidth_m *= 2.0;
            nbins_m = math.ceil(nbins_m/2);
            max_m = min_m + nbins_m*binwidth_m
            data_array0a = numpy.zeros(nbins_m)
            for i in range(len(data_array0a)):
                data_array0a[i] = data_array0[2*i] + data_array0[2*i+1]
            data_array0 = data_array0a
            data_array = numpy.log10(data_array0)
            bin_mag_cut_float = (mag_cuts[zbin]-min_m)/binwidth_m
            bin_mag_cut = math.floor(bin_mag_cut_float + 0.0001)
            # print >> sys.stderr, "bin of the mag cut = %f -> %d" %(bin_mag_cut_float, bin_mag_cut)
    
        # what is the slope between that bin and the one above?
        slope_oneup = (data_array[bin_mag_cut+1] - data_array[bin_mag_cut])/binwidth_m
        slope_onedown = (data_array[bin_mag_cut] - data_array[bin_mag_cut-1])/binwidth_m
        slope_twodown = (data_array[bin_mag_cut] - data_array[bin_mag_cut-2])/(2*binwidth_m)

        print >> sys.stderr, "slopes %f %f %f, n=%d" %(slope_oneup, slope_onedown, slope_twodown, 10.0**data_array[bin_mag_cut])

        # fit the data with a polynomial to try to get a less noisy slope
        # include one point after the cut
        # only include points with log(N)>=2
        # this parameterization seems poor.  maybe we can improve this.
        y_array = []
        x_array = []
        for i in range(len(data_array)):
            if( data_array[i] < 2.0 ):
                continue
            if( i > bin_mag_cut+1.5 ):
                continue
            # to decrease the effects of the parameterizations, don't fit too many points on the brighter side.
            if( i < bin_mag_cut - 1.0/binwidth_m ):
                continue
            x_array.append(i)
            y_array.append(data_array[i])
        fit_terms = numpy.polyfit(x_array, y_array, 3)
        fit_slope = (fit_terms[0]*3.0*bin_mag_cut_float*bin_mag_cut_float + fit_terms[1]*2.0*bin_mag_cut_float + fit_terms[2])/binwidth_m
    
        # smooth the curve?
        data_array_smooth = numpy.ones(len(data_array))
        data_array_smooth_lin = numpy.ones(len(data_array))
        data_array1 = data_array0.copy()
        old_slope_oneup_smooth = 99999.0;
        slope_oneup_smooth = 0.;
        slope_onedown_smooth = 0.;
        slope_twodown_smooth = 0.;
    
        # first, cut out any last bin that only holds a small number of objects.
        for i in reversed(range(len(data_array1))):
            if( data_array1[i] < 0.67*data_array1[i-1] ):
                data_array1[i] = 1;
                # print >> sys.stderr, "cutting out bin %d, mag %f" %(i, i*binwidth_m + min_m)
    
        while( abs(old_slope_oneup_smooth-slope_oneup_smooth) > 0.1*abs(old_slope_oneup_smooth) ):
            old_slope_oneup_smooth = slope_oneup_smooth
            for i in range( 1,len(data_array_smooth)-1,1 ):
                if( data_array1[i-1] > 1.5 and data_array1[i] > 1.5 and data_array1[i+1] > 1.5 ):
                    data_array_smooth_lin[i] = (0.5*data_array1[i] + 0.25*data_array1[i-1] + 0.25*data_array1[i+1])

            # print "DATA_SMOOTH:"
            # print data_array_smooth_lin

            data_array_smooth = numpy.log10(data_array_smooth_lin)
            slope_oneup_smooth = (data_array_smooth[bin_mag_cut+1] - data_array_smooth[bin_mag_cut])/binwidth_m
            slope_onedown_smooth = (data_array_smooth[bin_mag_cut] - data_array_smooth[bin_mag_cut-1])/binwidth_m
            slope_twodown_smooth = (data_array_smooth[bin_mag_cut] - data_array_smooth[bin_mag_cut-2])/(2*binwidth_m)
            print >> sys.stderr, "slopes smooth %f %f %f" %(slope_oneup_smooth, slope_onedown_smooth, slope_twodown_smooth)
            data_array1 = data_array_smooth_lin.copy()
            data_array_smooth_lin = numpy.ones(len(data_array1))
    
        # rebin this a bit wider and redo the slopes
        # this is no better than smoothing?
        # bin_mag_cut_float = (mag_cut-min_m)/(2*binwidth_m)
        # bin_mag_cut = math.floor(bin_mag_cut_float + 0.0001)
        # #print >> sys.stderr, "bin of the mag cut2 = %f -> %d" %(bin_mag_cut_float, bin_mag_cut)
        # data_array2 = numpy.zeros(len(data_array0)/2)
        # for i in range(len(data_array2)):
        #     data_array2[i] = numpy.log10(data_array1[2*i] + data_array1[2*i+1])
        # slope2_oneup = (data_array2[bin_mag_cut+1] - data_array2[bin_mag_cut])/(2*binwidth_m)
        # slope2_onedown = (data_array2[bin_mag_cut] - data_array2[bin_mag_cut-1])/(2*binwidth_m)
        # slope2_twodown = (data_array2[bin_mag_cut] - data_array2[bin_mag_cut-2])/(4*binwidth_m)
    
        #print >> sys.stderr, "slopes2 %f %f %f" %(slope2_oneup, slope2_onedown, slope2_twodown)

        print >> sys.stderr, "fit slope %f" %fit_slope
        # for i in range(len(data_array)):
        #     print "%d %f %f" %(i, data_array[i], fit_terms[0]*i*i*i + fit_terms[1]*i*i + fit_terms[2]*i + fit_terms[3])


        # for the actual result, let's use the fit slope, and use the other two to estimate an error.
        slope_array[zbin] = 2.5*fit_slope - 1
        slope_error_array[zbin] = 2.5*math.sqrt(( (fit_slope-slope_oneup)*(fit_slope-slope_oneup) + (fit_slope-slope_oneup_smooth)*(fit_slope-slope_oneup_smooth) )/2.0)
        print >> sys.stderr, "z bin %d best alpha = %f +/- %f" %(zbin, slope_array[zbin], slope_error_array[zbin])
    
    
    
    # now, write outputs...
    # N(z) and slopes in an hdf5 file.
    # f = tables.openFile('sample_info.hdf5', 'w')

    # photoz: a table
    f.createGroup('/', 'photoz')
    f.createGroup('/photoz', 'catalog')
    # photoz_data = numpy.zeros((2,len(z_phot)))
    # photoz_data[0:] = z_phot
    # photoz_data[1:] = z_spec
    # phzcatalog = f.createArray('/photoz/catalog', pop, photoz_data)
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
        tot = data[z].sum()
        if tot > 0:
            noise_array[z] = 1.0/tot

    noise = f.createArray('/noise', 'noise1', np.diag(noise_array))
    noise.setAttr('ftype0', json.dumps('counts'))
    noise.setAttr('pop0', json.dumps(pop))
    noise.setAttr('ftype1', json.dumps('counts'))
    noise.setAttr('pop1', json.dumps(pop))

    # slopes: an array of one slope per redshift bin
    f.createGroup('/', 'slopes')
    slopes = f.createArray('/slopes', 'slope1', slope_array)
    slopes.setAttr('ftype', json.dumps('counts'))
    slopes.setAttr('pop', json.dumps(pop))
    # f.close()
    
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
