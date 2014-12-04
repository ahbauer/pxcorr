#!/usr/bin/env python
# encoding: UTF8

import sys
import os
import numpy as np
import pdb
import json
import tables

def main():
    
    galaxy_pop_cov = 'mice_bench'
    galaxy_pop_data = 'benchmark'
    cov_zbins = [0]
    data_zbins = [0]
    galaxy_pop_cov = '\"' + galaxy_pop_cov + '\"'
    galaxy_pop_data = '\"' + galaxy_pop_data + '\"'
    
    
    if len(sys.argv) != 4:
        print 'Usage: replace_covariance good_data.5 good_covariance.h5 output.h5'
        print 'Will only replace the auto-covariance of the population {0}'.format(galaxy_pop)
        exit(1)
    
    assert(len(cov_zbins) == len(data_zbins)), 'The desired covariance and data redshift bins are not the same length.'
    
    datafile = tables.openFile(sys.argv[1])
    covfile = tables.openFile(sys.argv[2])
    
    datafile.copy_file(sys.argv[3], overwrite=True)
    datafile.close()
    outfile = tables.openFile(sys.argv[3],'a')
    
    for leaf in covfile.root.cov._f_walknodes('Leaf'):
        if leaf.attrs['pop0'] == galaxy_pop_cov and leaf.attrs['pop1'] == galaxy_pop_cov and leaf.attrs['pop2'] == galaxy_pop_cov and leaf.attrs['pop3'] == galaxy_pop_cov :
            print 'Found the covariance: {0}'.format(leaf)
            
            for leaf2 in outfile.root.cov._f_walknodes('Leaf'):
                if leaf2.attrs['pop0'] == galaxy_pop_data and leaf2.attrs['pop1'] == galaxy_pop_data and leaf2.attrs['pop2'] == galaxy_pop_data and leaf2.attrs['pop3'] == galaxy_pop_data :
                    print 'Found the data: {0}'.format(leaf2)
                    
                    for z in range(len(cov_zbins)):
                        cz = cov_zbins[z]
                        dz = data_zbins[z]
                        leaf2[:,:,dz,dz,dz,dz] = leaf[:,:,cz,cz,cz,cz]
                        print 'Copied covariance zbin {0} to data zbin {1}'.format(cz,dz)
    outfile.close()
    covfile.close()

if __name__ == '__main__':
    main()