#!/usr/bin/env python
# encoding: UTF8

# Module for reading in 2D correlations. The HDF5 format
# is defined in doc/corr_format.

import json
import os
import pdb
import numpy as np
import tables
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
sys.path.append("~/software/python/pausci")
from sci.lib.libh5attr import h5getattr, h5setattr

class pkg_read(object): 
    def __init__(self, file_name, root= '/', only_open=False):
        """Initializes the input of correlations, covariances and metadata
           from the data.
        """
   
        file_name = os.path.expanduser(file_name)
        self.h5file = tables.openFile(file_name, title = 'Angular correlations',
                                      rootUEP=root)

        if not only_open:
            self.read()

    def _test_meta(self, meta):
        """Test meta data consistency."""

        for pop_name in meta['pop']:
            msg_mean = 'Negative redshift in: ' + pop_name
            msg_width = 'Negative width in: ' + pop_name

            assert (0 <= meta[pop_name]['z_mean']).all(), msg_mean
            if 'z_width' in meta[pop_name]:
                assert (0 <= meta[pop_name]['z_width']).all(), msg_width
            if 'z_sigma' in meta[pop_name]:
                assert (0 <= meta[pop_name]['z_sigma']).all(), msg_width

    def _get_metapop(self, mpopObject):
        """Reads the metadata from the HDF5 file and store
           them in a dictionary.
        """
        
        metapop = {}

        attr_list = mpopObject.attrs._f_list('user')
        for attr in attr_list:
            metapop[attr] = h5getattr(mpopObject, attr)
            
        for attr in  ['z_mean', 'z_width']:
            msg = 'Missing metadata: '+attr
            assert hasattr(mpopObject.attrs, attr), msg
          
        metapop['nbins'] = len(metapop['z_mean']) 

        return metapop

    def _get_meta(self):
        """Reads the metadata from the HDF5 file and store
           them in a dictionary.
        """

        meta = {}
        mObject = self.h5file.getNode('/meta', 'meta')
        attr_list = mObject.attrs._f_list('user')

        for attr in ['fourier', 'pop']: #, 'ftype']:
            meta[attr] = h5getattr(mObject, attr)

        meta['ftype'] = ['counts'] # HACK

        key = 'l' if meta['fourier'] else 'ang'
        for part in ['mean', 'width']:
            from_key = '{0}_{1}'.format(key, part)
            to_key = 'u_{1}'.format(key, part)

            meta[to_key] = h5getattr(mObject, from_key)

        # TODO: Looks redundant.
        for attr in ['fourier']: # 'pop', 
            msg = 'Missing metadata: '+attr
            assert hasattr(mObject.attrs, attr), msg

        if hasattr(mObject.attrs, 'fid'):
            meta['fid'] = h5getattr(mObject, 'fid')

        for i, pop_name in enumerate(meta['pop']):
            test_key = '/meta/{}/meta'.format(pop_name)

            assert test_key in self.h5file, 'Missing group: {}'.format(test_key)
            #pdb.set_trace()
            metapopObject = self.h5file.getNode('/meta/'+pop_name, 'meta')
            metapop = self._get_metapop(metapopObject)
            meta[pop_name] = metapop
 
        self._test_meta(meta)

        return meta

    def _get_corr(self, meta):
        """Reads the correlations from the HDF5 file and store
           them in a dictionary.
        """

        corr = dict((u, {}) for u in meta['u_mean'])
        corr_grp = self.h5file.getNode('/', 'corr')
        for corr_name, corrObj in corr_grp._v_children.iteritems():

            y = lambda x: h5getattr(corrObj, x)
            A = ((y('ftype0'), y('pop0')), (y('ftype1'), y('pop1')))
            attr_corr_list = corrObj.attrs._f_list('user')
            corrArray = corrObj.read()
            
            aShape = (len(meta['u_mean']), meta[y('pop0')]['nbins'], meta[y('pop1')]['nbins']) 
            msg = 'Wrong dimension of the correlation: ',corr_name , 'Expected: '\
                  ,str(list(aShape)),'Shape input: ',str(list(corrArray.shape))
            
            assert aShape == corrArray.shape, msg
            
            for i,u in enumerate(meta['u_mean']):
                corr[u][A] = corrArray[i]
            
            # for i in range (0,len(meta['pop'])):
            #     pop_name = y('pop' + str(i))
            #     aShape = (len(meta['u_mean']), meta[pop_name]['nbins'], meta[pop_name]['nbins']) 
            #     msg = 'Wrong dimension of the correlation: ',corr_name , 'Expected: '\
            #           ,str(list(aShape)),'Shape input: ',str(list(corrArray.shape))
            #   
            #     assert aShape == corrArray.shape, msg
            # 
            #     for i,u in enumerate(meta['u_mean']):
            #         corr[u][A] = corrArray[i]
                    
        return corr


    def _get_cov(self, meta):
        """Reads the covariance matrices from the HDF5 file and store
           them in a dictionary.
        """    

        cov = dict((u, {}) for u in meta['u_mean'])
        for cov_name, covObj in self.h5file.getNode('/', 'cov')._v_children.iteritems():
            
            y = lambda x: h5getattr(covObj, x)
            A = ((y('ftype0'), y('pop0')), (y('ftype1'), y('pop1')), \
                 (y('ftype2'), y('pop2')), (y('ftype3'), y('pop3')))
           
            z = lambda i: meta[h5getattr(covObj, 'pop'+str(i))]['nbins'] 
            nbins = tuple([z(i) for i in range(4)])

            # There can either be 5 or 6 indexes.
            expected = lambda i: i*(len(meta['u_mean']),) + nbins

            covArray = covObj.read()
            msg_covdim = 'Wrong dimension of the covariance: {0}\nExpected: {1} or {2} \nInput: {3}'\
                         .format(cov_name, expected(1), expected(2), covArray.shape)

            assert covArray.shape in [expected(1), expected(2)], msg_covdim

            for i,u in enumerate(meta['u_mean']):
                cov[u][A] = covArray[i]

        return cov

    def _get_ngal(self, meta):
        """Reads the galaxy ditribution for each population
           from the HDF5 file and store them in a dictionary
        """

        # Note: This functionality has not been officially documented or
        #       even used. How we are going to specify n(z) is still an
        #       open question.
        ngal = {}
        for ngal_name, ngalObj in self.h5file.getNode('/', 'ngal')._v_children.iteritems():
            A = h5getattr(ngalObj, 'pop')
            ngal[A] = ngalObj.read()
            
        return ngal

    def _get_noise(self, meta):
        """Noise which needs to be added to the theoretical models of the correlations."""


        noise = {}
        for name, noiseObj in self.h5file.getNode('/noise')._v_children.iteritems():
            y = lambda x: h5getattr(noiseObj, x)

            ftype0, pop0 = y('ftype0'), y('pop0')
            ftype1, pop1 = y('ftype1'), y('pop1')
            A = ((ftype0, pop0), (ftype1, pop1))

            part = noiseObj.read()
            meta_shape = meta[pop0]['nbins'], meta[pop1]['nbins']

            msg = 'The noise for {0}.{1} and {2}.{3} is expected to be {4}.\nInput shape {5}.'
            assert meta_shape == part.shape, \
                   msg.format(ftype0, pop0, ftype1, pop1, meta_shape, part.shape)


            noise[A] = part

        return noise

    def _get_catalogs(self, meta):
        """Read in photo-z tables used to determine the n(z) for each bin."""

        msg_nottable = 'Galaxy catalogs should be tables.'
        msg_cols = 'Catalogs require the fields z_s and z_p'

        photoz = self.h5file.getNode('/photoz/catalogs')
        catalogs = {x:[] for x in meta['pop']}
        for cat_name, catObj in photoz._v_children.iteritems():
            # TODO; Check photo-z is specified for all the populations..

            assert isinstance(catObj, tables.table.Table), msg_nottable
            colset = set(catObj.cols._v_colnames)
            assert set(['z_s', 'z_p']).issubset(colset), msg_cols

            z_s = catObj.read(field='z_s')
            z_p = catObj.read(field='z_p')

            pop_name = h5getattr(catObj, 'pop')
            catalogs[pop_name].append({'z_s': z_s, 'z_p': z_p})

        return catalogs

    def _get_nz(self, meta):
        """Static n(z). Sometime people like to define those."""

        nz_grp = self.h5file.getNode('/photoz/nz')
        nz_info = {x:[] for x in meta['pop']}
        for nz_name, grp in nz_grp._v_children.iteritems():
            z_path = '/photoz/nz/{}/z'.format(nz_name)
            nz_path = '/photoz/nz/{}/val'.format(nz_name)

            z_node = self.h5file.getNode(z_path)
            z = z_node.read()
            nz = self.h5file.getNode(nz_path).read()

            assert z.ndim == 1, 'Redshifts array should be an 1D array.'
            assert nz.ndim == 2, 'Redshifts array should be an 2D array.'
            assert len(z) == nz.shape[1], 'Different length of the z and nz array.'

            # The solution here is not the optimal...
            popid = h5getattr(z_node, 'pop')
            nz_info[popid].append({'z': z, 'nz': nz})

        return nz_info

    def _get_photoz(self, meta):
        """Read in the photo-z informations which can be specified in different
           ways.
        """

        photoz = {}
        if '/photoz/catalogs' in self.h5file:
            photoz['catalogs'] = self._get_catalogs(meta)
        if '/photoz/nz' in self.h5file:
            photoz['nz'] = self._get_nz(meta)

        return photoz

    def _get_slopes(self, meta):
        """Slopes used for magnification."""

        slopes = {}
        for slope_table in self.h5file.getNode('/slopes')._v_children.itervalues():
            pop = h5getattr(slope_table, 'pop')
            assert pop in meta['pop'], 'Not a population:{}'.format(pop)

# TODO: Check that all the slopes are included....
#            assert hasattr(slope_table.cols, 'z'), 'Missing z information.'

            slopes[pop] = slope_table.read()

        return slopes

    def read(self):
        """Read in the data."""

        meta = self._get_meta()
        self.corr =  self._get_corr(meta)

        if '/cov' in self.h5file:
            self.cov = self._get_cov(meta)
            
        if '/ngal' in self.h5file:
            self.ngal = self._get_ngal(meta)
            
        # AHB cheating:
        # if '/noise' in self.h5file:
        #     self.noise = self._get_noise(meta)

        if '/photoz' in self.h5file:
            self.photoz = self._get_photoz(meta)

        if '/slopes' in self.h5file:
            self.slopes = self._get_slopes(meta)

 
        self.meta = meta
        self.h5file.close()

    def plot(self):
        """Make some plots of the correlations"""
        
        # open the plot file
        filename = "pxcorr_out.pdf"
        pp = PdfPages(filename)
        xlab = "theta (degrees)"
        ylab= "w(theta)"
        
        us = self.meta['u_mean']
        As = self.corr[self.corr.keys()[0]].keys()
        umin = 0.9*us[0] # - (us[1]-us[0])
        if umin <= 0.:
            umin = 0.003
        umax = us[-1] + (us[-1]-us[-2])
                    
                
        for A in As:
            AA = A + A
            nzbins0 = self.corr[us[0]][A].shape[0]
            nzbins1 = self.corr[us[0]][A].shape[1]
            for zbin1 in range(nzbins0):
                for zbin2 in range(nzbins1):
                    ttl = str(A) + ", z bins " + str(zbin1) + "x" + str(zbin2)
                    fig = plt.figure()
                    ws = []
                    for u in us:
                        ws.append(self.corr[u][A][:][zbin1][zbin2])
                    wmin = min(ws)
                    if wmin < 0:
                        wmin *= 1.5
                    else:
                        wmin *= 0.5
                    wmax = 1.5*max(ws)
                    ews = []
                    for i,u in enumerate(us):
                        ews.append(np.sqrt(self.cov[u][AA][i][zbin1][zbin2][zbin1][zbin2]))
                    ax0 = fig.add_subplot(1,1,1, xlim=[umin,umax], ylim=[wmin,wmax], xlabel=xlab, ylabel=ylab, title=ttl)
                    ax0.errorbar(us,ws,yerr=ews,fmt='.')
                    ax0.set_xscale('log')
                    pp.savefig()
        
        pp.close()
                    
    
    def print_corr(self):
        As = self.corr[self.corr.keys()[0]].keys()
        us = self.meta['u_mean']
        nzbins = self.corr[us[0]][As[0]].shape[0]
        for A in As:
            AA = A + A
            f = open("pxcorr_corr_{0}".format(A), 'w')
            for iu, u in enumerate(us):
                f.write( "{0} ".format(u) )
                for zbin1 in range(nzbins):
                    for zbin2 in range(nzbins):
                        f.write( "{0} {1}".format(self.corr[u][A][:][zbin1][zbin2], np.sqrt(self.cov[u][AA][iu][zbin1][zbin2][zbin1][zbin2])) )
                f.write( "\n" )
            f.close()
    
    def print_cov(self):
        As = self.corr[self.corr.keys()[0]].keys()
        us = self.meta['u_mean']
        nzbins = self.corr[us[0]][As[0]].shape[0]
        cov_matrix = np.zeros((nzbins*len(us), nzbins*len(us)))
        for A in As:
            AA = A + A
            for iu1, u1 in enumerate(us):
                for iu2, u2 in enumerate(us):
                    for zbin1 in range(nzbins):
                        for zbin2 in range(nzbins):
                            icov = iu1*nzbins + zbin1
                            jcov = iu2*nzbins + zbin2
                            cov_matrix[icov,jcov] = self.cov[u1][AA][iu2][zbin1][zbin2][zbin1][zbin2]
            np.savetxt("pxcorr_cov_{0}".format(A), cov_matrix)
    
    def print_simple(self):
        As = self.corr[self.corr.keys()[0]].keys()
        us = self.meta['u_mean']
        for A in As:
            AA = A + A
            nzbins0 = self.corr[us[0]][A].shape[0]
            nzbins1 = self.corr[us[0]][A].shape[1]
            for zbin1 in range(nzbins0):
                for zbin2 in range(nzbins1):
                    f = open("pxcorr_simple_{0}_z{1}z{2}".format(A,zbin1,zbin2), 'w')
                    # only print the diagonal part of the covariance
                    for iu1, u1 in enumerate(us):
                        correlation = self.corr[u1][A][:][zbin1][zbin2]
                        error = np.sqrt(self.cov[u1][AA][iu1][zbin1][zbin2][zbin1][zbin2])
                        f.write( "{0} {1} {2}\n".format(u1, correlation, error) )
                    f.close()

class meta_read(pkg_read, object): 
    def __init__(self, file_name, root = '/'):
        """Initializes the metadata
           from the data
        """

        self.h5file = tables.openFile(file_name, title = 'Metadata',
                                      rootUEP=root)

    def read(self):
        """Read in the meta part."""

        meta = self._get_meta()
        self.corr =  self._get_corr(meta)

        # Note: Need to be improved to handle group names without ending slash...
        if '/cov' in self.h5file:
            self.cov = self._get_cov(meta)
            
        if '/ngal' in self.h5file:
            self.ngal = self._get_ngal(meta)

        self.meta = meta
        self.h5file.close()


def main():
    filename = sys.argv[1]
    print "reading filename %s" %filename
    pkg = pkg_read(filename)
    pkg.plot()
    pkg.print_simple()


if __name__ == '__main__':
    main()
    
