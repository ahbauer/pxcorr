#!/usr/bin/env python
# encoding: UTF8

import numpy as np
import tables
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
from lib_mbe import h5getattr, h5setattr
import healpy
import partpix_map

def main():
    filename = sys.argv[1]
    print "Reading filename %s" %filename
    h5file = tables.openFile(filename)
    
    plot_filename = "pxcorr_jks.pdf"
    pp = PdfPages(plot_filename)
    
    # read in the metadata
    meta = {}
    mObject = h5file.getNode('/meta', 'meta')
    attr_list = mObject.attrs._f_list('user')

    for attr in ['fourier', 'pop']: #, 'ftype']:
        meta[attr] = h5getattr(mObject, attr)

    meta['ftype'] = ['counts'] # HACK

    key = 'l' if meta['fourier'] else 'ang'
    for part in ['mean', 'width']:
        from_key = '{0}_{1}'.format(key, part)
        to_key = 'u_{1}'.format(key, part)

        meta[to_key] = h5getattr(mObject, from_key)

    for i, pop_name in enumerate(meta['pop']):
        test_key = '/meta/{}/meta'.format(pop_name)

        assert test_key in h5file, 'Missing group: {}'.format(test_key)
        #pdb.set_trace()
        metapopObject = h5file.getNode('/meta/'+pop_name, 'meta')
        
        metapop = {}
        attr_list = metapopObject.attrs._f_list('user')
        for attr in attr_list:
            metapop[attr] = h5getattr(metapopObject, attr)
        for attr in  ['z_mean', 'z_width']:
            msg = 'Missing metadata: '+attr
            assert hasattr(metapopObject.attrs, attr), msg
          
        metapop['nbins'] = len(metapop['z_mean']) 
        
        meta[pop_name] = metapop
    
    u_index = 0
    # if len(meta['u_mean']) > 1:
    #     print "There's more than one angular bin, but we'll just look at the first one."
    #     print " JUST KIDDING, the 7th one!"
    #     u_index = 6
    
    # read in the covariance matrices
    # # print h5file.getNode('/', 'cov')._v_children.keys()
    # cov_children = h5file.getNode('/', 'cov')._v_children
    # for cov_name, covObj in cov_children.items():  # WEIRD fails if you do iteritems
    cov = dict((u, {}) for u in meta['u_mean'])
    cov_grp = h5file.getNode('/', 'cov')
    for cov_name, covObj in cov_grp._v_children.iteritems():
        
        # print cov_name

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
    
    
    # read in the correlation matrices
    # and the jackknife samples
    corr = dict((u, {}) for u in meta['u_mean'])
    corr_grp = h5file.getNode('/', 'corr')
    jks = dict((u, {}) for u in meta['u_mean'])
    jk_group = h5file.getNode('/', 'jk')
    for corr_name, corrObj in corr_grp._v_children.iteritems():

        y = lambda x: h5getattr(corrObj, x)
        A = ((y('ftype0'), y('pop0')), (y('ftype1'), y('pop1')))
        attr_corr_list = corrObj.attrs._f_list('user')
        corrArray = corrObj.read()
        
        if y('pop0') != y('pop1'):
            continue
        print "Looking at autocorrelations of {0}, {1} redshift bins".format(y('pop0'), meta[y('pop0')]['nbins'])
        
        aShape = (len(meta['u_mean']), meta[y('pop0')]['nbins'], meta[y('pop1')]['nbins']) 
        msg = 'Wrong dimension of the correlation: ',corr_name , 'Expected: '\
              ,str(list(aShape)),'Shape input: ',str(list(corrArray.shape))
        
        assert aShape == corrArray.shape, msg
        
        u_dict = dict()
        for i,u in enumerate(meta['u_mean']):
            corr[u][A] = corrArray[i]
            jks[u][A] = dict()
            u_dict[i] = u
    
        # read in the jackknife info
        msg = 'Correlation array name does not start with corr...'
        assert corr_name[0:4] == 'corr', msg
        corr_num = corr_name[4:]
        jk_subgroupname = 'jk' + corr_num
        jk_subgroup = h5file.getNode('/jk/', jk_subgroupname)
        for jk_name, jkTable in jk_subgroup._v_children.iteritems():
            jk_partpix = partpix_map.Partpix_Map()
            i = int(jkTable._f_getAttr('ang_index'))
            z0 = int(jkTable._f_getAttr('zbin0'))
            z1 = int(jkTable._f_getAttr('zbin1'))
            jk_partpix.read_from_hdf5table(jkTable)
            if i not in jks[u_dict[i]][A]:
                jks[u_dict[i]][A][i] = dict()
            if z0 not in jks[u_dict[i]][A][i]:
                jks[u_dict[i]][A][i][z0] = dict()
            jks[u_dict[i]][A][i][z0][z1] = jk_partpix
            
        for zbin in range(meta[y('pop0')]['nbins']):
            fig = plt.figure()
            ax0 = fig.add_subplot(1,1,1, title='Jackknives, pop {0} zbin {1} (Corr = {2:4.3e} +/- {3:4.3e})'.format(y('pop0'), zbin, corr[u_dict[u_index]][A][:][zbin][zbin], np.sqrt(cov[u_dict[u_index]][A+A][u_index][zbin][zbin][zbin][zbin])), xlabel='correlation at 1 degree')
            plt.hist(jks[u_dict[u_index]][A][u_index][zbin][zbin].partmap, bins=50)
            plt.axvline(corr[u_dict[u_index]][A][:][zbin][zbin], color='b', linestyle='dashed', linewidth=2)
            pp.savefig()
            # thetas, phis = healpy.pixelfunc.pix2ang(jks[u_dict[0]][A][0][zbin][zbin].Nside(), jks[u_dict[0]][A][0][zbin][zbin].pixel_mapping_arraytohigh)
            # ras = phis*180./3.1415926
            # decs = dec = 90. - thetas*180./3.1415926;
            # fig = plt.figure()
            # ax0 = fig.add_subplot(1,1,1, title='Jackknives, zbin {0}'.format(zbin),xlabel='RA',ylabel='Dec')
            # plt.hexbin(ras,decs,jks[u_dict[0]][A][0][zbin][zbin].partmap, gridsize=int(1.2*np.sqrt(len(jks[u_dict[0]][A][0][zbin][zbin].partmap))),extent=[60,90,-62,-42])
            # pp.savefig()
            
            
            plt.clf()
            hpmap = np.zeros(jks[u_dict[u_index]][A][u_index][zbin][zbin].Npix())
            for i in range(jks[u_dict[u_index]][A][u_index][zbin][zbin].Npartpix()):
                hpmap[jks[u_dict[u_index]][A][u_index][zbin][zbin].pixel_mapping_arraytohigh[i]] = jks[u_dict[u_index]][A][u_index][zbin][zbin].partmap[i]
            minval = np.min(jks[u_dict[u_index]][A][u_index][zbin][zbin].partmap)
            healpy.visufunc.cartview(hpmap,lonra=[60,90],latra=[-62,-42],min=minval,flip='geo')
            pp.savefig()
        
    pp.close()
    h5file.close()
        
if __name__ == '__main__':
    main()
    