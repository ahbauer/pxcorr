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

test_plots = False
galaxy_string = 'gold'

def main():
    filename = sys.argv[1]
    print "Reading filename %s" %filename
    h5file = tables.openFile(filename)
    
    
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
    
    
    # read in the correlation matrices
    corr = dict((u, {}) for u in meta['u_mean'])
    corr_grp = h5file.getNode('/', 'corr')
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
    
    
    # read in the covariance matrices
    cov = dict((u, {}) for u in meta['u_mean'])
    # print h5file.getNode('/', 'cov')._v_children.keys()
    cov_children = h5file.getNode('/', 'cov')._v_children
    for cov_name, covObj in cov_children.items():  # WEIRD fails if you do iteritems
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
    
    # test: make some plots
    if test_plots:
        filename = "pxcorr_out.pdf"
        pp = PdfPages(filename)
        xlab = "theta (degrees)"
        ylab= "w(theta)"
    
        us = meta['u_mean']
        As = corr[corr.keys()[0]].keys()
        umin = 0.9*us[0] # - (us[1]-us[0])
        if umin <= 0.:
            umin = 0.003
        umax = us[-1] + (us[-1]-us[-2])
        for A in As:
            AA = A + A
            nzbins0 = corr[us[0]][A].shape[0]
            nzbins1 = corr[us[0]][A].shape[1]
            for zbin1 in range(nzbins0):
                for zbin2 in range(nzbins1):
                    ttl = str(A) + ", z bins " + str(zbin1) + "x" + str(zbin2)
                    fig = plt.figure()
                    ws = []
                    for u in us:
                        ws.append(corr[u][A][:][zbin1][zbin2])
                    wmin = min(ws)
                    if wmin < 0:
                        wmin *= 1.5
                    else:
                        wmin *= 0.5
                    wmax = 1.5*max(ws)
                    ews = []
                    for i,u in enumerate(us):
                        ews.append(np.sqrt(cov[u][AA][i][zbin1][zbin2][zbin1][zbin2]))
                    ax0 = fig.add_subplot(1,1,1, xlim=[umin,umax], ylim=[wmin,wmax], xlabel=xlab, ylabel=ylab, title=ttl)
                    ax0.errorbar(us,ws,yerr=ews,fmt='.')
                    ax0.set_xscale('log')
                    pp.savefig()
        pp.close()
    
    # ok, find the galaxy correlations in this big file.
    galaxy_pop = 0
    n_galpops=0
    systematics_pops = []
    for pop in meta['pop']:
        if pop.find(galaxy_string) >= 0:
            print "Found the galaxy population %s" %pop
            galaxy_pop = pop
            n_galpops += 1
        else:
            systematics_pops.append(pop)
    if n_galpops != 1:
        raise AssertionError, "%d populations found with the galaxy string %s" %(n_galpops, galaxy_string)
    
    print "Systematics populations found:"
    print systematics_pops
    
    # now, find the galaxy counts autocorrelation errors
    A = (('counts',galaxy_pop),('counts',galaxy_pop))
    AA = A + A
    us = meta['u_mean']
    nzbins = corr[us[0]][A].shape[0]
    autocorr_errors = np.zeros((nzbins,len(us)))
    autocorr_vals = np.zeros((nzbins,len(us)))
    for zbin in range(nzbins):
        for i,u in enumerate(us):
            autocorr_errors[zbin,i] = np.sqrt(cov[u][AA][i][zbin][zbin][zbin][zbin])
            autocorr_vals[zbin,i] = corr[u][A][:][zbin][zbin]
    # open plot file
    filename = "pxcorr_eval_sys.pdf"
    pp = PdfPages(filename)
    xlab = "theta (degrees)"
    umin = 0.9*us[0] # - (us[1]-us[0])
    if umin <= 0.:
        umin = 0.003
    umax = us[-1] + (us[-1]-us[-2])
    
    # for each systematic, calculate the cross^2 over the auto
    for systematic in systematics_pops:
        # print systematic
        A_cross1 = (('counts',systematic),('counts',galaxy_pop))
        A_cross2 = (('counts',galaxy_pop),('counts',systematic))
        A_auto = (('counts',systematic),('counts',systematic))
        autocorr_sys = np.zeros(len(us))
        autocorr_sys_err = np.zeros(len(us))
        crosscorr_sys = np.zeros((nzbins,len(us)))
        crosscorr_err = np.zeros((nzbins,len(us)))
        for i,u in enumerate(us):
            autocorr_sys[i] = corr[u][A_auto][:][0][0]
            autocorr_sys_err[i] = np.sqrt(cov[u][A_auto+A_auto][i][0][0][0][0])
            for zbin in range(nzbins):
                try:
                    crosscorr_sys[zbin,i] = corr[u][A_cross1][:][0][zbin]
                    crosscorr_err[zbin,i] = np.sqrt(cov[u][A_cross1+A_cross1][i][0][zbin][0][zbin])
                except KeyError, IndexError:
                    crosscorr_sys[zbin,i] = corr[u][A_cross2][:][zbin][0]
                    crosscorr_err[zbin,i] = np.sqrt(cov[u][A_cross2+A_cross2][i][zbin][0][zbin][0])
        
        fig = plt.figure()
        for zbin in range(nzbins):
            z = meta[galaxy_pop]['z_mean'][zbin]
            
            correction = (crosscorr_sys[zbin,:]*crosscorr_sys[zbin,:])/autocorr_sys
            requirement = correction/autocorr_errors[zbin,:]
            
            ttl = 'galaxy autocorr z {0}'.format(z)
            ylab = 'w(theta)'
            wmin = min(autocorr_vals[zbin,:]-autocorr_errors[zbin,:]-correction)
            if wmin < 0:
                wmin *= 1.3
            else:
                wmin *= 0.7
            if wmin<0:
                wmin=0
            wmax = 1.2*(max(autocorr_vals[zbin,:]+autocorr_errors[zbin,:]+correction))
            # ax0 = fig.add_subplot(2,2,1, xlim=[umin,umax], ylim=[wmin,wmax], xlabel=xlab, ylabel=ylab, title=ttl)
            ax0 = fig.add_subplot(2,2,1, xlim=[umin,umax], xlabel=xlab, ylabel=ylab, title=ttl)
            ax0.errorbar(us,autocorr_vals[zbin,:],yerr=autocorr_errors[zbin,:],fmt='.-')
            ax0.plot(us,autocorr_vals[zbin,:]-correction,'--',color='green')
            ax0.set_xscale('log')
            ax0.set_yscale('log', nonposy='clip')
            
            ttl = '{0} autocorr'.format(systematic)
            ylab = 'w(theta)'
            wmin = min(autocorr_sys[:]) - max(autocorr_sys_err[:])
            if wmin < 0:
                wmin *= 1.2
            else:
                wmin *= 0.8
            wmax = 1.2*(max(autocorr_sys)+max(autocorr_sys_err))
            ax0 = fig.add_subplot(2,2,2, xlim=[umin,umax], ylim=[wmin,wmax], xlabel=xlab, ylabel=ylab, title=ttl)
            ax0.errorbar(us,autocorr_sys,yerr=autocorr_sys_err,fmt='.-')
            ax0.set_xscale('log')
            
            ttl = 'gal-{0} crosscorr z {1}'.format(systematic, z)
            ylab = 'w(theta)'
            wmin = min(crosscorr_sys[zbin,:])-max(crosscorr_err[zbin,:])
            if wmin < 0:
                wmin *= 1.2
            else:
                wmin *= 0.8
            wmax = 1.2*(max(crosscorr_sys[zbin,:])+max(crosscorr_err[zbin,:]))
            ax0 = fig.add_subplot(2,2,3, xlim=[umin,umax], ylim=[wmin,wmax], xlabel=xlab, ylabel=ylab, title=ttl)
            ax0.errorbar(us,crosscorr_sys[zbin,:],yerr=crosscorr_err[zbin,:],fmt='.-')
            ax0.set_xscale('log')
            
            ttl = '{0} requirement, gal z {1}'.format(systematic, z)
            ylab= "correction/error"
            wmin = min(requirement)
            if wmin < 0:
                wmin *= 1.2
            else:
                wmin *= 0.8
            wmax = 1.2*max(requirement)
            ax0 = fig.add_subplot(2,2,4, xlim=[umin,umax], ylim=[wmin,wmax], xlabel=xlab, ylabel=ylab, title=ttl)
            ax0.plot(us,requirement,'.-')
            ax0.set_xscale('log')
            
            fig.tight_layout()
            pp.savefig()
            plt.clf()
            
            # for i,u in enumerate(us):
            #     print "%f %f %f %f %f" %(correction[i], autocorr_errors[zbin,i], autocorr_vals[zbin,i], crosscorr_sys[zbin,i], autocorr_sys[i])
    
    pp.close()
    h5file.close()

if __name__ == '__main__':
    main()
