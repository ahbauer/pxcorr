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
from scipy import linalg
from scipy import optimize
import sys
from lib_mbe import h5getattr, h5setattr
import partpix_map


def is_significant(correction, correction_cov, chi2):
    significant = False
    
    if chi2 > 1.0:
        significant = True
    
    # for i in range(len(correction)):
        # if abs(correction[i]) > 1.*np.sqrt(correction_cov[i,i]):
            # significant = True
    
    return significant

def load_theory(theory_choice, z):
    if theory_choice == 'desdm':
        zstrings = {0.3:'0p3', 0.5:'0p5', 0.7:'0p7', 0.9:'0p9', 1.1:'1p1'}
        zindices = {0.3: '0', 0.5:'1', 0.7:'2', 0.9:'3', 1.1:'4'}
        if z not in zstrings.keys():
            print "Warning, no theory for redshift {0}".format(z)
        theory_filename = '/Users/bauer/surveys/DES/sva1/lss_cats/gold/theory/wtheta_Planck_N' + zindices[z] + '-DESiauto22.5_z' + zstrings[z] + '.dat'
        theory_index = 2
    elif theory_choice == 'tpz1':
        theory_filename = '/Users/bauer/surveys/DES/sva1/lss_cats/gold/theory/tpz/N_z_' + str(z-0.1) + '_' + str(z+0.1) + '_tpz.dat'
        theory_index = 1
    elif theory_choice == 'tpz2':
        theory_filename = '/Users/bauer/surveys/DES/sva1/lss_cats/gold/theory/tpz2/wtheta_Planck_zbin_' + str(z-0.1) + '-' + str(z+0.1) + '_bench_iauto.lt.22p5_tpz2_SG.lt.0p001.dat'
        theory_index = 2
    else:
        msg = 'Error: theory choice not understood: {0}'.format(theory_choice)
        raise Exception(msg)
    return theory_filename, theory_index
    
    
def fit_bias_real(corr, cov, theory_bias1, invert_cov=False):
    
    ndf = len(corr)-1
    
    cov_results = None
    if invert_cov:
        cov_results = linalg.inv(cov)
    else:
        cov_results = linalg.lu_factor(cov)
    
    result = optimize.fmin( fit_bias_real_guts, x0=1.0, args=(corr,cov_results,theory_bias1,invert_cov), disp=False )
    b = result[0]
    chi2 = fit_bias_real_guts( b, corr, cov_results, theory_bias1, invert_cov )
    
    # the error is when the chi2 goes up by 1
    eb_high = 0.
    for eb_high in np.arange(0.01,10.0,0.01):
        chi2_high = fit_bias_real_guts( b+eb_high, corr, cov_results, theory_bias1, invert_cov )
        if chi2_high > chi2+1:
            break
    eb_low = 0.
    for eb_low in np.arange(0.01,10.0,0.01):
        chi2_low = fit_bias_real_guts( b-eb_low, corr, cov_results, theory_bias1, invert_cov )
        if chi2_low > chi2+1:
            break
    
    return b, eb_low, eb_high, chi2, ndf


def fit_bias_real_guts(bias, corr, cov_results, theory_bias1, invert_cov=False):
    
    corr_minus_theory = corr - theory_bias1*bias*bias
    
    chi2=-1.0
    if invert_cov:
        chi2 = np.dot(corr_minus_theory.T,np.dot(cov_results,corr_minus_theory))
    else:
        x = linalg.lu_solve(cov_results, corr_minus_theory)
        chi2 = np.sum(corr_minus_theory*x)
    
    return chi2


def main():
    if len(sys.argv) < 2:
        print 'Usage: pxcorr_eval_sys.py hdf5_filename (nosys,allsys)'
        print '       nosys means ignore any systematics and just plot the galaxy autocorr and theory'
        print '       allsys means correct for all systematics, not just the significant ones'
        exit(1)
    
    
    galaxy_strings = ['bench', 'gold']
    do_sys = True
    do_all_sys = False
    ignore_systematic_cross = False
    diagonal_covariance = False
    
    fit_theory = True
    theory_choice = 'tpz2' # 'desdm' 'tpz1' 'tpz2'
    
    print 'Looking for galaxy samples with the strings: {0}'.format(galaxy_strings)
    if ignore_systematic_cross:
        print 'WARNING: Ignoring the cross-correlation between systematics!!'
    if diagonal_covariance:
        print 'WARNING: Only using the diagonal covariance!!'
    
    if fit_theory:
        print 'Using theory predictions for photo-z: {0}'.format(theory_choice)
    
    filename = sys.argv[1]
    print "Reading filename %s" %filename
    h5file = tables.openFile(filename)
    
    if len(sys.argv) > 2 and sys.argv[2] == 'nosys':
        print 'nosys: Ignoring systematics measurements'
        do_sys = False
    
    if len(sys.argv) > 2 and sys.argv[2] == 'allsys':
        print 'allsys: Correcting for all systematics'
        do_all_sys = True
        
    
    # read in the metadata
    # this code originally copied from martin eriksen
    # this block tells us what populations are in the file, and also the angular binning
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
    # and the jackknife samples
    corr = dict((u, {}) for u in meta['u_mean'])
    corr_grp = h5file.getNode('/', 'corr')
    jks = None
    jk_group = None
    
    
    if do_sys:
        jks = dict((u, {}) for u in meta['u_mean'])
        try:
            jk_group = h5file.getNode('/', 'jk')
        except tables.exceptions.NoSuchNodeError:
            print 'WARNING: No jk node in the hdf5 file, so not dealing with any systematics.'
            do_sys = False
            
    for corr_name, corrObj in corr_grp._v_children.iteritems():

        y = lambda x: h5getattr(corrObj, x)
        A = ((y('ftype0'), y('pop0')), (y('ftype1'), y('pop1')))
        attr_corr_list = corrObj.attrs._f_list('user')
        corrArray = corrObj.read()
        
        aShape = (len(meta['u_mean']), meta[y('pop0')]['nbins'], meta[y('pop1')]['nbins']) 
        msg = 'Wrong dimension of the correlation: ',corr_name , 'Expected: '\
              ,str(list(aShape)),'Shape input: ',str(list(corrArray.shape))
        
        assert aShape == corrArray.shape, msg
        
        u_dict = dict()
        for i,u in enumerate(meta['u_mean']):
            corr[u][A] = corrArray[i]
            if do_sys:
                jks[u][A] = dict()
                u_dict[i] = u
    
        # read in the jackknife info
        msg = 'Correlation array name does not start with corr...'
        assert corr_name[0:4] == 'corr', msg
        corr_num = corr_name[4:]
        if do_sys:
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
        
    
    # read in the covariance matrices
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
    
    
    # ok, find the galaxy correlations in this big dataset.
    galaxy_pop = 0
    n_galpops=0
    systematics_pops = []
    for pop in meta['pop']:
        are_gals = False
        for galaxy_string in galaxy_strings:
            if pop.find(galaxy_string) >= 0:
                print "Found the galaxy population %s" %pop
                galaxy_pop = pop
                n_galpops += 1
                are_gals = True
        if not are_gals:
            systematics_pops.append(pop)
    if n_galpops != 1:
        raise AssertionError, "{0} galaxy populations found with possibilities {1}" %(n_galpops, galaxy_strings)
    
    if len(systematics_pops) == 0:
        do_sys = False
    
    if do_sys:
        print "Systematics populations found:"
        print systematics_pops
    else:
        print "Found these systematics but am not using them:"
        print systematics_pops
        systematics_pops = []
    
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
    
    # let's print out the autocorr covariance while we're at it
    for zbin in range(nzbins):
        outfile = open('cov_z' + str(zbin), 'w')
        for i1,u1 in enumerate(us):
            for i2,u2 in enumerate(us):
                outfile.write('{0} {1} {2}\n'.format(i1,i2,cov[u1][AA][i2][zbin][zbin][zbin][zbin]))
    
    # for each systematic, calculate the cross^2 over the auto
    sys_auto_matrix = []
    sys_cross_matrix = []
    significant_jkmaps = []
    systematic_pops = []
    for zbin in range(nzbins):
        significant_jkmaps.append([])
        sys_auto_matrix.append([])
        sys_cross_matrix.append([])
        systematic_pops.append([])
    for systematic in systematics_pops:
    
        sys.stderr.write( ' Starting systematic {0}                              \r'.format(systematic))
        if systematic.find('magerr') >= 0:
            print 'Skipping!'
            continue
        
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
    
        for zbin in range(nzbins):
            z = meta[galaxy_pop]['z_mean'][zbin]
        
            correction = (crosscorr_sys[zbin,:]*crosscorr_sys[zbin,:])/autocorr_sys
            correction_cov = np.zeros((len(us),len(us)))
            requirement = correction/autocorr_errors[zbin,:]
            requirement_error = np.zeros(len(requirement))
        
            # let's calculate the errors on the requirement!
            # as we don't have errors on the autocorr_errors, we'll really calculate the errors on the correction.
            jk_array = []
            for i,u in enumerate(us):
                jkMap_autosys = jks[u][A_auto][i][0][0]
                jkMap_cross = None
                try:
                    jkMap_cross = jks[u][A_cross1][i][0][zbin]
                except KeyError, IndexError:
                    jkMap_cross = jks[u][A_cross2][i][zbin][0]
                jk_autosys_intersection, jk_cross_intersection = jkMap_autosys.intersection(jkMap_cross)
                if jk_autosys_intersection.npartpix < 0.5*jkMap_autosys.npartpix:
                    print "Warning, small jackknife overlap area for {0}! {1} vs {2} pixels".format(A_cross1,jk_autosys_intersection.npartpix, jkMap_autosys.npartpix)
                jk_array.append(np.zeros(jk_autosys_intersection.npartpix))
                for p in range(jk_autosys_intersection.npartpix):
                    # calculate the correction
                    jk_array[i][p] = jk_cross_intersection.partmap[p]*jk_cross_intersection.partmap[p]/jk_autosys_intersection.partmap[p]
                jk_mean = np.mean(jk_array[i])
                ecorr = 0.
                for p in range(jk_autosys_intersection.npartpix):
                    ecorr += (jk_array[i][p]-jk_mean)*(jk_array[i][p]-jk_mean)
                ecorr = (float(jk_autosys_intersection.npartpix-1.0))/float(jk_autosys_intersection.npartpix) * ecorr
                requirement_error[i] = np.sqrt(ecorr)/autocorr_errors[zbin,i]
        
            # do the whole covariance matrix
            for i,u in enumerate(us):
                imean = np.mean(jk_array[i])
                for j,v in enumerate(us):
                    msg = 'jk arrays are different lengths for different angular bins!: {0}:{1}, {2}:{3}'.format(i,len(jk_array[i]),j,len(jk_array[j]))
                    assert len(jk_array[i]) == len(jk_array[j]), msg
                    jmean = np.mean(jk_array[j])
                    for p in range(len(jk_array[i])):
                        correction_cov[i,j] += (jk_array[i][p]-imean)*(jk_array[j][p]-jmean)
                    correction_cov[i,j] *= (float(len(jk_array[i])-1.0))/float(len(jk_array[i]))
            correction_error = np.sqrt(np.diagonal(correction_cov))
            
            # calculate chi2 difference from zero
            cov_inv = np.linalg.inv(correction_cov)
            chi2 = 0.
            for i,u in enumerate(us):
                for j,v in enumerate(us):
                    chi2 += correction[i]*cov_inv[i,j]*correction[j]
            chi2 /= len(us)
            
            # is the systematic significant?
            significant = is_significant(correction, correction_cov, chi2)
            
            if do_all_sys:
                significant = True
        
            # if it's significant, add its contribution to the galaxy correction calculation
            # also, make some plots about this systematic's contribution.
            if significant:
                sys_auto_matrix[zbin].append(autocorr_sys)
                sys_cross_matrix[zbin].append(crosscorr_sys[zbin,:])
                significant_jkmaps[zbin].append((jk_autosys_intersection, jk_cross_intersection))
                systematic_pops[zbin].append(systematic)
                
                fig = plt.figure()
                ttl = 'galaxy autocorr z {0}: {1}'.format(z, significant)
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
            
                ttl = '{0} req., gal z {1}, chi2 {2:.2f}'.format(systematic, z, chi2)
                ylab= "alpha"
                wmin = min(correction)
                if wmin < 0:
                    wmin *= 1.2
                else:
                    wmin *= 0.8
                wmax = 1.2*max(correction)
                ax0 = fig.add_subplot(2,2,4, xlim=[umin,umax], ylim=[wmin,wmax], xlabel=xlab, ylabel=ylab, title=ttl)
                ax0.errorbar(us,correction,yerr=correction_error,fmt='.-')
                ax0.set_xscale('log')
        
                fig.tight_layout()
                pp.savefig()
                plt.close(fig)
                
                # write out info for future plots (paper)
                out_array = np.array([us, autocorr_vals[zbin,:],autocorr_errors[zbin,:],autocorr_sys,autocorr_sys_err,crosscorr_sys[zbin,:],crosscorr_err[zbin,:],correction,correction_error])
                np.savetxt('z{0}_{1}.txt'.format(zbin,systematic),out_array)
                
    
    # ok, now calculate the total correction to the galaxy autocorrelation
    alphas = []
    for zbin in range(nzbins):
        galaxy_corrected = np.zeros(len(us))
        corrected_errs = np.zeros(len(us))
        corrected_cov = np.zeros((len(us),len(us)))
        z = meta[galaxy_pop]['z_mean'][zbin]
        if do_sys:
            print "zbin {0}: {1} significant systematics used: {2}".format(zbin,len(systematic_pops[zbin]), systematic_pops[zbin])
            # crap, we need to solve this for each angular bin!
            alphas_per_angle = []
            for i,u in enumerate(us):
                sys_auto_matrix2 = np.zeros((len(sys_auto_matrix[zbin]),len(sys_auto_matrix[zbin])))
                sys_cross_vector2 = np.zeros(len(sys_auto_matrix[zbin]))
                for s1 in range(len(sys_auto_matrix[zbin])):
                    for s2 in range(len(sys_auto_matrix[zbin])):
                        A_auto1 = (('counts',systematic_pops[zbin][s1]),('counts',systematic_pops[zbin][s2]))
                        A_auto2 = (('counts',systematic_pops[zbin][s2]),('counts',systematic_pops[zbin][s1]))
                        A_cross1 = (('counts',systematic_pops[zbin][s2]),('counts',galaxy_pop))
                        A_cross2 = (('counts',galaxy_pop),('counts',systematic_pops[zbin][s2]))
                        try:
                            sys_auto_matrix2[s1,s2] = corr[u][A_auto1][:][0][0]
                        except KeyError, IndexError:
                            sys_auto_matrix2[s1,s2] = corr[u][A_auto2][:][0][0]
                        try:
                            sys_cross_vector2[s2] = corr[u][A_cross1][:][0][zbin]
                        except KeyError, IndexError:
                            sys_cross_vector2[s2] = corr[u][A_cross2][:][zbin][0]
                if ignore_systematic_cross:
                    sys_auto_matrix2 = np.identity(len(sys_auto_matrix[zbin]))*sys_auto_matrix2
                result = np.linalg.solve(sys_auto_matrix2, sys_cross_vector2)
                alphas_per_angle.append(result)
            alphas_per_angle = np.vstack(alphas_per_angle)
            alphas.append(np.median(alphas_per_angle,axis=1))
            
            # now apply the correction for each jackknife subsample
            jk_corrected = []
            jk_corrected_means = np.zeros(len(us))
            A_auto = (('counts',galaxy_pop),('counts',galaxy_pop))
            for i,u in enumerate(us):
                jk_corrected.append([])
                jkGalAutoMap = jks[u][A_auto][i][zbin][zbin]
                for k in range(jkGalAutoMap.npartpix):
                    correction = 0.
                    for s1 in range(len(sys_auto_matrix[zbin])):
                        A_cross1 = (('counts',galaxy_pop),('counts',systematic_pops[zbin][s1]))
                        A_cross2 = (('counts',systematic_pops[zbin][s1]),('counts',galaxy_pop))
                        jkval = None
                        try:
                            jkCrossMap = jks[u][A_cross1][i][zbin][0]
                        except KeyError, IndexError:
                            jkCrossMap = jks[u][A_cross2][i][0][zbin]
            
                        if jkGalAutoMap.pixel_mapping_arraytohigh[k] not in jkCrossMap.pixel_mapping_arraytohigh:
                            continue
        
                        for s2 in range(len(sys_auto_matrix[zbin])):
                            A_auto1 = (('counts',systematic_pops[zbin][s1]),('counts',systematic_pops[zbin][s2]))
                            A_auto2 = (('counts',systematic_pops[zbin][s2]),('counts',systematic_pops[zbin][s1]))
                            jkAutoMap = None
                            try:
                                jkAutoMap = jks[u][A_auto1][i][0][0]
                            except KeyError, IndexError:
                                jkAutoMap = jks[u][A_auto2][i][0][0]
                
                            if jkGalAutoMap.pixel_mapping_arraytohigh[k] not in jkAutoMap.pixel_mapping_arraytohigh:
                                continue
                
                            if (not ignore_systematic_cross) or (s1 == s2):
                                correction += alphas_per_angle[i][s1]*alphas_per_angle[i][s2]*jkAutoMap.partmap[k]
        
                        correction -= 2.0 * alphas_per_angle[i][s1]*jkCrossMap.partmap[k]
        
                    jk_corrected[i].append(jkGalAutoMap.partmap[k] + correction)
                jk_corrected_means[i] = np.mean(jk_corrected[i])
    
                # and do the overall value
                correction = 0.
                for s1 in range(len(sys_auto_matrix[zbin])):
                    A_cross1 = (('counts',galaxy_pop),('counts',systematic_pops[zbin][s1]))
                    A_cross2 = (('counts',systematic_pops[zbin][s1]),('counts',galaxy_pop))
                    crosscorrelation = None
                    try:
                        crosscorrelation = corr[u][A_cross1][:][zbin][0]
                    except KeyError, IndexError:
                        crosscorrelation = corr[u][A_cross2][:][0][zbin]
        
                    for s2 in range(len(sys_auto_matrix[zbin])):
                        A_auto1 = (('counts',systematic_pops[zbin][s1]),('counts',systematic_pops[zbin][s2]))
                        A_auto2 = (('counts',systematic_pops[zbin][s2]),('counts',systematic_pops[zbin][s1]))
                        autocorrelation = None
                        try:
                            autocorrelation = corr[u][A_auto1][:][0][0]
                        except KeyError, IndexError:
                            autocorrelation = corr[u][A_auto2][:][0][0]
            
                        if (not ignore_systematic_cross) or (s1 == s2):
                            correction += alphas_per_angle[i][s1]*alphas_per_angle[i][s2]*autocorrelation
                    
                    correction -= 2.0 * alphas_per_angle[i][s1]*crosscorrelation
                A_gal = (('counts',galaxy_pop),('counts',galaxy_pop))
                galaxy_corrected[i] = corr[u][A_gal][:][zbin][zbin] + correction

            # calculate new covariance matrix (we're not actually saving this for now)
            for i,u in enumerate(us):
                imean = np.mean(jk_corrected[i])
                for j,v in enumerate(us):
                    jmean = np.mean(jk_corrected[j])
                    msg = 'jk corrected arrays are different lengths for different angular bins!: {0}:{1}, {2}:{3}'.format(i,len(jk_corrected[i]),j,len(jk_corrected[j]))
                    assert len(jk_corrected[i]) == len(jk_corrected[j]), msg
                    for p in range(len(jk_corrected[i])):
                        corrected_cov[i,j] += (jk_corrected[i][p]-imean)*(jk_corrected[j][p]-jmean)
                    corrected_cov[i,j] *= (float(len(jk_corrected[i])-1.0))/float(len(jk_corrected[i]))
                corrected_errs[i] = np.sqrt(corrected_cov[i,i])

            # write out the corrected covariance matrix 
            outfile = open('cov_z' + str(zbin) + '_sys', 'w')
            for i1,u1 in enumerate(us):
                for i2,u2 in enumerate(us):
                    outfile.write('{0} {1} {2}\n'.format(i1,i2,corrected_cov[i1,i2]))
                
    
        else:
            galaxy_corrected = autocorr_vals[zbin,:]
            corrected_errs = autocorr_errors[zbin,:]
            for i,u in enumerate(us):
                for j,v in enumerate(us):
                    corrected_cov[i,j] = cov[u][AA][j][zbin][zbin][zbin][zbin]
        
        # make some plots!
        fig = plt.figure()
        ttl = 'galaxy autocorr z {0} with real correction'.format(z)
        ylab = 'w(theta)'
        ax0 = fig.add_subplot(1,1,1, xlim=[umin,umax], ylim=[1.e-5,0.1],xlabel=xlab, ylabel=ylab, title=ttl)
        vals = autocorr_vals[zbin,:]
        pos_indices = vals>=0.
        neg_indices = np.logical_not(pos_indices)
        ax0.errorbar(us[pos_indices],autocorr_vals[zbin,:][pos_indices],yerr=autocorr_errors[zbin,:][pos_indices],fmt='.',color='blue')
        if not pos_indices.all():
            ax0.errorbar(us[neg_indices],-1.*autocorr_vals[zbin,:][neg_indices],yerr=autocorr_errors[zbin,:][neg_indices],fmt='.',color='cyan')
        pos_indices = galaxy_corrected>=0.
        neg_indices = np.logical_not(pos_indices)
        ax0.errorbar(us[pos_indices],galaxy_corrected[pos_indices],yerr=corrected_errs[pos_indices],fmt='.',color='red')
        if not pos_indices.all():
            ax0.errorbar(us[neg_indices],-1.*galaxy_corrected[neg_indices],yerr=corrected_errs[neg_indices],fmt='.',color='magenta')
    
        if fit_theory:
            min_index = 0
            theory_filename = None
            theory_index = None
            theory_filename, theory_index = load_theory(theory_choice, z)
            
            theory_file = open(theory_filename)
            if os.path.isfile(theory_filename):
                theory = np.loadtxt(theory_file,unpack=True)
                theory_ws = np.array(np.interp( us, theory[0], theory[theory_index], right=0. ))
                valid_indices = np.sum(us<theory[0][-1])
                
                if diagonal_covariance:
                    print 'USING DIAGONAL COVARIANCE ONLY'
                    corrected_cov *= np.identity(len(corrected_cov[0]))
                
                bias, bias_error_low, bias_error_high, chi2, ndf = fit_bias_real( galaxy_corrected[min_index:valid_indices], corrected_cov[min_index:valid_indices,min_index:valid_indices], theory_ws[min_index:valid_indices] )
                print "bias = {0:3.2f} + {1:3.2f} - {2:3.2f}, {3}/{4} chi2/ndf".format(bias, bias_error_high, bias_error_low, chi2, ndf)
                # plot it on the plot
                ax0.plot(us[min_index:valid_indices],theory_ws[min_index:valid_indices]*bias*bias, label='bias={0:3.2f}+/-{1:3.2f}, chi2/ndf={2:4.1f}/{3}'.format(bias,0.5*(bias_error_low+bias_error_high),chi2,ndf))
                low_bias2 = bias - bias_error_low
                low_bias2 *= low_bias2
                high_bias2 = bias + bias_error_high
                high_bias2 *= high_bias2
                theory_ws_low = low_bias2*theory_ws
                theory_ws_high = high_bias2*theory_ws
                # print us
                # print galaxy_corrected
                # print corrected_errs
                ax0.fill_between( us[min_index:valid_indices], theory_ws_low[min_index:valid_indices], theory_ws_high[min_index:valid_indices], alpha=0.2 )
                ax0.legend()
            else:
                print "Warning, file expected but not found for redshift {0}".format(z)
                
        ax0.set_xscale('log')
        ax0.set_yscale('log', nonposy='clip')
        pp.savefig()
        plt.close(fig)
        
        # write some info out to make pretty paper plots
        out_array = np.array([us, autocorr_vals[zbin,:],autocorr_errors[zbin,:],galaxy_corrected,corrected_errs])
        np.savetxt('z{0}_total.txt'.format(zbin,systematic),out_array)
    
    pp.close()
    h5file.close()

if __name__ == '__main__':
    main()
