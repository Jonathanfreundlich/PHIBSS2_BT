################################################################################
# RUN BULGE-DICK FITS ON ALL GALAXIES OR ON SPECIFIC CASES, METHOD USED IN FREUNDLICH ET AL. 2019, A&A, 622, 105

import os
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits #import pyfits

#from shutil import copyfile
from scipy.interpolate import interp2d
#import pyfits
#import pickle
#import logging
#import os.path
#import copy
#import shutil

script_directory='.'
os.chdir(script_directory)

execfile('./BT_satellites.py')
execfile('./BT_pixel.py')
execfile('./BT_profile.py')
execfile('./BT_colormap.py')
execfile('./BT_header.py')
execfile('./BT_print.py')
execfile('./BT_noise.py')
execfile('./BT_aux.py')
execfile('./BT_figure.py')
execfile('./BT_sigma.py')

ID_list=['L14CO008'] # LIST OF IDS IF MULTIPLE IDS

################################################################################
# STEP 1: DEFINE FUNCTION TO GET B/T GIVEN Rb/Rd AND nsersic FROM SINGLE SERSIC FIT TO DEFINE THE INITIAL GUESSES
# USES THE RESULT FROM execfile('Sersic_ideal.py')
   
[BT,RR,n,dn,R,dR] = np.loadtxt('composite_poisson_psfmean.dat').reshape((6, 10, 11))

BT1=BT.reshape((110))
RR1=RR.reshape((110))
n1=n.reshape((110))
R1=R.reshape((110))

# INTERPOLATION
interp_BT = interp2d(RR1,n1,BT1,kind='linear')
interp_R = interp2d(RR1,BT1,1./R1)

################################################################################
# STEP 2: SINGLE SERSIC, BULGE AND DISK FITS TO THE IMAGES

EXPTIME=1
GAIN=1
NCOMBINE=1
log_galfit=[]
sigma_file='sigma.fits'
psf_type='ACS1_F814W_G2V'

# Create uniform sigma matrix
do_create_sigma=True
if do_create_sigma:
    for ID in ID_list:
        galpath='./EXAMPLE'
        uniform_sigma(galpath+'/'+ID.lower()+'-acs-I-r.fits',galpath+'/sigma.fits')
    
# Carry out single Sersic, single bulge and single disk fits
do_sersic=True
if do_sersic:
    for ID in ID_list:
        galpath='./EXAMPLE'
        copyfile('./'+psf_type+'.fits',galpath+'/'+psf_type+'.fits')
        create_sersic_constraints(galpath,ID)
        create_feedme(galpath,'.',ID,psf_type=psf_type,gtype='sersic',sigma_file=sigma_file)
        launch_galfit(galpath,'sersic',ID,plot_image=True,save_image=True,redo=True,figsize=(16,4),fontsize=12)
        create_feedme(galpath,'.',ID,psf_type=psf_type,gtype='bulge',sigma_file=sigma_file)
        launch_galfit(galpath,'bulge',ID,plot_image=True,save_image=True,redo=True,figsize=(16,4),fontsize=12)
        create_feedme(galpath,'.',ID,psf_type=psf_type,gtype='disk',sigma_file=sigma_file)
        launch_galfit(galpath,'disk',ID,plot_image=True,save_image=True,redo=True,figsize=(16,4),fontsize=12)
        create_composite_constraints(galpath,with_bulge_constraints=True)
	
################################################################################
# STEP 3: TWO COMPONENT FITS

# List all fit types, including with initial conditions starting from the single bulge and single disk fits 
# (composite_disk* and composite_bulge*)
#logging.info(' ')
do_composite=True
redo_composite=True
fit_types=['composite%02.0f'%(RRi*10) for RRi in linspace(0.1,1,10)]
fit_types.append('disk')
fit_types.append('bulge')
for RRi in linspace(0.1,1,10):
    fit_types.append('composite_disk10%02.0f'%(RRi*10))
    fit_types.append('composite_disk50%02.0f'%(RRi*10))
    fit_types.append('composite_bulge50%02.0f'%(RRi*10))
    fit_types.append('composite_bulge90%02.0f'%(RRi*10))
    
# Define and carry out the bulge+disk composite fits
if do_composite:
    create_composite_constraints(galpath,with_bulge_constraints=True)
    for (i,ID) in enumerate(ID_list):
            print ID
            print 'PSF type: ',psf_type
            #logging.info('\n')
            #logging.info('Start ID=%s \n'%ID)
            variable_sky=True
            variable_AR=True
            create_sersic_constraints(galpath,ID,nmax=4)
            create_feedme(galpath,'.',ID,psf_type=psf_type,gtype='sersic',variable=True,variable_satellite=False,sigma_file='sigma.fits',variable_sky=variable_sky,variable_AR=variable_AR)
            launch_galfit(galpath,'sersic',ID,plot_image=True,save_image=True,redo=True,figsize=(16,4),fontsize=16)
            launch_composite(galpath,ID,redo=redo_composite,plot_image=True,save_image=True,figsize=(16,4),fontsize=16,psf_type=psf_type,sigma_file=sigma_file,logfile='./run_galfit_all.log',variable_sky=variable_sky,variable_AR=variable_AR)
            [[typebest,chiredbest,BTbest,Rdbest,Rbbest,mdbest,mbbest,ARdbest,ARbbest,PAdbest,PAbbest],[fit_types,chired,BT,Rd,Rb,md,mb,ARd,ARb,PAd,PAb]]=get_bestfit(galpath,ID,ftypes=fit_types,typechi='CHI2NU')
            #logging.info('Finished ID=%s \n'%ID)
            #logging.info('typebest=%s, BTbest=%.2f \n'%(typebest,BTbest))
            plt.close('all')
        
do_composite_others=False
if do_composite_others:
    #log_galfit.append('composite_disk:')
    for ID in phibss_project:
        for RRi in linspace(0.1,1,10):
            launch_composite(galpath,ID,RRi=RRi,redo=redo_composite,plot_image=True,save_image=True,figsize=(16,4),fontsize=16,psf_type=psf_type,sigma_file=sigma_file,reference='disk')
        #log_galfit.append('   '+ID)
        
    #log_galfit.append('composite_bulge:')
    for ID in phibss_project:
        for RRi in linspace(0.1,1,10):
            launch_composite(galpath,ID,RRi=RRi,redo=redo_composite,plot_image=True,save_image=True,figsize=(16,4),fontsize=16,psf_type=psf_type,sigma_file=sigma_file,reference='bulge')
        #log_galfit.append('   '+ID)                

        
################################################################################