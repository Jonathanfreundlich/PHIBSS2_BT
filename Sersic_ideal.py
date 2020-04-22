# -*- coding: utf-8 -*-
# SINGLE SERSIC FITS ON NOISE-FREE IDEAL DISK+BULGE SYSTEMS
# (PRELIMINARY STEP FOR THE BULGE+DISK DECOMPOSITION CARRIED OUT IN FREUNDLICH ET AL. 2019, A&A, 622, 105)
# -- For each Rbulge/Rdisk and B/T ratios, compute the Sersic fit parameters

import os
import re
import glob
from astropy.io import fits as pyfits

# Sersic parameters bn
b1=1.67835
b4=7.669

galpath='./MOCK'
os.chdir(galpath)

execfile('../Sersic_aux.py')

#######################################################

feedme_ideal='ideal.feedme'
feedme_sersic='sersic.feedme'
sigma_file='none' 
galfit_output='sersic.fits'
directories_are_set=True
psf_file='none'
savefile=True

Rd=20.
md=0.

noise_mean=0 #ELECTRONS
noise_rms=1e-3 # ELECTRONS
EXPTIME=1
GAIN=1
NCOMBINE=1

BT=linspace(0,1,11)
RR=linspace(0.1,1,10)

# Bulge magnitude mb if md = 0
mb=-2.5*log10(BT/(1.-BT))

BT_array=zeros((size(RR),size(BT)))
RR_array=zeros((size(RR),size(BT)))
n_array=zeros((size(RR),size(BT)))
dn_array=zeros((size(RR),size(BT)))
R_array=zeros((size(RR),size(BT)))
dR_array=zeros((size(RR),size(BT)))

for i in range(size(RR)):
    for j in range(size(BT)):
        print ' '
        print 'Rb/Rd = ', RR[i], 'B/T = ', BT[j],'...'
        BT_array[i,j]=BT[j]
        RR_array[i,j]=RR[i]

        Rb=RR[i]*Rd
        
        if BT[j]==1:
            mb=0.
        else: 
            mb=-2.5*log10(BT[j]/(1.-BT[j]))
        
        # CREATE IDEAL GALAXY
        print '  - create ideal galaxy'
        os.remove('ideal_input.fits')
        replace_header('moffat_psf.fits','ideal_input.fits',['EXPTIME','GAIN','NCOMBINE'],[EXPTIME,GAIN,NCOMBINE])
        create_feedme(feedme_ideal)
        add_parameters(feedme_ideal,galfit_input='ideal_input.fits',galfit_output='ideal.fits',psf_file=psf_file,galfit_sigma=sigma_file)
        if BT[j]<1:
            add_sersic(feedme_ideal,magnitude=md,Re=Rd,nsersic=1.)
        if BT[j]>0:
            add_sersic(feedme_ideal,magnitude=mb,Re=Rb,nsersic=4.)
        
        add_sky(feedme_ideal,variable=False)
        
        os.system('galfit'+' -outsig '+feedme_ideal) 
        
        # SINGLE SERSIC FIT
        print '  - single sersic fit'
        os.remove('ideal_input.fits')
        replace_header('ideal.fits','ideal_input.fits',['EXPTIME','GAIN','NCOMBINE'],[EXPTIME,GAIN,NCOMBINE],iframe=2)
        create_feedme(feedme_sersic)
        add_parameters(feedme_sersic,galfit_input='ideal_input.fits',galfit_output='sersic.fits',psf_file=psf_file,galfit_sigma=sigma_file)
        add_sersic(feedme_sersic,variable=True,magnitude=0.,Re=Rd,nsersic=1.)
        add_sky(feedme_sersic,variable=True)
     
        os.system('galfit -skyped %f -skyrms %f '%(noise_mean,noise_rms)+'-outsig '+feedme_sersic) 
        
        # Retrieve nsersic
        print '  - retrieve nsersic and Rtot'
        output_list=pyfits.open(galfit_output)
        model = output_list[2]
        nsersic=float(re.findall("\d+\.\d+",model.header['1_N'])[0])
        dnsersic=float(re.findall("\d+\.\d+",model.header['1_N'])[1])
        Rsersic=float(re.findall("\d+\.\d+",model.header['1_RE'])[0])/Rd
        dRsersic=float(re.findall("\d+\.\d+",model.header['1_RE'])[1])/Rd
        
        n_array[i,j]=nsersic
        dn_array[i,j]=dnsersic
        R_array[i,j]=Rsersic
        dR_array[i,j]=dRsersic
        
        print '  - summary: '
        print '       Rb/Rd      = ', RR[i]
        print '       B/T        = ', BT[j]
        print '       nsersic    = ', nsersic, ' +/- ', dnsersic
        print '       Rsersic/Rd = ', Rsersic, ' +/- ', dRsersic
        print ' '

remove_logs(galpath)

data=array([BT_array,RR_array,n_array,dn_array,R_array,dR_array])
slice_name=['B/T','Rb/Rd','nsersic','dnsersic','Rsersic/Rd','dRsersic/Rd']

if savefile:
    if psf_file=='none':
        fname='composite_%s'%psf_file
    else:
        fname='composite_%s'%psf_file[:-5]
    if sigma_file<>'none':
        fname=fname+'_'+sigma_file.split('/')[-1][:-5]
    filename='../%s.dat'%fname
    with file(filename, 'w') as outfile:
        outfile.write('# Array shape: {0}\n'.format(data.shape))
        outfile.write('# slices: BT, Rb/Rd, nsersic, dnsersic, Rsersic/Rd, dRsersic/Rd\n')
        for (slice_2d,name) in zip(data,slice_name):
            outfile.write('# Slice %s \n'%name)
            savetxt(outfile, slice_2d)
        
# TO READ THE DATA (FOR EXAMPLE)
#[BT,RR,n,dn] = np.loadtxt(galpath+'/composite_poisson_psfmean.dat').reshape((4, 10, 11))
