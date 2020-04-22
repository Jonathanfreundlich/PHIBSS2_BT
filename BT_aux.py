# -*- coding: utf-8 -*-
# AUXILIARY FUNCTIONS FOR BULGE DISK DECOMPOSITION

print 'Running BT_aux.py'

from shutil import copyfile
from scipy.special import gamma
from tempfile import mkstemp
from shutil import move
from os import remove, close

################################################################################
# REPLACE PATTERN OR LINE IN FILE

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def replace_line(file_path,pattern,subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:  
                if pattern in line:
                    new_file.write(subst+'\n')
                else:
                    new_file.write(line)
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def find_line(file_path,pattern):
    output=[]
    with open(file_path) as old_file:
        for line in old_file:
            if pattern in line:
                output.append(line)
    return output
    
################################################################################

def replace_header(oldfile,newfile,keys,vals,iframe=0):
    old_file=pyfits.open(oldfile)[iframe]
    new_header=old_file.header
    for (key,val) in zip(keys,vals):
        new_header[key]=(val)
    
    new_hdu=pyfits.PrimaryHDU(old_file.data,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    new_hdulist.writeto(newfile,clobber=True)

################################################################################

def redress_image(oldfile,newfile,technique='minimum',iframe=0,add_val=0.):
    if technique=='minimum':
        old_file=pyfits.open(oldfile)[iframe]
        old_data=old_file.data
        new_data=old_file.data-array_nonan(old_data.reshape(size(old_data))).min()
        new_header=old_file.header
        new_header['COMMENT']='This image was redressed by J. Freundlich'
        new_header['COMMENT']='Namely by substracting %f'%old_file.data.min()
        new_header['COMMENT']='In order to be all positive'
        new_hdu=pyfits.PrimaryHDU(new_data,new_header)
        new_hdulist = pyfits.HDUList([new_hdu])
        new_hdulist.writeto(newfile,clobber=True)
    elif technique=='none':
        old_file=pyfits.open(oldfile)[iframe]
        old_header=old_file.header
        old_data=old_file.data
        new_hdu=pyfits.PrimaryHDU(old_data,old_header)
        new_hdulist = pyfits.HDUList([new_hdu])
        new_hdulist.writeto(newfile,clobber=True)
    elif technique=='add':
        old_file=pyfits.open(oldfile)[iframe]
        old_data=old_file.data
        new_data=old_file.data+add_val
        new_header=old_file.header
        new_header['COMMENT']='This image was redressed by J. Freundlich'
        new_header['COMMENT']='Namely by adding %f'%add_val
        new_hdu=pyfits.PrimaryHDU(new_data,new_header)
        new_hdulist = pyfits.HDUList([new_hdu])
        new_hdulist.writeto(newfile,clobber=True)    	
         
################################################################################
# CREATE A SPECIFIC DIRECTORY FOR THE GALFIT FIT AND COPY THE ACS IMAGE, THE PSF AND THE SIGMA MATRIX INSIDE
# -- the old path structure is specific to PHIBSS here
# -- the new path is of the form galpath/IDsmall where IDsmall=ID.lower()

def create_galfit_directory(galpath,oldpath,EGN=[],redress=False,prefix='',subdir=''):
    # NOTE : EGN = [EXPTIME,GAIN,NCOMBINE]
    for i in range(size(ID_list)):
        ID=ID_list[i]
        print ID
        IDsmall=ID.lower()
        if not os.path.exists(galpath):
            os.makedirs(galpath)
        
        if not os.path.exists(galpath+'/'+IDsmall):
            os.makedirs(galpath+'/'+IDsmall)
        
        oldfile=oldpath+'/'+prefix+IDsmall+'/'+subdir+IDsmall+'-acs-I-r.fits'
        newfile=galpath+'/'+IDsmall+'/'+IDsmall+'-acs-I-r.fits'
        copyfile(oldfile,newfile)
        
        if size(EGN)==3:
            [EXPTIME,GAIN,NCOMBINE] = EGN
            replace_header(newfile,newfile,['EXPTIME','GAIN','NCOMBINE'],[EXPTIME,GAIN,NCOMBINE])

        if redress:
            redress_image(newfile,newfile)
        
def copy_sigma(galpath,oldpath):
    for i in range(size(ID_list)):
        ID=ID_list[i]
        print ID
        IDsmall=ID.lower()
        try:
            copyfile(oldpath+'/'+IDsmall+'/'+IDsmall+'_sigma.fits',galpath+'/'+IDsmall+'/'+IDsmall+'_sigma.fits')
        except:
            print '--------- No sigma file for %s: uniform sigma instead'%IDsmall
            uniform_sigma(oldpath+'/'+IDsmall+'/'+'sigma.fits',galpath+'/'+IDsmall+'/'+IDsmall+'_sigma.fits')
            
def copy_psf(galpath,psf_type='ACS1_F814W_G2V'):
    psf_file='/Users/jonathanf/Desktop/PHIBSS2/galfit/PSF/%s/%s.fits'%(psf_type,psf_type)
    for i in range(size(ID_list)):
        ID=ID_list[i]
        print ID
        IDsmall=ID.lower()
        copyfile(psf_file,galpath+'/'+IDsmall+'/'+psf_type+'.fits')
 
################################################################################
# CREATE TXT FILES FOR GALFIT (CONSTRAINTS, FEEDME)
# -- the new path is of the form galpath/IDsmall where IDsmall=ID.lower()

def create_sersic_constraints(galpath,ID=[],nmax=8):
    if ID==[]: 
        for i in range(size(ID_list)):
            ID=ID_list[i]
            create_sersic_constraints(galpath+'/'+ID.lower(),ID)
    else:
        print ID
        with open(galpath+'/'+'sersic.constraints','w') as f:
            f.write('# Component/    parameter   constraint  Comment \n')
            f.write('# operation     (see below)   range \n')
            f.write('  \n')
            f.write('    1              n        0.2 to %i     # Soft constraint  \n'%nmax)
            f.write('  \n')

def create_composite_constraints(galpath,with_bulge_constraints=False):
    for i in range(size(ID_list)):
        ID=ID_list[i]
        print ID
        IDsmall=ID.lower()
        with open(galpath+'/'+'composite.constraints','w') as f:
            f.write('# Component/    parameter   constraint  Comment \n')
            f.write('# operation     (see below)   range \n')
            f.write('  \n')
            f.write('    1-2            x        -2 2     # Soft constraint  \n')
            f.write('    1-2            y        -2 2     # Soft constraint  \n')     
            if with_bulge_constraints:
                f.write('    2_1          PA         offset     # Hard constraint  \n') 
                f.write('    2/1          b/a         1 100     # Soft constraint  \n') 
            f.write('  \n')

def create_feedme(galpath,oldpath,ID=[],galfit_input=[],psf_type='ACS1_F814W_G2V',gtype='sersic',sigma_file='none',magnitude=-5.,RR='',interp=True,variable=True,variable_satellite=False,variable_sky=True,variable_AR=True):
    if ID==[]:
        for i in range(size(ID_list)):
            ID=ID_list[i]
            create_feedme(galpath,oldpath,ID=ID,psf_type=psf_type,gtype=gtype,magnitude=magnitude,RR=RR,sigma_file=sigma_file,variable_sky=variable_sky,variable_AR=variable_AR)
    else:
        print ID
        IDsmall=ID.lower()
        
        if variable_sky:variable_sky=variable
        if variable_AR:variable_AR=variable
        
        if isinstance(RR,float):
            strRR='%02.0f'%(RR*10)
        else:
            strRR=''
        
        if galfit_input==[]:
            galfit_input=IDsmall+'-acs-I-r.fits'
            
        new_file=galpath+'/'+gtype+strRR+'.feedme'
        old_file=oldpath+'/'+'sersic'+'.feedme'
        
        try: 
            region_line=find_line(old_file,'H) ')[0]
            center_line=find_line(old_file,'1) ')[0]
        except:
            print '--------- No sersic.feedme file: bulge_disk_I.feedme instead'
            old_file=oldpath+'/'+IDsmall+'/'+'bulge_disk_I.feedme'
            region_line=find_line(old_file,'H) ')[0]
            center_line=find_line(old_file,'1) ')[0]
                    
        region=[int(i) for i in region_line.split()[1:5]]
        if center_line.split()[0]=='#':
            center=[float(i) for i in center_line.split()[2:4]]
        else:
            center=[float(i) for i in center_line.split()[1:3]]
        
        #os.chdir(galpath)
        open(new_file,'w')
        psf_file=psf_type+'.fits'
        print psf_file
        image_size=[region[1]-region[0],region[3]-region[2]]
        psf_size=[pyfits.open(psf_file)[0].header['NAXIS1'],pyfits.open(psf_file)[0].header['NAXIS2']]
        conv_box=[max(image_size[0],psf_size[0]),max(image_size[1],psf_size[1])]
        constraint_file=gtype+'.constraints'
        if gtype in ['composite_disk10','composite_bulge90','composite_disk50','composite_bulge50']:
            constraint_file='composite.constraints'
        add_parameters(new_file,galfit_input=galfit_input,galfit_output=gtype+strRR+'.fits',galfit_sigma=sigma_file,psf_file=psf_file,constraint_file=constraint_file,region=region,conv_box=conv_box)
        
        if gtype<>'sersic':
            output_list=pyfits.open(galpath+'/sersic.fits')
            sersic_model = output_list[2]
            csersic=[float(sersic_model.header['1_XC'].split(' +/- ')[0].replace('[','').replace(']','')),float(sersic_model.header['1_YC'].split(' +/- ')[0].replace('[','').replace(']',''))]
            msersic=float(sersic_model.header['1_MAG'].split(' +/- ')[0].replace('[','').replace(']',''))
            Rsersic=float(sersic_model.header['1_RE'].split(' +/- ')[0].replace('[','').replace(']',''))
            nsersic=float(sersic_model.header['1_N'].split(' +/- ')[0].replace('[','').replace(']',''))
            ARsersic=float(sersic_model.header['1_AR'].split(' +/- ')[0].replace('[','').replace(']',''))
            PAsersic=float(sersic_model.header['1_PA'].split(' +/- ')[0].replace('[','').replace(']',''))
            for N_SKY in ('2_SKY','3_SKY','4_SKY','5_SKY'):
    			try:
    				skysersic=float(sersic_model.header[N_SKY].split(' +/- ')[0].replace('[','').replace(']',''))
    				break
    			except:
    				continue
    	
        if gtype=='sersic':
            add_sersic(new_file,position=center,magnitude=magnitude,Re=10.,nsersic=1.,BA=1.,PA=0.,variable=variable,variable_AR=variable_AR)
        
        elif gtype=='bulge':
            add_bulge(new_file,position=center,magnitude=msersic,Re=Rsersic,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
        
        elif gtype=='disk':
            add_disk(new_file,position=center,magnitude=msersic,Re=Rsersic,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
        
        elif gtype=='composite':
            if interp:
                BT=max(0.,interp_BT(RR,nsersic)[0])
                BT=min(BT,1.)
                Rd=max(0.01,interp_R(RR,BT)[0]*Rsersic)
                Rb=RR*Rd
            else:
                BT=BTval
                
            
            if BT==0:
                BT=0.01
            elif BT==1:
                BT=0.99
                
            md=msersic-2.5*log10(1.-BT)
            mb=msersic-2.5*log10(BT)
            
            add_disk(new_file,position=csersic,magnitude=md,Re=Rd,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            add_bulge(new_file,position=csersic,magnitude=mb,Re=Rb,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            
        elif gtype=='composite_disk10':
            output_list=pyfits.open(galpath+'/'+IDsmall+'/disk.fits')
            disk_model = output_list[2]
            cdisk=[float(disk_model.header['1_XC'].split(' +/- ')[0].replace('[','').replace(']','')),float(disk_model.header['1_YC'].split(' +/- ')[0].replace('[','').replace(']',''))]
            mdisk=float(disk_model.header['1_MAG'].split(' +/- ')[0].replace('[','').replace(']',''))
            Rdisk=float(disk_model.header['1_RE'].split(' +/- ')[0].replace('[','').replace(']',''))
            ndisk=float(disk_model.header['1_N'].split(' +/- ')[0].replace('[','').replace(']',''))
            ARdisk=float(disk_model.header['1_AR'].split(' +/- ')[0].replace('[','').replace(']',''))
            PAdisk=float(disk_model.header['1_PA'].split(' +/- ')[0].replace('[','').replace(']',''))
            
            Rd=Rdisk
            Rb=RR*Rd
            BT=0.1
            md=mdisk-2.5*log10(1.-BT)
            mb=mdisk-2.5*log10(BT)
            
            add_disk(new_file,position=cdisk,magnitude=md,Re=Rd,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            add_bulge(new_file,position=cdisk,magnitude=mb,Re=Rb,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            
        elif gtype=='composite_bulge90':
            output_list=pyfits.open(galpath+'/'+IDsmall+'/bulge.fits')
            disk_model = output_list[2]
            cbulge=[float(disk_model.header['1_XC'].split(' +/- ')[0].replace('[','').replace(']','')),float(disk_model.header['1_YC'].split(' +/- ')[0].replace('[','').replace(']',''))]
            mbulge=float(disk_model.header['1_MAG'].split(' +/- ')[0].replace('[','').replace(']',''))
            Rbulge=float(disk_model.header['1_RE'].split(' +/- ')[0].replace('[','').replace(']',''))
            nbulge=float(disk_model.header['1_N'].split(' +/- ')[0].replace('[','').replace(']',''))
            ARbulge=float(disk_model.header['1_AR'].split(' +/- ')[0].replace('[','').replace(']',''))
            PAbulge=float(disk_model.header['1_PA'].split(' +/- ')[0].replace('[','').replace(']',''))
            
            Rb=Rbulge
            Rd=Rb/RR
            BT=0.9
            md=mbulge-2.5*log10(1.-BT)
            mb=mbulge-2.5*log10(BT)
            
            add_disk(new_file,position=cbulge,magnitude=md,Re=Rd,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            add_bulge(new_file,position=cbulge,magnitude=mb,Re=Rb,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
        
        elif gtype=='composite_disk50':
            output_list=pyfits.open(galpath+'/'+IDsmall+'/disk.fits')
            disk_model = output_list[2]
            cdisk=[float(disk_model.header['1_XC'].split(' +/- ')[0].replace('[','').replace(']','')),float(disk_model.header['1_YC'].split(' +/- ')[0].replace('[','').replace(']',''))]
            mdisk=float(disk_model.header['1_MAG'].split(' +/- ')[0].replace('[','').replace(']',''))
            Rdisk=float(disk_model.header['1_RE'].split(' +/- ')[0].replace('[','').replace(']',''))
            ndisk=float(disk_model.header['1_N'].split(' +/- ')[0].replace('[','').replace(']',''))
            ARdisk=float(disk_model.header['1_AR'].split(' +/- ')[0].replace('[','').replace(']',''))
            PAdisk=float(disk_model.header['1_PA'].split(' +/- ')[0].replace('[','').replace(']',''))
            
            Rd=Rdisk
            Rb=RR*Rd
            BT=0.5
            md=mdisk-2.5*log10(1.-BT)
            mb=mdisk-2.5*log10(BT)
            
            add_disk(new_file,position=cdisk,magnitude=md,Re=Rd,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            add_bulge(new_file,position=cdisk,magnitude=mb,Re=Rb,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            
        elif gtype=='composite_bulge50':
            output_list=pyfits.open(galpath+'/'+IDsmall+'/bulge.fits')
            disk_model = output_list[2]
            cbulge=[float(disk_model.header['1_XC'].split(' +/- ')[0].replace('[','').replace(']','')),float(disk_model.header['1_YC'].split(' +/- ')[0].replace('[','').replace(']',''))]
            mbulge=float(disk_model.header['1_MAG'].split(' +/- ')[0].replace('[','').replace(']',''))
            Rbulge=float(disk_model.header['1_RE'].split(' +/- ')[0].replace('[','').replace(']',''))
            nbulge=float(disk_model.header['1_N'].split(' +/- ')[0].replace('[','').replace(']',''))
            ARbulge=float(disk_model.header['1_AR'].split(' +/- ')[0].replace('[','').replace(']',''))
            PAbulge=float(disk_model.header['1_PA'].split(' +/- ')[0].replace('[','').replace(']',''))
            
            Rb=Rbulge
            Rd=Rb/RR
            BT=0.5
            md=mbulge-2.5*log10(1.-BT)
            mb=mbulge-2.5*log10(BT)
            
            add_disk(new_file,position=cbulge,magnitude=md,Re=Rd,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
            add_bulge(new_file,position=cbulge,magnitude=mb,Re=Rb,BA=ARsersic,PA=PAsersic,variable=variable,variable_AR=variable_AR)
                                                    
        add_satellite(new_file,ID,variable=variable,all_variable=variable_satellite)
        
        if gtype in ['composite','composite_disk','composite_bulge']:
            add_sky(new_file,variable=variable_sky,sky=skysersic)
        else:
            add_sky(new_file,variable=variable_sky)

################################################################################
# LAUNCH GALFIT PROCEDURE
          
def launch_galfit(galpath,gtype='sersic',ID=[],redo=False,plot_image=False,save_image=False,figsize=(24,6),spaces=[0.1,0.05,0.21,0.84],fontsize=20,logfile=[]):
    if ID==[]:
        for i in range(size(ID_list)):
            ID=ID_list[i]
            launch_galfit(galpath,gtype,ID,redo,plot_image,save_image,figsize,spaces,fontsize)
    else:
        print ID
        IDsmall=ID.lower()
        feedme_file=gtype+'.feedme'
        currentdir=os.getcwd()
        os.chdir(galpath) 

        if not os.path.exists(gtype+'.fits') or redo: 
            os.system('galfit '+feedme_file)
            #if logfile<>[]:
            #    logging.basicConfig(filename='/Users/jonathanf/Desktop/PHIBSS2/galfit/Lang2014/python/run_galfit_all.log',level=logging.DEBUG,format='%(message)s')
            #    logging.info('ID=%s: run galfit %s \n'%(ID,gtype))
        else:
            print '%s.fits exists and we do not create it again'%gtype
        
        if plot_image:
            plot_figure('.',ID,gtype,save_image=save_image,figsize=figsize,spaces=spaces,fontsize=fontsize,plot_BT=True)
        
        os.chdir(currentdir)
                                                                                                                                
def launch_composite(galpath,ID=[],RRi=[],redo=False,plot_image=False,save_image=False,figsize=(24,6),spaces=[0.1,0.05,0.21,0.84],fontsize=20,psf_type='ACS1_F814W_G2V',sigma_file='none',reference='sersic',logfile=[],variable_sky=True,variable_AR=True):
    if ID==[]:
        for ID in ID_list:
            print ID
            launch_composite(galpath,ID,RRi=RRi,redo=redo,plot_image=plot_image,save_image=save_image,figsize=figsize,spaces=spaces,fontsize=fontsize,psf_type=psf_type,sigma_file=sigma_file,reference=reference,variable_sky=variable_sky,variable_AR=variable_AR)
    else:
        if RRi==[]:
            for RRi in linspace(0.1,1,10):
                launch_composite(galpath,ID,RRi=RRi,redo=redo,plot_image=plot_image,save_image=save_image,figsize=figsize,spaces=spaces,fontsize=fontsize,psf_type=psf_type,sigma_file=sigma_file,reference=reference,variable_sky=variable_sky,variable_AR=variable_AR)
        else:
            strRRi='%02.0f'%(RRi*10)
            gtype='composite'
            if reference=='disk10': gtype+='_disk10'
            if reference=='bulge90': gtype+='_bulge90'
            if reference=='disk50': gtype+='_disk50'
            if reference=='bulge50': gtype+='_bulge50'
            print gtype
            if not os.path.exists(galpath+'/'+gtype+strRRi+'.fits') or redo:
                create_feedme(galpath,galpath,ID=ID,psf_type=psf_type,gtype=gtype,RR=RRi,sigma_file=sigma_file,variable_sky=variable_sky,variable_AR=variable_AR)
                try:
                    launch_galfit(galpath,gtype+strRRi,ID,plot_image=plot_image,save_image=save_image,redo=redo,figsize=figsize,spaces=spaces,fontsize=fontsize,logfile=logfile)
                except:
                    print 'It seems the fitting crashed'
                plt.close('all')
            else:
                print 'File '+gtype+strRRi+'.fits '+'already exists'
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
################################################################################ 
            
def add_parameters(feedme_file,galfit_input='moffat_psf.fits',galfit_output='model.fits',galfit_sigma='none',mask_file='none',psf_file='moffat_psf.fits',constraint_file='none',region=[0,264,0,264],conv_box=[264,264]):
    with open(feedme_file,'a') as f:
        f.write('===============================================================================\n')
        f.write('# IMAGE and GALFIT CONTROL PARAMETERS \n')
        f.write('A) %s                  # Input data image (FITS file) \n'%galfit_input)
        f.write('B) %s                  # Output data image block \n'%galfit_output)
        f.write('C) %s                  # Sigma image name (made from data if blank or "none") \n'%galfit_sigma)
        f.write('D) %s                  # Input PSF image and (optional) diffusion kernel \n'%psf_file)
        f.write('E) 1                   # PSF fine sampling factor relative to data \n')
        f.write('F) %s                # Bad pixel mask (FITS image or ASCII coord list) \n'%mask_file)
        f.write('G) %s                  # File with parameter constraints (ASCII file) \n'%constraint_file)
        f.write('H) %i %i %i %i         # Image region to fit (xmin xmax ymin ymax) \n'%(region[0],region[1],region[2],region[3]))
        f.write('I) %i %i               # Size of the convolution box (x y) \n'%(conv_box[0],conv_box[1]))
        f.write('J) 0                   # Magnitude photometric zeropoint \n')
        f.write('K) 0.038  0.038        # Plate scale (dx dy)    [arcsec per pixel] \n')
        f.write('O) regular             # Display type (regular, curses, both) \n')
        f.write('P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps \n')

def add_sersic(feedme_file,position=[132,132],magnitude=0.,Re=1.,nsersic=1.,BA=0.,PA=0.,variable=False,variable_AR=False,only_mag_variable=False):
    if only_mag_variable:
        variable=False
        variable_AR=False
    with open(feedme_file,'a') as f:
        f.write('\n')
        f.write('# Object: sersic \n')
        f.write(' 0) sersic                 #  Component type\n')
        f.write(' 1) %f %f %i %i              #  position x, y \n'%(position[0],position[1],int(variable),int(variable)))
        if only_mag_variable:
            f.write(' 3) %f %i                  #  Integrated magnitude \n'%(magnitude,1))
        else:
            f.write(' 3) %f %i                  #  Integrated magnitude \n'%(magnitude,int(variable)))
        f.write(' 4) %f %i                  #  R_e (half-light radius) [pix] \n'%(Re,int(variable)))
        f.write(' 5) %f %i                  #  Sersic index n (de Vaucouleurs n=4)\n'%(nsersic,int(variable)))
        f.write(' 6) 0.0000      0          #     ----- \n')
        f.write(' 7) 0.0000      0          #     ----- \n')
        f.write(' 8) 0.0000      0          #     ----- \n')
        f.write(' 9) %f %i                   #  axis ratio (b/a) \n'%(BA,int(variable_AR)))
        f.write(' 10) %f %i                  #  position angle (PA) [deg: Up=0, Left=90] \n'%(PA,int(variable_AR)))
        f.write(' Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n')

def add_bulge(feedme_file,position=[132,132],magnitude=0.,Re=1.,BA=0.,PA=0.,variable=False,variable_AR=False):
    with open(feedme_file,'a') as f:
        f.write('\n')
        f.write('# Object: bulge \n')
        f.write(' 0) sersic                 #  Component type\n')
        f.write(' 1) %f %f %i %i              #  position x, y \n'%(position[0],position[1],int(variable),int(variable)))
        f.write(' 3) %f %i                  #  Integrated magnitude \n'%(magnitude,int(variable)))
        f.write(' 4) %f %i                  #  R_e (half-light radius) [pix] \n'%(Re,int(variable)))
        f.write(' 5) 4 0                    #  Sersic index n (de Vaucouleurs n=4)\n')
        f.write(' 6) 0.0000      0          #     ----- \n')
        f.write(' 7) 0.0000      0          #     ----- \n')
        f.write(' 8) 0.0000      0          #     ----- \n')
        f.write(' 9) %f %i                   #  axis ratio (b/a) \n'%(BA,int(variable_AR)))
        f.write(' 10) %f %i                  #  position angle (PA) [deg: Up=0, Left=90] \n'%(PA,int(variable_AR)))
        f.write(' Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n')

def add_disk(feedme_file,position=[132,132],magnitude=0.,Re=1.,BA=0.,PA=0.,variable=False,variable_AR=False):
    with open(feedme_file,'a') as f:
        f.write('\n')
        f.write('# Object: disk \n')
        f.write(' 0) sersic                 #  Component type\n')
        f.write(' 1) %f %f %i %i              #  position x, y \n'%(position[0],position[1],int(variable),int(variable)))
        f.write(' 3) %f %i                  #  Integrated magnitude \n'%(magnitude,int(variable)))
        f.write(' 4) %f %i                  #  R_e (half-light radius) [pix] \n'%(Re,int(variable)))
        f.write(' 5) 1 0                    #  Sersic index n (de Vaucouleurs n=4)\n')
        f.write(' 6) 0.0000      0          #     ----- \n')
        f.write(' 7) 0.0000      0          #     ----- \n')
        f.write(' 8) 0.0000      0          #     ----- \n')
        f.write(' 9) %f %i                   #  axis ratio (b/a) \n'%(BA,int(variable_AR)))
        f.write(' 10) %f %i                  #  position angle (PA) [deg: Up=0, Left=90] \n'%(PA,int(variable_AR)))
        f.write(' Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n')
                                      
def add_sky(feedme_file,variable=False,sky=0.):
    with open(feedme_file,'a') as f:
        f.write('\n')
        f.write('# Object: sky \n')
        f.write(' 0) sky                    #  object type \n')
        f.write(' 1) %f %i                  # sky background at center of fitting region [ADUs] \n'%(sky,int(variable)))
        f.write(' 2) 0.0000      0          # dsky/dx (sky gradient in x) \n')
        f.write(' 3) 0.0000      0          # dsky/dy (sky gradient in y) \n')
        f.write(' Z) 0                      # output option (0 = resid., 1 = Do not subtract) \n')

################################################################################
# GET PARAMETERS FROM GALFIT MODEL

def print_rn(filename,scale_kpc):
    output_list=pyfits.open(filename)
    image = output_list[1]
    model = output_list[2]
    try:
            Rs,dRs=re.findall("\d+\.\d+",model.header['1_RE']) # pix
    except:
            Rs=re.findall("\d+\.\d+",model.header['1_RE']) # pix
            dRs=0.
    R12_sec,dR12_sec=get_r12('custom',image,model,ncomponent=1,kpcscale=1,return_error=True)# in arcsec
    R12_kpc,dR12_kpc=get_r12('custom',image,model,ncomponent=1,kpcscale=scale_kpc,return_error=True)# in kpc
    nsersic,dnsersic=get_nsersic(model,ncomponent=1,return_error=True)
    #
    print 'R12     = ', Rs,' +/- ', dRs, ' pixels'
    print 'R12     = ', R12_sec,' +/- ', dR12_sec, ' arcsec'
    print 'R12     = ', R12_kpc,' +/- ', dR12_kpc, ' kpc'
    print 'nsersic = ', nsersic,' +/- ', dnsersic

def get_nsersic(model,ncomponent=1,return_error=False):
    if (model.header['COMP_%i'%ncomponent]=='sersic'):
        try:
            n,dn=re.findall("\d+\.\d+",model.header['%i_N'%ncomponent])
        except:
            n=re.findall("\d+\.\d+",model.header['%i_N'%ncomponent]) # pix
            dn=0.
        if return_error:
            return n,dn   
        else:
            return n  
    else: 
        print 'This is not a Sersic model'
     
def get_sersic(galpath,ID,typechi='CHI2NU'):
    IDsmall=ID.lower()
    output_list=pyfits.open(galpath+'/'+IDsmall+'/'+'sersic.fits')
    sersic_model = output_list[2]
    msersic=float(sersic_model.header['1_MAG'].split(' +/- ')[0])
    Rsersic=float(sersic_model.header['1_RE'].split(' +/- ')[0])
    nsersic=float(sersic_model.header['1_N'].split(' +/- ')[0])
    ARsersic=float(sersic_model.header['1_AR'].split(' +/- ')[0].replace('[','').replace(']',''))
    PAsersic=float(sersic_model.header['1_PA'].split(' +/- ')[0].replace('[','').replace(']',''))
    chisersic=sersic_model.header[typechi]
    for N_SKY in ('2_SKY','3_SKY','4_SKY','5_SKY'):
    	try:
    		skysersic=float(sersic_model.header[N_SKY].split(' +/- ')[0])
    		break
    	except:
    		continue
    		
    return [chisersic,nsersic,Rsersic,msersic,ARsersic,PAsersic]
          
def get_r12(ID,image,model,ncomponent=1,kpcscale=[],return_error=False): # kpc
    try:
        pixel_deg=def_pixel_scale_new(image) #deg
        pixel_sec=pixel_deg[0]*3600. # arcsec
        print 'get_r12 pixel_sec=%.2f'%pixel_sec
    except:
        print 'get_r12 WARNING: WE USE pixel_sec=0.03'
        pixel_sec=nan #0.03 removed in case the telescope is not the HST
    try:
        scale_kpc=phibss_scale[ID.upper()] 
    except:
        print 'get_r12 WARNING: NO PHIBSS_SCALE DEFINED FOR %s'%ID.upper()
        if kpcscale<>[]:
            print 'get_r12 WARNING: WE USE scale_kpc = %.2f FROM USER'%kpcscale
            scale_kpc=kpcscale
        else:
            print 'get_r12 WARNING: scale_kpc = 1'
            scale_kpc=1
    if (model.header['COMP_%i'%ncomponent]=='expdisk'):
        try:
            Rs,dRs=re.findall("\d+\.\d+",model.header['%i_RS'%ncomponent]) # pix
            Rs_sec=float(Rs)*pixel_sec # arcsec
            Rs_kpc=Rs_sec*scale_kpc # kpc
            dRs_sec=float(dRs)*pixel_sec # arcsec
            dRs_kpc=dRs_sec*scale_kpc # kpc    
        except:
            Rs=re.findall("\d+\.\d+",model.header['%i_RS'%ncomponent]) # pix
            Rs_sec=float(Rs[0])*pixel_sec # arcsec
            Rs_kpc=Rs_sec*scale_kpc # kpc
            dRs_sec=0.
            dRs_kpc=0.
        R0=1.678*Rs_kpc
        dR0=1.678*dRs_kpc
    if (model.header['COMP_%i'%ncomponent]=='sersic'):
        try:
            Re,dRe=re.findall("\d+\.\d+",model.header['%i_RE'%ncomponent]) # pix
            Re_sec=float(Re)*pixel_sec # arcsec
            dRe_sec=float(dRe)*pixel_sec # arcsec
            Re_kpc=Re_sec*scale_kpc # kpc
            dRe_kpc=dRe_sec*scale_kpc # kpc
        except:
            Re=re.findall("\d+\.\d+",model.header['%i_RE'%ncomponent]) # pix
            Re_sec=float(Re[0])*pixel_sec # arcsec
            Re_kpc=Re_sec*scale_kpc # kpc
            dRe_sec=0.
            dRe_kpc=0.
        R0=Re_kpc
        dR0=dRe_kpc
    if return_error:
        return R0,dR0
    else:
        return R0

def get_scale(galpath,ID): # kpc/pixel
    galfit_output=galpath+'/'+ID.lower()+'/sersic.fits'
    output_list=pyfits.open(galfit_output)
    image = output_list[1]
    pixel_deg=def_pixel_scale_new(image)[0] #deg
    pixel_sec=pixel_deg*3600. # arcsec
    scale_kpc=phibss_scale[ID] # kpc/arcsec
    return pixel_sec*scale_kpc # kpc/pixel   
            
def get_BT(galpath,ID,RR=''):
    if isinstance(RR,float):
        strRR='%02.0f'%(RR*10)
    else:
        strRR=''
            
    IDsmall=ID.lower()
    output_list=pyfits.open(galpath+'/'+IDsmall+'/composite'+strRR+'.fits')
    model = output_list[2]
    md=float(model.header['1_MAG'].split(' +/- ')[0])
    mb=float(model.header['2_MAG'].split(' +/- ')[0])
    BT=1./(1.+10**((mb-md)/2.5))
    return BT


def Sigma_e(mtot,n,Re,mzero=0.,texp=1.):
    b1=1.67835
    b4=7.669
    if n==1:b=b1
    elif n==4:b=b4
    else:
        print 'ERROR: n should be 1 or 4'
        b=nan
    Ftot=10**(-(mtot-mzero)/2.5)*texp
    Se=Ftot/(2.*pi*Re**2*exp(b)*n*b**(-2.*n)*gamma(2.*n))
    return Se

def Sigma_sersic(mtot,n,Re,mzero=0.,texp=1.):
    Se=Sigma_e(mtot,n,Re,mzero,texp)
    return Sersic1D(amplitude=Se, r_eff=Re, n=n)
    
################################################################################
# COPY ALL IMAGES TO NEW DIRECTORY AND GET BEST FIT

def copy_pdfs(galpath,newdirectory,gtype='sersic',dirname='',typechi='CHI2NU'):
    if not os.path.exists(newdirectory+'/'+'MODELS_'+gtype.upper()):
        os.makedirs(newdirectory+'/'+'MODELS_'+gtype.upper())
    if gtype=='all':
        copy_pdfs(galpath,newdirectory,gtype='sersic',dirname='all')
        copy_pdfs(galpath,newdirectory,gtype='bulge',dirname='all')
        copy_pdfs(galpath,newdirectory,gtype='disk',dirname='all')
        copy_pdfs(galpath,newdirectory,gtype='composite',dirname='all')
        copy_pdfs(galpath,newdirectory,gtype='best',dirname='all')
    else:
        dirname=gtype
        print ' '
        print 'Copy %s files to %s'%(gtype,galpath)+'/'+'MODELS_'+dirname.upper()
        if gtype=='best-redo':
            for ID in ID_list:
                print ID
                IDsmall=ID.lower()
                fit_types=['composite%02.0f'%(RRi*10) for RRi in linspace(0.1,1,10)]
                fit_types.append('disk')
                fit_types.append('bulge')
                for RRi in linspace(0.1,1,10):
                    fit_types.append('composite_disk10%02.0f'%(RRi*10))
                    fit_types.append('composite_disk50%02.0f'%(RRi*10))
                    fit_types.append('composite_bulge50%02.0f'%(RRi*10))
                    fit_types.append('composite_bulge90%02.0f'%(RRi*10))
                [[typebest,chiredbest,BTbest,Rdbest,Rbbest,mdbest,mbbest,ARdbest,ARbbest,PAdbest,PAbbest],[fit_types,chired,BT,Rd,Rb,md,mb,ARd,ARb,PAd,PAb]]=get_bestfit(galpath,ID,ftypes=fit_types,typechi=typechi)
                copyfile(galpath+'/'+IDsmall+'/'+ID+'_'+typebest+'_model.pdf', newdirectory+'/'+'MODELS_'+dirname.upper()+'/'+ID+'_'+typebest+'_model.pdf')
        elif gtype=='best':
            for ID in ID_list:
                print ID
                IDsmall=ID.lower()
                typebest=galfit_type[ID]
                copyfile(galpath+'/'+IDsmall+'/'+ID+'_'+typebest+'_model.pdf', newdirectory+'/'+'MODELS_'+dirname.upper()+'/'+ID+'_'+typebest+'_model.pdf')
        else: 
            for i in range(size(ID_list)):
                ID=ID_list[i]
                print ID
                IDsmall=ID.lower()
                if gtype=='composite':
                    RR=linspace(0.1,1,10)
                    for RRi in RR:
                        strRRi='%02.0f'%(RRi*10)
                        try:
                            copyfile(galpath+'/'+IDsmall+'/'+ID+'_'+gtype+strRRi+'_model.pdf', newdirectory+'/'+'MODELS_'+dirname.upper()+'/'+ID+'_'+gtype+strRRi+'_model.pdf')
                        except:
                            print 'No '+ID+'_'+gtype+strRRi+'_model.pdf'+' found' 
                else:   
                    try:   
                        copyfile(galpath+'/'+IDsmall+'/'+ID+'_'+gtype+'_model.pdf', newdirectory+'/'+'MODELS_'+dirname.upper()+'/'+ID+'_'+gtype+'_model.pdf')
                    except:
                        print 'No '+ID+'_'+gtype+'_model.pdf'+' found' 
                        
def get_bestfit(galpath,ID,with_composite_disk=False,with_composite_bulge=False,ftypes=[],typechi='CHI2NU',size_condition=True,get_uncertainties=False):
    #
    if ftypes==[]:
        RR=linspace(0.1,1,10)
        ftypes=['composite%02.0f'%(RRi*10) for RRi in RR]
        if with_composite_disk:
            for RRi in RR:
                ftypes.append('composite_disk%02.0f'%(RRi*10))
        if with_composite_bulge:
            for RRi in RR:
                ftypes.append('composite_bulge%02.0f'%(RRi*10))
        ftypes.append('disk')
        ftypes.append('bulge')
    #
    IDsmall=ID.lower()
    chired=nan*ones(size(ftypes))
    chi=nan*ones(size(ftypes))
    Rd=nan*ones(size(ftypes))
    Rb=nan*ones(size(ftypes))
    md=nan*ones(size(ftypes))
    mb=nan*ones(size(ftypes))
    BT=nan*ones(size(ftypes))
    ARb=nan*ones(size(ftypes))
    ARd=nan*ones(size(ftypes))
    PAb=nan*ones(size(ftypes))
    PAd=nan*ones(size(ftypes))
    if get_uncertainties:
        dRd=nan*ones(size(ftypes))
        dRb=nan*ones(size(ftypes))
        dmd=nan*ones(size(ftypes))
        dmb=nan*ones(size(ftypes))
        dBT=nan*ones(size(ftypes))
        dARb=nan*ones(size(ftypes))
        dARd=nan*ones(size(ftypes))
        dPAb=nan*ones(size(ftypes))
        dPAd=nan*ones(size(ftypes))
    for (i,fit_type) in zip(range(size(fit_types)),ftypes):
        #print fit_type
        try:
            model = pyfits.open(galpath+'/'+fit_type+'.fits')[2]
            chiredi=model.header['CHI2NU']
            chii=model.header['CHISQ']
            if fit_type=='bulge':
                Rbi=float(model.header['1_RE'].split(' +/- ')[0].replace('*',''))
                Rdi=nan
                mbi=float(model.header['1_MAG'].split(' +/- ')[0].replace('*',''))
                mdi=nan
                BTi=1.
                ARbi=float(model.header['1_AR'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                ARdi=nan
                PAbi=float(model.header['1_PA'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                PAdi=nan
                if get_uncertainties:
                    dRbi=float(model.header['1_RE'].split(' +/- ')[1].replace('*',''))
                    dRdi=nan
                    dmbi=float(model.header['1_MAG'].split(' +/- ')[1].replace('*',''))
                    dmdi=nan
                    dBTi=0.
                    dARbi=float(model.header['1_AR'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dARdi=nan
                    dPAbi=float(model.header['1_PA'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dPAdi=nan

            elif fit_type=='disk':
                Rdi=float(model.header['1_RE'].split(' +/- ')[0].replace('*',''))
                Rbi=nan
                mdi=float(model.header['1_MAG'].split(' +/- ')[0].replace('*',''))
                mbi=nan  
                BTi=0.  
                ARdi=float(model.header['1_AR'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                ARbi=nan
                PAdi=float(model.header['1_PA'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                PAbi=nan
                if get_uncertainties:
                    dRdi=float(model.header['1_RE'].split(' +/- ')[1].replace('*',''))
                    dRbi=nan
                    dmdi=float(model.header['1_MAG'].split(' +/- ')[1].replace('*',''))
                    dmbi=nan
                    dBTi=0.
                    dARdi=float(model.header['1_AR'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dARbi=nan
                    dPAdi=float(model.header['1_PA'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dPAbi=nan
                                    
            elif fit_type=='sersic':
                Rdi=float(model.header['1_RE'].split(' +/- ')[0].replace('*',''))
                Rbi=nan
                mdi=float(model.header['1_MAG'].split(' +/- ')[0].replace('*',''))
                mbi=nan  
                BTi=float(model.header['1_N'].split(' +/- ')[0].replace('*',''))
                ARdi=float(model.header['1_AR'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                ARbi=nan
                PAdi=float(model.header['1_PA'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                PAbi=nan
                if get_uncertainties:
                    dRdi=float(model.header['1_RE'].split(' +/- ')[1].replace('*',''))
                    dRbi=nan
                    dmdi=float(model.header['1_MAG'].split(' +/- ')[1].replace('*',''))
                    dmbi=nan
                    dBTi=float(model.header['1_N'].split(' +/- ')[1].replace('*',''))
                    dARdi=float(model.header['1_AR'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dARbi=nan
                    dPAdi=float(model.header['1_PA'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dPAbi=nan                
            else:            
                Rdi=float(model.header['1_RE'].split(' +/- ')[0].replace('*',''))
                Rbi=float(model.header['2_RE'].split(' +/- ')[0].replace('*',''))
                mdi=float(model.header['1_MAG'].split(' +/- ')[0].replace('*',''))
                mbi=float(model.header['2_MAG'].split(' +/- ')[0].replace('*',''))
                BTi=1./(1.+10**((mbi-mdi)/2.5))
                ARdi=float(model.header['1_AR'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                ARbi=float(model.header['2_AR'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))
                PAdi=float(model.header['1_PA'].split(' +/- ')[0].replace('*','').replace('[','').replace(']',''))  
                PAbi=float(model.header['2_PA'].split(' +/- ')[0].replace('*','').replace('{','').replace('}',''))
                if get_uncertainties:
                    dRdi=float(model.header['1_RE'].split(' +/- ')[1].replace('*',''))
                    dRbi=float(model.header['2_RE'].split(' +/- ')[1].replace('*',''))
                    dmdi=float(model.header['1_MAG'].split(' +/- ')[1].replace('*',''))
                    dmbi=float(model.header['2_MAG'].split(' +/- ')[1].replace('*',''))
                    dBTi=BTi**2*log(10)*10**((mbi-mdi)/2.5)*sqrt(dmdi**2+dmbi**2)
                    dARdi=float(model.header['1_AR'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dARbi=float(model.header['2_AR'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dPAdi=float(model.header['1_PA'].split(' +/- ')[1].replace('*','').replace('[','').replace(']',''))
                    dPAbi=nan
                
                if Rbi<0.1:
                    chiredi=nan
                    chii=nan
                if size_condition and Rbi>Rdi:
                    chiredi=nan
                    chii=nan                                

        except:
            Rbi=nan
            Rdi=nan
            mbi=nan
            mdi=nan
            BTi=nan
            ARdi=nan
            ARbi=nan
            PAdi=nan
            PAbi=nan
            chiredi=nan
            chii=nan
            if get_uncertainties:
                dRdi=nan
                dRbi=nan
                dmdi=nan
                dmbi=nan
                dBTi=nan
                dARdi=nan
                dARbi=nan
                dPAdi=nan
                dPAbi=nan

        Rb[i]=Rbi
        Rd[i]=Rdi
        mb[i]=mbi
        md[i]=mdi
        BT[i]=BTi
        ARd[i]=ARdi
        ARb[i]=ARbi
        PAd[i]=PAdi
        PAb[i]=PAbi
        chired[i]=chiredi
        chi[i]=chii
        if get_uncertainties:
            dRb[i]=dRbi
            dRd[i]=dRdi
            dmb[i]=dmbi
            dmd[i]=dmdi
            dBT[i]=dBTi
            dARd[i]=dARdi
            dARb[i]=dARbi
            dPAd[i]=dPAdi
            dPAb[i]=dPAbi

    if typechi=='CHI2NU':
        ibest=where(chired==min(array_nonan(chired)))[0][0] 
    elif typechi=='CHISQ':
        ibest=where(chi==min(array_nonan(chi)))[0][0] 
    chiredbest=chired[ibest]    
    Rdbest=Rd[ibest]
    Rbbest=Rb[ibest]
    mdbest=md[ibest]
    mbbest=mb[ibest]
    BTbest=BT[ibest]
    ARdbest=ARd[ibest]
    ARbbest=ARb[ibest]
    PAdbest=PAd[ibest]
    PAbbest=PAb[ibest]
    if get_uncertainties:
        dRdbest=dRd[ibest]
        dRbbest=dRb[ibest]
        dmdbest=dmd[ibest]
        dmbbest=dmb[ibest]
        dBTbest=dBT[ibest]
        dARdbest=dARd[ibest]
        dARbbest=dARb[ibest]
        dPAdbest=dPAd[ibest]
        dPAbbest=dPAb[ibest]
    
    typebest=ftypes[ibest]
    
    if get_uncertainties:
        return [typebest,chiredbest,BTbest,dBTbest,Rdbest,dRdbest,Rbbest,dRbbest,mdbest,dmdbest,mbbest,dmbbest,ARdbest,dARdbest,ARbbest,dARbbest,PAdbest,dPAdbest,PAbbest,dPAbbest]
    else:
        return [[typebest,chiredbest,BTbest,Rdbest,Rbbest,mdbest,mbbest,ARdbest,ARbbest,PAdbest,PAbbest],[ftypes,chired,BT,Rd,Rb,md,mb,ARd,ARb,PAd,PAb]]

################################################################################
# CREATE TEX FILE FOR FIGURE

def create_figtex(filename,add_spectra=False):
    lines_of_text=[]
    for i in range(size(ID_list)):
        ID=ID_list[i]
        typebest=galfit_type[ID]
        if i % 7 == 0.:
            lines_of_text.append(r"\begin{figure*}[h!]")
            if i <> 0.:
                lines_of_text.append(r"\ContinuedFloat")
            lines_of_text.append(r"\flushleft")
            lines_of_text.append(r"\includegraphics[height=0.175\textwidth,clip,trim=0 0 4cm 0]{MODELS_BEST/%s_%s_model.pdf}" %(ID,typebest))
            if add_spectra:
                lines_of_text.append(r"\hfill\includegraphics[height=0.175\textwidth,clip,trim=0 0 32cm 0]{/Users/jonathanf/Desktop/PHIBSS2/ARTICLE/python/plots/spectra/%s-spectra.pdf}"%ID)
            lines_of_text.append(r"\\")
        else:
            lines_of_text.append(r"\includegraphics[height=0.1623\textwidth,clip,clip,trim=0 0 4cm 1.1cm]{MODELS_BEST/%s_%s_model.pdf}" %(ID,typebest))
            if add_spectra:
                lines_of_text.append(r"\hfill\includegraphics[height=0.1623\textwidth,clip,trim=0 0 32cm 1.1cm]{/Users/jonathanf/Desktop/PHIBSS2/ARTICLE/python/plots/spectra/%s-spectra.pdf}"%ID)
            lines_of_text.append(r"\\")
            if i % 7 == 6.:
                if i==6:
                    lines_of_text.append(r"\label{fig:images}")
                    lines_of_text.append(r"\caption{\captiontext}")
                lines_of_text.append(r"\end{figure*}")
                lines_of_text.append(r" ")

    lines_of_text.append(r"\end{figure*}")
    lines_of_text.append(r" ")
           

    f = open(filename, 'w')
    for i in lines_of_text:
        f.write(i+ '\n')
    f.close()

################################################################################
