# WRITE SIGMA FROM MODEL

print 'Running BT_sigma.py'
 
def sigma_from_wht(weightimage,sigma_output='sigma.fits'):
    '''
    Create sigma image from HST weight image (fits file)
    The weight map should be equal to the expected inverse variance (i.e., 1/RMS^2) per pixel. 
    '''
    weight=pyfits.open(weightimage)[0]
    data=weight.data
    header=weight.header
    sigma=1./sqrt(data)
    sigma[where(sigma==nan)]=1e10*ones_like(sigma[where(sigma==nan)])
    header['COMMENT']=('Sigma from HST weight (1/sqrt(weight)), created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,header)
    new_hdulist = pyfits.HDUList([new_hdu])
    if os.path.exists(sigma_output): os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output   
    
def uniform_sigma(image_file,sigma_output='sigma.fits',ncomponent=0):
    '''
    Create uniform sigma matrix
    image_file corresponds to the ACS fits image
    '''
    fitsfile=pyfits.open(image_file)[ncomponent]
    data=fitsfile.data
    header=fitsfile.header
    if size(data.shape)==4:
        image_data=data[0,0,:,:]
    elif size(data.shape)==2:
        image_data=data     
    #
    sigma=ones_like(image_data)
    new_header=header
    new_header['COMMENT']=('Uniform sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output,clobber=True)
    print 'A new %s was created or modified'%sigma_output   

def sigma_poisson(image_file,sigma_output='sigma.fits',nthreshold=3,conversion_factor=1):
    '''
    Create sigma matrix corresponding to Poisson noise
    
    '''
    image_full=pyfits.open(image_file)[0]
    data=image_full.data
    header=image_full.header
    data=data[0,0,:,:]
    mean_sky,sigma_sky=get_noise(data,savefigure=False)  
    rms_sky=sqrt(mean_sky**2+sigma_sky**2)
    threshold=nthreshold*rms_sky
    print 'sky_rms = ', rms_sky, ', nthreshold = ', nthreshold
    #
    sigma=ones_like(data)*rms_sky**2
    sigma[where(data>threshold)]+=(data[where(data>threshold)]-threshold)/conversion_factor
    sigma=sqrt(sigma)  
    #
    new_header=header
    new_header['COMMENT']=('Non uniform sigma, created by J. Freundlich')
    new_header['COMMENT']=('sigma^2 = sky_rms^2 below %i*sky_rms'%nthreshold)
    new_header['COMMENT']=('sigma^2 = sky_rms^2 + (signal-%i*sky_rms) above %i*sky_rms '%(nthreshold,nthreshold))
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+galpath+'/'+ID.lower()+'/'+sigma_output)
    new_hdulist.writeto(galpath+'/'+ID.lower()+'/'+sigma_output,clobber=True)
    

def sky_sigma(image_file,sigma_output,ncomponent=0,nthreshold=3.):
    '''
    rms_sky, mean_sky, sigma_sky
    if the sky is Gaussian, rms_sky=sqrt(mean_sky^2+sigma_sky^2)
    
    if data>mean_sky+nthreshold*:
        sigma^2 = (data-mean_sky)+mean_sky+sigma_sky^2
    if data<rms_sky:
        sigma^2 = mean_sky+sigma_sky^2
    '''
    fitsfile=pyfits.open(image_file)[ncomponent]
    data=fitsfile.data
    header=fitsfile.header
    if size(data.shape)==4:
        image_data=data[0,0,:,:]
    elif size(data.shape)==2:
        image_data=data     
    mean_sky,sigma_sky=get_noise(image_data,savefigure=False)    
    data_min=image_data.min()
    rms_sky=sqrt(mean_sky**2+sigma_sky**2)
    #
    sigma2=ones_like(image_data)*(mean_sky+sigma_sky**2)
    sigma2[where(image_data>mean_sky)]+=(image_data[where(image_data>mean_sky)]-mean_sky)
    sigma=sqrt(sigma2)
    #
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('sigma^2=(data-rms_sky)+rms_sky+sigma_sky^2 above rms_sky')
    new_header['COMMENT']=('sigma^2=rms_sky+sigma_sky^2 below')
    new_header['COMMENT']=('Created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output   
    
    
            
def redress_sigma(sigma_file,sigma_output,minval):
    sigma_data=pyfits.open(sigma_file)[0].data
    sigma_data[where(sigma_data<minval)]=minval*ones_like(sigma_data[where(sigma_data<minval)])
    sigma_data[where(sigma_data==1e10)]=minval*ones_like(sigma_data[where(sigma_data==1e10)])
    #
    new_header=pyfits.open(sigma_file)[0].header
    new_header['COMMENT']=('Sigma redressed with a minimum of %f'%minval)
    new_header['COMMENT']=('Created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma_data,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output       
    

    #from scipy.ndimage.filters import *
    #smoothdata=gaussian_filter(data,nsigma)
    
    
def sigma_circular(image_file,sigma_output,center,radius):
    from astropy.io import fits as pyfits
    image=pyfits.open(image_file)[0]
    y,x = np.indices((image.data[0,0,:,:].shape))
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    sigma=ones(image_data.shape)
    sigma[where(r<radius)]=1.E10
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('Circular sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output    

def sigma_outskirts(image_file,sigma_output,center,radius):
    from astropy.io import fits as pyfits
    image=pyfits.open(image_file)[0]
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    y,x = np.indices((image_data.shape))
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    sigma=ones(image_data.shape)
    sigma[where(r<radius)]=1.E10
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('Circular sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output    
    

def sigma_central(image_file,sigma_output,center,radius):
    '''
    image_file=galfit_output
    center=center in galfit_output (not in the initial acs image)
    if neeeded, do shift_center(galfit_output,center_acs)
    '''
    from astropy.io import fits as pyfits
    image=pyfits.open(image_file)[0]
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    y,x = np.indices((image_data.shape))
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    sigma=ones(image_data.shape)
    sigma[where(r>radius)]=1.E10
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('Circular sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output  


def sigma_gn032(image_file,sigma_output):
    from astropy.io import fits as pyfits
    image=pyfits.open(image_file)[0]
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    y,x = np.indices((image_data.shape))
    sigma=ones(image_data.shape)
    # remove spiral
    A1=(180.,322.)#(44.,200-78.)
    A2=(310.,124.)#(102.,200-177.)
    sigma[where(y<A1[1]+(A2[1]-A1[1])/(A2[0]-A1[0])*(x-A1[0]))]=1.E10
    # remove other galaxy
    #center=(181,70)
    #radius=10.
    #r = np.sqrt((x-center[0])**2+(y-center[1])**2)    
    #sigma[where(r<radius)]=1.E10
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('XC54 sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output 


def sigma_noring(image_file,sigma_output,center,radius1,radius2):
    from astropy.io import fits as pyfits
    image=pyfits.open(image_file)[0]
    y,x = np.indices((image.data[0,0,:,:].shape))
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    sigma=ones(image_data.shape)
    sigma[where((r>radius1) & (r<radius2))]=1.E10
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('Circular sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output  

def view_sigma(sigma_file):
    # PLOT SIGMA
    sigma_list=pyfits.open(sigma_file)
    sigma=sigma_list[0]
    #
    figure('sigma')
    clf()
    imshow(flipud(sigma.data))
    ax=gca()
    #ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)
    title('sigma')
    colorbar()

def sigma_xc54(image_file,sigma_output):
    from astropy.io import fits as pyfits
    image=pyfits.open(image_file)[0]
    sigma=ones(image.data[0,0,:,:].shape)
    y,x = np.indices((image.data[0,0,:,:].shape))
    # remove spiral
    A1=(122.7,188.7)
    A2=(197.,271.4)
    A3=(363.2,356.1)
    sigma[where(y>A1[1]+(A2[1]-A1[1])/(A2[0]-A1[0])*(x-A1[0]))]=1.E10
    sigma[where(y>A2[1]+(A3[1]-A2[1])/(A3[0]-A2[0])*(x-A2[0]))]=1.E10
    # remove other galaxy
    center=(241.4,162.1)
    radius=57.
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)    
    sigma[where(r<radius)]=1.E10
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('XC54 sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output 

def sigma_circle(image_file,sigma_output,center,radius):
    from astropy.io import fits as pyfits
    image=pyfits.open(image_file)[0]
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    y,x = np.indices((image_data.shape))
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    sigma=ones(image_data.shape)
    sigma[where(r>radius)]=1.E10
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('Circular sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+sigma_output)
    new_hdulist.writeto(sigma_output)
    print 'A new %s was created or modified'%sigma_output 


def view_masked_image(image_file,galaxy,vmin=1.E10,vmax=1.E10):
    # PLOT SIGMA
    image_list=pyfits.open(image_file)
    image_origin=image_list[0]
    #
    galfit_sigma='%s_sigma.fits'%project[galaxy].lower()
    sigma_list=pyfits.open(galfit_sigma)
    sigma=sigma_list[0]
    #
    galfit_output='%s_model.fits'%project[galaxy].lower()
    output_list=pyfits.open(galfit_output)
    image = output_list[1]
    model = output_list[2]
    limits=re.findall("\d+\d+",model.header['FITSECT'])
    #
    if vmin==1.E10:
        vmin=image.data.min()
    if vmax==1.E10:
        vmax=image.data.max()
    #
    figure('sigma')
    clf()
    if size(image_origin.data.shape)==4:
        image_data=image_origin.data[0,0,:,:]
    elif size(image_origin.data.shape)==2:
        image_data=image_origin.data 
    fig=imshow(flipud(image_data[0,0,:,:]/sigma.data),vmin=vmin,vmax=vmax)
    #contour(flipud(sigma.data),levels=[1])
    axis([int(limits[0]),int(limits[1]),int(limits[3]),int(limits[2])])
    axis('off')
    title('sigma')
    #colorbar()   

def view_masked(galfit_output,galfit_sigma,vmin=1.E10,vmax=1.E10):
    # PLOT SIGMA
    image_list=pyfits.open(galfit_output)
    image_origin=image_list[0]
    #
    sigma_list=pyfits.open(galfit_sigma)
    sigma=sigma_list[0]
    #
    output_list=pyfits.open(galfit_output)
    image = output_list[1]
    model = output_list[2]
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    limits=re.findall("\d+\d+",model.header['FITSECT'])
    #
    if vmin==1.E10:
        vmin=image.data.min()
    if vmax==1.E10:
        vmax=image.data.max()
    #
    figure('sigma')
    clf()
    fig=imshow(flipud(image_data/sigma.data),vmin=vmin,vmax=vmax)
    #contour(flipud(sigma.data),levels=[1])
    #axis([int(limits[0]),int(limits[1]),int(limits[3]),int(limits[2])])
    #axis('off')
    title('sigma')
    #colorbar()  

def view_masked_acs(acs_file,galfit_sigma,vmin=1.E10,vmax=1.E10):
    # PLOT SIGMA
    image_list=pyfits.open(acs_file)
    image=image_list[0]
    #
    sigma_list=pyfits.open(galfit_sigma)
    sigma=sigma_list[0]
    #
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    limits=re.findall("\d+\d+",model.header['FITSECT'])
    #
    if vmin==1.E10:
        vmin=image.data.min()
    if vmax==1.E10:
        vmax=image.data.max()
    #
    figure('sigma')
    clf()
    fig=imshow(flipud(image_data/sigma.data),vmin=vmin,vmax=vmax)
    #contour(flipud(sigma.data),levels=[1])
    #axis([int(limits[0]),int(limits[1]),int(limits[3]),int(limits[2])])
    #axis('off')
    title('sigma')
    #colorbar()  

def sigma_md59(acs_file,galfit_sigma):
    from astropy.io import fits as pyfits
    image=pyfits.open(acs_file)[0]
    if size(image.data.shape)==4:
        image_data=image.data[0,0,:,:]
    elif size(image.data.shape)==2:
        image_data=image.data 
    sigma=ones(image_data.shape)
    y,x = np.indices((image_data.shape))
    # remove spiral
    A1=(10,24.)
    A2=(100.,64)
    sigma[where(y<A1[1]+(A2[1]-A1[1])/(A2[0]-A1[0])*(x-A1[0]))]=1.E10
    #
    new_header=pyfits.open('sigma.fits')[0].header
    new_header['COMMENT']=('XC54 sigma, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+galfit_sigma)
    new_hdulist.writeto(galfit_sigma)
    print 'A new %s was created or modified'%galfit_sigma 


def view_mask_full(acs_file,galfit_sigma,vmin=1.E10,vmax=1.E10):
    # PLOT SIGMA
    image_list=pyfits.open(acs_file)
    image_origin=image_list[0]
    if size(image_origin.data.shape)==4:
        image_data=image_origin.data[0,0,:,:]
    elif size(image_origin.data.shape)==2:
        image_data=image_origin.data 
    #
    sigma_list=pyfits.open(galfit_sigma)
    sigma=sigma_list[0]
    #
    #
    if vmin==1.E10:
        vmin=image_data.min()
    if vmax==1.E10:
        vmax=image_data.max()
    #
    figure('sigma')
    clf()
    fig=imshow(flipud(image_data/sigma.data),vmin=vmin,vmax=vmax)
    #contour(flipud(sigma.data),levels=[1])
    #axis([int(limits[0]),int(limits[1]),int(limits[3]),int(limits[2])])
    axis('off')
    title('sigma')
    #colorbar()   


def remove_nan(acs_file):
    import math
    image_list=pyfits.open(acs_file)
    image_origin=image_list[0]
    if size(image_origin.data.shape)==4:
        image_data=image_origin.data[0,0,:,:]
    elif size(image_origin.data.shape)==2:
        image_data=image_origin.data 
    #
    image_output=image_data
    for i in range(image_output.shape[0]):
        for j in range(image_output.shape[1]):
            if math.isnan(image_output[i,j]):
                image_output[i,j]=-1.E10
    #
    new_file=acs_file[:-5]+'_nonan.fits'
    new_header=pyfits.open(acs_file)[0].header
    new_header['COMMENT']=('NaNs were replaced by 1.E10, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(image_output,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    os.system('rm '+new_file)
    new_hdulist.writeto(new_file)
    print 'A new %s was created or modified'%new_file
