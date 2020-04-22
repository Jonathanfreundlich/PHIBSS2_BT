# AUXILIARY FUNCTIONS FOR THE SINGLE SERSIC FITS ON NOISE-FREE IDEAL DISK+BULGE SYSTEMS (Sersic_ideal.py)

print 'Running Sersic_aux.py'

def create_feedme(feedme_file):
    open(feedme_file,'w')
    
def add_parameters(feedme_file,galfit_input='moffat_psf.fits',galfit_output='model.fits',galfit_sigma='none',psf_file='moffat_psf.fits',constraint_file='none',region=[0,264,0,264],conv_box=[264,264]):
    with open(feedme_file,'a') as f:
        f.write('===============================================================================\n')
        f.write('# IMAGE and GALFIT CONTROL PARAMETERS \n')
        f.write('A) %s                  # Input data image (FITS file) \n'%galfit_input)
        f.write('B) %s                  # Output data image block \n'%galfit_output)
        f.write('C) %s                  # Sigma image name (made from data if blank or "none") \n'%galfit_sigma)
        f.write('D) %s                  # Input PSF image and (optional) diffusion kernel \n'%psf_file)
        f.write('E) 1                   # PSF fine sampling factor relative to data \n')
        f.write('F) none                # Bad pixel mask (FITS image or ASCII coord list) \n')
        f.write('G) %s                  # File with parameter constraints (ASCII file) \n'%constraint_file)
        f.write('H) %i %i %i %i         # Image region to fit (xmin xmax ymin ymax) \n'%(region[0],region[1],region[2],region[3]))
        f.write('I) %i %i               # Size of the convolution box (x y) \n'%(conv_box[0],conv_box[1]))
        f.write('J) 0                   # Magnitude photometric zeropoint \n')
        f.write('K) 0.038  0.038        # Plate scale (dx dy)    [arcsec per pixel] \n')
        f.write('O) regular             # Display type (regular, curses, both) \n')
        f.write('P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps \n')

def add_sersic(feedme_file,position=[132,132],magnitude=0.,Re=1.,nsersic=1.,variable=False):
    with open(feedme_file,'a') as f:
        f.write('\n')
        f.write('# Object: sersic \n')
        f.write(' 0) sersic                 #  Component type\n')
        f.write(' 1) %f %f 0 0              #  position x, y \n'%(position[0],position[1]))
        f.write(' 3) %f %i                  #  Integrated magnitude \n'%(magnitude,int(variable)))
        f.write(' 4) %f %i                  #  R_e (half-light radius) [pix] \n'%(Re,int(variable)))
        f.write(' 5) %f %i                  #  Sersic index n (de Vaucouleurs n=4)\n'%(nsersic,int(variable)))
        f.write(' 6) 0.0000      0          #     ----- \n')
        f.write(' 7) 0.0000      0          #     ----- \n')
        f.write(' 8) 0.0000      0          #     ----- \n')
        f.write(' 9) 1. 0                   #  axis ratio (b/a) \n')
        f.write(' 10) 0. 0                  #  position angle (PA) [deg: Up=0, Left=90] \n')
        f.write(' Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n')
        
def add_sky(feedme_file,variable=False):
    with open(feedme_file,'a') as f:
        f.write('\n')
        f.write('# Object: sky \n')
        f.write(' 0) sky                    #  object type \n')
        f.write(' 1) 0. %i                  # sky background at center of fitting region [ADUs] \n'%int(variable))
        f.write(' 2) 0.0000      0          # dsky/dx (sky gradient in x) \n')
        f.write(' 3) 0.0000      0          # dsky/dy (sky gradient in y) \n')
        f.write(' Z) 0                      # output option (0 = resid., 1 = Do not subtract) \n')

def replace_header(oldfile,newfile,keys,vals,iframe=0):
    old_file=pyfits.open(oldfile)[iframe]
    new_header=old_file.header
    for (key,val) in zip(keys,vals):new_header[key]=(val)
    new_hdu=pyfits.PrimaryHDU(old_file.data,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    new_hdulist.writeto(newfile)
    
def create_uniform_sigma(sigma_file):
    uniform_file='/Users/freundlich/Desktop/PROJETS/PHIBSS2/galfit/Lang2014/models/uniform_sigma.fits'
    dirty_sigma=pyfits.open(sigma_file)[0]
    uniform_sigma=ones(dirty_sigma.data.shape)
    new_header=dirty_sigma.header
    new_header['COMMENT']=('This is a uniform sigma image, created by J. Freundlich')
    new_hdu = pyfits.PrimaryHDU(uniform_sigma,new_header)
    new_hdulist = pyfits.HDUList([new_hdu])
    new_hdulist.writeto(uniform_file)

def remove_logs(galpath):
    for filename in glob.glob('galfit.*'):
        os.remove(filename)
    
def plot_model(galpath,galfit_output):
    # LOAD GALFIT MODEL
    output_list=pyfits.open(galpath+'/'+galfit_output)
    image = output_list[1]
    model = output_list[2]
    residual = output_list[3]
    #
    vmin=1e-4#5.5e-2
    vmax=vmax=image.data.max()
    #
    space_down=0.1
    space_left=0.05
    width=0.21
    height=0.84
	#
    fig, (ax1, ax2, ax3) = subplots(nrows=1, ncols=3,figsize=(24,6))#, sharey=True,sharex=True)
    #
    # IMAGE
    im=ax1.imshow(flipud(image.data),norm=LogNorm(vmin=vmin,vmax=vmax)) 
    ax1.set_title(r'HST map')#,fontsize=fontsize)
    ax1.set_position([space_left,space_down, width, height])
	#
    # MODEL
    ax2.imshow(flipud(model.data),norm=LogNorm(vmin=vmin,vmax=vmax))
    ax2.set_title(r'Model')#,fontsize=fontsize)
    ax2.set_position([space_left+width,space_down, width, height])
	#
    # RESIDUAL
    ax3.imshow(flipud(residual.data),norm=LogNorm(vmin=vmin,vmax=vmax))
    ax3.set_title(r'Residuals')#,fontsize=fontsize)
    ax3.set_position([space_left+2*width,space_down, width, height])
	#
    cbar_ax = fig.add_axes([space_left+3*width,space_down, 0.02, height])
    fig.colorbar(im, cax=cbar_ax)
    cbar_ax.set_ylabel(r'arbitrary')#,fontsize=fontsize)#, rotation=270)     
	#
    show()