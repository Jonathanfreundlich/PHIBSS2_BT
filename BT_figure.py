# PLOT FIGURE AS IN FREUNDLICH ET AL. 2019 (Fig. A1)

print 'Running BT_figure.py'

from matplotlib.colors import LogNorm

def plot_figure(galpath,ID=[],gtype=[],view_center=True,view_orientation=True,view_ellipse=True,view_r12=False,view_scale=True,view_satellite=True,cross_sec=0.5,ccolor='red',vmin=1e-4,midpoint=0.6,fontsize=20,num_figure=[],save_image=False,save_path='',figsize=(24,6),spaces=[0.14,0.05,0.195,0.78],plot_BT=False,resize_xyl=[],indicate_BT=False,indicate_chi=False,indicate_galname=False):
    #view_center=True         # View center as a cross
    #view_orientation=True    # View orientation
    #view_ellipse=True        # View size ellipse
    #view_r12=False            # View the size on the profile
    #view_scale=True          # View length scale
    #cross_sec=0.5 # Size of the cross on the figure in arcsec
    #ccolor='red' # Color of the comments on the figure
    #vmin=1e-4#5.5e-2
    #midpoint=0.6 # For the Grey color scale  
    #fontsize=20
    
    [space_down,space_left,width,height]=spaces

    if ID==[]:
        for i in range(size(phibss_project)):
            ID=phibss_project[i]
            plot_figure(galpath,ID,gtype,view_center,view_orientation,view_ellipse,view_r12,view_scale,view_satellite,cross_sec,ccolor,vmin,midpoint,fontsize,num_figure,save_image,figsize=figsize,plot_BT=plot_BT,resize_xyl=resize_xyl,indicate_BT=indicate_BT,indicate_chi=indicate_chi,indicate_galname=indicate_galname)
            plt.close('all')
    else:
        if gtype==[]:
            for gtype in ['sersic','bulge','disk']+['composite%02.0f'%i for i in linspace(0.1,1,10)*10]:
                plot_figure(galpath,ID,gtype,view_center,view_orientation,view_ellipse,view_r12,view_scale,view_satellite,cross_sec,ccolor,vmin,midpoint,fontsize,num_figure,save_image,figsize=figsize,plot_BT=plot_BT,resize_xyl=resize_xyl,indicate_BT=indicate_BT,indicate_chi=indicate_chi,indicate_galname=indicate_galname)
                plt.close('all')
        else:
            IDsmall=ID.lower()
            i=where(array(ID_list)==ID)[0][0]
            galfit_output=galpath+'/'+gtype+'.fits'
            
            # LOAD GALFIT MODEL
            output_list=pyfits.open(galfit_output)
            image = output_list[1]
            model = output_list[2]
            residual = output_list[3]
            
            # OBTAIN THE PIXEL SCALE FROM IMAGE
            pixel_deg=def_pixel_scale(image,ID) #deg
            pixel_sec=pixel_deg*3600. # arcsec
            try:
                scale_kpc=phibss_scale[ID.upper()] 
            except:
                print 'WARNING: NO PHIBSS_SCALE DEFINED FOR %s'%ID.upper()
                scale_kpc=1
                    
            # GET THE CENTER OF THE GALAXY FROM MODEL
            ncomponent=1
            header=get_header_parameters(galfit_output,ncomponent)
            center=header['center']
            #print center
            PA=header['PA']
            AR=header['AR']
            PA_rad=header['PA_rad']
            RMAX=header['RMAX']
            
            # OFFSET
            offset=median(image.data)
            
            # DETERMINE RADIAL PROFILES
            if size(image.data.shape)==4:
                image_data=image.data[0,0,:,:]-offset
                radial_image=radial_profile(image_data,center)
            elif size(image.data.shape)==2:
                image_data=image.data-offset
                radial_image=radial_profile(image_data,center)
            else:
                print 'ERROR: image.data.shape = ', image.data.shape
            model_data=model.data-offset
            radial_model=radial_profile(model_data,center)
            radius=range(size(radial_image))
            r12=get_r12(ID,image,model,1)/(pixel_sec*scale_kpc) # in pixels
        
            try:
                r_positive=min(where((radial_model<0.) | (radial_image<0.))[0])-1
            except:
                r_positive=RMAX
            
            size_cross=cross_sec/pixel_sec # pixels
            vmax=image_data.max()
    
            # DEFINE COLOR SCALE
            cmap = matplotlib.cm.get_cmap('Greys')
            shifted_cmap=shift_colormap(cmap,start=0.,midpoint=midpoint,stop=1.,name='shifted')
        
            if view_orientation:
                xx=linspace(0,image.header['NAXIS1'])
                ya=center[1]-1./tan(PA_rad)*(xx-center[0])
                ya=image.header['NAXIS2']-ya
                yb=center[1]+tan(PA_rad)*(xx-center[0])
                yb=image.header['NAXIS2']-yb
            
            if resize_xyl<>[]:
                if isnan(resize_xyl[0]): resize_xyl[0]=center[0]
                if isnan(resize_xyl[1]): resize_xyl[1]=image.header['NAXIS2']-center[1] 
                resize_xyl[2]= min(resize_xyl[2],image.header['NAXIS2'],image.header['NAXIS1'],2*resize_xyl[0],2*resize_xyl[1],2*(image.header['NAXIS1']-resize_xyl[0]),2*(image.header['NAXIS2']-resize_xyl[1]))    
                if isnan(resize_xyl[3]): resize_xyl[3]=6.*r12     
               
            if view_scale:
                scale_dist=0.05
                scale_x=scale_dist*image.header['NAXIS1']
                scale_length=10./(scale_kpc*pixel_sec)/image.header['NAXIS2']
                scale_y1=image.header['NAXIS2']*(1.-scale_dist)
                scale_y2=image.header['NAXIS2']*(1.-scale_dist-scale_length)
                
                if resize_xyl<>[]:
                    scale_x=scale_dist*resize_xyl[2]+(resize_xyl[0]-resize_xyl[2]/2)
                    scale_length=10./(scale_kpc*pixel_sec)/resize_xyl[2]
                    scale_y1=resize_xyl[2]*(1.-scale_dist)+(resize_xyl[1]-resize_xyl[2]/2)
                    scale_y2=resize_xyl[2]*(1.-scale_dist-scale_length)+(resize_xyl[1]-resize_xyl[2]/2)
            
            fig, (ax1, ax2, ax3) = subplots(nrows=1, ncols=3,figsize=figsize, sharey=True,sharex=True)

            # PLOT 1: IMAGE
            data=image_data.copy()
            data[numpy.isnan(data)]=0.
            data[where(data<=0.)]=vmin
            ax1.set_position([space_left,space_down, width, height])
            im=ax1.imshow(flipud(data),cmap=shifted_cmap,norm=LogNorm(vmin=vmin,vmax=vmax)) 
            ax1.set_title(r'$\rm HST$ $\rm map$',fontsize=fontsize,y=1.02)
            if view_center:
                ax1.annotate('', xy=(center[0],image.header['NAXIS2']-center[1]-size_cross), xycoords='data',xytext=(center[0],image.header['NAXIS2']-center[1]+size_cross), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                ax1.annotate('', xy=(center[0]-size_cross,image.header['NAXIS2']-center[1]), xycoords='data',xytext=(center[0]+size_cross,image.header['NAXIS2']-center[1]), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                #ax1.plot(center[0],image.header['NAXIS2']-center[1],'or')
                ax1.axis([0,image.header['NAXIS1'],image.header['NAXIS2'],0])
            if view_orientation:
                ax1.plot(xx,ya,'--',color=ccolor)
                ax1.plot(xx,yb,'--',color=ccolor)
                
            if view_ellipse:
                ####
                theta=linspace(0,2*pi,1000)
                x_ellipse=r12*cos(theta)
                y_ellipse=AR*r12*sin(theta)
                angle=PA_rad+pi/2
                x_ellipse_rot=x_ellipse*cos(angle)-y_ellipse*sin(angle)
                y_ellipse_rot=x_ellipse*sin(angle)+y_ellipse*cos(angle)
                #
                x_ellipse=center[0]+x_ellipse
                y_ellipse=center[1]+y_ellipse
                x_ellipse_rot=center[0]+x_ellipse_rot
                y_ellipse_rot=center[1]+y_ellipse_rot
                y_ellipse=image.header['NAXIS2']-y_ellipse
                y_ellipse_rot=image.header['NAXIS2']-y_ellipse_rot
                #
                ax1.plot(x_ellipse_rot,y_ellipse_rot,color='r')
                ####
                # IF THERE IS A SECOND COMPONENT
                try: 
                    r12_b=get_r12(ID,image,model,2)/(pixel_sec*scale_kpc)
                    AR_b=get_header_parameters(galfit_output,2)['AR']
                    PA_rad_b=get_header_parameters(galfit_output,2)['PA_rad']
                    center_b=get_header_parameters(galfit_output,2)['center']
                    angle_b=PA_rad_b+pi/2
                    x_ellipse_b=r12_b*cos(theta)
                    y_ellipse_b=AR_b*r12_b*sin(theta)
                    x_ellipse_rot_b=x_ellipse_b*cos(angle_b)-y_ellipse_b*sin(angle_b)
                    y_ellipse_rot_b=x_ellipse_b*sin(angle_b)+y_ellipse_b*cos(angle_b)
                    x_ellipse_b=center_b[0]+x_ellipse_b
                    y_ellipse_b=center_b[1]+y_ellipse_b
                    x_ellipse_rot_b=center_b[0]+x_ellipse_rot_b
                    y_ellipse_rot_b=center_b[1]+y_ellipse_rot_b
                    y_ellipse_b=image.header['NAXIS2']-y_ellipse_b
                    y_ellipse_rot_b=image.header['NAXIS2']-y_ellipse_rot_b
                    #
                    ax1.plot(x_ellipse_rot_b,y_ellipse_rot_b,color='b',linestyle='--')
                    ax2.plot(x_ellipse_rot_b,y_ellipse_rot_b,color='b',linestyle='--')
                    ax3.plot(x_ellipse_rot_b,y_ellipse_rot_b,color='b',linestyle='--')
                    #print x_ellipse_rot_b,y_ellipse_rot_b
                except:
                    print 'No second component'
                    
            if view_satellite:
                if gtype in ['sersic','disk','bulge']:
                    nsat=2
                else:
                    nsat=3
                
                print "satellite"
                while model.header['COMP_%i'%nsat][:4]<>'sky':
                	rsat=get_r12(ID,image,model,nsat)/(pixel_sec*scale_kpc) # in pixels
                	hsat=get_header_parameters(galfit_output,nsat)
                	csat=hsat['center']
                	PAsat=hsat['PA']
                	ARsat=hsat['AR']
                	PAsat_rad=hsat['PA_rad']
                	RMAXsat=hsat['RMAX']
                	####
                	thetasat=linspace(0,2*pi,1000)
                	x_ellipsesat=rsat*cos(thetasat)
                	y_ellipsesat=ARsat*rsat*sin(thetasat)
                	anglesat=PAsat_rad+pi/2
                	x_ellipse_rotsat=x_ellipsesat*cos(anglesat)-y_ellipsesat*sin(anglesat)
                	y_ellipse_rotsat=x_ellipsesat*sin(anglesat)+y_ellipsesat*cos(anglesat)
                	#
                	x_ellipsesat=csat[0]+x_ellipsesat
                	y_ellipsesat=csat[1]+y_ellipsesat
                	x_ellipse_rotsat=csat[0]+x_ellipse_rotsat
                	y_ellipse_rotsat=csat[1]+y_ellipse_rotsat
                	y_ellipsesat=image.header['NAXIS2']-y_ellipsesat
                	y_ellipse_rotsat=image.header['NAXIS2']-y_ellipse_rotsat 
                	ax1.plot(x_ellipse_rotsat,y_ellipse_rotsat,color='g')
                	nsat=nsat+1
            if view_scale:
                ax1.axvline(scale_x,ymin=scale_dist,ymax=scale_dist+scale_length, color='k',linewidth=2)
                ax1.axhline(y=scale_y1, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
                ax1.axhline(y=scale_y2, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
            
            ax1.axis([0,image.header['NAXIS1'],image.header['NAXIS2'],0])
            if resize_xyl<>[]:
                ax1.axis([resize_xyl[0]-resize_xyl[2]/2,resize_xyl[0]+resize_xyl[2]/2,resize_xyl[1]+resize_xyl[2]/2,resize_xyl[1]-resize_xyl[2]/2])
            
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            
            if indicate_galname:
                ax1.text(-0.14, 0.,r'#%i   %s'%(i+1,ID),horizontalalignment='right',verticalalignment='bottom',rotation='vertical',transform=ax1.transAxes, fontsize=30)
                ax1.text(-0.05, 0.,r'%s %s'%(phibss_field[ID],phibss_ID1[ID]),horizontalalignment='right',verticalalignment='bottom',rotation='vertical',transform=ax1.transAxes, fontsize=30)
            else:
                ax1.text(-0.05, 0.5,ID,horizontalalignment='right',verticalalignment='center',rotation='vertical',transform=ax1.transAxes, fontsize=30)


            if plot_BT:
                if gtype=='sersic':
                    BT=nan
                elif gtype=='bulge':
                    BT=1
                elif gtype=='disk':
                    BT=0
                else:
                    md=float(model.header['1_MAG'].split(' +/- ')[0].replace('*',''))
                    mb=float(model.header['2_MAG'].split(' +/- ')[0].replace('*',''))
                    BT=1./(1.+10**((mb-md)/2.5))
                
                if indicate_BT:
                    ax1.text(0.05,0.92,'B/T = %.2f'%BT,horizontalalignment='left',verticalalignment='center',transform=ax1.transAxes,fontsize=fontsize)
                if indicate_chi:
                    ax1.text(0.05,0.88,'chi2red = %.2E'%model.header['CHI2NU'],horizontalalignment='left',verticalalignment='center',transform=ax1.transAxes,fontsize=fontsize)
                                                       
            # PLOT 2: MODEL
            data=model_data.copy()
            data[numpy.isnan(data)]=0.
            data[where(data<=0.)]=abs(data[where(data<>0.)]).min()
        
            ax2.imshow(flipud(data),cmap=shifted_cmap,norm=LogNorm(vmin=vmin,vmax=vmax))
            ax2.set_position([space_left+width,space_down, width, height])
            ax2.set_title(r'$\rm Model$',fontsize=fontsize,y=1.02)
            if view_center:
                ax2.annotate('', xy=(center[0],image.header['NAXIS2']-center[1]-size_cross), xycoords='data',xytext=(center[0],image.header['NAXIS2']-center[1]+size_cross), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                ax2.annotate('', xy=(center[0]-size_cross,image.header['NAXIS2']-center[1]), xycoords='data',xytext=(center[0]+size_cross,image.header['NAXIS2']-center[1]), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
            if view_orientation:
                ax2.plot(xx,ya,'--',color=ccolor)
                ax2.plot(xx,yb,'--',color=ccolor)
            if view_ellipse:
                ax2.plot(x_ellipse_rot,y_ellipse_rot,color='r')
            if view_satellite:
                if gtype in ['sersic','disk','bulge']:
                    nsat=2
                else:
                    nsat=3
                
                print "satellite"
                while model.header['COMP_%i'%nsat][:4]<>'sky':
                	rsat=get_r12(ID,image,model,nsat)/(pixel_sec*scale_kpc) # in pixels
                	hsat=get_header_parameters(galfit_output,nsat)
                	csat=hsat['center']
                	PAsat=hsat['PA']
                	ARsat=hsat['AR']
                	PAsat_rad=hsat['PA_rad']
                	RMAXsat=hsat['RMAX']
                	####
                	thetasat=linspace(0,2*pi,1000)
                	x_ellipsesat=rsat*cos(thetasat)
                	y_ellipsesat=ARsat*rsat*sin(thetasat)
                	anglesat=PAsat_rad+pi/2
                	x_ellipse_rotsat=x_ellipsesat*cos(anglesat)-y_ellipsesat*sin(anglesat)
                	y_ellipse_rotsat=x_ellipsesat*sin(anglesat)+y_ellipsesat*cos(anglesat)
                	#
                	x_ellipsesat=csat[0]+x_ellipsesat
                	y_ellipsesat=csat[1]+y_ellipsesat
                	x_ellipse_rotsat=csat[0]+x_ellipse_rotsat
                	y_ellipse_rotsat=csat[1]+y_ellipse_rotsat
                	y_ellipsesat=image.header['NAXIS2']-y_ellipsesat
                	y_ellipse_rotsat=image.header['NAXIS2']-y_ellipse_rotsat 
                	ax2.plot(x_ellipse_rotsat,y_ellipse_rotsat,color='g')
                	nsat=nsat+1
            if view_scale:
                ax2.axvline(scale_x,ymin=scale_dist,ymax=scale_dist+scale_length, color='k',linewidth=2)
                ax2.axhline(y=scale_y1, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
                ax2.axhline(y=scale_y2, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
            ax2.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)

            # PLOT 3: RESIDUAL
            data=residual.data.copy()
            data[numpy.isnan(data)]=0.
            data[where(data<=0.)]=vmin
            
            ax3.imshow(flipud(data),cmap=shifted_cmap,norm=LogNorm(vmin=vmin,vmax=vmax))
            ax3.set_position([space_left+2*width,space_down, width, height])
            ax3.set_title(r'$\rm Residual$',fontsize=fontsize,y=1.02)
            if view_center:
                ax3.annotate('', xy=(center[0],image.header['NAXIS2']-center[1]-size_cross), xycoords='data',xytext=(center[0],image.header['NAXIS2']-center[1]+size_cross), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                ax3.annotate('', xy=(center[0]-size_cross,image.header['NAXIS2']-center[1]), xycoords='data',xytext=(center[0]+size_cross,image.header['NAXIS2']-center[1]), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
            if view_orientation:
                ax3.plot(xx,ya,'--',color=ccolor)
                ax3.plot(xx,yb,'--',color=ccolor)
            if view_ellipse:
                ax3.plot(x_ellipse_rot,y_ellipse_rot,color='r')
            if view_satellite:
                if gtype in ['sersic','disk','bulge']:
                    nsat=2
                else:
                    nsat=3
                
                print "satellite"
                while model.header['COMP_%i'%nsat][:4]<>'sky':
                	rsat=get_r12(ID,image,model,nsat)/(pixel_sec*scale_kpc) # in pixels
                	hsat=get_header_parameters(galfit_output,nsat)
                	csat=hsat['center']
                	PAsat=hsat['PA']
                	ARsat=hsat['AR']
                	PAsat_rad=hsat['PA_rad']
                	RMAXsat=hsat['RMAX']
                	####
                	thetasat=linspace(0,2*pi,1000)
                	x_ellipsesat=rsat*cos(thetasat)
                	y_ellipsesat=ARsat*rsat*sin(thetasat)
                	anglesat=PAsat_rad+pi/2
                	x_ellipse_rotsat=x_ellipsesat*cos(anglesat)-y_ellipsesat*sin(anglesat)
                	y_ellipse_rotsat=x_ellipsesat*sin(anglesat)+y_ellipsesat*cos(anglesat)
                	#
                	x_ellipsesat=csat[0]+x_ellipsesat
                	y_ellipsesat=csat[1]+y_ellipsesat
                	x_ellipse_rotsat=csat[0]+x_ellipse_rotsat
                	y_ellipse_rotsat=csat[1]+y_ellipse_rotsat
                	y_ellipsesat=image.header['NAXIS2']-y_ellipsesat
                	y_ellipse_rotsat=image.header['NAXIS2']-y_ellipse_rotsat 
                	ax3.plot(x_ellipse_rotsat,y_ellipse_rotsat,color='g')
                	nsat=nsat+1
            if view_scale:
                ax3.axvline(scale_x,ymin=scale_dist,ymax=scale_dist+scale_length, color='k',linewidth=2)
                ax3.axhline(y=scale_y1, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
                ax3.axhline(y=scale_y2, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
            ax3.xaxis.set_visible(False)
            ax3.yaxis.set_visible(False)
            tick_params(axis='both', labelsize=20)
            
            # COLOR    BAR   
            cbar_ax = fig.add_axes([space_left+3*width,space_down, 0.02, height])
            fig.colorbar(im, cax=cbar_ax)
            cbar_ax.set_ylabel(r'$\rm arbitrary$',fontsize=fontsize)#, rotation=270)     
            tick_params(labelsize=20)
            
            # LIGHT PROFILE
            profile_ax=fig.add_axes([space_left+3*width+0.1,space_down, width, height])
            profile_ax.plot(array(radius)*pixel_sec*scale_kpc,radial_image,'r',linewidth=2,label=r'$\rm data$')
            profile_ax.plot(array(radius)*pixel_sec*scale_kpc,radial_model,'b--',linewidth=2,label=r'$\rm model$')
            if view_r12:
                i_r12=interp(r12,array(radius)*pixel_sec*scale_kpc,radial_model)
                profile_ax.axvline(r12,ymin=0.,ymax=(log10(i_r12)-log10(vmin))/(log10(vmax)-log10(vmin)), color='k')
                profile_ax.axhline(y=i_r12, xmin=0, xmax=r12/(RMAX*pixel_sec*scale_kpc), color='k')
            profile_ax.set_xlabel(r'$r$ $\rm [kpc]$',fontsize=fontsize)
            if resize_xyl<>[]:
                profile_ax.axis([0,min(RMAX,resize_xyl[3]/pixel_sec/scale_kpc,r_positive)*pixel_sec*scale_kpc,vmin,vmax])
            else:
                profile_ax.axis([0,min(RMAX,r12*6,r_positive)*pixel_sec*scale_kpc,vmin,vmax])
            profile_ax.set_yscale('log')
            tick_params(labelsize=20)
            profile_ax.legend(loc='upper right',fontsize=fontsize-4)
            profile_ax.set_title(r'$\rm Light$ $\rm profile$',fontsize=fontsize,y=1.02)
                    
            if save_image:
                print 'saving file '+ID+'_'+gtype+'_model.pdf'
                if save_path=='':
                    savefig(galpath+'/'+ID+'_'+gtype+'_model.pdf')
                else:
                    savefig(save_path+'/'+ID+'_'+gtype+'_model.pdf')

def plot_fig_model(galfit_output,gtype=[],view_center=True,view_orientation=True,view_ellipse=False,view_r12=False,view_scale=False,view_satellite=False,cross_sec=0.5,ccolor='red',vmin=1e-4,midpoint=0.6,fontsize=20,num_figure=[],save_image=False,savefile='model.pdf',figsize=(24,6),spaces=[0.1,0.05,0.21,0.84],plot_BT=False,kpcscale=[],scale_length=10):
    '''  
    view_center=True         # View center as a cross
    view_orientation=True    # View orientation
    view_ellipse=True        # View size ellipse
    view_r12=False            # View the size on the profile
    view_scale=True          # View length scale
    cross_sec=0.5 # Size of the cross on the figure in arcsec
    ccolor='red' # Color of the comments on the figure
    vmin=1e-4#5.5e-2
	midpoint=0.6 # For the Grey color scale  
    fontsize=20
    '''
    [space_down,space_left,width,height]=spaces
    if gtype==[]:
            for gtype in ['sersic','bulge','disk']+['composite%02.0f'%i for i in linspace(0.1,1,10)*10]:
                plot_figure(galfit_output,gtype,view_center,view_orientation,view_ellipse,view_r12,view_scale,view_satellite,cross_sec,ccolor,vmin,midpoint,fontsize,num_figure,save_image,figsize=figsize,plot_BT=plot_BT)
                plt.close('all')
    else:
            # LOAD GALFIT MODEL
            output_list=pyfits.open(galfit_output)
            image = output_list[1]
            model = output_list[2]
            residual = output_list[3]
            
            # OBTAIN THE PIXEL SCALE FROM IMAGE
            pixel_deg=def_pixel_scale_new(image)[0] #deg
            pixel_sec=pixel_deg*3600. # arcsec
            if kpcscale<>[]:
                scale_kpc=kpcscale
            else:
                print 'WARNING: kpc_scale =1'
                scale_kpc=1
                    
            # GET THE CENTER OF THE GALAXY FROM MODEL
            ncomponent=1
            header=get_header_parameters(galfit_output,ncomponent)
            center=header['center']
            PA=header['PA']
            AR=header['AR']
            PA_rad=header['PA_rad']
            RMAX=header['RMAX']
            
            # OFFSET
            offset=median(image.data)
            
            # DETERMINE RADIAL PROFILES
            if size(image.data.shape)==4:
                image_data=image.data[0,0,:,:]-offset
                radial_image=radial_profile(image_data,center)
            elif size(image.data.shape)==2:
                image_data=image.data-offset
                radial_image=radial_profile(image_data,center)
            else:
                print 'ERROR: image.data.shape = ', image.data.shape
            model_data=model.data-offset
            radial_model=radial_profile(model_data,center)
            radius=range(size(radial_image))
        
            try:
                r_positive=min(where((radial_model<0.) | (radial_image<0.))[0])-1
            except:
                r_positive=RMAX
            
            size_cross=cross_sec/pixel_sec # pixels
            vmax=image_data.max()
    
            # DEFINE COLOR SCALE
            cmap = matplotlib.cm.get_cmap('Greys')
            shifted_cmap=shift_colormap(cmap,start=0.,midpoint=midpoint,stop=1.,name='shifted')
        
            if view_orientation:
                xx=linspace(0,image.header['NAXIS1'])
                ya=center[1]-1./tan(PA_rad)*(xx-center[0])
                ya=image.header['NAXIS2']-ya
                yb=center[1]+tan(PA_rad)*(xx-center[0])
                yb=image.header['NAXIS2']-yb
            
            r12=get_r12('custom',image,model,1,scale_kpc)/(pixel_sec*scale_kpc) # in pixels
            print 'r12 = ', r12, ' pixels'
            if view_ellipse:
                ####
                theta=linspace(0,2*pi,1000)
                x_ellipse=r12*cos(theta)
                y_ellipse=AR*r12*sin(theta)
                angle=PA_rad+pi/2
                x_ellipse_rot=x_ellipse*cos(angle)-y_ellipse*sin(angle)
                y_ellipse_rot=x_ellipse*sin(angle)+y_ellipse*cos(angle)
                #
                x_ellipse=center[0]+x_ellipse
                y_ellipse=center[1]+y_ellipse
                x_ellipse_rot=center[0]+x_ellipse_rot
                y_ellipse_rot=center[1]+y_ellipse_rot
                y_ellipse=image.header['NAXIS2']-y_ellipse
                y_ellipse_rot=image.header['NAXIS2']-y_ellipse_rot
        
            if ID not in ['XA53','XE53','L14EG014']:
                view_satellite=False # as there is no satellite
            if view_satellite:
                if gtype in ['sersic','disk','bulge']:
                    nsat=2
                else:
                    nsat=3
                
                rsat=get_r12(ID,image,model,nsat,scale_kpc)/(pixel_sec*scale_kpc) # in pixels
                hsat=get_header_parameters(galfit_output,nsat)
                csat=hsat['center']
                PAsat=hsat['PA']
                ARsat=hsat['AR']
                PAsat_rad=hsat['PA_rad']
                RMAXsat=hsat['RMAX']
                ####
                thetasat=linspace(0,2*pi,1000)
                x_ellipsesat=rsat*cos(thetasat)
                y_ellipsesat=ARsat*rsat*sin(thetasat)
                anglesat=PAsat_rad+pi/2
                x_ellipse_rotsat=x_ellipsesat*cos(anglesat)-y_ellipsesat*sin(anglesat)
                y_ellipse_rotsat=x_ellipsesat*sin(anglesat)+y_ellipsesat*cos(anglesat)
                #
                x_ellipsesat=csat[0]+x_ellipsesat
                y_ellipsesat=csat[1]+y_ellipsesat
                x_ellipse_rotsat=csat[0]+x_ellipse_rotsat
                y_ellipse_rotsat=csat[1]+y_ellipse_rotsat
                y_ellipsesat=image.header['NAXIS2']-y_ellipsesat
                y_ellipse_rotsat=image.header['NAXIS2']-y_ellipse_rotsat           
               
            if view_scale:
                #scale_length=10.#kpc
                scale_dist=0.05
                scale_x=scale_dist*image.header['NAXIS1']
                scale_length=scale_length/(scale_kpc*pixel_sec)/image.header['NAXIS2']
                scale_y1=image.header['NAXIS2']*(1.-scale_dist)
                scale_y2=image.header['NAXIS2']*(1.-scale_dist-scale_length)
            
            if num_figure<>[]: 
                fig=figure(num_figure)
            else: 
                fig=figure()
            fig, (ax1, ax2, ax3) = subplots(nrows=1, ncols=3,figsize=figsize, sharey=True,sharex=True)

            # PLOT 1: IMAGE
            data=image_data.copy()
            data[numpy.isnan(data)]=0.
            data[where(data<=0.)]=vmin
            ax1.set_position([space_left,space_down, width, height])
            im=ax1.imshow(flipud(data),cmap=shifted_cmap,norm=LogNorm(vmin=vmin,vmax=vmax)) 
            if view_center:
                ax1.annotate('', xy=(center[0],image.header['NAXIS2']-center[1]-size_cross), xycoords='data',xytext=(center[0],image.header['NAXIS2']-center[1]+size_cross), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                ax1.annotate('', xy=(center[0]-size_cross,image.header['NAXIS2']-center[1]), xycoords='data',xytext=(center[0]+size_cross,image.header['NAXIS2']-center[1]), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                ax1.axis([0,image.header['NAXIS1'],image.header['NAXIS2'],0])
            if view_orientation:
                ax1.plot(xx,ya,'--',color=ccolor)
                ax1.plot(xx,yb,'--',color=ccolor)
                ax1.axis([0,image.header['NAXIS1'],image.header['NAXIS2'],0])
            if view_ellipse:
                ax1.plot(x_ellipse_rot,y_ellipse_rot,color='r')
            if view_satellite:
                ax1.plot(x_ellipse_rotsat,y_ellipse_rotsat,color='g')
            if view_scale:
                ax1.axvline(scale_x,ymin=scale_dist,ymax=scale_dist+scale_length, color='k',linewidth=2)
                ax1.axhline(y=scale_y1, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
                ax1.axhline(y=scale_y2, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
            
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.text(0.03,0.92,r'$\rm HST map$',fontsize=fontsize+4,fontweight='bold',transform=ax1.transAxes)
            
            if plot_BT:
                if gtype=='sersic':
                    BT=nan
                elif gtype=='bulge':
                    BT=1
                elif gtype=='disk':
                    BT=0
                else:
                    md=float(model.header['1_MAG'].split(' +/- ')[0])
                    mb=float(model.header['2_MAG'].split(' +/- ')[0])
                    BT=1./(1.+10**((mb-md)/2.5))
                
                ax1.text(0.05,0.95,'B/T = %.2f'%BT,horizontalalignment='left',verticalalignment='center',transform=ax2.transAxes,fontsize=fontsize+4)
                  
            # PLOT 2: MODEL
            data=model_data.copy()
            data[numpy.isnan(data)]=0.
            data[where(data<=0.)]=abs(data[where(data<>0.)]).min()
        
            ax2.imshow(flipud(data),cmap=shifted_cmap,norm=LogNorm(vmin=vmin,vmax=vmax))
            ax2.set_position([space_left+width,space_down, width, height])
            if view_center:
                ax2.annotate('', xy=(center[0],image.header['NAXIS2']-center[1]-size_cross), xycoords='data',xytext=(center[0],image.header['NAXIS2']-center[1]+size_cross), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                ax2.annotate('', xy=(center[0]-size_cross,image.header['NAXIS2']-center[1]), xycoords='data',xytext=(center[0]+size_cross,image.header['NAXIS2']-center[1]), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
            if view_orientation:
                ax2.plot(xx,ya,'--',color=ccolor)
                ax2.plot(xx,yb,'--',color=ccolor)
            if view_ellipse:
                ax2.plot(x_ellipse_rot,y_ellipse_rot,color='r')
            if view_satellite:
                ax2.plot(x_ellipse_rotsat,y_ellipse_rotsat,color='g')
            if view_scale:
                ax2.axvline(scale_x,ymin=scale_dist,ymax=scale_dist+scale_length, color='k',linewidth=2)
                ax2.axhline(y=scale_y1, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
                ax2.axhline(y=scale_y2, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
            ax2.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)
            ax2.text(0.03,0.92,r'$\rm Model$',fontsize=fontsize+4,fontweight='bold',transform=ax2.transAxes)
      
            # PLOT 3: RESIDUAL
            data=residual.data.copy()
            data[numpy.isnan(data)]=0.
            data[where(data<=0.)]=vmin
            
            ax3.imshow(flipud(data),cmap=shifted_cmap,norm=LogNorm(vmin=vmin,vmax=vmax))
            ax3.set_position([space_left+2*width,space_down, width, height])
            if view_center:
                ax3.annotate('', xy=(center[0],image.header['NAXIS2']-center[1]-size_cross), xycoords='data',xytext=(center[0],image.header['NAXIS2']-center[1]+size_cross), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
                ax3.annotate('', xy=(center[0]-size_cross,image.header['NAXIS2']-center[1]), xycoords='data',xytext=(center[0]+size_cross,image.header['NAXIS2']-center[1]), textcoords='data',arrowprops=dict(arrowstyle='-',color=ccolor,linewidth=0.5))
            if view_orientation:
                ax3.plot(xx,ya,'--',color=ccolor)
                ax3.plot(xx,yb,'--',color=ccolor)
            if view_ellipse:
                ax3.plot(x_ellipse_rot,y_ellipse_rot,color='r')
            if view_satellite:
                ax3.plot(x_ellipse_rotsat,y_ellipse_rotsat,color='g')
            if view_scale:
                ax3.axvline(scale_x,ymin=scale_dist,ymax=scale_dist+scale_length, color='k',linewidth=2)
                ax3.axhline(y=scale_y1, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
                ax3.axhline(y=scale_y2, xmin=scale_dist-0.01, xmax=scale_dist+0.01, color='k',linewidth=2)
            ax3.xaxis.set_visible(False)
            ax3.yaxis.set_visible(False)
            ax3.text(0.03,0.92,r'$\rm Residuals$',fontsize=fontsize+4,fontweight='bold',transform=ax3.transAxes)
            
            # COLOR    BAR   
            cbar_ax = fig.add_axes([space_left+3*width,space_down, 0.02, height])
            fig.colorbar(im, cax=cbar_ax)
            cbar_ax.set_ylabel(r'arbitrary',fontsize=fontsize)#, rotation=270)     
            
            # LIGHT PROFILE
            profile_ax=fig.add_axes([space_left+3*width+0.08,space_down, width, height])
            profile_ax.plot(array(radius)*pixel_sec*scale_kpc,radial_image,'r',linewidth=2,label=r'data')
            profile_ax.plot(array(radius)*pixel_sec*scale_kpc,radial_model,'b--',linewidth=2,label=r'model')
            if view_r12:
                i_r12=interp(r12,array(radius)*pixel_sec*scale_kpc,radial_model)
                profile_ax.axvline(r12,ymin=0.,ymax=(log10(i_r12)-log10(vmin))/(log10(vmax)-log10(vmin)), color='k')
                profile_ax.axhline(y=i_r12, xmin=0, xmax=r12/(RMAX*pixel_sec*scale_kpc), color='k')
            profile_ax.set_xlabel(r'$r$ $\rm [kpc]$',fontsize=fontsize+4)
            profile_ax.axis([0,min(RMAX,r12*6,r_positive)*pixel_sec*scale_kpc,vmin,vmax])
            profile_ax.set_yscale('log')
            profile_ax.legend(loc='upper right',fontsize=fontsize)
            profile_ax.text(0.03,0.92,r'$\rm Light$ $\rm profile$',fontsize=fontsize+4,fontweight='bold',transform=profile_ax.transAxes)    
                    
            if save_image:
                savefig(savefile)
            