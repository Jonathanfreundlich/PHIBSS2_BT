# PRINT MODEL PARAMETERS

print 'Running BT_print.py'

def print_parameters(model,pixel_sec,scale_kpc):
    print ' '
    ncomponent=1
    while ('COMP_%i'%ncomponent in model.header):
        #print '  # Component number: %i'%ncomponent
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
            try:
                AR,dAR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR))*180./pi #deg 
                dinclination=float(dAR)/sqrt(1.-float(AR)**2)*180/pi #deg
            except:
                AR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR[0]))*180./pi #deg 
                dinclination=0.
            #
            print '# Component type:                             %s'%model.header['COMP_%i'%ncomponent]
            print '  X Position [pixel]:                         %s'%model.header['%i_XC'%ncomponent]
            print '  Y Position [pixel]:                         %s'%model.header['%i_YC'%ncomponent]
            print '  Integrated magnitude [mag]:                 %s'%model.header['%i_MAG'%ncomponent]
            print '  Effective radius Re [pixel]:                %s'%model.header['%i_RE'%ncomponent]
            print '                                              %f +/- %f arcsec'%(Re_sec,dRe_sec)
            print '                                              %f +/- %f kpc'%(Re_kpc,dRe_kpc)
            print '  Sersic index:                               %s'%model.header['%i_N'%ncomponent]
            print '  Axis ratio (B/A):                           %s'%model.header['%i_AR'%ncomponent]
            print '  Inclination i [deg]                         %s +/- %s'%(inclination,dinclination)
            print '  Position angle PA [degrees: Up=0, Left=90]: %s'%model.header['%i_PA'%ncomponent]
            print ' '
        elif (model.header['COMP_%i'%ncomponent]=='sky'):
            print '# Component type:                             %s'%model.header['COMP_%i'%ncomponent]
            print '  Sky background [ADUs]:                      %s'%model.header['%i_SKY'%ncomponent]
            print ' '
        elif (model.header['COMP_%i'%ncomponent]=='expdisk'):
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
            try:
                AR,dAR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR))*180./pi #deg 
                dinclination=float(dAR)/sqrt(1.-float(AR)**2)*180/pi #deg
            except:
                AR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR[0]))*180./pi #deg  
                dinclination=0.
            #
            print '# Component type:                             %s'%model.header['COMP_%i'%ncomponent]
            print '  X Position [pixel]:                         %s'%model.header['%i_XC'%ncomponent]
            print '  Y Position [pixel]:                         %s'%model.header['%i_YC'%ncomponent]
            print '  Integrated magnitude [mag]:                 %s'%model.header['%i_MAG'%ncomponent]
            print '  Effective radius Rs [pixel]:                %s'%model.header['%i_RS'%ncomponent]
            print '                                              %f +/- %f arcsec'%(Rs_sec,dRs_sec)
            print '                                              %f +/- %f kpc'%(Rs_kpc,dRs_kpc)
            print '  Axis ratio (B/A):                           %s'%model.header['%i_AR'%ncomponent]
            print '  Inclination i [deg]                         %s +/- %s'%(inclination,dinclination)
            print '  Position angle PA [degrees: Up=0, Left=90]: %s'%model.header['%i_PA'%ncomponent]
            print '  '
        else:
            print '  Component type:                             %s'%model.header['COMP_%i'%ncomponent]
        ncomponent=ncomponent+1

def condensed_parameters(model,pixel_sec,scale_kpc):
    print ' '
    print '# Summary'
    ncomponent=1
    while ('COMP_%i'%ncomponent in model.header):
        #print '# Component number: %i'%ncomponent
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
            try:
                AR,dAR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR))*180./pi #deg 
                dinclination=float(dAR)/sqrt(1.-float(AR)**2)*180/pi #deg
            except:
                AR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR[0]))*180./pi #deg 
                dinclination=0.
            try:
                PA,dPA=re.findall("\d+\.\d+",model.header['%i_PA'%ncomponent]) # deg
                PA=float(PA)
                dPA=float(dPA)
            except:
                PA=re.findall("\d+\.\d+",model.header['%i_PA'%ncomponent]) # deg
                PA=float(PA[0])
                dPA=0.
            print '  Component type:                             %s'%model.header['COMP_%i'%ncomponent]
            print '  i=%.2f+/-%.2f deg, PA=%.2f+/-%.2f deg, Re=%.3f+/-%.3f arcsec'%(float(inclination),float(dinclination),float(PA),float(dPA),Re_sec,dRe_sec)
            print ' '
        elif (model.header['COMP_%i'%ncomponent]=='expdisk'):
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
            try:
                AR,dAR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR))*180./pi #deg 
                dinclination=float(dAR)/sqrt(1.-float(AR)**2)*180/pi #deg
            except:
                AR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
                inclination=arccos(float(AR[0]))*180./pi #deg  
                dinclination=0.
            try:
                PA,dPA=re.findall("\d+\.\d+",model.header['%i_PA'%ncomponent]) # deg
                PA=float(PA)
                dPA=float(dPA)
            except:
                PA=re.findall("\d+\.\d+",model.header['%i_PA'%ncomponent]) # deg
                PA=float(PA[0])
                dPA=0.
            #
            print '  Component type:                             %s'%model.header['COMP_%i'%ncomponent]
            print '  i=%.2f+/-%.2f deg, PA=%.2f+/-%.2f deg, Re=%.3f+/-%.3f arcsec'%(float(inclination),float(dinclination),float(PA),float(dPA),Rs_sec,dRs_sec)
            print ' '            
        else:
            pass #print '  Component type:                             %s'%model.header['COMP_%i'%ncomponent]
        ncomponent=ncomponent+1
