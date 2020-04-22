# Find pixel scale (in degrees) from fits image with header

print 'Running BT_pixel.py'

def def_pixel_scale_new(image):
    try:
        CD1_1=image.header['CD1_1']
        CD1_2=image.header['CD1_2']
        CD2_2=image.header['CD2_2']
        CD2_1=image.header['CD2_1']
        CDELT1 = sqrt(CD1_1**2+CD2_1**2)
        CDELT2 = sqrt(CD1_2**2+CD2_2**2)
    except:
        print 'WARNING: NO SCALE FOUND IN IMAGE, WE TOOK 0.03 ARCSEC'
        CDELT1= 0.03/3600. 
        CDELT2= 0.03/3600. 
    return CDELT1,CDELT2


def def_pixel_scale(image,project):
    try:
        CDELT1=image.header['CDELT1']
        CDELT2=image.header['CDELT2']
        if log10(max(abs(CDELT1),abs(CDELT2)))>-2.:            
            try:
                CD1_1=image.header['CD1_1']
                CD1_2=image.header['CD1_2']
                CD2_2=image.header['CD2_2']
                CD2_1=image.header['CD2_1']
                if (abs(CD1_1)==abs(CD1_1)):
                    return abs(CD1_1)
                else:
                    print 'WARNING: THE SCALES CD1_1 AND CD2_2 ARE DIFFERENT'
                    print '         we take pixel_scale=sqrt(CD1_1^2+CD2_2^2)'
                    print 'pixel_scale = %g deg = %g arcsec'%(sqrt(CD1_1**2+CD2_2**2),sqrt(CD1_1**2+CD2_2**2)*3600.)
                    return sqrt(CD1_1**2+CD2_2**2)
            except: 
                print 'WARNING: NO SCALE FOUND IN IMAGE FOR %s, WE TOOK 0.03 ARCSEC'%project
                return 0.03/3600.   
        elif (abs(CDELT1)==abs(CDELT2)):
            return abs(CDELT1)
        else:
            print 'WARNING: THE SCALES CDELT1 AND CDELT2 ARE DIFFERENT'
            print '         we take pixel_scale=sqrt(CDELT1^2+CDELT2^2)'
            print 'pixel_scale = %g deg = %g arcsec'%(sqrt(CDELT1**2+CDELT2**2),sqrt(CDELT1**2+CDELT2**2)*3600.)
            return sqrt(CDELT1**2+CDELT2**2)
    except:
        try:
            CD1_1=image.header['CD1_1']
            CD2_2=image.header['CD2_2']
            if (abs(CD1_1)==abs(CD1_1)):
                return abs(CD1_1)
            else:
                print 'WARNING: THE SCALES CD1_1 AND CD2_2 ARE DIFFERENT'
                print '         we take pixel_scale=sqrt(CD1_1^2+CD2_2^2)'
                print 'pixel_scale = %g deg = %g arcsec'%(sqrt(CD1_1**2+CD2_2**2),sqrt(CD1_1**2+CD2_2**2)*3600.)
                return sqrt(CD1_1**2+CD2_2**2)
        except:
            print 'WARNING: NO SCALE FOUND IN IMAGE FOR %s, WE TOOK 0.03 ARCSEC'%project
            return 0.03/3600.

