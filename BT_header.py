# GET HEADER PARAMETERS

'''
DEFINES    XC, DXC FROM HEADER
           YC, DXC
           PA, dPA
           AR, dAR
           XMIN, XMAX, YMIN, YMAX, RMAX
           size_image

REQUIRES image, model
'''

print 'Running BT_header.py'

def init_header_parameters():
    header=dict()
    header['XC']=0.
    header['DXC']=0.
    header['YC']=0.
    header['DYC']=0.
    header['PA']=0.   
    header['PA_rad']=0.
    header['AR']=0.
    header['center']=0.
    header['RMAX']=0.
    return header
                   
def get_header_parameters(galfit_output,ncomponent):
    output_list=pyfits.open(galfit_output)
    image = output_list[1]
    model = output_list[2]
    header=dict()
    try:
        XC,DXC=re.findall("\d+\.\d+",model.header['%i_XC'%ncomponent]) # pix
        header['XC']=float(XC)
        header['DXC']=float(DXC)
    except:
        XC=re.findall("\d+\.\d+",model.header['%i_XC'%ncomponent]) # pix
        XC=XC[0]
        header['XC']=float(XC)
        header['DXC']=0. 
    try:
        YC,DYC=re.findall("\d+\.\d+",model.header['%i_YC'%ncomponent]) # pix
        header['YC']=float(YC)
        header['DYC']=float(DYC)
    except:
        YC=re.findall("\d+\.\d+",model.header['%i_YC'%ncomponent]) # pix
        YC=YC[0]
        header['YC']=float(YC)
        header['DYC']=0. 
    try:
        PA,dPA=model.header['%i_PA'%ncomponent].split(' +/- ')
        #re.findall("\d+\.\d+",model.header['%i_PA'%ncomponent]) # deg
        PA=float(PA)
        header['PA']=PA
    except:
        PA=model.header['%i_PA'%ncomponent][1:-1]
        #re.findall("\d+\.\d+",model.header['%i_PA'%ncomponent]) # deg
        PA=float(PA)
        header['PA']=PA
    PA_rad=PA*pi/180. # rad
    header['PA_rad']=PA_rad
    try:
        AR,dAR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
        AR=float(AR)
        header['AR']=AR
    except:
        AR=re.findall("\d+\.\d+",model.header['%i_AR'%ncomponent])
        AR=float(AR[0])
        header['AR']=AR
    try:
        XMIN=int(re.findall("\d+", image.header['OBJECT'])[-4])
        XMAX=int(re.findall("\d+", image.header['OBJECT'])[-3])
        YMIN=int(re.findall("\d+", image.header['OBJECT'])[-2])
        YMAX=int(re.findall("\d+", image.header['OBJECT'])[-1])
        center=array([float(XC)-XMIN,float(YC)-YMIN])
        header['center']=center
    except:
        center=array([float(XC),float(YC)])
        header['center']=center
    size_image= (image.header['NAXIS1'],image.header['NAXIS2'])
    RMAXX=max(center[0],center[1],size_image[0]-center[0],size_image[1]-center[1])
    header['RMAX']=RMAXX
    return header
    
def shift_center(galfit_output,center):
    output_list=pyfits.open(galfit_output)
    image = output_list[1]
    try:
        XMIN=int(re.findall("\d+", image.header['OBJECT'])[-4])
        XMAX=int(re.findall("\d+", image.header['OBJECT'])[-3])
        YMIN=int(re.findall("\d+", image.header['OBJECT'])[-2])
        YMAX=int(re.findall("\d+", image.header['OBJECT'])[-1])  
        return array([center[0]-XMIN,center[1]-YMIN])
    except:
        print 'ERROR: DID NOT FIND THE HEADER OF %s'%galfit_output     

