# TRACE THE RADIAL DENSITY PROFILE

print 'Running BT_profile.py'

def radial_profile(data, center):
    y,x = np.indices((data.shape)) # first determine radii of all pixels
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    ind = np.argsort(r.flat) # get sorted indices
    sr = r.flat[ind] # sorted radii
    sim = data.flat[ind] # image values sorted by radii    
    ri = sr.astype(np.int32) # integer part of radii (bin size = 1)
    # determining distance between changes
    deltar = ri[1:] - ri[:-1] # assume all radii represented
    rind = np.where(deltar)[0] # location of changed radius
    nr = rind[1:] - rind[:-1] # number in radius bin
    csim = np.cumsum(sim, dtype=np.float64) # cumulative sum to figure out sums for each radii bin
    tbin = csim[rind[1:]] - csim[rind[:-1]] # sum for image values in radius bins
    radialprofile = tbin/nr # the answer
    return radialprofile

def deprojected_profile(data,center,PA_rad,AR):
    y,x = np.indices((data.shape)) # first determine position of all pixels
    x=x-center[0] # take center as origin
    y=y-center[1]
    xp=x*cos(PA_rad)+y*sin(PA_rad) # rotate by the position angle PA
    yp=-x*sin(PA_rad)+y*cos(PA_rad)
    rp=np.sqrt((xp/AR)**2+yp**2) # determine the deprojected radius
    ind = np.argsort(rp.flat) # get sorted indices
    sr = rp.flat[ind] # sorted radii
    sim = data.flat[ind] # image values sorted by radii    
    rpi = sr.astype(np.int32) # integer part of radii (bin size = 1)
    # determining distance between changes
    deltarp = rpi[1:] - rpi[:-1] # assume all radii represented
    rpind = np.where(deltarp)[0] # location of changed radius
    nrp = rpind[1:] - rpind[:-1] # number in radius bin
    csim = np.cumsum(sim, dtype=np.float64) # cumulative sum to figure out sums for each radii bin
    tbin = csim[rpind[1:]] - csim[rpind[:-1]] # sum for image values in radius bins
    radialprofile = tbin/nrp # the answer
    return radialprofile  

def center_coord(model,center):
    import re
    limits=re.findall("\d+\d+",model.header['FITSECT'])
    xc=int(limits[0])+center[0]
    yc=int(limits[2])+center[1]
    return xc, yc
    
    


