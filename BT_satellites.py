# CUSTOM FUNCTION TO ADD SATELLITES IN THE IMAGE IF NEEDED

'''
Satellites are fitted as Sersic fits
- The magnitudes and half light radii of the satellites can be fixed (variable=False) or variable (variable=True)
- The positions, Sersic indices, axis ratios and position angles of the satellites can be fixed (all_variable=False) or variable (all_variable=True)
'''

print 'Running BT_satellites.py'

def add_satellite(feedme_file,ID,variable=True,all_variable=True):
    with open(feedme_file,'a') as f:
        if ID=='XA53':
            print '--------- A satellite was added'
            f.write('\n')
            f.write('# Object: satellite \n')
            f.write(' 0) sersic                 #  Component type\n')
            f.write(' 1) 229.64 242.95 %i %i  #  position x, y \n'%(int(all_variable),int(all_variable)))
            f.write(' 3) -3.72     %i          #  Integrated magnitude \n'%int(variable))
            f.write(' 4) 4.22      %i          #  R_e (half-light radius) [pix] \n'%int(variable))
            f.write(' 5) 1.07      %i          #  Sersic index n (de Vaucouleurs n=4)\n'%int(all_variable))
            f.write(' 6) 0.0000      0          #     ----- \n')
            f.write(' 7) 0.0000      0          #     ----- \n')
            f.write(' 8) 0.0000      0          #     ----- \n')
            f.write(' 9) 0.80        %i           #  axis ratio (b/a) \n'%int(all_variable))
            f.write(' 10) 63.14      %i             #  position angle (PA) [deg: Up=0, Left=90] \n'%int(all_variable))
            f.write(' Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n')
       
        elif ID=='L14CO007': 
            print '--------- Two satellites were added'
            f.write('\n')
            f.write('# Object: satellite \n')
            f.write(' 0) sersic                 #  Component type\n')
            f.write(' 1) 177.75 319.09 %i %i  #  position x, y \n'%(int(all_variable),int(all_variable)))
            f.write(' 3) -3.30     %i          #  Integrated magnitude \n'%int(variable))
            f.write(' 4) 14.64     %i          #  R_e (half-light radius) [pix] \n'%int(variable))
            f.write(' 5) 0.87      %i          #  Sersic index n (de Vaucouleurs n=4)\n'%int(all_variable))
            f.write(' 6) 0.0000      0          #     ----- \n')
            f.write(' 7) 0.0000      0          #     ----- \n')
            f.write(' 8) 0.0000      0          #     ----- \n')
            f.write(' 9) 0.70      %i          #  axis ratio (b/a) \n'%int(all_variable))
            f.write(' 10) -20.16   %i         #  position angle (PA) [deg: Up=0, Left=90] \n'%int(all_variable))
            f.write(' Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n')
            f.write('\n')
            f.write('# Object: satellite \n')
            f.write(' 0) sersic                 #  Component type\n')
            f.write(' 1) 312.03 212.39 %i %i  #  position x, y \n'%(int(all_variable),int(all_variable)))
            f.write(' 3) -3.83     %i          #  Integrated magnitude \n'%int(variable))
            f.write(' 4) 17.22     %i          #  R_e (half-light radius) [pix] \n'%int(variable))
            f.write(' 5) 1.34      %i          #  Sersic index n (de Vaucouleurs n=4)\n'%int(all_variable))
            f.write(' 6) 0.0000      0          #     ----- \n')
            f.write(' 7) 0.0000      0          #     ----- \n')
            f.write(' 8) 0.0000      0          #     ----- \n')
            f.write(' 9) 0.56     %i          #  axis ratio (b/a) \n'%int(all_variable))
            f.write(' 10) 39.08    %i          #  position angle (PA) [deg: Up=0, Left=90] \n'%int(all_variable))
            f.write(' Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n')
