import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from constants import rsun

def rotation_law(h):
    mass = h.star_mass[0]
    vsurf = 1.5e7
    if mass < 1.2:
        vsurf *= (mass/1.2)**2
    omega = vsurf / (rsun*10**h.log_R)
    return omega

def tri_area(xs,ys):
  arr = np.ones((3,3))
  arr[0] = xs
  arr[1] = ys
  area = 0.5 * np.linalg.det(arr)
  return area

def tri_max_length(xs,ys):
    side_lengths = list(((xs[i]-xs[i-1])**2+(ys[i]-ys[i-1])**2)**0.5 for i in range(3))
    return max(side_lengths)

def find_zams(logl,loglh,model):
    zams=1
    while (loglh[zams] < 1.0*logl[zams]): 
     zams=zams+1
    return zams; 

def find_h(dh,center_h1,model):
    zams=1
    while (center_h1[zams] > (center_h1[1] - dh)): 
     zams=zams+1
    return zams; 

def find_mams(center_h1,model):
    mams=1
    while (center_h1[mams] > 0.5 * center_h1[1]): 
     mams=mams+1
    return mams; 

def find_tams(center_h1,model):
    tams=1
    while (center_h1[tams] > 0.05): 
     tams=tams+1
    return tams;    

def find_max(a,b,c,d):
    z= [0] * len(a)
    for i in range(0, len(a)):
      z[i]=max(a[i],b[i],c[i],d[i])   
    return z;

def find_mid_ms(model,star_age,zams,tams):
    mid_ms=1
    age_ms=(star_age[tams]-star_age[zams])/2
    while (star_age[mid_ms] < age_ms): 
     mid_ms=mid_ms+1
    return mid_ms;  

def find_frac_ms(model,star_age,zams,tams,frac):
    frac_ms=1
    age_frac_ms=(star_age[tams]-star_age[zams])*frac
    while ((star_age[frac_ms] - star_age[zams]) < age_frac_ms): 
     frac_ms=frac_ms+1
    return frac_ms;   

def concat(a,b):
    return a+'_'+b

def CustomCmap(from_rgb,to_rgb):

    # from color r,g,b
    r1,g1,b1 = from_rgb

    # to color r,g,b
    r2,g2,b2 = to_rgb

    cdict = {'red': ((0, r1, r1),
                   (1, r2, r2)),
           'green': ((0, g1, g1),
                    (1, g2, g2)),
           'blue': ((0, b1, b1),
                   (1, b2, b2))}

    cmap = LinearSegmentedColormap('custom_cmap', cdict)
    return cmap