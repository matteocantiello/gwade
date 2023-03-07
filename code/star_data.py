import numpy as np

# Stars
fi = open('../data/roAp_stars_2.csv','r')
stars = []
for i,line in enumerate(fi):
    if i > 0:
        s = line.split(',')
        if len(s[6]) > 0 and len(s[2]) > 0:
            logL = float(s[6])
            logT = float(s[2])
            stars.append([logT,logL])
stars = np.array(stars)