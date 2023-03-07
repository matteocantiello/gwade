import mesa_reader as mr
import numpy as np
from pylab import *
from math import log10, pi
import pickle

from scipy import interpolate
from scipy.interpolate import griddata

from numpy import loadtxt
import pandas as pd 
 
from constants import *
from plot_settings import *
from functions import *
from getters import *
from star_data import stars

prefix,DIR,mods,hs,pf = pickle.load(open('parsed.data','rb'))

N = 1000
#r_grid = 1 - 10**np.linspace(-5,0,num=nR,endpoint=True)
tau_grid = 10**np.linspace(0,5,num=N,endpoint=True)

mS = []

B = []

for m,p in zip(*(mods,pf)):
	if m < 1.0 or m > 2.0:
		continue

	mS.append(m)

	tau = np.array(p.tau)

	B.append(np.interp(tau_grid, tau, p.B_shutoff))

B = np.array(B).T
B[B < 1e2] = 0

fig = plt.figure()
ax = plt.subplot(111)
cntr = plt.pcolormesh(mS, tau_grid, np.log10(B))
cbar = fig.colorbar(cntr, ax=ax)
cbar.ax.set_ylabel(r'$\log B_{\rm shutoff}/\mathrm{G}$')
plt.xlabel(r'$M/M_\odot$')
plt.ylabel(r'$\tau$')
plt.yscale('log')
plt.ylim([1e5,1])
plt.savefig('../figures/m_tau.png')

