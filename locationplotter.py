import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

plt.rcParams['xtick.minor.visible'], plt.rcParams['xtick.top'] = True,True 
plt.rcParams['ytick.minor.visible'], plt.rcParams['ytick.right'] = True,True 
plt.rcParams['xtick.direction'], plt.rcParams['ytick.direction'] = 'in','in'
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["text.usetex"] = False
plt.rcParams["text.latex.preamble"] = r"\usepackage{txfonts}"
plt.rcParams['font.size'] = 14
#%%

""" Importing Data"""
data = pd.read_csv('Final_new.csv',  na_values='Nan')
data = data.drop(columns = ['Index','Finding chart name'])
data = data.to_numpy()

#Masking some data
gaia_confirm = data[:,1]
virac_confirm = data[:,3]
masker = np.where((gaia_confirm > 0) + (virac_confirm > 0))
dist = data[:,34]/1000

bensby17 = ascii.read('tablea2.dat', readme = 'ReadMe.txt')
glon = bensby17[9][:]
glat = bensby17[10][:]
glon1 = np.array(glon[masker])
glat1 = np.array(glat[masker])
# Importing contour data
data1 = np.genfromtxt('cobebulgecontour.dat')
arr = np.hsplit(data1,2)

#Converting data to a better format
l = len(data1)
for i in np.arange(len(data1)):
    x = arr[0][np.arange(2401)]
    y = arr[1][np.arange(2401)]


#%%

""" Plot with no function of distance"""

#Function for our circle
def circle(x, y, r):
     theta = np.arange(0, 360+1, 2) * np.pi / 180
     x1 = x + r*np.cos(theta)
     y1 = y + r*np.sin(theta)
     plt.plot(x1, y1, color='k', ls=':', linewidth=0.6)
     return
 
plt.axis('scaled')
for i in range(1, 7): circle(0, 0, i*2)
plt.ylim(-7,7)
plt.xlim(9,-9)
plt.scatter(x, y, linewidths=(0.01), c='k', ls=':', s=1)

color = dist[masker]
plt.scatter(glon, glat,s=25, color = 'red', edgecolors='black', alpha=1)

#plt.colorbar(label="Distance [kpc]", orientation="vertical")
plt.xlabel("Galactic longitude [deg]")
plt.ylabel("Galactic latitude [deg]")

#%%

""" Plot with a function of distance"""

def circle(x, y, r):
     theta = np.arange(0, 360+1, 2) * np.pi / 180
     x1 = x + r*np.cos(theta)
     y1 = y + r*np.sin(theta)
     plt.plot(x1, y1, color='k', ls=':', linewidth=0.6)
     return
plt.axis('scaled')
for i in range(1, 7): circle(0, 0, i*2)
plt.ylim(-7,7)
plt.xlim(9,-9)
plt.scatter(x, y, linewidths=(0.01), c='k', ls=':', s=1)

mask_main = np.where(gaia_confirm > 0)
color = dist[mask_main]
plt.scatter(glon, glat, alpha=0.3, s = 20 , edgecolors='k', color='white')
plt.scatter(glon[mask_main],glat[mask_main],s=25,c=color, cmap="magma", edgecolors='black', linewidth=0.4, vmax=9)
plt.colorbar(label="Distance [kpc]", orientation="vertical", pad=0.015)

plt.xlabel("Galactic longitude [deg]")
plt.ylabel("Galactic latitude [deg]")

