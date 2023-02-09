from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from PyAstronomy import pyasl
from galpy.orbit import Orbit
from galpy.potential import LogarithmicHaloPotential
from galpy.potential import MWPotential2014
from astropy.io import ascii
from tqdm import tqdm

plt.rcParams['xtick.minor.visible'], plt.rcParams['xtick.top'] = True, True
plt.rcParams['ytick.minor.visible'], plt.rcParams['ytick.right'] = True, True
plt.rcParams['xtick.direction'], plt.rcParams['ytick.direction'] = 'in', 'in'
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["text.usetex"] = False
plt.rcParams["text.latex.preamble"] = r"\usepackage{txfonts}"
plt.rcParams['font.size'] = 14
# %%

""" Importing Data """

bensby17 = ascii.read('tablea2.dat', readme='ReadMe.txt')
bensby21 = ascii.read('table2.dat', readme='ReadMe(1).txt')
ider = bensby17['Num']
data_new = np.genfromtxt('gaia_output.txt', skip_header=True, delimiter=",")
data_virac = np.genfromtxt('virac_output.txt', skip_header=True, delimiter=",")

data = pd.read_csv('Final_new.csv',  na_values='Nan')
data = data.drop(columns=['Index', 'Finding chart name'])
data = data.to_numpy()

HRV, peri_gaia, apo_gaia, ecc_gaia, z_max_gaia, energy_gaia = data[:, 8], data[:,18], data[:, 19], data[:, 20], data[:, 21], data[:, 22]
age, gaia_geodist =  data[:,29], data[:, 34]
coord_ra, coord_dec, pm_ra, pm_dec, hrv = data_new[:,3], data_new[:, 4], data_new[:, 1], data_new[:, 2], data_new[:, 5]
pmra_err,pmdec_err = data_new[:,6], data_new[:,7]
# %%
""" Test Run """
# Initialize our empty lists

lister = np.empty([90, 9])
lister[:] = np.nan

for i in np.arange(0, 90):
    idx_coord = i
    if np.isnan(pm_ra)[idx_coord]:
        continue
    print(idx_coord)

    # Gaia data
    ra = coord_ra[idx_coord]
    dec = coord_dec[idx_coord]
    pmra = pm_ra[idx_coord]
    pmdec = pm_dec[idx_coord]
    d = gaia_geodist[idx_coord]
    rv = hrv[idx_coord]
    
    # Calculate our galactic space velocities
    gal_vel = pyasl.gal_uvw(ra, dec, pmra, pmdec, d, rv)    
    U,V,W = gal_vel[0], gal_vel[1], gal_vel[2]

    # Galpy Potentials
    lp = LogarithmicHaloPotential(normalize=1.)
    mw = MWPotential2014
    # Galpy orbital integrations
    o = Orbit([ra, dec, d/1000, U, V, W], radec=True, uvw=True, ro=8., vo=220.)
    ts = np.linspace(0, 1, 10)*u.Gyr
    o.integrate(ts, mw)
    o.plot('k.')
    #orbits.plot()

    print(o.rap(), o.rperi(), o.e(), o.zmax(), o.E())

    lister[i, :] = [ider[i],U,V,W, o.rap(), o.rperi(), o.e(), o.zmax(), o.E()]

    #ascii.write(lister, 'gaia_parameters_mp.txt', overwrite=True, format='csv')

#Saving our results
    name = [idx_coord]
    save_results_to = '/Users/calum/Desktop/Master_Thesis/Findings/bulge_findingcharts/Orbit_integration_gaia/R_Z/ name={}'.format(
        name) + '.png'
    #plt.savefig(save_results_to)
    

#%%    

# Cleaning some data
dister = np.genfromtxt('b_j.tsv', delimiter="\t")

split_inds = np.where(np.diff(dister[:,0]) != 0)[0] + 1
split_vir = np.split(dister, split_inds)

num_cols = split_vir[0].shape[1]
result = np.zeros((len(split_vir), num_cols))
                  
for i, arr in enumerate(split_vir):
    min_ind = np.argmin(arr[:,2])
    result[i] = arr[min_ind]
#np.savetxt("dister.txt", result, delimiter=";", fmt="%-13.8f")

# Adding (upper+lower)/2 is an assumption we make to have normal distribution of our distances

bailer = np.genfromtxt('dister.txt', delimiter=";")
lower, upper = bailer[:,3], bailer[:,4]
approxer = (lower + upper)/2
sigmer = (upper-lower)/2
#(upper - lower)/2 for a pseudo uncertainty

#%%

"""Bootstrapping Run"""

# Number of samples
sampleno = 100

# Empty list that stores: [star, sample number, [U, V, W, rap, peri, ecc, zmax, E]]
storer = np.zeros((90, sampleno, 3 + 5))
ts = np.linspace(0, 1, 10)*u.Gyr

#Randomly draw our samples
for star in tqdm(np.arange(0, 90)):    
    pmra_samp = np.random.normal(pm_ra[star], pmra_err[star], size=sampleno)
    pmdec_samp = np.random.normal(pm_dec[star], pmdec_err[star], size=sampleno)
    dist_samp = np.random.normal(gaia_geodist[star], sigmer[star], size=sampleno)

#Generate our galactic space velocities U,V,W    
    for sample in tqdm(np.arange(sampleno)):
        pmra, pmdec, d = pmra_samp[sample], pmdec_samp[sample], dist_samp[sample]
        U, V, W = pyasl.gal_uvw(coord_ra[star], coord_dec[star], pmra, pmdec, d, hrv[star])
        storer[star, sample, 0:3] = [U, V, W]
#Integrate our orbits       
        o = Orbit([coord_ra[star], coord_dec[star], d/1000, U, V, W], radec=True, uvw=True, ro=8., vo=220.)
        o.integrate(ts, mw)
        storer[star, sample, 3:] = [o.rap(), o.rperi(), o.e(), o.zmax(), o.E()]
        
# verify the mean
verifier = abs(pmdec_err[star]-np.mean(pmra_samp))

# Save 
for i in range(90):    
    np.savetxt(f"gaia_error/gaia_error_1gyr_{i+1}.txt", storer[i], delimiter=';', fmt='%1.4f', overwrite=False)

#%%

#Empty list to store mean and stds

list_mean = np.zeros([90,5])
list_std = np.zeros([90,5])

#Calculate mean and stds
for i in range(0,90):
    gaia_uncert = np.genfromtxt(f'gaia_error/gaia_error_1gyr_{i+1}.txt', delimiter=";")
    apo, peri, ecc, zmax, E = gaia_uncert[:,3], gaia_uncert[:,4], gaia_uncert[:,5], gaia_uncert[:,6], gaia_uncert[:,7]
    list_mean[i] = np.array([np.mean(apo), np.mean(peri), np.mean(ecc), np.mean(zmax), np.mean(E)])
    list_std[i] = np.array([np.std(apo), np.std(peri), np.std(ecc), np.std(zmax), np.std(E)])

#np.savetxt('gaia_mean_std.txt', list_mean[i], list_std[i] ,delimiter=';')
ascii.write(list_mean ,'gaia_mean_1gyr.txt', overwrite = False, format = 'csv')
ascii.write(list_std ,'gaia_std_1gyr.txt', overwrite = False, format = 'csv')
