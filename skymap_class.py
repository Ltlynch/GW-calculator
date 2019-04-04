
# coding: utf-8

# In[16]:


get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import healpy as hp
import pandas as pd
import scipy.constants as const
from astropy.io import fits
from astropy.cosmology import WMAP5
from astropy import units as u
from astropy.coordinates import SkyCoord
from mpmath import *
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d


# In[21]:


class Skymap:
    def __init__(self, ulfile, pixelfile):
        assert(type(ulfile) == str and type(pixelfile) == str)
        assert(os.path.isfile(ulfile) and os.path.isfile(pixelfile))
        
        ul = np.loadtxt(ulfile)
        pixel = np.genfromtxt(pixelfile, names = True)
        
        self.theta = pixel['theta']
        self.phi = pixel['phi']
        
        self.map_dictionary = {}
        self.map_dictionary[95] = self.create_map(ul[:,1])
        self.map_dictionary[75] = self.create_map(ul[:,2])
        self.map_dictionary[50] = self.create_map(ul[:,3])
        self.map_dictionary[25] = self.create_map(ul[:,4])
        self.map_dictionary[1] = self.create_map(ul[:,5])
    
    def create_map(self, ulcolumn):
        npix = len(self.theta)
        nside = hp.npix2nside(npix)
        indices = hp.ang2pix(nside, self.theta, self.phi)
        skymap = np.zeros(npix, dtype=np.float)
        skymap[indices] += ulcolumn[indices]
        return skymap
    
    def interpolate(self, conf, theta, phi):
        if conf == 'all':
            maps = []
            for key in sorted(self.map_dictionary.keys()):
                maps.append(hp.pixelfunc.get_interp_val(self.map_dictionary[key], theta, phi))
            return maps
        else:
            return hp.pixelfunc.get_interp_val(self.map_dictionary[conf], theta, phi)
    
    def plot(self, confidence):
        hp.mollview(np.log10(self.map_dictionary[confidence]), title= 'Upper Limit Skymap')
        
skymap3nHz = Skymap('uls_smoothed_3nHz.txt', 'skymap_pixels.txt')
skymap10nHz = Skymap('uls_smoothed_10nHz.txt', 'skymap_pixels.txt')
skymap31nHz = Skymap('uls_smoothed_31nHz.txt','skymap_pixels.txt')
skymap100nHz = Skymap('uls_smoothed_100nHz.txt','skymap_pixels.txt')
skymap318nHz = Skymap('uls_smoothed_318nHz.txt','skymap_pixels.txt')

#skymap3nHz.plot(95)

input_theta = .3
input_phi = .4
input_freq = 5
input_strain = 5e-15

#print(skymap3nHz.interpolate('all', input_theta, input_phi))
#print(skymap10nHz.interpolate('all', input_theta, input_phi))
#print(skymap31nHz.interpolate('all', input_theta, input_phi))
#print(skymap100nHz.interpolate('all', input_theta, input_phi))
#print(skymap318nHz.interpolate('all', input_theta, input_phi))


# In[26]:


interp3nHz= np.array(skymap3nHz.interpolate('all', input_theta, input_phi))
interp10nHz = np.array(skymap10nHz.interpolate('all', input_theta, input_phi))
interp31nHz = np.array(skymap31nHz.interpolate('all', input_theta, input_phi))
interp100nHz = np.array(skymap100nHz.interpolate('all', input_theta, input_phi))
interp318nHz = np.array(skymap318nHz.interpolate('all', input_theta, input_phi))

combined_columns = np.array([interp3nHz, interp10nHz, interp31nHz, interp100nHz, interp318nHz])

#print(interp3nHz)
#print(interp10nHz)

#f = np.concatenate((interp1, interp2), axis=1)

freq_array = np.array([3, 10, 31, 100, 318])
#print(interp1)
#print(freq_array)

row1 = np.array(combined_columns[:,0])
row25 = np.array(combined_columns[:,1])
row50 = np.array(combined_columns[:,2])
row75 = np.array(combined_columns[:,3])
row95 = np.array(combined_columns[:,4])

combined_rows = np.array([row1, row25, row50, row75, row95])

#freq_interp = interp1d(freq_array, combined)
freq_interp1 = interp1d(freq_array, combined_rows[0], kind = 'cubic')
freq_interp2 = interp1d(freq_array, combined_rows[1], kind = 'cubic')
freq_interp3 = interp1d(freq_array, combined_rows[2], kind = 'cubic')
freq_interp4 = interp1d(freq_array, combined_rows[3], kind = 'cubic')
freq_interp5 = interp1d(freq_array, combined_rows[4], kind = 'cubic')

# = np.array([freq_interp1(input_freq), freq_interp2(input_freq), freq_interp3(input_freq), freq_interp4(input_freq), freq_interp5(input_freq)])


confidence_array = [1, 25, 50, 75, 95]

new_freq_rows = np.array([freq_interp1(input_freq), freq_interp2(input_freq), freq_interp3(input_freq), freq_interp4(input_freq), freq_interp5(input_freq)])

#print(new_freq_rows)

conf_interp = interp1d(new_freq_rows, confidence_array, kind='cubic')

print(float(conf_interp(input_strain)))

#confidence_interp = interp1d(new_freq_array,)
#print(freq_interp2(input_freq))
#print(freq_interp(input_freq))
#freq_interp2(input_freq)
#freq_interp3(input_freq)

#print(full)
#print(full[0][0])

