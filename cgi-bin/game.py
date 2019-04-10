
# coding: utf-8

# # Calculations for the GravCalc game
# * Step 1: write function that takes a black hole binary position, and a pulsar position, calculates the Hellings and Downs factor, and returns an amplitude.  Only the earth term. This factor will multiply the residual amplitude done in the next function.  This function really just handles the angular part.
# ```
# def calc_hdamp(gwra, gwdec, pra, pdec):
# ```
# * Step 2: write function that calculates the strain amplitude for two black holes with total mass, period, and distance
# We can set e=0, mass_ratio=1, and i=0 in your routine
# THIS IS LIAM'S CALCULATION (that you probably already have somewhere)
# ```
# def calc_resid(totalmass, period, distance):
# ```
# * Step 3: write function that takes an amplitude, a period, and a cadence, and creates a sine wave at that amplitude and period, with a randomized phase
# ```
# def create_sine(amp, period, cadence):
# ```
# The "amp" that we'll put in there will be a produce of the amps from Step 1 and 2.
# 
# * Step 4: write a function that takes a list of black hole positions, and periods, and adds all the sine waves together, and fits for a period and period derivative and subtracts it
# ```
# def sum_BHBs(ras, decs, periods, totalmasses, distances, cadence, lengthofdata):
# ```

# In[11]:


# Importing libraries and setting up the length and cadence of data
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# Set number of residual points
# 10 years of data, with residuals every 3 weeks
cadence = 21  # number of days between data points
lengthofdata = 10*365.25 # days 
cpc = 9.71561e-9 # speed of light in pc per second


# In[2]:


def calc_hdamp(gwra, gwdec, pra, pdec):
    # By "hellings and downs" I just mean that this function handles
    # the angular part of the amplitude - the part that's just based
    # on the angle between the pulsar and the GW.
    # We do the dot product in spherical coordinates
    # RA is like phi, and dec is like 90 - theta
    gphi = -(np.pi/2 - gwdec) # We use the opposite of the GW position
    # because we want the propagation direction
    pphi = np.pi/2 - pdec
    gthe = gwra + np.pi # Again, the opposite
    pthe = pra
    dotprod = np.sin(gphi)*np.sin(pphi)*np.cos(gthe - pthe) + np.cos(gphi)*np.cos(pphi)
    # oneplus is 1 + k\dot n where k is the propagation direction
    # of the GW, and n is the direction of the pulsar
    oneplus = 1 + dotprod
    # We need to account for exactly aligned which will be zero
    if oneplus > .01:  # Not sure what the right tolerance is
        amp = 1/oneplus
    else:
        amp = 0.
    return amp
    


# In[3]:


def calc_resid(total_m, orbital_p, distance):
    # This is in nanoseconds
    # It expects mass in solar masses
    # orbital_period in years
    # distance in Mpc
    orbital_e = 0.0
    orbital_i= 0.0
    m_ratio = 1.0
    result_rs = 190 * (((total_m)/10e9)**(5./3.)) * ((orbital_p)**(1./3.)) * (100/distance) * (m_ratio/((1 + m_ratio))**2) * np.sqrt(1 - (orbital_e)**2) * (1 + np.cos(orbital_i**2))/2
    return result_rs


# In[4]:


def create_sine(amp, period, cadence, lengthofdata):
    '''
    function that takes an amplitude, a GW period, and a cadence, 
    and creates a sine wave at that amplitude and period, with 
    a randomized phase Cadence is given in days
    '''
    # Create time series
    x = np.arange(0,lengthofdata, cadence) # in days
    # Create randomized phase
    phase = np.random.uniform(0,2*np.pi)
    # Make sine wave
    y = amp*np.sin(2*np.pi*x/period + phase)
    return [x,y]


# In[5]:


def sum_BHBs(ras, decs, periods, totalmasses, distances, cadence, lengthofdata):
    # Set the pulsar position to whatever is in the middle
    # of the plot
    pra = np.pi
    pdec = 0.0
    resids = np.zeros(np.size(np.arange(0,lengthofdata, cadence))) # This should 
        # at some point. It'll work until we change the line that looks
        # like it in create_sine
    for i in range(np.size(ras)):
        # Looping through the BHBs
        hdamp=calc_hdamp(ras[i], decs[i], pra, pdec)
        # these positions are all expected in radians
        resid_amp=calc_resid(totalmasses[i], periods[i], distances[i])
        x,res=create_sine(hdamp*resid_amp, periods[i]*365.25, cadence, lengthofdata)
        resids += res
    return x, resids


# In[6]:


# David, these arrays would be populated by user input, in this
# case the user moving around the binary black holes on the game
# board

# These first two (ras and decs) are created by the position that the user places
# the black hole. They are the longitude and latitude on the sky.
ras = np.array([.1, .2, .5, .7])
decs = np.array([.2, .4, .6, .8])  # in radians

# These next three (period, totalmass, distance) are input by 
# which little bubbles the user chooses to drag
periods = np.array([.4, 1, 2, 16])  #these are in years
totalmasses = np.array([3e8, 1e9,3e9, 5e9]) #in solarmasses
distances = np.array([50,100, 200, 50]) #in Mpc


# In[7]:



def func(x, a, b, c):
    return a * x**2 + b*x + c

def fit_p(x, res):
    # fit a polynomial to the residuals
    popt, pcov = curve_fit(func, x, res)
    return (res-func(x, *popt))
    


# In[8]:


time,res = sum_BHBs(ras, decs, periods, totalmasses, distances, cadence, lengthofdata)
subtracted_res=fit_p(time,res)


# In[10]:


# This is sort of what we would like the plots to look like
plt.plot(time, res)
plt.ylabel("Residual in ns")
plt.xlabel("Time in days")
plt.title("Residuals prior to pulsar model fitting")
plt.show() 
plt.plot(time, subtracted_res)
plt.ylabel("Residual in ns")
plt.xlabel("Time in days")
plt.title("Residuals after pulsar model fitting")
plt.show()

