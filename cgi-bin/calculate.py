#!/usr/bin/env python36
import sys
import traceback
import cgi
import cgitb
import json
import numpy as np
import scipy.constants as const
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy.cosmology import WMAP5
import os
from scipy.interpolate import interp1d
#from .skymap_class import Skymap
from uuid import uuid4

sys.stderr = sys.stdout
cgitb.enable()

try:

    form = cgi.FieldStorage()

    print("Content-Type: application/json")
    print("Access-Control-Allow-Origin: *")
    print("")

    action = form.getfirst("action","calculate")

    def yrs_to_s(years):
        return years*365*24*3600

    def s_to_yrs(seconds):
        return seconds/(365*24*3600)

    def pc_to_m(pc):
        return pc*3.0856776e16

    def m_to_pc(m):
        return m/3.085677e16

    def kg_to_M0(kg):
        return kg/1.98847e30

    def M0_to_kg(M0):
        return M0*1.98847e30

    if action == 'redshiftcalc':

        H0 = 70
        distance = float(form.getfirst("distance",""))
        redshift = round((distance*H0)/(const.c*1e-3), 4) # const.c is in meters, converting to km


        print(json.dumps({'redshift':redshift}))

    elif action == 'distancecalc':

        redshift = float(form.getfirst("redshift",""))
        #distance = WMAP5.comoving_distance(redshift)	
        H0 = 70
        distance = round((redshift*const.c*1e-3)/H0 ,2)

        print(json.dumps({'distance':distance}))

    elif action == 'orbitalsepcalc':

        orbital_p = float(form.getfirst("orbital_p",""))
        total_m = float(form.getfirst("total_m",""))
        orbital_s = round(m_to_pc(2*np.cbrt((const.G * M0_to_kg(total_m) * yrs_to_s(orbital_p)**2)/(4*(const.pi)**2))), 2)

        print(json.dumps({'orbital_s':orbital_s}))

    elif action == 'orbitalpercalc':

        orbital_s = float(form.getfirst("orbital_s",""))
        total_m = float(form.getfirst("total_m",""))
        orbital_p = s_to_yrs(np.sqrt((4 * (const.pi**2) * (pc_to_m(orbital_s)/2)**3)/(const.G * M0_to_kg(total_m))))

        print(json.dumps({'orbital_p':orbital_p}))

    elif action == 'totalmasscalc':	
        
        orbital_p = float(form.getfirst("orbital_p",""))
        orbital_s = float(form.getfirst("orbital_s",""))
        total_m = round(float(kg_to_M0((4*(const.pi**2)*(pc_to_m(orbital_s/2))**3)/(const.G * (yrs_to_s(orbital_p))**2))), 4)

        print(json.dumps({'total_m':total_m}))

    else:
        # calculate
        m1 = float(form.getfirst("m1", ""))
        m2 = float(form.getfirst("m2", ""))
        total_m = float(form.getfirst("total_m",""))
        m_ratio = float(form.getfirst("m_ratio",""))
        distance = float(form.getfirst("distance", ""))
        redshift = float(form.getfirst("redshift", ""))
        orbital_p = float(form.getfirst("orbital_p",""))
        orbital_s = float(form.getfirst("orbital_s",""))
        orbital_i = float(form.getfirst("orbital_i",""))
        orbital_e = float(form.getfirst("orbital_e",""))
        #TODO figure out how to parse the RAJ string and convert it into phi
        position_ra = form.getfirst("position_ra","")
        position_dec = form.getfirst("position_dec","")
        
        #hours, minutes = position_dec.split('.')
        #dec_hours = '+{0}d{1}m'.format(hours, minutes)
        
        #coords = SkyCoord(position_ra, dec_hours, frame = 'icrs')
        
        RA_hours = position_ra.replace('h',',').replace('m',',').replace('s', ',').split(',')
        input_phi = float(position_dec)
        input_theta = float(RA_hours[0]) + float(RA_hours[1])/60 + float(RA_hours[2])/3600
        
        result_rs = 190 * (((total_m)/10e9)**(5./3.)) * ((orbital_p)**(1./3.)) * (100/distance) * (m_ratio/((1 + m_ratio))**2) * np.sqrt(1 - (orbital_e)**2) * (1 + np.cos(orbital_i**2))/2

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
         
        skymap3nHz = Skymap('../cw_maps_database/uls_smoothed_3nHz.txt', '../cw_maps_database/skymap_pixels.txt')
        skymap10nHz = Skymap('../cw_maps_database/uls_smoothed_10nHz.txt', '../cw_maps_database/skymap_pixels.txt')
        skymap31nHz = Skymap('../cw_maps_database/uls_smoothed_31nHz.txt','../cw_maps_database/skymap_pixels.txt')
        skymap100nHz = Skymap('../cw_maps_database/uls_smoothed_100nHz.txt','../cw_maps_database/skymap_pixels.txt')
        skymap318nHz = Skymap('../cw_maps_database/uls_smoothed_318nHz.txt','../cw_maps_database/skymap_pixels.txt')

        interp3nHz= np.array(skymap3nHz.interpolate('all', input_theta, input_phi))
        interp10nHz = np.array(skymap10nHz.interpolate('all', input_theta, input_phi))
        interp31nHz = np.array(skymap31nHz.interpolate('all', input_theta, input_phi))
        interp100nHz = np.array(skymap100nHz.interpolate('all', input_theta, input_phi))
        interp318nHz = np.array(skymap318nHz.interpolate('all', input_theta, input_phi))

        combined_columns = np.array([interp3nHz, interp10nHz, interp31nHz, interp100nHz, interp318nHz])

        #print(interp3nHz)
        #print(interp10nHz)

        #f = np.concatenate((interp1, interp2), axis=1)

        freq_array = np.array([3e-9, 10e-9, 31e-9, 100e-9, 318e-9])
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


        confidence_array = np.array([1, 25, 50, 75, 95])
        orbital_f = 1/(orbital_p*365*24*3600)
        new_freq_rows = np.array([freq_interp1(orbital_f), freq_interp2(orbital_f), freq_interp3(orbital_f), freq_interp4(orbital_f), freq_interp5(orbital_f)])

        #print(new_freq_rows)

        conf_interp = interp1d(new_freq_rows, confidence_array, kind='cubic')

        result_ds = round(float(conf_interp(result_rs*1e-9*orbital_f)), 4)
        PNG_name = uuid4().hex
        print(json.dumps({'result_rs':result_rs,'result_ds':result_ds,'PNG_name':PNG_name}))
except:
    print("Content-Type: text/html")
    print("")  
    print("\n\n<PRE>")
    traceback.print_exc()
