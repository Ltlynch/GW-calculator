#!/usr/bin/env python36
import sys
import traceback
import cgi
import cgitb;
import json
import numpy as np
import scipy.constants as const
import healpy as hp
from astropy.cosmology import WMAP5

sys.stderr = sys.stdout
cgitb.enable()

try:
    #from astropy import units as u

    form = cgi.FieldStorage()

    print("Content-Type: application/json")
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
        redshift = (distance*H0)/const.c

        print(json.dumps({'redshift':redshift}))

    elif action == 'distancecalc':
    	
        redshift = float(form.getfirst("redshift",""))
        distance = WMAP5.comoving_distance(redshift)	

        print(json.dumps({'distance':distance}))

    elif action == 'orbitalsepcalc':

        orbital_p = float(form.getfirst("orbital_p",""))
        total_m = float(form.getfirst("total_m",""))
        orbital_s = s_to_yrs(np.sqrt((4 * (const.pi**2) * (pc_to_m(orbital_s)/2)**3)/(const.G * M0_to_kg(total_m))))

        print(json.dumps({'orbital_s':orbital_s}))

    elif action == 'orbitalpercalc':

        orbital_s = float(form.getfirst("orbital_s",""))
        total_m = float(form.getfirst("total_m",""))
        orbital_p = 2*np.cbrt((const.G * M0_to_kg(total_m) * yrs_to_s(orbital_p)**2)/(4*(const.pi)**2))

        print(json.dumps({'orbital_p':orbital_p}))

    elif action == 'totalmasscalc':	
        
        orbital_p = float(form.getfirst("orbital_p",""))
        orbital_s = float(form.getfirst("orbital_s",""))
        total_m = kg_to_M0((4*(const.pi**2)*(pc_to_m(orbital_s/2))**3)/(const.G * (yrs_to_s(orbital_p))**2))

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
        position_ra = float(form.getfirst("position_ra",""))
        position_dec = float(form.getfirst("position_dec",""))

        result_rs = 190 * (((total_m)/10e9)**(5./3.)) * ((orbital_p)**(1./3.)) * (100/distance)

        data = np.genfromtxt('11yr_skymap_v4.txt', names=True)	    
        nside = 8
        npix = hp.nside2pix(nside)
        indices = hp.ang2pix(nside, data['theta'], data['phi'])
        skymap = np.zeros(npix, dtype=np.float)
        skymap[indices] += data['ul'][indices]
        ul_value = hp.pixelfunc.get_interp_val(skymap, theta, phi)


        result_ds = 56.3  # This last thing is not quite done yet. The last thing to do it to connect the ul_value above, with the residual signature in order to produce a detection probability.

        print(json.dumps({'result_rs':result_rs,'result_ds':result_ds}))
except:
    print("Content-Type: text/html")
    print("")  
    print("\n\n<PRE>")
    traceback.print_exc()
