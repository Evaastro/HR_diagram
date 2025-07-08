#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  2 14:06:10 2025

@author: fink
"""
import numpy as np
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import dustmaps
from dustmaps.config import config
import dustmaps.gaia_tge
from dust_extinction.parameter_averages import G23


plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'figure.dpi': 350
})

"""This is the query I use to create the file results gaia
SELECT TOP 1000000
parallax, bp_rp, phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,a_g_val, ra, dec
FROM gaiadr2.gaia_source
WHERE visibility_periods_used>8
AND astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))
AND parallax_over_error>10
AND phot_g_mean_flux_over_error>50
AND phot_rp_mean_flux_over_error>20
AND phot_bp_mean_flux_over_error>20
AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
AND 1000/parallax <= 200
ORDER BY random_index
"""

chemin="results_gaia.csv"

def mag_abs(mag,dist,A):
    """
    Computation of the absolute magnitude

    Args:
        mag (float): apparent magnitude
        dist (float): distance in parallax (mas)
    Returns:
        float: absolute magnitude
    """
    return mag+5+5*np.log10(dist/1000)-A

def Absorption(lam, coord):
    """
    Computation of the interstellar absorption in function of lambda

    Args:
        lam (float): wavelenght in amstrong units
        coord (astropy.coordinates): coordinates
    Returns:
        float: interstellar absorption in function of lambda
    """
    gaia_tge=dustmaps.gaia_tge.GaiaTGEQuery()
    E_BV=gaia_tge(coord)
    ext_model=G23(Rv=3.1)
    return ext_model(lam * u.AA) * E_BV

def HR_plot(path=chemin,only_HR=True):
    """
    Plot of the Hertzsprung-Russell diagram

    Args:
        path (str): path to the file.cvs after Gaia query
        only_HR (bool): if True display the figure
    Returns:
        None: plot of the Hertzsprung-Russell diagram
    """
    data=pd.read_csv(path)
    plt.figure(figsize=(5,6))
    plt.hexbin(data['bp_rp'],mag_abs(data['phot_g_mean_mag'],data['parallax'],0), cmap='hot',gridsize=400, mincnt=0.1)
    plt.gca().invert_yaxis()
    plt.xlabel("$G_{BP}-G_{RP}$")
    plt.ylabel("$M_G$")
    plt.tight_layout()
    if only_HR:
        plt.show()


def add_star(ra=76.85812,dec=53.69585,radius=1* u.arcsec):
    """
    On the Hertzsprung-Russell diagram is the addition of the stars

    Args:
        ra (float): right ascension
        dec (float): declinaison
        radius (float): radius for the cone search
    Returns:
        None: display the Hertzsprung-Russell diagram with the added star
    """
    coord=SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    search=Gaia.cone_search_async(coord, radius=radius)
    result=search.get_results()
    HR_plot(only_HR=False)
    if np.ma.is_masked(result['parallax'][0]):
        return print("No Gaia distances found !")
    if len(result)>0:
        A_G=result['ag_gspphot'][0]
        Ebpminrp=result['ebpminrp_gspphot'][0]
        dist=1000/result['distance_gspphot'][0]
        if np.ma.is_masked(dist):
            dist=result['parallax'][0]
        coord=SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree),\
                       distance=1/dist*u.pc, frame='icrs')
        if np.ma.is_masked(A_G):
            lambda_G=5850.88 #G band
            A_G=Absorption(lambda_G, coord)
        if np.ma.is_masked(Ebpminrp):
            lambda_BP=5041.61 #BP band
            lambda_RP=7690.74 #RP band
            A_BP=Absorption(lambda_BP, coord)
            A_RP=Absorption(lambda_RP, coord)
            Ebpminrp=A_BP-A_RP
        color=result['bp_rp'][0]-Ebpminrp
        mag_g_GSP=mag_abs(result['phot_g_mean_mag'][0],1000/result['distance_gspphot'][0],A_G)
        mag_g_gaia=mag_abs(result['phot_g_mean_mag'][0],result['parallax'][0],A_G)
        mag_gz=mag_abs(result['phot_g_mean_mag'][0],result['parallax'][0],0)
        mag_g=mag_abs(result['phot_g_mean_mag'][0],dist,A_G)
        plt.scatter(color, mag_g, marker='*', s=500,color='#3FE93F',edgecolor='black', linewidth=0.7)
    else:
        print("No Gaia object it is found at theses coordinates.")
    plt.show()

print("\nProgramme terminer : HR function \n")
