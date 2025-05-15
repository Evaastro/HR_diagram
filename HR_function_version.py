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
from dustmaps.config import config
config.reset()
import dustmaps.gaia_tge
import dustmaps
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

chemin="/Users/fink/Desktop/eva/Code/HR_diagram/results_gaia.csv"


def mag_abs(mag,dist,A):
    """This function compute the absoluted magnitude"""
    return mag+5+5*np.log10(dist/1000)+A 

def Absorption(lam, coord):
    gaia_tge=dustmaps.gaia_tge.GaiaTGEQuery()
    E_BV=gaia_tge(coord)
    ext_model = G23(Rv=3.1)
    return ext_model(lam * u.AA) * E_BV

def HR_plot(path=chemin,only_HR=True):
    data=pd.read_csv(path)
    plt.figure(figsize=(4,6))
    plt.hexbin(data['bp_rp'],mag_abs(data['phot_g_mean_mag'],data['parallax'],0), cmap='hot',gridsize=400, mincnt=0.1)
    plt.gca().invert_yaxis()
    plt.xlabel("$G_{BP}-G_{RP}$")
    plt.ylabel("$M_G$")
    plt.tight_layout()
    if only_HR:
        plt.show()

#HR_plot()

def add_star(ra=256.5229102004,dec=-26.5805651308,radius=1* u.arcsec): 
    coord=SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    search=Gaia.cone_search_async(coord, radius=radius)
    result=search.get_results()
    #HR_plot(only_HR=False)
    if len(result)>0:
        for row in result:
            print(f"source_id: {row['source_id']}")
            print(f"G_mag: {row['phot_g_mean_mag']}")
            print(f"BP_mag: {row['phot_bp_mean_mag']}")
            print(f"RP_mag: {row['phot_rp_mean_mag']}")
            print("----")
        A_G=result['ag_gspphot'][0]
        Ebpminrp=result['ebpminrp_gspphot'][0]
        if np.ma.is_masked(A_G): 
            print("Do not have any Gaia Absorption")
            lambda_G=5850.88 #G band 
            print("AG ici",A_G)
            A_G=Absorption(lambda_G, coord)
            print("AG ici apres comp ",A_G)
        if np.ma.is_masked(Ebpminrp):
            print("Do not have the exctinction E(BP-RP)")
            lambda_BP=5041.61 #BP band 
            lambda_RP=7690.74 #RP band
            A_BP=Absorption(lambda_BP, coord)
            A_RP=Absorption(lambda_RP, coord)
            Ebpminrp=A_BP-A_RP
        color=result['bp_rp'][0]+Ebpminrp
        print("absorption coefficient ",A_G)
        mag_g=mag_abs(result['phot_g_mean_mag'][0],result['parallax'][0],A_G)
        mag_gz=mag_abs(result['phot_g_mean_mag'][0],result['parallax'][0],0)
        print("mag_g",mag_g)
        print("mag_gz",mag_gz)
        plt.scatter(color, mag_g, marker='*', s=500)
    else:
        print("Aucun objet Gaia trouvé à ces coordonnées.")
    #plt.show()

print("\nDefault object")
add_star()
#print("\nthis object do not have any absorption in gaia of the G mag")
#add_star(ra=316.390485,dec=51.057441)

#for the plot of the M-dwarf uncatalogued in the 10 fields
"""HR_plot(only_HR=False)
print("\nthis is for a M-dwarf 505213300002885")
add_star(ra=46.01105,dec=14.27167)
print("\nThis is for M-dwarf 505214400002794")
add_star(ra=45.09710,dec=14.29321)
print("\nThis is for M-dwarf 250212300004457")
add_star(ra=34.42344,dec=-23.82088)
print("\nThis is for M-dwarf 293212400005130")
add_star(ra=351.50534,dec=-23.81982)
plt.savefig("/Users/fink/Desktop/eva/Code/plots/HR_with_Mdwarf.pdf")
plt.show()"""


"""This is not working this is the M-dwarf in the 1000 trees"""
HR_plot(only_HR=False)
print("\nThis is for M-dwarf 807211100054997")
add_star(ra=19.81237,dec=63.72292)
print("This is for M-dwarf 795207300005922")
add_star(ra=244.13867,dec=53.32193)
print("This is for M-dwarf 807212400034410")
add_star(ra=15.11587,dec=62.23060)
plt.savefig("/Users/fink/Desktop/eva/Code/plots/HR_with_Mdwarf_1000trees.pdf")
plt.show()

#print("\nbinary of p=10.7383")
#add_star(ra=42.54547,dec=14.76320)
#print("\nbinary of p=9.0082")
#add_star(ra=90.52052, dec=30.15200)
#print("\nOne of the variable ZTF ID 384201200022110")
#add_star(ra=282.56958,dec=-12.15728)

print("\nProgramme terminer : HR function version \n")