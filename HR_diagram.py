# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np 
import matplotlib.pyplot as plt
import astropy.io 
from astroquery.gaia import Gaia
from dustmaps.config import config
config.reset()
import dustmaps.gaia_tge
import time 

start_time = time.time()
plt.close('all')

def mag_abs(mag,dist):
    """This function compute the absoluted magnitude"""
    return mag+5+5*np.log10(dist/1000)

def A(EBV,Gbp_0_Grp_0,c1,c2,c3,c4,c5,c6,c7):
    """This function compute the absorption"""
    A0=3.1*EBV
    k=c1+(c2*(Gbp_0_Grp_0))+(c3*(Gbp_0_Grp_0)**2)+(c4*(Gbp_0_Grp_0)**3)+(c5*A0)+(c6*A0**2)+(c7*(Gbp_0_Grp_0)*A0)
    return k*A0

def E(Ax,Ay):
    return Ax-Ay

def A_G(EBV,phot_bp_mean_mag,phot_rp_mean_mag,dist):
    Gbp_first=mag_abs(phot_bp_mean_mag,dist)
    Grp_first=mag_abs(phot_rp_mean_mag,dist)
    Abp_first=A(EBV,Gbp_first-Grp_first,1.1517,-0.0871,-0.0333,0.0173,-0.0230,0.0006,0.0043)
    Arp_first=A(EBV,Gbp_first-Grp_first,0.6104,-0.0170,-0.0026,-0.0017,-0.0078,0.00005,0.0006)
    E_bp_rp=E(Abp_first,Arp_first)
    Gbp_Grp=(Gbp_first-Grp_first)-E_bp_rp
    return A(EBV,Gbp_Grp,0.9761,-0.1704,0.0086,0.0011,-0.0438,0.0013,0.0099)

#%%
query= """SELECT 
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
ORDER BY random_index
"""

#AND a_g_val<0.045 
#AND 1000/parallax <= 200
#AND e_bp_min_rp_val < 0.015
#This is the same criterion with the justification 
#every citation comme from Gaia Data Release 2 Observational Hertzsprung-Russell diagrams of 2018
#SELECT 
#TOP 4276690
#do not change the shape of the diagram

#parallax, bp_rp, phot_g_mean_mag

#FROM gaiadr2.gaia_source
#the figure is done with the DR2

#WHERE visibility_periods_used>8
#Form the article "This removes strong outliers, in particular at the faint end of the local HRD (Arenou et al. 2018). It also leads to more incompleteness, but this is not an issue for this paper."

#AND astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))
#artefact removal 

#AND parallax_over_error>10
#I think it is for the parallax precision. In the artcile is written "We built the Gaia HRDs by simply estimating the absolute
#Gaia magnitude in the G band for individual stars using MG =
#G + 5 + 5 log10(ϖ/1000), with ϖthe parallax in milliarcseconds
#(plus the extinction, see next section). This is valid only when
#the relative precision on the parallax is lower than about 20%
#(Luri et al. 2018). We aim here to examine the fine structures in
#the HRD revealed by Gaia and therefore adopt a 10% relative
#precision criterion, which corresponds to an uncertainty on MG
#smaller than 0.22 mag: parallax_over_error>10."

#AND phot_g_mean_flux_over_error>50
#AND phot_rp_mean_flux_over_error>20
#AND phot_bp_mean_flux_over_error>20
#this is for the removal of variable stars (article argument)

#AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
#AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
#this two criterions is to remove artefact of the HR diagram p3 section 2.1

#AND e_bp_min_rp_val < 0.015 
#this is explain by the article p7 section 4 "E(B−V) < 0.015 mag extinction criteria, as these sources mostly lie within the local bubble"

#ORDER BY random_index
#this do not come from the article but me. Before all the criterion I got the WD part of the diagram

job=Gaia.launch_job_async(query)
results=job.get_results()
print(results.columns)

#%%
from astropy.coordinates import SkyCoord
import astropy.units as u
import dustmaps

gaia=dustmaps.gaia_tge.GaiaTGEQuery()

mag_abs_g=mag_abs(results['phot_g_mean_mag'],results['parallax'])

ra=results['ra'].data
dec=results['dec'].data
coords=SkyCoord(ra=ra*u.degree, dec=dec*u.degree,distance=1/results['parallax']*u.pc, frame='icrs')
E=gaia.query(coords)

#%%

filtre=np.where(results['a_g_val']<=A_G(0.015,results['phot_bp_mean_mag'],results['phot_rp_mean_mag'],results['parallax']))
filtre_E=np.where(E<0.015)

plt.figure()
#plt.hexbin(results['bp_rp'][filtre],mag_abs_g[filtre], cmap='hot',gridsize=50, mincnt=1)
plt.hexbin(results['bp_rp'][filtre_E],mag_abs_g[filtre_E], cmap='hot',gridsize=50, mincnt=1)
#plt.hexbin(results['bp_rp'],mag_abs_g, cmap='hot',gridsize=50, mincnt=1)
plt.gca().invert_yaxis()
plt.xlabel("diff_BP_RP")
plt.ylabel("mag_abs_g")
#plt.xlim(None,5)
#plt.savefig("/Users/fink/Desktop/eva/Code/plots/HR_diag_as_article_gaia_2018.pdf")
plt.show()

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed Time: {elapsed_time} seconds")

print(query)
print("Programme terminer : HR digram\n")
