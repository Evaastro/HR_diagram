# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np 
import matplotlib.pyplot as plt
import astropy.io 
from astroquery.gaia import Gaia
import time 

start_time = time.time()
plt.close('all')

#TOP 4 276 690
#200 000 000
#phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL 
#TOP 4276690
query= """SELECT 
TOP 4276690
parallax, bp_rp, phot_g_mean_mag
FROM gaiadr2.gaia_source
WHERE visibility_periods_used>8
AND astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))
AND parallax_over_error>10
AND phot_g_mean_flux_over_error>50
AND phot_rp_mean_flux_over_error>20
AND phot_bp_mean_flux_over_error>20
AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
AND e_bp_min_rp_val < 0.015
ORDER BY random_index
"""

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

def mag_abs(mag,dist=results['parallax']):
    return mag+5+5*np.log10(dist/1000)

mag_abs_g=mag_abs(results['phot_g_mean_mag'])

#%%

plt.figure()
plt.hexbin(results['bp_rp'],mag_abs_g, cmap='hot',gridsize=50, mincnt=1)
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
