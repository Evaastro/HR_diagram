# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np 
import matplotlib.pyplot as plt
import astropy.io 
from astroquery.gaia import Gaia

plt.close('all')

#TOP 4 276 690
#ra, dec,source_id, 
query= """SELECT 
TOP 4750
parallax, phot_bp_mean_mag, phot_rp_mean_mag, phot_g_mean_mag
FROM gaiadr2.gaia_source
WHERE phot_bp_mean_mag IS NOT NULL 
AND phot_rp_mean_mag IS NOT NULL
AND parallax >1
AND visibility_periods_used>8
AND astrometric_excess_noise<1
AND parallax_over_error>10
AND phot_g_mean_flux_over_error>50
AND phot_rp_mean_flux_over_error>20
AND phot_bp_mean_flux_over_error>20
ORDER BY random_index
"""
#ORDER BY random_index

job=Gaia.launch_job_async(query)
results=job.get_results()
print(results.columns)

def mag_abs(mag,dist=results['parallax']):
    parsec=1/dist   #dist is in parallax
    A=0.            #absorption
    #return mag-(5*np.log10(parsec))+5-A
    return mag+5+5*np.log10(dist/1000)

mag_abs_g=mag_abs(results['phot_g_mean_mag'],results['parallax'])
diff_BP_RP=(mag_abs(results['phot_bp_mean_mag'])-mag_abs(results['phot_rp_mean_mag']))

#%%

plt.figure()
plt.hexbin(diff_BP_RP,mag_abs_g, cmap='hot',gridsize=50, mincnt=1)
plt.gca().invert_yaxis()
plt.xlabel("diff_BP_RP")
plt.ylabel("mag_abs_g")
plt.xlim(None,5)
plt.savefig("/Users/fink/Desktop/eva/Code/plots/HR_diag_as_article_gaia_2018.pdf")
plt.show()

print("\nProgramme terminer : HR digram\n")
