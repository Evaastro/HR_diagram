# Hertzsprung-Russell diagram

This repository contains code allowing reproducibility of a observational Hertzsprung-Russell diagram (HR diagram) using the data of Gaia data release 2 (DR2).
You can also use the code to place a star (Gaia DR3) on this HR diagram with the position of the star in equatorial coordinated. 

# Gaia query 

You need to do a Gaia query on [Gaia archive](https://gea.esac.esa.int/archive/) to construct results_gaia.csv file. 
```
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
```
A new version of this file will be provided from Gaia DR3. 

# Installation 

You can install the other dependencies using pip:
```
python3 -m pip install -r requirements.txt
```
