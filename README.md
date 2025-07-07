# Hertzsprung-Russell Diagram

This repository contains code that allows the reproducibility of an observational Hertzsprung-Russell diagram (HR diagram) using data from Gaia Data Release 2 (DR2).\
You can also use the code to place a star (from Gaia DR3) on this HR diagram using the star's position in equatorial coordinates.

---

## Gaia Query

You need to run a Gaia query on the [Gaia Archive](https://gea.esac.esa.int/archive/) to generate the `results_gaia.csv` file:

```sql
SELECT TOP 1000000
  parallax, bp_rp, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, a_g_val, ra, dec
FROM gaiadr2.gaia_source
WHERE visibility_periods_used > 8
  AND astrometric_chi2_al / (astrometric_n_good_obs_al - 5) < 1.44 * GREATEST(1, EXP(-0.4 * (phot_g_mean_mag - 19.5)))
  AND parallax_over_error > 10
  AND phot_g_mean_flux_over_error > 50
  AND phot_rp_mean_flux_over_error > 20
  AND phot_bp_mean_flux_over_error > 20
  AND phot_bp_rp_excess_factor > 1.0 + 0.015 * POWER(phot_bp_mean_mag - phot_rp_mean_mag, 2)
  AND phot_bp_rp_excess_factor < 1.3 + 0.06 * POWER(phot_bp_mean_mag - phot_rp_mean_mag, 2)
  AND 1000 / parallax <= 200
ORDER BY random_index
```

A new version of this file can be generated using Gaia DR3.

---

## Installation

You can install the required dependencies using:

```bash
python3 -m pip install -r requirements.txt
```

> ⚠️ Some difficulties may be encountered when installing `dustmaps`. Please refer to the [official installation instructions](https://dustmaps.readthedocs.io/en/latest/installation.html).
