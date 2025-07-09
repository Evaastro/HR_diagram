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
Do not forget to download the data in the format csv when generating the data. The format by default of the Gaia query is VOTable. It can be changed before downloading of the data.

A new version of this file can be generated using Gaia DR3.

---

## Installation

To install this version of HR_diagram, you should follow these steps:

Get a local copy of the code with:
```
git clone https://github.com/Evaastro/HR_diagram.git
```

and then place yourself into the HR_diagram directory:
```
cd HR_diagram
```

You can install the required dependencies using:

```bash
python3 -m pip install -r requirements.txt
```

This project also requires the `dustmaps` module, and especially the **Gaia TGE** map from  
[Delchambre et al. (2022)](https://doi.org/10.1051/0004-6361/202243423).

 ⚠️ The installation of `dustmaps` is not always simple.  
Please follow the official documentation:  
[https://dustmaps.readthedocs.io/en/latest/installation.html](https://dustmaps.readthedocs.io/en/latest/installation.html)


## How to Use the Code (Step-by-Step)

The main code is `HR_function.py`.  
It contains four functions:

- `mag_abs()` Computes the absolute magnitude.
- `Absorption()` Computes interstellar absorption depending on wavelength.
- `HR_plot()` Plots the HR diagram.
- `add_star(ra=..., dec=...)` Adds a star from Gaia DR3 using coordinates.

### 1. Run the script

Start by changing the line 44 of the HR_function.py by the path to your data file you generated. 

```python
chemin="results_gaia.csv"
```

Import the functions:

```python
from HR_function import HR_plot, add_star
```

### 2. Plot the HR Diagram

Run:

```python
HR_plot()
```

This will display and save the HR diagram.
📷 Example:

<img src="/HR_empty.png" alt="Example HR Diagram" width="300"/>

### 3. Add a Star (Optional)

To add a specific star from Gaia DR3, use:

```python
add_star(ra=76.85812, dec=53.69585)
```

📷 Example:

<img src="/add_star_exemple.png" alt="Example star HR Diagram" width="300"/>

Replace the values with your star's coordinates.

## License

This project is open-source and available under the MIT License.
