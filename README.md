# energy-forecasting-tools

This repository is a place for energy forecasting tools. It currently just includes solar-related tools, but may expand to include wind and (maybe one day) electric load.

Look at the notebook [solar_example.ipynb](solar_example.ipynb) for a quick example.

## Suggested environment setup:
Using miniforge:
```
conda create --name energy_forecasting_env python=3.12 -y
conda activate energy_forecasting_env
conda install -c conda-forge herbie-data -y
pip install pvlib ipykernel
```

## References
This project uses several Python packages, including pvlib, an open-source solar PV modeling package [1, 2], and Herbie [3, 4], a package for accessing weather forecast data from NOAA. `pv_model.py` (with the `model_pv_power()` function used here) comes from [5] which leverages some functions from [6].

<img src="images/pvlib_powered_logo_horiz.png" width="200"/>


[1] Anderson, K., Hansen, C., Holmgren, W., Jensen, A., Mikofski, M., and Driesse, A. “pvlib python: 2023 project update.” Journal of Open Source Software, 8(92), 5994, (2023). [DOI: 10.21105/joss.05994](http://dx.doi.org/10.21105/joss.05994).

[2] https://github.com/pvlib/pvlib-python

[3] Blaylock, B. K. (2025). Herbie: Retrieve Numerical Weather Prediction Model Data (Version 2025.3.1) [Computer software]. https://doi.org/10.5281/zenodo.4567540

[4] https://github.com/blaylockbk/Herbie

[5] https://github.com/williamhobbs/pv-system-model

[6] Hobbs, W., Anderson, K., Mikofski, M., and Ghiz, M. "An approach to modeling linear and non-linear self-shading losses with pvlib." 2024 PV Performance Modeling Collaborative (PVPMC). https://github.com/williamhobbs/2024_pvpmc_self_shade 