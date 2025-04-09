Forecaster notes:

- speeding up ensemble:
  - GEFS: only get every 6 hours. because it's rolling avg, no daily energy info is lost
  - only get irradiance - assume temp and wind speed can be ignored
    - this should allow fastherbie to be used efficiently
  - version history: https://vlab.noaa.gov/web/emc/gefs
  - GEFSv12 is based on GFSv15.1 (https://journals.ametsoc.org/view/journals/mwre/150/3/MWR-D-21-0245.1.xml)
- GFS (and GEFS) irradiance issues:
  - bug in solar zenith angle was discovered in Sep 2018 and fixed in GFSv15.1. It caused an overestimation of solar energy [https://www.emc.ncep.noaa.gov/emc/docs/FV3GFS_OD_Briefs_10-01-18_4-1-2019.pdf, pg 92, 115, accessed via https://www.emc.ncep.noaa.gov/emc/pages/numerical_forecast_systems/gfs.php, link to "Details on changes in GFS V15")
    - there *may* have been additional radiation improvements in GFSv16, implemented March 2021 (https://www.emc.ncep.noaa.gov/emc/pages/numerical_forecast_systems/gfs/implementations.php, "Updated radiation package...")
  - GEFSv12 is based on GFSv15.1 (https://journals.ametsoc.org/view/journals/mwre/150/3/MWR-D-21-0245.1.xml)
  - version history: https://vlab.noaa.gov/web/emc/gfs-docs
- ECMWF IFS
  - 'oper', high-resolution forecast, atmospheric fields, only available at 00z and 12z runs, https://herbie.readthedocs.io/en/stable/gallery/ecmwf_models/ecmwf.html#Products
    - step=0h to 144h by 3h, 144h to 240h by 6h 
  - 'scda', short cut-off high-resolution forecast, atmospheric fields, only available at 06z and 18z runs
    - step=0h to 90h by 3h
	- NOTE: available out to 144h starting 2024-11-12 06:00
- NBM: 
  - odd schedule: very different fxx list for different cycles based on https://www.nco.ncep.noaa.gov/pmb/products/blend/
    - see https://vlab.noaa.gov/web/mdl/nbm-data-availability-v4.2
- wind resource forecast notes:
  - NBM has 10, 30, 80 m wind speed and direction (not u and v, WIND and WDIR), e.g., https://www.nco.ncep.noaa.gov/pmb/products/blend/conus/00/blend.t00z.core.f012.co.grib2.shtml
  - GFS has 10, 20, 30, 40, 50, 80, 100 m wind
    - 10, 80, 100 m temp
	- 80 m pres, specific humidity 
	- surface pressure, 2m temp
	- e.g., https://www.nco.ncep.noaa.gov/pmb/products/gfs/gfs.t00z.pgrb2.0p25.f003.shtml
  - HRRR has 10 and 80 m wind (UGRD and VGRD)
    - surface pressure, 2m temp
    - https://www.nco.ncep.noaa.gov/pmb/products/hrrr/hrrr.t00z.wrfsfcf02.grib2.shtml
  - GEFS
    - 80 m temp, pres, wind (U and V)
	- 100 m temp, wind (U and V)
	- https://www.nco.ncep.noaa.gov/pmb/products/gens/gec00.t00z.pgrb2b.0p50.f003.shtml
	- Note: these are in "pgrb2b", secondary parameters, not "pgrb2a"
  - ECMWF IFS has 10 and 100 m wind (U and V)
    - https://www.ecmwf.int/en/forecasts/datasets/open-data
  
sdswrf
sdswrf_prev: sdswrf from previous interval
int_len: interval length, hours
hr: hour of day (UTC)
hour_of_mixed_period: hours since initialization
irr_int_ave: interval-average irradiance (using end of interval, aka right labeled, timestamps)

irr_int_ave = (hour_of_mixed_period * sdswrf - (hour_of_mixed_period - int_len)*sdswrf_prev)/int_len

Check units:

((h * w*m^-2) - (h * w*m^-2)) / h = w*m^-2 [GOOD!]


=(E3*B3-(E3-F3)*B2)/F3