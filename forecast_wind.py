import pandas as pd
import xarray as xr
from herbie import Herbie
import time
import warnings


def _model_input_formatter(init_date, run_length, lead_time_to_start=0,
                           model='gfs', resource_type='solar'):
    """
    Helper function to format model-specific inputs for Herbie.

    In the case where the user selects an invalid intitialization date, or
    combination of init date and lead time, it tries to update the init date
    and lead time to match a valid init date for the selected model, but this
    hasn't been fully tested.

    Parameters
    ----------
    init_date : pandas-parsable datetime
        Targetted initialization datetime.

    run_length : int
        Length of the forecast in hours.

    lead_time_to_start : int
        Number of hours from the init_date to the first interval in the
        forecast.

    model : {'gfs', 'ifs', 'aifs', 'hrrr'}
        Forecast model name, case insensitive. Default is 'gfs'.

    resource_type : {'solar, 'wind'}
        Resrouce type. Default is 'solar'.

    Returns
    -------
    date : pandas-parsable datetime
        initialization date, rounded down to the last valid date for the given
        model if needed.

    fxx_range : int or list of ints
        fxx (lead time) values.

    product : string
        model product, e.g., 'pgrb2.0p25' for 'gfs'

    search_str : string
        wgrib2-style search string for Herbie to select variables of interest.
    """

    if model == 'gfs':
        # GFS:
        # 0 to 120 by 1, 123 to 384 by 3
        # runs every 6 hours starting at 00z
        update_freq = '6h'
        # round down to last actual initialization time
        date = init_date.floor(update_freq)

        # offset in hours between selected init_date and fcast run
        init_offset = int((init_date - date).total_seconds()/3600)
        lead_time_to_start = lead_time_to_start + init_offset

        # maximum forecast horizon, update with new lead time
        fxx_max = run_length + lead_time_to_start

        # set forecast lead times
        if lead_time_to_start <= 120 and fxx_max > 120:
            fxx_max = round(fxx_max/3)*3
            fxx_range = [*range(lead_time_to_start, 120+1, 1),
                         *range(123, fxx_max + 1, 3)]
        elif lead_time_to_start > 120:
            fxx_max = round(fxx_max/3)*3
            lead_time_to_start = round(lead_time_to_start/3)*3
            fxx_range = range(lead_time_to_start, fxx_max, 3)
        else:
            fxx_range = range(lead_time_to_start, fxx_max, 1)

        # Herbie inputs
        product = 'pgrb2.0p25'
        if resource_type == 'solar':
            search_str = 'DSWRF|:TMP:2 m above|[UV]GRD:10 m above'
        elif resource_type == 'wind':
            search_str = (
                '[UV]GRD:10 m above|[UV]GRD:80 m above|'
                '[UV]GRD:100 m above|:TMP:2 m above|PRES:surface'
            )

    elif model == 'ifs':
        # From https://www.ecmwf.int/en/forecasts/datasets/open-data
        # For times 00z &12z: 0 to 144 by 3, 150 to 240 by 6.
        # For times 06z & 18z: 0 to 90 by 3.
        # From:
        # https://confluence.ecmwf.int/display/DAC/ECMWF+open+data%3A+real-time+forecasts+from+IFS+and+AIFS
        # Product "oper" runs 00z, 12z, 0h to 144h by 3h, 144h to 240h by 6h
        # Product "scda" runs 06z, 18z, 0h to 90h by 3h
        # **BUT**, see https://github.com/blaylockbk/Herbie/discussions/421
        # Starting 2024-11-12 06:00, 'scda' runs to 144h by 3h

        # round to last 6 hours to start
        date = init_date.floor('6h')
        init_offset = int((init_date - date).total_seconds()/3600)
        lead_time_to_start = lead_time_to_start + init_offset
        fxx_max = run_length + lead_time_to_start

        # pick init time based on forecast max lead time:
        # check if 'scda' product is ideal
        if init_date.hour == 6 or init_date.hour == 18:
            if init_date >= pd.to_datetime('2024-11-12 06:00'):
                scda_fxx_max = 144
            else:
                scda_fxx_max = 90
            if fxx_max > scda_fxx_max:  # forecast beyond scda
                update_freq = '12h'  # must use 'oper' runs
                warnings.warn(
                    ("You have specified an init_date which would have mapped "
                     "to a 06z or 18z. Those runs the IFS 'scda' product, and "
                     "'scda' only goes out 144 hours (90h prior to 2024-11-12)"
                     ". You will get forecasts from the 'oper' run 6 hours "
                     "earlier, instead."))
            else:
                update_freq = '6h'  # can use 'oper' or 'scda'
        else:
            update_freq = '6h'  # can use 'oper' or 'scda'
        # round down to last actual initialization time
        date = init_date.floor(update_freq)

        # offset in hours between selected init_date and fcast run
        init_offset = int((init_date - date).total_seconds()/3600)
        lead_time_to_start = lead_time_to_start + init_offset
        if lead_time_to_start > 141:
            run_length = max(run_length, 6)  # make sure it's long enough
        fxx_max = run_length + lead_time_to_start  # update this

        # set forecast intervals
        if lead_time_to_start <= 144 and fxx_max > 144:
            lead_time_to_start = round(lead_time_to_start/3)*3
            fxx_max = round(fxx_max/6)*6
            # make sure it goes to at least the next interval
            fxx_max = max(fxx_max, 150)
            fxx_range = [*range(lead_time_to_start, 145, 3),
                         *range(150, fxx_max + 1, 6)]
        elif lead_time_to_start > 144:
            lead_time_to_start = round(lead_time_to_start/6)*6
            fxx_max = round(fxx_max/6)*6
            fxx_range = range(lead_time_to_start, fxx_max + 1, 6)
        else:
            lead_time_to_start = round(lead_time_to_start/3)*3
            fxx_max = round(fxx_max/3)*3
            fxx_range = range(lead_time_to_start, fxx_max + 1, 3)

        # Herbie inputs
        if date.hour == 6 or date.hour == 18:
            product = 'scda'
        else:
            product = 'oper'

        if resource_type == 'solar':
            search_str = ':ssrd|10[uv]|2t:sfc'
        elif resource_type == 'wind':
            search_str = ':10[uv]|:100[uv]|:2t:sfc|:sp:'

    elif model == 'aifs':
        # From https://www.ecmwf.int/en/forecasts/datasets/set-ix,
        # https://www.ecmwf.int/en/forecasts/dataset/set-x
        # 4 forecast runs per day (00/06/12/18)
        # 6 hourly steps to 360 (15 days)

        # round to last 6 hours to start
        date = init_date.floor('6h')
        init_offset = int((init_date - date).total_seconds()/3600)
        lead_time_to_start = lead_time_to_start + init_offset
        fxx_max = run_length + lead_time_to_start

        update_freq = '6h'
        # round down to last actual initialization time
        date = init_date.floor(update_freq)

        # offset in hours between selected init_date and fcast run
        init_offset = int((init_date - date).total_seconds()/3600)
        lead_time_to_start = lead_time_to_start + init_offset
        if lead_time_to_start > 141:
            run_length = max(run_length, 6)  # make sure it's long enough
        fxx_max = run_length + lead_time_to_start  # update this

        # set forecast intervals
        fxx_range = range(lead_time_to_start, fxx_max + 1, 6)

        # Herbie inputs
        product = 'oper'  # deterministic

        if resource_type == 'solar':
            search_str = ':ssrd|10[uv]|2t:sfc'
        elif resource_type == 'wind':
            search_str = ':10[uv]|:100[uv]|:2t:sfc|:sp:'

    elif model == 'hrrr':
        # maximum forecast horizon
        fxx_max = run_length + lead_time_to_start
        product = 'sfc'

        if resource_type == 'solar':
            search_str = 'DSWRF|:TMP:2 m above|[UV]GRD:10 m above'
        elif resource_type == 'wind':
            search_str = (
                '[UV]GRD:10 m above|[UV]GRD:80 m above|'
                ':TMP:2 m above|PRES:surface'
            )

        update_freq = '1h'

        # round down to last actual initialization time
        date = init_date.floor(update_freq)

        fxx_range = range(lead_time_to_start, fxx_max, 1)

    return date, fxx_range, product, search_str


def get_wind_forecast(latitude, longitude, init_date, run_length,
                      lead_time_to_start=0, model='gfs', attempts=2,
                      hrrr_hour_middle=True, hrrr_coursen_window=None):
    """
    Get a wind resource forecast for one or several sites from one of several
    NWPs. This function uses Herbie [1]_ and pvlib [2]_.

    Parameters
    ----------
    latitude : float or list of floats
        Latitude in decimal degrees. Positive north of equator, negative
        to south.

    longitude : float or list of floats
        Longitude in decimal degrees. Positive east of prime meridian,
        negative to west.

    init_date : pandas-parsable datetime
        Model initialization datetime.

    run_length : int
        Length of the forecast in hours - number of hours forecasted

    lead_time_to_start : int, optional
        Number of hours between init_date (initialization) and
        the first forecasted interval. NOAA GFS data goes out
        384 hours, so run_length + lead_time_to_start must be less
        than or equal to 384.

    model : string, default 'gfs'
        Forecast model. Default is NOAA GFS ('gfs'), but can also be
        ECMWF IFS ('ifs') or NOAA HRRR ('hrrr')

    attempts : int, optional
        Number of times to try getting forecast data. The function will pause
        for n^2 minutes after each n attempt, e.g., 1 min after the first
        attempt, 4 minutes after the second, etc.

    hrrr_hour_middle : bool, default True
        If model is 'hrrr', setting this False keeps the forecast at the
        native instantaneous top-of-hour format. True (default) shifts
        the forecast to middle of the hour, more closely representing an
        integrated hourly forecast that is centered in the middle of the
        hour.

    hrrr_coursen_window : int or None, default None
        If model is 'hrrr', optional setting that is the x and y window size
        for coarsening the xarray dataset, effectively applying spatial
        smoothing to the HRRR model. The HRRR has a native resolution of
        about 3 km, so a value of 10 results in approx. 30 x 30 km grid.

    Returns
    -------
    data : pandas.DataFrane
        timeseries forecasted wind resource data

    References
    ----------

    .. [1] `Blaylock, B. K. (YEAR). Herbie: Retrieve Numerical Weather
       Prediction Model Data (Version 20xx.x.x) [Computer software].
       <https://doi.org/10.5281/zenodo.4567540>`_
    .. [2] `Anderson, K., et al. “pvlib python: 2023 project update.” Journal
       of Open Source Software, 8(92), 5994, (2023).
       <http://dx.doi.org/10.21105/joss.05994>`_
    """

    # # set clear sky model. could be an input variable at some point
    # model_cs = 'haurwitz'

    # variable formatting
    # if lat, lon are single values, convert to lists for pickpoints later
    if type(latitude) is float or type(latitude) is int:
        latitude = [latitude]
        longitude = [longitude]
    # convert init_date to datetime
    init_date = pd.to_datetime(init_date)

    # get model-specific Herbie inputs
    date, fxx_range, product, search_str = _model_input_formatter(
        init_date, run_length, lead_time_to_start, model, resource_type='wind')

    i = []
    for fxx in fxx_range:
        # get solar, 10m wind, and 2m temp data
        # try n times based loosely on
        # https://thingspython.wordpress.com/2021/12/05/how-to-try-something-n-times-in-python/
        for attempts_remaining in reversed(range(attempts)):
            attempt_num = attempts - attempts_remaining
            try:
                if attempt_num == 1:
                    # try downloading
                    ds = Herbie(
                        date,
                        model=model,
                        product=product,
                        fxx=fxx
                        ).xarray(search_str)
                else:
                    # after first attempt, set overwrite=True to overwrite
                    # partial files
                    ds = Herbie(
                        date,
                        model=model,
                        product=product,
                        fxx=fxx
                        ).xarray(search_str, overwrite=True)
            except Exception:
                if attempts_remaining:
                    print('attempt ' + str(attempt_num) + ' failed, pause for '
                          + str((attempt_num)**2) + ' min')
                    time.sleep(60*(attempt_num)**2)
            else:
                break
        else:
            raise ValueError('download failed, ran out of attempts')

        # merge - override avoids hight conflict between 2m temp and 10m wind
        ds = xr.merge(ds, compat='override')
        # calculate wind speed from u and v components
        ds = ds.herbie.with_wind('both')

        if model == 'hrrr' and hrrr_coursen_window is not None:
            ds = ds.coarsen(x=hrrr_coursen_window,
                            y=hrrr_coursen_window,
                            boundary='trim').mean()

        # use pick_points for single point or list of points
        i.append(
            ds.herbie.pick_points(
                pd.DataFrame(
                    {
                        "latitude": latitude,
                        "longitude": longitude,
                    }
                )
            )
        )
    ts = xr.concat(i, dim="valid_time")  # concatenate

    # convert to dataframe, convert names and units
    if model == 'gfs':
        df_temp = ts.to_dataframe()[
            ['si10', 'ws', 'si100', 'wdir10', 'wdir', 'wdir100', 't2m', 'sp']
            ]
        df_temp['t2m'] = df_temp['t2m'] - 273.15
        df_temp.rename(columns={
            'si10': 'wind_speed_10m',
            'ws': 'wind_speed_80m',
            'si100': 'wind_speed_100m',
            'wdir10': 'wind_direction_10m',
            'wdir': 'wind_direction_80m',
            'wdir100': 'wind_direction_100m',
            't2m': 'temp_air_2m',
            'sp': 'pressure_0m',
            }, inplace=True)
    elif model == 'hrrr':
        df_temp = ts.to_dataframe()[
            ['si10', 'ws', 'wdir10', 'wdir', 't2m', 'sp']
            ]
        df_temp['t2m'] = df_temp['t2m'] - 273.15
        df_temp.rename(columns={
            'si10': 'wind_speed_10m',
            'ws': 'wind_speed_80m',
            'wdir10': 'wind_direction_10m',
            'wdir': 'wind_direction_80m',
            't2m': 'temp_air_2m',
            'sp': 'pressure_0m',
            }, inplace=True)
    elif model == 'ifs' or model == 'aifs':
        df_temp = ts.to_dataframe()[
            ['si10', 'si100', 'wdir10', 'wdir100', 't2m', 'sp']
            ]
        df_temp['t2m'] = df_temp['t2m'] - 273.15
        df_temp.rename(columns={
            'si10': 'wind_speed_10m',
            'si100': 'wind_speed_100m',
            'wdir10': 'wind_direction_10m',
            'wdir100': 'wind_direction_100m',
            't2m': 'temp_air_2m',
            'sp': 'pressure_0m',
            }, inplace=True)

    # work through sites
    dfs = {}  # empty list of dataframes
    if type(latitude) is float or type(latitude) is int:
        num_sites = 1
    else:
        num_sites = len(latitude)

    for j in range(num_sites):
        df = df_temp[df_temp.index.get_level_values('point') == j]
        df = df.droplevel('point')

        # 60min version of data, centered at bottom of the hour
        # 1min interpolation, then 60min mean
        df_60min = (
            df
            .resample('1min')
            .interpolate()
            .resample('60min').mean()
        )
        dfs[j] = df_60min

    # concatenate creating multiindex with keys of the list of point numbers
    # assigned to 'point', reorder indices, and sort by valid_time
    df_60min = (
        pd.concat(dfs, keys=list(range(num_sites)), names=['point'])
        .reorder_levels(["valid_time", "point"])
        .sort_index(level='valid_time')
    )

    # set "point" index as a column
    df_60min = df_60min.reset_index().set_index('valid_time')

    # drop unneeded columns if they exist
    # df_60min = df_60min.drop(['t2m', 'sdswrf'], axis=1, errors='ignore')

    return df_60min
