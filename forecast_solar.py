import numpy as np
import pandas as pd
import xarray as xr
from herbie import Herbie, FastHerbie
import pvlib
import time
import warnings


def _model_input_formatter(init_date, run_length, lead_time_to_start=0,
                           model='gfs'):
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

    model : {'gfs', 'ifs', 'hrrr'}
        Forecast model name, case insensitive. Default is 'gfs'.

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
        search_str = 'DSWRF|:TMP:2 m above|[UV]GRD:10 m above'

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
        search_str = ':ssrd|10[uv]|2t:sfc'

    elif model == 'hrrr':
        # maximum forecast horizon
        fxx_max = run_length + lead_time_to_start
        product = 'sfc'
        search_str = 'DSWRF|:TMP:2 m above|[UV]GRD:10 m above'
        update_freq = '1h'

        # round down to last actual initialization time
        date = init_date.floor(update_freq)

        fxx_range = range(lead_time_to_start, fxx_max, 1)

    return date, fxx_range, product, search_str


def get_solar_forecast(latitude, longitude, init_date, run_length,
                       lead_time_to_start=0, model='gfs', attempts=2,
                       hrrr_hour_middle=True, hrrr_coursen_window=None):
    """
    Get a solar resource forecast for one or several sites from one of several
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
        timeseries forecasted solar resource data

    References
    ----------

    .. [1] `Blaylock, B. K. (YEAR). Herbie: Retrieve Numerical Weather
       Prediction Model Data (Version 20xx.x.x) [Computer software].
       <https://doi.org/10.5281/zenodo.4567540>`_
    .. [2] `Anderson, K., et al. “pvlib python: 2023 project update.” Journal
       of Open Source Software, 8(92), 5994, (2023).
       <http://dx.doi.org/10.21105/joss.05994>`_
    """

    # set clear sky model. could be an input variable at some point
    model_cs = 'haurwitz'

    # variable formatting
    # if lat, lon are single values, convert to lists for pickpoints later
    if type(latitude) is float or type(latitude) is int:
        latitude = [latitude]
        longitude = [longitude]
    # convert init_date to datetime
    init_date = pd.to_datetime(init_date)

    # get model-specific Herbie inputs
    date, fxx_range, product, search_str = _model_input_formatter(
        init_date, run_length, lead_time_to_start, model)

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
        ds = ds.herbie.with_wind('speed')

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
    # rename 'ssrd' to 'sdswrf' in ifs
    if model == 'ifs':
        ts = ts.rename({'ssrd': 'sdswrf'})
    # convert to dataframe
    df_temp = ts.to_dataframe()[['sdswrf', 't2m', 'si10']]
    # add timezone
    df_temp = df_temp.tz_localize('UTC', level='valid_time')
    # rename wind speed
    df_temp = df_temp.rename(columns={'si10': 'wind_speed'})
    # convert air temperature units
    df_temp['temp_air'] = df_temp['t2m'] - 273.15

    # work through sites
    dfs = {}  # empty list of dataframes
    if type(latitude) is float or type(latitude) is int:
        num_sites = 1
    else:
        num_sites = len(latitude)

    for j in range(num_sites):
        df = df_temp[df_temp.index.get_level_values('point') == j]
        df = df.droplevel('point')

        loc = pvlib.location.Location(
            latitude=latitude[j],
            longitude=longitude[j],
            tz=df.index.tz
            )

        if model == 'gfs':
            # for gfs ghi: we have to "unmix" the rolling average irradiance
            # that resets every 6 hours
            mixed = df[['sdswrf']].copy()
            mixed['hour'] = mixed.index.hour
            mixed['hour'] = mixed.index.hour
            mixed['hour_of_mixed_period'] = ((mixed['hour'] - 1) % 6) + 1
            mixed['sdswrf_prev'] = mixed['sdswrf'].shift(
                periods=1,
                fill_value=0
                )
            mixed['int_len'] = mixed.index.diff().total_seconds().values / 3600

            # set the first interval length:
            if lead_time_to_start >= 120:
                mixed.loc[mixed.index[0], 'int_len'] = 1
            else:
                mixed.loc[mixed.index[0], 'int_len'] = 3
            unmixed = ((mixed['hour_of_mixed_period'] * mixed['sdswrf']
                        - (mixed['hour_of_mixed_period'] - mixed['int_len'])
                        * mixed['sdswrf_prev']) / mixed['int_len'])
            df['ghi'] = unmixed
        elif model == 'ifs':
            # for ifs ghi: cumulative J/m^s to average W/m^2 over the interval
            # ending at the valid time. calculate difference in measurement
            # over diff in time to get avg J/s/m^2 = W/m^2
            df['ghi'] = df['sdswrf'].diff() / df.index.diff().seconds.values

        elif model == 'hrrr':
            df['ghi'] = df['sdswrf']

        if model == 'gfs' or model == 'ifs':
            # make 1min interval clear sky data covering our time range
            times = pd.date_range(
                start=df.index[0],
                end=df.index[-1],
                freq='1min',
                tz='UTC')

            cs = loc.get_clearsky(times, model=model_cs)

            # calculate average CS ghi over the intervals from the forecast
            # based on list comprehension example in
            # https://stackoverflow.com/a/55724134/27574852
            ghi = cs['ghi']
            dates = df.index
            ghi_clear = [
                ghi.loc[(ghi.index > dates[i]) & (ghi.index <= dates[i+1])]
                .mean() for i in range(len(dates) - 1)
                ]

            # write to df and calculate clear sky index of ghi
            df['ghi_clear'] = [np.nan] + ghi_clear
            df['ghi_csi'] = df['ghi'] / df['ghi_clear']

            # avoid divide by zero issues
            df.loc[df['ghi'] == 0, 'ghi_csi'] = 0

            # 60min version of data, centered at bottom of the hour
            # 1min interpolation, then 60min mean
            df_60min = (
                df[['temp_air', 'wind_speed']]
                .resample('1min')
                .interpolate()
                .resample('60min').mean()
            )
            # make timestamps center-labeled for instantaneous pvlib modeling
            # later
            df_60min.index = df_60min.index + pd.Timedelta('30min')
            # drop last row, since we don't have data for the last full hour
            # (just an instantaneous end point)
            df_60min = df_60min.iloc[:-1]
            # "backfill" ghi csi
            # merge based on nearest index from 60min version looking forward
            # in 3hr version
            df_60min = pd.merge_asof(
                left=df_60min,
                right=df.ghi_csi,
                on='valid_time',
                direction='forward'
            ).set_index('valid_time')

            # make 60min interval clear sky, centered at bottom of the hour
            times = pd.date_range(
                start=df.index[0]+pd.Timedelta('30m'),
                end=df.index[-1]-pd.Timedelta('30m'),
                freq='60min',
                tz='UTC')
            cs = loc.get_clearsky(times, model=model_cs)

            # calculate ghi from clear sky and backfilled forecasted clear sky
            # index
            df_60min['ghi'] = cs['ghi'] * df_60min['ghi_csi']

            # dni and dhi using pvlib erbs. could also DIRINT or erbs-driesse
            sp = loc.get_solarposition(times)
            out_erbs = pvlib.irradiance.erbs(
                df_60min.ghi,
                sp.zenith,
                df_60min.index,
            )
            df_60min['dni'] = out_erbs.dni
            df_60min['dhi'] = out_erbs.dhi

            # add clearsky ghi
            df_60min['ghi_clear'] = df_60min['ghi'] / df_60min['ghi_csi']

            dfs[j] = df_60min

        elif model == 'hrrr':
            if hrrr_hour_middle is True:
                # clear sky index
                times = df.index
                cs = loc.get_clearsky(times, model=model_cs)
                df['csi'] = df['ghi'] / cs['ghi']
                # avoid divide by zero issues
                df.loc[df['ghi'] == 0, 'csi'] = 0

                # make 1min interval clear sky data covering our time range
                times = pd.date_range(
                    start=df.index[0],
                    end=df.index[-1],
                    freq='1min',
                    tz='UTC')

                cs = loc.get_clearsky(times, model=model_cs)
                # calculate 1min interpolated temp_air, wind_speed, csi
                df_01min = (
                    df[['temp_air', 'wind_speed', 'csi']]
                    .resample('1min')
                    .interpolate()
                )
                # add ghi_clear
                df_01min['ghi_clear'] = cs['ghi']
                # calculate hour averages centered labelled at bottom of the
                # hour
                df_60min = df_01min.resample('1h').mean()
                df_60min.index = df_60min.index + pd.Timedelta('30min')
                # calculate new ghi
                df_60min['ghi'] = df_60min['csi'] * df_60min['ghi_clear']

            else:
                df_60min = df.copy()

            # dni and dhi using pvlib erbs. could also DIRINT or erbs-driesse
            sp = loc.get_solarposition(df_60min.index)
            out_erbs = pvlib.irradiance.erbs(
                df_60min.ghi,
                sp.zenith,
                df_60min.index,
            )
            df_60min['dni'] = out_erbs.dni
            df_60min['dhi'] = out_erbs.dhi

            # add clearsky ghi
            cs = loc.get_clearsky(df_60min.index, model=model_cs)
            df_60min['ghi_clear'] = cs['ghi']

            dfs[j] = df_60min.copy()

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
    df_60min = df_60min.drop(['t2m', 'sdswrf'], axis=1, errors='ignore')

    return df_60min


def get_solar_forecast_fast(latitude, longitude, init_date, run_length,
                            lead_time_to_start=0, model='gfs', attempts=2,
                            hrrr_hour_middle=True, hrrr_coursen_window=None):
    """
    Get a solar resource forecast for one or several sites from one of several
    NWPs. This function uses Herbie [1]_ and pvlib [2]_. This version
    uses FastHerbie and may be about 15% faster. It currently only works
    with a single init_date, not a list of dates like FastHerbie can use.

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
        timeseries forecasted solar resource data

    References
    ----------

    .. [1] `Blaylock, B. K. (YEAR). Herbie: Retrieve Numerical Weather
        Prediction Model Data (Version 20xx.x.x) [Computer software].
        <https://doi.org/10.5281/zenodo.4567540>`_
    .. [2] `Anderson, K., et al. “pvlib python: 2023 project update.” Journal
        of Open Source Software, 8(92), 5994, (2023).
        <http://dx.doi.org/10.21105/joss.05994>`_
    """

    # set clear sky model. could be an input variable at some point
    model_cs = 'haurwitz'

    # variable formatting
    # if lat, lon are single values, convert to lists for pickpoints later
    if type(latitude) is float or type(latitude) is int:
        latitude = [latitude]
        longitude = [longitude]
    # convert init_date to datetime
    init_date = pd.to_datetime(init_date)

    # get model-specific Herbie inputs
    date, fxx_range, product, search_str = _model_input_formatter(
        init_date, run_length, lead_time_to_start, model)

    delimiter = '|'
    search_string_list = search_str.split(delimiter)

    i = []
    ds_dict = {}
    FH = FastHerbie([date], model=model, product=product, fxx=fxx_range)
    for j in range(0, len(search_string_list)):
        # get solar, 10m wind, and 2m temp data
        # try n times based loosely on
        # https://thingspython.wordpress.com/2021/12/05/how-to-try-something-n-times-in-python/
        for attempts_remaining in reversed(range(attempts)):
            attempt_num = attempts - attempts_remaining
            try:
                if attempt_num == 1:
                    # try downloading
                    ds_dict[j] = FH.xarray(search_string_list[j],
                                           remove_grib=True)
                else:
                    # after first attempt, set overwrite=True to overwrite
                    # partial files
                    ds_dict[j] = FH.xarray(search_string_list[j],
                                           remove_grib=True,
                                           overwrite=True)
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
        ds = xr.merge(ds_dict.values(), compat='override')
        # calculate wind speed from u and v components
        ds = ds.herbie.with_wind('speed')

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
    # convert to dataframe
    # rename 'ssrd' to 'sdswrf' in ifs
    if model == 'ifs':
        df_temp = i[-1].to_dataframe()[['valid_time', 'ssrd', 't2m', 'si10']]
        df_temp = df_temp.rename(columns={'ssrd': 'sdswrf'})
    else:
        df_temp = i[-1].to_dataframe()[['valid_time', 'sdswrf', 't2m', 'si10']]

    # make 'valid_time' an index with 'point', drop 'step'
    df_temp = (df_temp.reset_index().set_index(['valid_time', 'point'])
               .drop('step', axis=1))

    # add timezone
    df_temp = df_temp.tz_localize('UTC', level='valid_time')
    # rename wind speed
    df_temp = df_temp.rename(columns={'si10': 'wind_speed'})
    # convert air temperature units
    df_temp['temp_air'] = df_temp['t2m'] - 273.15

    # work through sites
    dfs = {}  # empty list of dataframes
    if type(latitude) is float or type(latitude) is int:
        num_sites = 1
    else:
        num_sites = len(latitude)

    for j in range(num_sites):
        df = df_temp[df_temp.index.get_level_values('point') == j]
        df = df.droplevel('point')

        loc = pvlib.location.Location(
            latitude=latitude[j],
            longitude=longitude[j],
            tz=df.index.tz
            )

        if model == 'gfs':
            # for gfs ghi: we have to "unmix" the rolling average irradiance
            # that resets every 6 hours
            mixed = df[['sdswrf']].copy()
            mixed['hour'] = mixed.index.hour
            mixed['hour'] = mixed.index.hour
            mixed['hour_of_mixed_period'] = ((mixed['hour'] - 1) % 6) + 1
            mixed['sdswrf_prev'] = mixed['sdswrf'].shift(
                periods=1,
                fill_value=0
                )
            mixed['int_len'] = mixed.index.diff().total_seconds().values / 3600

            # set the first interval length:
            if lead_time_to_start >= 120:
                mixed.loc[mixed.index[0], 'int_len'] = 1
            else:
                mixed.loc[mixed.index[0], 'int_len'] = 3
            unmixed = ((mixed['hour_of_mixed_period'] * mixed['sdswrf']
                        - (mixed['hour_of_mixed_period'] - mixed['int_len'])
                        * mixed['sdswrf_prev']) / mixed['int_len'])
            df['ghi'] = unmixed
        elif model == 'ifs':
            # for ifs ghi: cumulative J/m^s to average W/m^2 over the interval
            # ending at the valid time. calculate difference in measurement
            # over diff in time to get avg J/s/m^2 = W/m^2
            df['ghi'] = df['sdswrf'].diff() / df.index.diff().seconds.values

        elif model == 'hrrr':
            df['ghi'] = df['sdswrf']

        if model == 'gfs' or model == 'ifs':
            # make 1min interval clear sky data covering our time range
            times = pd.date_range(
                start=df.index[0],
                end=df.index[-1],
                freq='1min',
                tz='UTC')

            cs = loc.get_clearsky(times, model=model_cs)

            # calculate average CS ghi over the intervals from the forecast
            # based on list comprehension example in
            # https://stackoverflow.com/a/55724134/27574852
            ghi = cs['ghi']
            dates = df.index
            ghi_clear = [
                ghi.loc[(ghi.index > dates[i]) & (ghi.index <= dates[i+1])]
                .mean() for i in range(len(dates) - 1)
                ]

            # write to df and calculate clear sky index of ghi
            df['ghi_clear'] = [np.nan] + ghi_clear
            df['ghi_csi'] = df['ghi'] / df['ghi_clear']

            # avoid divide by zero issues
            df.loc[df['ghi'] == 0, 'ghi_csi'] = 0

            # 60min version of data, centered at bottom of the hour
            # 1min interpolation, then 60min mean
            df_60min = (
                df[['temp_air', 'wind_speed']]
                .resample('1min')
                .interpolate()
                .resample('60min').mean()
            )
            # make timestamps center-labeled for instantaneous pvlib modeling
            # later
            df_60min.index = df_60min.index + pd.Timedelta('30min')
            # drop last row, since we don't have data for the last full hour
            # (just an instantaneous end point)
            df_60min = df_60min.iloc[:-1]
            # "backfill" ghi csi
            # merge based on nearest index from 60min version looking forward
            # in 3hr version
            df_60min = pd.merge_asof(
                left=df_60min,
                right=df.ghi_csi,
                on='valid_time',
                direction='forward'
            ).set_index('valid_time')

            # make 60min interval clear sky, centered at bottom of the hour
            times = pd.date_range(
                start=df.index[0]+pd.Timedelta('30m'),
                end=df.index[-1]-pd.Timedelta('30m'),
                freq='60min',
                tz='UTC')
            cs = loc.get_clearsky(times, model=model_cs)

            # calculate ghi from clear sky and backfilled forecasted clear sky
            # index
            df_60min['ghi'] = cs['ghi'] * df_60min['ghi_csi']

            # dni and dhi using pvlib erbs. could also DIRINT or erbs-driesse
            sp = loc.get_solarposition(times)
            out_erbs = pvlib.irradiance.erbs(
                df_60min.ghi,
                sp.zenith,
                df_60min.index,
            )
            df_60min['dni'] = out_erbs.dni
            df_60min['dhi'] = out_erbs.dhi

            # add clearsky ghi
            df_60min['ghi_clear'] = df_60min['ghi'] / df_60min['ghi_csi']

            dfs[j] = df_60min

        elif model == 'hrrr':
            if hrrr_hour_middle is True:
                # clear sky index
                times = df.index
                cs = loc.get_clearsky(times, model=model_cs)
                df['csi'] = df['ghi'] / cs['ghi']
                # avoid divide by zero issues
                df.loc[df['ghi'] == 0, 'csi'] = 0

                # make 1min interval clear sky data covering our time range
                times = pd.date_range(
                    start=df.index[0],
                    end=df.index[-1],
                    freq='1min',
                    tz='UTC')

                cs = loc.get_clearsky(times, model=model_cs)
                # calculate 1min interpolated temp_air, wind_speed, csi
                df_01min = (
                    df[['temp_air', 'wind_speed', 'csi']]
                    .resample('1min')
                    .interpolate()
                )
                # add ghi_clear
                df_01min['ghi_clear'] = cs['ghi']
                # calculate hour averages centered labelled at bottom of the
                # hour
                df_60min = df_01min.resample('1h').mean()
                df_60min.index = df_60min.index + pd.Timedelta('30min')
                # calculate new ghi
                df_60min['ghi'] = df_60min['csi'] * df_60min['ghi_clear']

            else:
                df_60min = df.copy()

            # dni and dhi using pvlib erbs. could also DIRINT or erbs-driesse
            sp = loc.get_solarposition(df_60min.index)
            out_erbs = pvlib.irradiance.erbs(
                df_60min.ghi,
                sp.zenith,
                df_60min.index,
            )
            df_60min['dni'] = out_erbs.dni
            df_60min['dhi'] = out_erbs.dhi

            # add clearsky ghi
            cs = loc.get_clearsky(df_60min.index, model=model_cs)
            df_60min['ghi_clear'] = cs['ghi']

            dfs[j] = df_60min.copy()

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
    df_60min = df_60min.drop(['t2m', 'sdswrf'], axis=1, errors='ignore')

    return df_60min


def get_solar_forecast_ensemble(latitude, longitude, init_date, run_length,
                                lead_time_to_start=0, model='ifs',
                                num_members=10, attempts=2):
    """
    Get an ensemble of solar resource forecasts for one or several sites.
    This function uses Herbie [1]_ and pvlib [2]_. This function
    uses FastHerbie. It currently only works with a single init_date,
    not a list of dates like FastHerbie can use. Temperature data comes
    from the ensemble mean, and wind speed is currently just a filler value
    of 2 m/s to save time.

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
        the first forecasted interval.

    model : string, default 'ifs'
        Forecast model. Default and only option is ECMWF IFS ('ifs'). NOAA
        GEFS may be added in the future.

    num_members : int, default 10
        Number of ensemble members to get. IFS has 50 members.

    attempts : int, optional
        Number of times to try getting forecast data. The function will pause
        for n^2 minutes after each n attempt, e.g., 1 min after the first
        attempt, 4 minutes after the second, etc.

    Returns
    -------
    data : pandas.DataFrane
        timeseries forecasted solar resource data

    References
    ----------

    .. [1] `Blaylock, B. K. (YEAR). Herbie: Retrieve Numerical Weather
        Prediction Model Data (Version 20xx.x.x) [Computer software].
        <https://doi.org/10.5281/zenodo.4567540>`_
    .. [2] `Anderson, K., et al. “pvlib python: 2023 project update.” Journal
        of Open Source Software, 8(92), 5994, (2023).
        <http://dx.doi.org/10.21105/joss.05994>`_
    """

    # set clear sky model. could be an input variable at some point
    model_cs = 'haurwitz'

    # check model
    if model.casefold() != ('ifs').casefold():
        raise ValueError('model must be ifs, you entered ' + model)

    # variable formatting
    # if lat, lon are single values, convert to lists for pickpoints later
    if type(latitude) is float or type(latitude) is int:
        latitude = [latitude]
        longitude = [longitude]
    # convert init_date to datetime
    init_date = pd.to_datetime(init_date)

    num_sites = len(latitude)

    # get model-specific Herbie inputs, except product and search string,
    # which are unique for the ensemble
    init_date, fxx_range, _, _ = _model_input_formatter(
        init_date, run_length, lead_time_to_start, model)

    dfs = []

    # set clear sky model. could be an input variable at some point
    model_cs = 'haurwitz'

    # loop through IFS ensemble members and get GHI data
    for number in range(1, num_members+1):
        search_str = ':ssrd:sfc:' + str(number) + ':'
        # try n times based loosely on
        # https://thingspython.wordpress.com/2021/12/05/how-to-try-something-n-times-in-python/
        for attempts_remaining in reversed(range(attempts)):
            attempt_num = attempts - attempts_remaining
            try:
                if attempt_num == 1:
                    # try downloading
                    ds = FastHerbie(DATES=[init_date],
                                    model='ifs',
                                    product='enfo',
                                    fxx=fxx_range).xarray(search_str)
                else:
                    # after first attempt, set overwrite=True to overwrite
                    # partial files
                    ds = FastHerbie(DATES=[init_date],
                                    model='ifs',
                                    product='enfo',
                                    fxx=fxx_range).xarray(search_str,
                                                          overwrite=True)
            except Exception:
                if attempts_remaining:
                    print('attempt ' + str(attempt_num) + ' failed, pause for '
                          + str((attempt_num)**2) + ' min')
                    time.sleep(60*(attempt_num)**2)
            else:
                break
        else:
            raise ValueError('download failed, ran out of attempts')

        # use pick_points for single point or list of points
        ds2 = ds.herbie.pick_points(pd.DataFrame({
                        "latitude": latitude,
                        "longitude": longitude,
                        }))
        # convert to dataframe
        df_temp = (ds2
                   .to_dataframe()
                   .reset_index()
                   .set_index('valid_time')[['point', 'ssrd']])
        # add timezone
        df_temp = df_temp.tz_localize('UTC', level='valid_time')
        # rename ssrd
        df_temp = df_temp.rename(columns={'ssrd': 'sdswrf'})

        # work through sites (points)
        if type(latitude) is float or type(latitude) is int:
            num_sites = 1
        else:
            num_sites = len(latitude)
        for point in range(num_sites):
            df = df_temp[df_temp['point'] == point].copy()

            loc = pvlib.location.Location(
                latitude=latitude[point],
                longitude=longitude[point],
                tz=df.index.tz
                )

            # convert cumulative J/m^s to average W/m^2
            df['ghi'] = df['sdswrf'].diff() / df.index.diff().seconds.values

            # make 1min interval clear sky data covering our time range
            times = pd.date_range(
                start=df.index[0],
                end=df.index[-1],
                freq='1min',
                tz='UTC')
            cs = loc.get_clearsky(times, model=model_cs)

            # calculate average CS ghi over the intervals from the forecast
            # based on list comprehension example in
            # https://stackoverflow.com/a/55724134/27574852
            ghi = cs['ghi']
            dates = df.index
            ghi_clear = [
                ghi.loc[(ghi.index > dates[i]) & (ghi.index <= dates[i+1])]
                .mean() for i in range(len(dates) - 1)
                ]

            # write to df and calculate clear sky index of ghi
            df['ghi_clear'] = [np.nan] + ghi_clear
            df['ghi_csi'] = df['ghi'] / df['ghi_clear']

            # avoid divide by zero issues
            df.loc[df['ghi'] == 0, 'ghi_csi'] = 0

            # make a dummy column
            df['dummy'] = 0

            # 60min version of data, centered at bottom of the hour
            # 1min interpolation, then 60min mean
            df_60min = (
                df['dummy']
                .resample('1min')
                .interpolate()
                .resample('60min').mean()
            )
            # make timestamps center-labeled for instantaneous pvlib modeling
            # later
            df_60min.index = df_60min.index + pd.Timedelta('30min')
            # drop last row, since we don't have data for the last full hour
            # (just an instantaneous end point)
            df_60min = df_60min.iloc[:-1]
            # "backfill" ghi csi
            # merge based on nearest index from 60min version looking forward
            # in 3hr version
            df_60min = pd.merge_asof(
                left=df_60min,
                right=df.ghi_csi,
                on='valid_time',
                direction='forward'
            ).set_index('valid_time')

            # make 60min interval clear sky, centered at bottom of the hour
            times = pd.date_range(
                start=df.index[0]+pd.Timedelta('30m'),
                end=df.index[-1]-pd.Timedelta('30m'),
                freq='60min',
                tz='UTC')
            cs = loc.get_clearsky(times, model=model_cs)

            # calculate ghi from clear sky and backfilled forecasted clear sky
            # index
            df_60min['ghi'] = cs['ghi'] * df_60min['ghi_csi']

            # dni and dhi using pvlib erbs. could also DIRINT or erbs-driesse
            sp = loc.get_solarposition(times)
            out_erbs = pvlib.irradiance.erbs(
                df_60min.ghi,
                sp.zenith,
                df_60min.index,
            )
            df_60min['dni'] = out_erbs.dni
            df_60min['dhi'] = out_erbs.dhi

            # add clearsky ghi
            df_60min['ghi_clear'] = df_60min['ghi'] / df_60min['ghi_csi']

            # add member number and point, drop dummy column
            df_60min['member'] = number
            df_60min['point'] = point
            df_60min = df_60min.drop(columns=['dummy'])

            # append
            dfs.append(df_60min)

    # convert to dataframe
    df_60min_irr = pd.concat(dfs)

    # get deterministic temp_air
    search_str = ':2t:sfc:g:0001:od:cf:enfo'

    # try n times based loosely on
    # https://thingspython.wordpress.com/2021/12/05/how-to-try-something-n-times-in-python/
    for attempts_remaining in reversed(range(attempts)):
        attempt_num = attempts - attempts_remaining
        try:
            if attempt_num == 1:
                # try downloading
                ds = FastHerbie(DATES=[init_date],
                                model='ifs',
                                product='enfo',
                                fxx=fxx_range).xarray(search_str)
            else:
                # after first attempt, set overwrite=True to overwrite
                # partial files
                ds = FastHerbie(DATES=[init_date],
                                model='ifs',
                                product='enfo',
                                fxx=fxx_range).xarray(search_str,
                                                      overwrite=True)
        except Exception:
            if attempts_remaining:
                print('attempt ' + str(attempt_num) + ' failed, pause for '
                      + str((attempt_num)**2) + ' min')
                time.sleep(60*(attempt_num)**2)
        else:
            break
    else:
        raise ValueError('download failed, ran out of attempts')

    # use pick_points for single point or list of points
    ds2 = ds.herbie.pick_points(pd.DataFrame({
                    "latitude": latitude,
                    "longitude": longitude,
                    }))

    # convert to dataframe
    df_temp = (ds2
               .to_dataframe()
               .reset_index()
               .set_index('valid_time')[['point', 't2m']])
    # add timezone
    df_temp = df_temp.tz_localize('UTC', level='valid_time')

    # convert air temperature units
    df_temp['temp_air'] = df_temp['t2m'] - 273.15

    dfs_temp_air = []
    # work through sites (points)
    if type(latitude) is float or type(latitude) is int:
        num_sites = 1
    else:
        num_sites = len(latitude)
    for point in range(num_sites):
        df = df_temp[df_temp['point'] == point].copy()

        # 60min version of data, centered at bottom of the hour
        # 1min interpolation, then 60min mean
        df_60min_temp_air = (
            df[['temp_air']]
            .resample('1min')
            .interpolate()
            .resample('60min').mean()
        )

        # make timestamps center-labeled for instantaneous pvlib modeling
        # later
        df_60min_temp_air.index = df_60min_temp_air.index + \
            pd.Timedelta('30min')
        # drop last row, since we don't have data for the last full hour
        # (just an instantaneous end point)
        df_60min_temp_air = df_60min_temp_air.iloc[:-1]

        # drop unneeded columns if they exist
        df_60min_temp_air = df_60min_temp_air.drop(['t2m'],
                                                   axis=1,
                                                   errors='ignore')

        # add member number and point, drop dummy column
        # df_60min_temp_air['member'] = pd.NA
        df_60min_temp_air['point'] = point

        # append
        dfs_temp_air.append(df_60min_temp_air)

    # concat
    df_60min_temp_air = pd.concat(dfs_temp_air)

    # final merge
    df_60min = pd.merge(df_60min_irr,
                        df_60min_temp_air,
                        on=['valid_time', 'point'])

    # add generic wind
    df_60min['wind_speed'] = 2

    return df_60min
