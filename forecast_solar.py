import numpy as np
import pandas as pd
import xarray as xr
from herbie import Herbie
import pvlib
import time


def get_solar_forecast(latitude, longitude, init_date, length_hours,
                       model='gfs', lead_time_hours=0, attempts=2):
    """
    Get a solar resource forecast for a single site from one of several
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

    length_hours : int
        Length of the forecast in hours - number of hours forecasted

    model : String, default 'gfs'
        Forecast model. Default is NOAA GFS ('gfs'), but can also be
        ECMWF IFS ('ifs').

    lead_time_hours : int, optional
        Number of hours between init_date (initialization) and
        the first forecasted interval. NOAA GFS data goes out
        384 hours, so length_hours + lead_time_hours must be less
        than or equal to 384.

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

    # variable formatting
    # if lat, lon are single values, convert to lists for pickpoints later
    if type(latitude) is float or type(latitude) is int:
        latitude = [latitude]
        longitude = [longitude]
    # convert init_date to datetime
    init_date = pd.to_datetime(init_date)

    if model == 'gfs':
        # GFS:
        # 0 to 120 by 1, 123 to 384 by 3
        # runs every 6 hours starting at 00z
        update_freq = '6h'

        # set forecast lead times
        if length_hours + lead_time_hours > 120:
            fxx_range = [*range(lead_time_hours, 120+1, 1),
                         *range(123, length_hours + lead_time_hours + 1, 3)]
        else:
            fxx_range = range(lead_time_hours,
                              length_hours + lead_time_hours, 1)

        # Herbie inputs
        product = 'pgrb2.0p25'
        search_str = 'DSWRF|:TMP:2 m above|[UV]GRD:10 m above'

    elif model == 'ifs':
        # From https://www.ecmwf.int/en/forecasts/datasets/open-data
        # For times 00z &12z: 0 to 144 by 3, 150 to 240 by 6.
        # For times 06z & 18z: 0 to 90 by 3.
        # pick init time based on forecast max lead time:
        if length_hours + lead_time_hours > 90:  # forecast beyond 90 h ahead
            update_freq = '12h'
        else:
            update_freq = '6h'

        # set forecast intervals
        if length_hours + lead_time_hours > 144:
            fxx_range = [*range(lead_time_hours, 145, 3),
                         *range(150, length_hours + lead_time_hours+1, 6)]
        else:
            fxx_range = range(lead_time_hours,
                              length_hours + lead_time_hours + 1, 3)

        # Herbie inputs
        product = 'oper'
        search_str = ':ssrd|10[uv]|2t:sfc'

    # round down to last actual initialization time
    date = init_date.floor(update_freq)

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
            if lead_time_hours >= 120:
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

        # make 5min interval clear sky data covering our time range
        times = pd.date_range(
            start=df.index[0],
            end=df.index[-1],
            freq='5min',
            tz='UTC')

        loc = pvlib.location.Location(
            latitude=latitude[j],
            longitude=longitude[j],
            tz=times.tz
            )
        cs = loc.get_clearsky(times, model=model_cs)

        # calculate average CS ghi over the intervals from the forecast
        # based on list comprehension example in
        # https://stackoverflow.com/a/55724134/27574852
        ghi = cs['ghi']
        dates = df.index
        ghi_clear = [
            ghi.loc[(ghi.index > dates[i]) & (ghi.index <= dates[i+1])].mean()
            for i in range(len(dates) - 1)
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
        # make timestamps center-labeled for instantaneous pvlib modeling later
        df_60min.index = df_60min.index + pd.Timedelta('30min')
        # drop last row, since we don't have data for the last full hour (just
        # an instantaneous end point)
        df_60min = df_60min.iloc[:-1]
        # "backfill" ghi csi
        # merge based on nearest index from 60min version looking forward in
        # 3hr version
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

    # concatenate creating multiindex with keys of the list of point numbers
    # assigned to 'point', reorder indices, and sort by valid_time
    df_60min = (
        pd.concat(dfs, keys=list(range(num_sites)), names=['point'])
        .reorder_levels(["valid_time", "point"])
        .sort_index(level='valid_time')
    )

    # set "point" index as a column
    df_60min = df_60min.reset_index().set_index('valid_time')

    return df_60min
