
import cmocean
from matplotlib import pyplot as plt
import seaborn as sns

def cmocean_to_plotly(cmap, pl_entries):
    h = 1.0/(pl_entries-1)
    pl_colorscale = []

    for k in range(pl_entries):
        C = map(np.uint8, np.array(cmap(k*h)[:3])*255)
        pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])

    return pl_colorscale



from calendar import day_abbr
from datetime import datetime, timedelta
from pathlib import Path
from time import time
from tracemalloc import start
import ctd
from matplotlib import dates
import matplotlib.dates as mdates
import gsw
import matplotlib.pyplot as plt
import matplotlib.colors
import pandas as pd
from geopy import Point
import numpy as np
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
# plt.style.use('ggplot')
# pd.options.plotti÷ng.backend = "matplotlib"

# WIND
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.calc import wind_components, wind_direction, wind_speed
from metpy.units import units








#%% cmd shift space is docu

path_parent = Path.cwd().parent.parent
figpath = path_parent.joinpath(Path('Figures'))
fname = "/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/Data/Moorings/20190612_SBE_GF10/20180828_SBE2989_5m.cnv" 

fpath = str(path_parent)+"/Data/Moorings/20190612_SBE_GF10/"
fpath_GF13 = str(path_parent.joinpath("Data","Moorings","20190802_SBE_GF13")) +"/"
fpath_GF10 = fpath


def open_cnv(fname, remove_5m=True):
    ''' Open cnv file and export dataframe down and metadata and cast'''
    cast = ctd.from_cnv(fname)#
    cast = rename_variables(cast, remove_5m)
    down, up = cast.split()
    metadata = cast._metadata
    return down,metadata,cast

def rename_variables(df,remove_5m):
    '''renames the variables as in cnv to shorter and equal names  AND REMOVES OUTLIERS'''
    df = df.rename(columns={'potemp090C':'temp', 'sal00':'sal', 
                        'density00':'dens', 'timeJV2':'time','timeJ':'time', 
                        'sigma-t00':"sigma-dens", 'sigma-È00':"sigma-dens", "sigma-�00":"sigma-dens"})
    df = df.reset_index()
    df = remove_above_zero(df)
    if remove_5m: 
        df=remove_outliers_5m(df)
        df=remove_outliers_5m(df)  #twice for the best result

    return df 
def remove_above_zero(cast):
    '''remove all vales lower than 0.1 dbar'''
    rows_before = len(cast.index)
    start = (cast["Pressure [dbar]"]>0.15).idxmax()
    cast = cast.iloc[start:]
    cast = cast[cast["Pressure [dbar]"]>0.15].reset_index()
    print(f"Nr. rows below water surface (dbar>0.15):{rows_before-len(cast.index)}")
    return cast

def remove_outliers_5m(cast, nr_per_hour=6):
    '''Removes outliers based on temp and sal by finding a x times std over 4 weeks
    nr_per_hour  = nr obs per hour, default every 10 min, is 6 times per hour, 
    with built in safety so you don't do it for other depth than 5 m'''
    if cast["Pressure [dbar]"].mean()<10:
        df5 = cast[['temp', 'sal']]
        window = nr_per_hour*24*28 #4 weeks
        max_z =3  # max std dev
        z_scores_alt = (df5-df5.rolling(window, center=True, min_periods=window//2).mean())/df5.rolling(window, center=True, min_periods=window//2).std()
        abs_z_scores = z_scores_alt.abs()
        filtered_entries = (abs_z_scores < max_z).all(axis=1)
        new_df = df5[filtered_entries]
        print(f"{len(df5) - len(new_df)} values removed, {(len(df5) - len(new_df))/len(new_df):.3f} %")
        return cast[filtered_entries]
    else: return cast

def remove_outliers(cast, nr_per_hour=6):
    '''Removes outliers based on temp and sal by finding a x times std over 4 weeks
    nr_per_hour  = nr obs per hour, default every 10 min, is 6 times per hour'''
    df5 = cast[['temp', 'sal']]
    window = nr_per_hour*24*28 #4 weeks
    max_z =3  # max std dev
    z_scores_alt = (df5-df5.rolling(window, center=True, min_periods=window//2).mean())/df5.rolling(window, center=True, min_periods=window//2).std()
    abs_z_scores = z_scores_alt.abs()
    filtered_entries = (abs_z_scores < max_z).all(axis=1)
    new_df = df5[filtered_entries]
    print(f"{len(df5) - len(new_df)} values removed, {(len(df5) - len(new_df))/len(new_df):.3f} %")
    return cast[filtered_entries]

def time_to_date(fname, cast, start_time='Jan 1 2018'):
    '''Changes nr of Julian days to datetime object
    Insert dataframe with at least one column 'time'  and get adjusted dataframe back
    '''
    
    dt_start_time = datetime.strptime(start_time, '%b %d %Y')
    cast['timedelta'] = cast.time.apply(lambda x: timedelta(x))
    cast['date'] = cast['timedelta'] + dt_start_time


def to_gsw(df, lon, lat):
    """Converts temperature and salinity to absolute salinity and conservative temperature"""
    df["SA"] = gsw.SA_from_SP(df["sal"].to_numpy(), df["Pressure [dbar]"].to_numpy(), lon, lat)
    df["CT"] = gsw.CT_from_t(df["SA"].to_numpy(), df["temp"].to_numpy(), df["Pressure [dbar]"].to_numpy())
    df['depSM'] = gsw.z_from_p(df['Pressure [dbar]'].to_numpy(), lat)
    return df


def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))


def TSdiagram(ax,Svals,Tvals,levels, alpha=1, fontsize=8): #,month):
    """
    Make a TS diagram background: sigma0 contours as function of SA and CT
    """
    Tg,Sg = np.meshgrid((Tvals),(Svals))
    sigma = gsw.sigma0(Sg,Tg)
    
    CS = ax.contour(Sg, Tg, sigma, linestyles='dashed', colors='grey', zorder=1, levels=levels, alpha=alpha)
    plt.clabel(CS, fontsize=fontsize, inline=1, fmt='%0.2f')



# Generate some data for the contour plot
salmin, salmax = 20, 34.5
tempmin, tempmax = -1, 4

def TS_plotly(salmin, salmax, tempmin, tempmax, Nlevels=10):
    x = np.linspace(salmin, salmax, 100)
    y = np.linspace(tempmin, tempmax, 100)
    SA = gsw.SA_from_SP(x,100, -50,64)
    CT = gsw.CT_from_pt(SA, y)
    X, Y = np.meshgrid(SA, CT)
    Z = gsw.density.sigma0(X, Y)

    fig = go.Figure(data =
    go.Contour(
        x=x,
        y=y,
        z=Z,
        contours_coloring='fill',
        colorscale=['rgba(255, 255, 255, 0)', 'rgba(255, 255, 255, 0)'],  # all levels the same color (white) and fully transparent
        hoverinfo='none',
        line=dict(width=1, color='rgba(0, 0, 0, 0.5)'), # contour line color    
        showscale=False, # don't show a color scale for this trace
        contours=dict(
            showlabels = True, 
            size=(np.max(Z)-np.min(Z))/Nlevels,
        ),
    )
)
    
    return fig



# %%
#============================================================================== 
# WIND
#==============================================================================

# find great circle distance between two coordinates
def great_circle_distance(lon1, lat1, lon2, lat2):
    '''Find great circle distance between two coordinates
    Takes longitude lon1, latitude lat1, longitude lon2 and latitude lat2 as input
    Returns the great circle distance between the two coordinates in kilometers'''
    # Convert to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    # Find distance
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km


# add grid lines for every degree of latitude and longitude and add labels
def add_grid_labels(ax):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

# add weather station to plot   
def add_weather_station(ax, lon, lat, name, color='red'):
    '''Add weather station to plot
    Takes an axis ax, longitude lon, latitude lat, name name and color color as input
    Returns an axis ax'''
    ax.plot(lon, lat,'o', color=color, markersize=5, transform=ccrs.PlateCarree())
    if hasattr(lon, "__len__"):
        if len(name)==1:
            name = np.repeat(name, len(lon))
        for i in range(len(lon)):
            ax.text(lon[i], lat[i], name[i], color=color, transform=ccrs.PlateCarree())
    else:
        ax.text(lon, lat, name, color=color, transform=ccrs.PlateCarree())
    return ax
def find_nearest_point(ds, lon, lat):
    '''Find nearest point in dataset
    Takes a dataset ds, longitude lon and latitude lat as input
    Returns the nearest point in dataset ds'''
    # Find nearest point
    lon_idx = (np.abs(ds.longitude - lon)).argmin()
    lat_idx = (np.abs(ds.latitude - lat)).argmin()
    # Extract nearest point
    nearest_point = ds.isel(latitude=lat_idx, longitude=lon_idx)
    return nearest_point

def find_nearest_point_CARRA_grid(ds, lon, lat):
    '''
    FIND CLOSEST GRID POINT TO NUUK DMI STATION FOR CARRA AND RETURN DS
    takes ds: xarray dataset of CARRA
    lon: longitude of point
    lat: latitude of point
    for dataset with latidue and longitude as variables
    '''
    ds['distance'] = great_circle_distance(ds.longitude, ds.latitude, lon2=lon, lat2= lat)
    return ds.where(ds.distance==ds.distance.min(), drop=True).squeeze()

# Make title in  the middle 
def fix_title_loc(fig):
    fig.update_layout(title_x=0.5)

def square_plot(fig):
    fig.update_yaxes(
    scaleanchor = "x",
    scaleratio = 1,
  )
    

def compute_wind_comp(ds, speed='si10', direction='wdir10'):
    if isinstance(ds, pd.DataFrame):
        df_units = pandas_dataframe_to_unit_arrays(w, {speed: units('m/s'), direction: units.deg})
        ds['u'], ds['v'] = wind_components(df_units[speed], df_units[direction])
    else:
        ds['u'], ds['v'] = wind_components(ds[speed]*units('m/s'), ds[direction]*units.deg)
    return ds

# =============================================================================
# # CREATE LAGGED CORRELATION MATRIX PLOT
# =============================================================================
def df_derived_by_shift(df,lag=0,NON_DER=[]):
    df = df.copy()
    if not lag:
        return df
    cols ={}
    for i in range(1,lag+1):
        for x in list(df.columns):
            if x not in NON_DER:
                if not x in cols:
                    cols[x] = ['{}_{}'.format(x, i)]
                else:
                    cols[x].append('{}_{}'.format(x, i))
    for k,v in cols.items():
        columns = v
        dfn = pd.DataFrame(data=None, columns=columns, index=df.index)    
        i = 1
        for c in columns:
            dfn[c] = df[k].shift(periods=i)
            i+=1
        df = pd.concat([df, dfn], axis=1)
    return df

def plot_correlation_matrix_lagged(df, lag, columns_not_correlated=[]):
    ''' Plots correlation matrix of a dataframe with lagged columns
    Takes a dataframe df (two columns) and calculates correlations between columns and lagged columns 
    lag is lag step
      columns_not_correlated is a list of columns in df that should not be correlated '''
    df_new = df_derived_by_shift(df, lag, columns_not_correlated)
    plt.figure(figsize=(15,10))
    plt.title(u'6 hours', y=1.05, size=16)
    mask = np.zeros_like(df_new.corr())
    mask[np.triu_indices_from(mask)] = True

    svm = sns.heatmap(df_new.corr(), mask=mask, linewidths=0.1, 
            square=True, cmap= plt.cm.RdBu_r, linecolor='white', annot=True, vmin=-1, vmax=1)
    return svm  



# =============================================================================
# DATES AND TIMES
# =============================================================================


def np_to_datetime(date):
    """
    Converts a numpy datetime64 object to a python datetime object 
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    """
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    return datetime.utcfromtimestamp(timestamp)

from scipy.signal import butter,filtfilt# Filter requirements.
def butter_lowpass_filter(data, cutoff, fs, order):
    '''Low pass filter
    data = data to be filtered
    cutoff = cutoff frequency, Hz
    fs = sample rate, Hz 
    order = order of the filter, # sin wave can be approx represented as quadratic, so 2
    returns filtered data
    '''
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y
