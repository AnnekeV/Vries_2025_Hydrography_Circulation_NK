
import pandas as pd
import gsw
import xarray as xr
import numpy as np
from paths import *
from bathymetry_along_fjord import df_bath, df_stat


df_all = pd.read_csv(file_monthly_ctd_all_years_all_stations) # Load the data via pahts
df_all = df_all.dropna(subset=['Distance']) # Dropping when not in main transect
df_all = df_all.rename(columns={'Pressure [dbar]':'Pressure'})
df_all['Depth'] = gsw.z_from_p(df_all.Pressure.values, df_all.Latitude.values).round(1)  # give accurate depth
df_all['Plot_date'] = df_all['Date'].copy()
df_all['SA'] = gsw.SA_from_SP(df_all['Salinity [PSU]'].values, df_all['Pressure'].values, lon=-51.4, lat=64)
df_all['CT'] = gsw.CT_from_pt(df_all['SA'].values, df_all['Potential temperature [°C]'].values)
df_all['Depth'] =df_all['Pressure']*-1
df_all['Sigma_dens'] = gsw.density.rho(SA=df_all['SA'], CT=df_all['CT'], p=df_all['Pressure'])-1000


def plot_individual_transect(fig, df, x, z, variable_fill, fill_levels, stations=None, variable_line=None, line_levels = [1026, 1029, 0.5], colorscale='Plasma', clabel=" ", row=1, col=1):
    '''
    colorscale = colorscale for background contour. default is plasma
    clabel = title colorbar, default = ""
    '''
    fig.add_trace(go.Contour(
        x=df[x], 
        y=df[z], 
        z=df[variable_fill],
        connectgaps=False,
        line_width=0.1,
        colorscale=colorscale,
        contours=dict( start=fill_levels[0],end=fill_levels[1], size=fill_levels[2]),             
        colorbar = dict(title=clabel, titleside='right',),
        # hoverinfo=[df[variable_fill], df['Date']],
        line_smoothing=1,
        ), row, col)
    if variable_line is not None:
        fig.add_trace(go.Contour(
            x=df[x], 
            y=df[z], 
            z=df[variable_line],
            connectgaps=False,
            contours_coloring='lines', 
            # line_width=1,
            colorscale=[[0,'black'], [1, 'black']],
            contours=dict(showlabels=True, start=line_levels[0],end=line_levels[1], size=line_levels[2]),
            showscale=False,
            ), row=row, col=col)
    if stations is not None:
        fig.add_trace(go.Scatter(x=stations["Loc"], 
                                 y=np.repeat(df[z].min(), len(stations["Loc"])), 
                                 mode="markers+text", 
                                 marker_symbol="triangle-down", marker_color="black", marker_line_color="white", marker_line_width=1), row, col)
        for i in range(len(stations['Loc'])):
            fig.add_annotation(x=stations["Loc"][i], 
                           y=df[z].min(), 
                            text=stations["Names"].iloc[i],
                showarrow=False,
                yshift=10, row=row,col=col)
        
    # fig.update_layout(hovermode="y unified")
    return fig

def plot_transect_from_big_file(date, variable_fill, plot_locations=True, levels_CT=[1.0,3.4,0.2], fig=None, variable_line="Sigma_dens", line_levels= [27, 30, 0.25],) :

    if fig == None:
        fig = make_subplots(rows=1, cols=1)
    dfdate = df_all[(df_all.Plot_date==date)]
    df = dfdate[(df_all.Pressure>100)]

    station_locations=None
    if plot_locations: 
        station_locations, idx = np.unique(df.Distance, return_index=True)
        station_names = df['St.'].iloc[idx]
        stations = {"Loc":station_locations, "Names":station_names}
    

    if variable_fill == 'CT':
        fig = plot_individual_transect(fig, df, 'Distance', "Pressure", 'CT', variable_line= variable_line, fill_levels=levels_CT, colorscale='Plasma', stations=stations, clabel= "Conservative Temperature [°C]", line_levels=line_levels)
    if variable_fill == 'SA':
        fig = plot_individual_transect(fig, df, 'Distance', "Pressure", 'SA', variable_line= variable_line, fill_levels=[33.2,33.70, 0.05], line_levels = line_levels, colorscale='Viridis', stations=stations, clabel= "Absolute Salinity")
    fig.update_layout(title=date, title_x=0.5)
    fig.update_yaxes(title="Pressure [dbar]")
    fig.update_xaxes(title="Distance [km]")
    return fig, dfdate


