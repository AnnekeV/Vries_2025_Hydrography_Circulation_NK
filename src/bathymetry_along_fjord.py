#%%
from pathlib import Path
from pyproj import Proj
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import cmocean
import pandas as pd
import plotly.graph_objects as go
from plotly.graph_objs.layout import YAxis,XAxis,Margin

cmap = cmocean.cm.deep_r

path_parent = Path.cwd().parent.parent
figpath = (Path.cwd().parent.parent).joinpath(Path('Figures'))
folder = path_parent.joinpath("Data","Other","Bathymetry")
bath_file = "Bathymetry GF1_GF19.csv"
station_path = path_parent.joinpath('Data',"Stations GHF.csv")
f_microcat_depth = path_parent.joinpath(Path('Processing', 'intermediate_files', 'mooring_mean_pressure_levels_microcat.csv'))

df_bath = pd.read_csv(f"{folder}/{bath_file}")
df_stat_all = pd.read_csv(station_path)
df_mc = pd.read_csv(f"{f_microcat_depth}") 


# %%

def generate_wavy_line(x1, x2, amplitude=1, wavelength=5, num_points=500):
    x = np.linspace(x1, x2, num_points)
    y = amplitude * np.sin(2 * np.pi * x / wavelength)
    return x, y

# Example usage:
x1, x2 = 0, 30  # Define the range for the line
x, y = generate_wavy_line(x1, x2, amplitude=1, wavelength=5, num_points=500)

# Create a Plotly trace for the wavy line
trace = go.Scatter(x=x, y=y, mode='lines')

# Create a layout for the plot
layout = go.Layout(
    title='Wavy Line',
    xaxis=dict(title='X'),
    yaxis=dict(title='Y'),
)






# %%

#drop final two stations
df_stat = df_stat_all.set_index('Station').loc[:"GF17"].reset_index()
# keep three stations
df_stat = df_stat_all.set_index('Station').loc[['GF3', 'GF5', 'GF10']].reset_index()




layout = go.Layout(
    title="Bathymetry along main fjord branch",
    xaxis=XAxis(
        title="Distance [km]"
    ),
    xaxis2 = XAxis(
        overlaying= 'x', 
        tickvals = df_stat["Distance"],
        ticktext = df_stat["Station"],
        side= 'top',
        # showgrid=False,
        # minor=dict(ticklen=6, tickcolor="black", showgrid=True),
        showticklabels = True,
    ),
    width=400, height=400)

fig_bath = go.Figure(layout=layout)

import scipy.signal as signal
# Design the Butterworth filter
b, a = signal.butter(2, 1/4, btype='low', analog=False)
# Apply the filter to the data
df_bath["Depth_filtered"] = signal.filtfilt(b, a,df_bath["Depth"].values)



def plot_bathymetry(fig, add_mooring_GF10=False):
    fig.add_trace(go.Scatter(
            x=[df_bath["Distance [km]"].min(), df_bath["Distance [km]"].max(), df_bath["Distance [km]"].max(), df_bath["Distance [km]"].min()],
            y=[-650, -650,-10,-10],
            fill='toself',
            mode='lines',
            line=dict(color="bisque", width=.2),
            fillcolor='bisque',))
            
        # Create a Plotly trace for the wavy line
    x, y = generate_wavy_line(df_bath["Distance [km]"].min(), df_bath["Distance [km]"].max(), amplitude=5, wavelength=10, num_points=500)
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', line=dict(color="lightblue", width=2)))

    fig.add_trace(
        go.Scatter(
            x=df_bath["Distance [km]"],
            y=df_bath["Depth_filtered"],
            name="Profile from BedMachine",
            # mode="lines",
            fill='tonexty',
            line=dict(color='peru', width=2),
            fillcolor='azure',
        ))

    fig.add_trace(
        go.Scatter(
            x=df_stat["Distance"],
            # y=df_stat["DepthStation"]*-1,
            y=np.zeros(len(df_stat["DepthStation"]))-610,
            marker=dict(size=1, color = 'rgba(0,0,0,0)'),
            text=df_stat["Station"],
            mode="markers",
            name = "Station location at depth",
            textposition="bottom center",
            # xaxis = 'x2',

        ))

         


    # fig.update_xaxes( minor=dict(
    #         tickvals = df_stat["Distance"],
            # ticktext = df_stat["Station"])
    #         # ticklen=6, tickcolor="black", showgrid=True)
    #         )
    if add_mooring_GF10:
        fig.add_trace(
        go.Scatter(
            x=np.repeat(94,len(df_mc)),
            y=df_mc.depth*-1,
            name="Microcat mooring GF10"
        ))

    return fig


if __name__ == '__main__':
    plot_bathymetry(fig_bath)
    fig_bath.update_layout(showlegend=False, )
    # fig.write_image(f"{figpath}/bathymetry.svg")
    fig_bath.show()

    plot_bathymetry(fig_bath, add_mooring_GF10=False)
    fig_bath.update_layout(showlegend=False, template = "simple_white" ,
                            title=""
                        )
    # set x range and xticks every 25 km
    fig_bath.update_xaxes(range=[-2, 190], dtick=50)
    # set y range and yticks every 100 m
    fig_bath.update_yaxes(range=[-660, 5], dtick=200)
    fig_bath.update_xaxes(visible=False)
    fig_bath.update_yaxes(visible=False)

    # fig_bath.write_image(f"{figpath}/bathymetry along fjord/bathymetry2_background_no_axes.svg")
    fig_bath.show()


# %%

# Open depth microcat




# %%
