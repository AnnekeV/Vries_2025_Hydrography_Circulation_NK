
from pathlib import Path

path_parent = Path('/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/Vries_2025_Hydrography_Circulation_NK')
figpath = path_parent.joinpath(Path('outputs'))
path_outputs = figpath
datapath = path_parent.joinpath(Path('data'))


# Used in the velocity plotting notebook

# get 

# ADCP mooringas at GF10 300 and 75 kHz
f75 = str( datapath /  'raw' / 'ADCP75.mat')
f300 = str(datapath / 'raw' / 'ADCP300.mat')

# MONTHLY CTD'S AND MOORING MICROCAT
file_monthly_CTD = "../data/processed/monthly_18_19_gf10.csv"
file_mooring_microcat = "../data/processed/mooring_gf10.csv"
file_monthly_ctd_all_years_all_stations = "../data/processed/monthly_all_years_all_stations.csv"

file_processed_niaqornaa_weather_station = "../data/processed/selected_niaqornaa_weatherstation.nc"


# BATHYMETRY
file_bathymetry_along_fjord  = f"{datapath}/raw/Bathymetry GF1_GF19.csv"

# TIDAL DATA as determined from the ADCP
file_slack_tides = "../data/processed/positive_slack_tide_moments.csv"


# Diffusion and vertical heat flux
file_vertical_heat_flux = "../data/processed/VerticalHeatFluxfromdiffusionprofiles.csv"