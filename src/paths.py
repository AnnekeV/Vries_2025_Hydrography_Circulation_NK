
from pathlib import Path

path_parent = Path.cwd().parent
figpath = path_parent.joinpath(Path('outputs'))
datapath = path_parent.joinpath(Path('data'))


# Used in the velocity plotting notebook

# ADCP mooringas at GF10 300 and 75 kHz
f75 = '../data/raw/ADCP75.mat'
f300 = '../data/raw/ADCP300.mat'

# MONTHLY CTD'S AND MOORING MICROCAT
file_monthly_CTD = "../data/processed/monthly_18_19_gf10.csv"
file_mooring_microcat = "../data/processed/mooring_gf10.csv"

file_processed_niaqornaa_weather_station = "../data/processed/selected_niaqornaa_weatherstation.nc"

