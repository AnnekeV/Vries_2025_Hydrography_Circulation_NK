
import pandas as pd
import gsw
from paths import file_monthly_ctd_all_years_all_stations


df_all = pd.read_csv(file_monthly_ctd_all_years_all_stations) # Load the data via pahts
df_all = df_all.dropna(subset=['Distance']) # Dropping when not in main transect
df_all = df_all.rename(columns={'Pressure [dbar]':'Pressure'})
df_all['Depth'] = gsw.z_from_p(df_all.Pressure.values, df_all.Latitude.values).round(1)  # give accurate depth
df_all['Plot_date'] = df_all['Date'].copy()
df_all['SA'] = gsw.SA_from_SP(df_all['Salinity [PSU]'].values, df_all['Pressure'].values, lon=-51.4, lat=64)
df_all['CT'] = gsw.CT_from_pt(df_all['SA'].values, df_all['Potential temperature [°C]'].values)
df_all['Depth'] =df_all['Pressure']*-1
df_all['Sigma_dens'] = gsw.density.rho(SA=df_all['SA'], CT=df_all['CT'], p=df_all['Pressure'])-1000

