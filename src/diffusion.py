#%%

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from functions import open_cnv, fpath_GF10,figpath, TSdiagram
# Diffusion for mooring at 540 m
from sklearn.linear_model import LinearRegression
import gsw
pd.options.plotting.backend = "matplotlib"
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from cycler import cycler
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.pyplot import cm
import plotly.graph_objects as go

stylesheet =  'seaborn-ticks',
plt.style.use(stylesheet)
import warnings
warnings.filterwarnings("ignore")
#%%


def to_gsw(df):
    '''Converts cast data into absolute salinity (SA), Conservative temperature (CT), Depth in m relative to sealevel and sigma density (Sigma_dens) [density-1000]
    needs a dataframe with input sal, temp, Pressure [dbar] '''
    df['SA'] = gsw.SA_from_SP(df['sal'].values, df['Pressure [dbar]'].values, lon=-51.4, lat=64)
    df['CT'] = gsw.CT_from_pt(df['SA'].values, df['temp'].values)
    df['Depth'] = gsw.z_from_p(df['Pressure [dbar]'].to_numpy(), lat=64.18)
    df['Sigma_dens'] = gsw.density.rho(SA=df['SA'], CT=df['CT'], p=df['Pressure [dbar]'])-1000
    return df






def make_XYZ(cast_diffusion, variable):
    'returns XYZ based on cast diffuision and specified variable'
    X = cast_diffusion["time"].values.reshape(-1, 1)  # values converts it into a numpy array
    Y = cast_diffusion[variable].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
    Z = cast_diffusion["Pressure [dbar]"].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
    dt_mooring = pd.to_timedelta(cast_diffusion['time'], 'days') + pd.to_datetime('2018-01-01')
    return X,Y,Z, dt_mooring


def lin_reg(X, Y, dt_mooring ,printcoef=False):
    ''' Linear regression for Array X and Y
    Creates plot also'''
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)  # make predictions
    # plt.scatter(dt_mooring, Y, s=5)
    # plt.plot(dt_mooring, Y_pred, color='blue')

    r_sq = linear_regressor.score(X,Y)
    if printcoef:
        print(f"coefficient of determination: {r_sq}")
        print(f"intercept: {linear_regressor.intercept_}")
        print(f"slope: {linear_regressor.coef_}")
    slope_d_day = np.squeeze(linear_regressor.coef_)
    slope_dt = np.squeeze(slope_d_day)/(24*3600)
    if printcoef: print(f"slope: {slope_dt} kg/m3/s")
    return slope_d_day, Y_pred, slope_dt, r_sq

def lin_reg_confidence(x, y, alpha=0.05):
    """
    Perform simple linear regression on x and y and return the slope, intercept, and confidence intervals for
    these parameters.
    
    Parameters:
    x (array-like): array of independent variables
    y (array-like): array of dependent variables
    alpha (float): significance level for confidence intervals
    alpha = 0.05 means 95 % of data is within these values
    
    Returns:
    slope (float): slope of the regression line
    intercept (float): intercept of the regression line
    slope_ci (tuple): tuple of lower and upper bounds for the slope at the specified alpha level
    intercept_ci (tuple): tuple of lower and upper bounds for the intercept at the specified alpha level
    """
    x = sm.add_constant(x) # add constant term for intercept
    model = sm.OLS(y, x).fit() # fit linear regression model
    intercept,slope = model.params # extract slope and intercept
    intercept_ci,slope_ci = model.conf_int(alpha=alpha) # compute confidence intervals
    standard_error = model.bse[1] # standard error of the slope
    return slope, intercept, slope_ci, intercept_ci, model.rsquared, standard_error
one_day = 24*60*60 # s

def plot_confidence_interval(X,Y, ax, dt_mooring, alpha=0.05):
    slope,  intercept, slope_conf, in_conf , Rsq = lin_reg_confidence(X,Y,alpha=alpha)
    slope_conf_sec =slope_conf/one_day
    x = np.squeeze(X)
    ax.plot(dt_mooring, intercept+slope*x)
    ax.fill_between(dt_mooring, in_conf[0]+slope_conf[0]*x, in_conf[1]+slope_conf[1]*x, label = f"{slope_conf_sec[0]:.2e}, {slope_conf_sec[1]:.2e}", alpha=0.5)


#%% ===============================
# DIFFUSION OVER TIME IN MOORING DATA
idx = pd.IndexSlice
plt.style.use('seaborn-ticks')


# IMPORT DATA
down, metadata,cast_540 = open_cnv(f"{fpath_GF10}20190612_SBE5968_540m.cnv" , False)
down, metadata,cast_330 = open_cnv(f"{fpath_GF10}20190612_SBE5969_330m.cnv" , False)

# fig, [ax_dens, ax_temp, ax_sal] = plt.subplots(3, sharex=True)
ylabelfontsize= 'large'
plt.rc('legend',fontsize='small') # using a named size
# plt.rc('title',fontsize='x-large') # using a named size


#DEFINE DATES
# date_start = '2018-07-01'
date_start = '2018-07-01'
date_end_regime1= '2018-10-01'
date_start_regime2 ='2018-12-01'
date_final = '2019-04-01'

df_dt = pd.DataFrame(index=[330,540],columns=['summer', 'winter'],  dtype='float64')
df_dt_con = pd.DataFrame(columns=["Variable", "Depth", "Period", "type", "value", ])

list_rsq, list_index, list_slope = [], [], []

both_casts = {}

linestype_depth = ["-", "--"]
count= 0

colors = ['gainsboro', 'lavender']
unit_drho_dt = r"$kg m^{-3} s^{-1}$"

def dt_confidence_construct_dataframe(df_dt_con, depth, X, Y, slope_dt, period, variable):
    slope,  intercept, slope_conf, in_conf , Rsq ,bse= lin_reg_confidence(X,Y,alpha=0.05)
    slope_conf_sec =slope_conf/one_day
    df_dt_con.loc[len(df_dt_con)]=[variable,depth, period, 'lower', slope_conf_sec[0]]
    df_dt_con.loc[len(df_dt_con)]=[variable, depth,period, 'middle', slope_dt]
    df_dt_con.loc[len(df_dt_con)]=[variable, depth,period, 'upper', slope_conf_sec[1]]
    return Rsq

fig, [ax_dens, ax_temp, ax_sal] = plt.subplots(3, sharex=True, figsize=(12,8))

for i, depth in zip([cast_330, cast_540],[330,540]):

    cast_diffusion = to_gsw(i.copy())
    cast_diffusion['date'] = pd.to_timedelta(cast_diffusion['time'], 'days') + pd.to_datetime('2018-01-01')
    cast_diffusion = cast_diffusion.set_index('date')#['2019-04-01':'2019-06-22'].reset_index()

    both_casts[depth] = cast_diffusion

    # DENSITY
    #all
    variable='density'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion.reset_index(),variable="Sigma_dens")
    slope_d_day, Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    slope_dt = np.squeeze(slope_d_day)/(24*3600)
    print(f"density total slope: {slope_dt} kg/m3/s, R2: {r_sq:.2f}")
    # slope_dt_dens = slope_dt

    if count ==1:
        ax_dens_left =ax_dens
        ax_dens = ax_dens.twinx()
    ax_dens.scatter(dt_mooring, Y, color=colors[count], s=1,zorder=2)
    # ax_dens.plot(dt_mooring, Y_pred, label=f'slope: {slope_dt:.2e} kg/m3/s', color='blue')

    #until nov
    period='summer'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion[date_start:date_end_regime1].reset_index(),variable="Sigma_dens")
    slope_d_day, Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    ax_dens.plot(dt_mooring, Y_pred, linestype_depth[count], label=f'{depth} m' , zorder=3,
                #  color='purple'
                 )
    #  annotate the slope
    slope, intercept, slope_ci, intercept_ci, rsquared, bse = lin_reg_confidence(X,Y,alpha=0.05)
    ax_dens.annotate(f'{slope_dt:.2e} {unit_drho_dt} {bse:.0f}', xy=(dt_mooring[0], Y_pred[0]), xytext=(dt_mooring[0], Y_pred[0]+0.05), fontsize=10, zorder=100, bbox = dict(boxstyle="round", fc=(1, 0.8, 1, 0.5)), arrowprops=dict(arrowstyle="->", ))
    df_dt.loc[depth, 'summer'] = slope_dt
    list_rsq.append(r_sq), list_index.append( [depth, period, variable]), list_slope.append(slope_dt)
    # plot_confidence_interval(X,Y,ax_dens,dt_mooring)
    Rsq = dt_confidence_construct_dataframe( df_dt_con, depth, X, Y, slope_dt, period, variable)
    print(f'{depth}m in {period}: {slope_dt:.2e} kg/m3/s, R2: {r_sq:.2f}')
    print(r_sq,Rsq)


    #from dec to april
    period='winter'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion[date_start_regime2:date_final].reset_index(),variable="Sigma_dens")
    slope_d_day, Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    ax_dens.plot(dt_mooring, Y_pred, linestype_depth[count] ,
                #  color='green'
                zorder=100,
                 )
    ax_dens.annotate(f'{slope_dt:.2e} {unit_drho_dt}', xy=(dt_mooring[0], Y_pred[0]), xytext=(dt_mooring[0], Y_pred[0]+0.05), fontsize=10, zorder=100, bbox = dict(boxstyle="round", fc=(1, 0.8, 1, 0.5)), arrowprops=dict(arrowstyle="->", ))
    df_dt.loc[depth, 'winter'] = slope_dt
    list_rsq.append(r_sq), list_index.append( [depth, period, variable]), list_slope.append(slope_dt)
    Rsq = dt_confidence_construct_dataframe(df_dt_con, depth, X, Y, slope_dt, period, variable)
    print(r_sq,Rsq)


    ax_dens.set_ylabel(
        r"$\sigma$ at " + str(depth) + "m\n" + r"[$kg m^{-3}$]"  ,
        fontsize=ylabelfontsize)
    start_ax_dens_y = Y_pred.mean().round()+0.2-0.1*count
    ax_dens.set_yticks( np.arange(start_ax_dens_y, start_ax_dens_y+0.47, 0.1))
    ax_dens.set_ylim(start_ax_dens_y, start_ax_dens_y+0.4)

    # TEMPERATURE
    variable='temperature'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion.reset_index(),variable="CT")
    slope_d_day,Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    slope_dt = np.squeeze(slope_d_day)/(24*3600)
    print(f"slope: {slope_dt:.2e} kg/m3/s")
    # ax_temp.plot(dt_mooring,Y_pred, color="blue", label=f"slope: {slope_dt:.2e} kg/m3/s")
    ax_temp.scatter(dt_mooring, Y, s=2, color=colors[count])

    #until nov
    period='summer'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion[date_start:date_end_regime1].reset_index(),variable="CT")
    slope_d_day, Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    ax_temp.plot(dt_mooring, Y_pred, linestype_depth[count], label=f'{depth}m, slope: {slope_dt:.2e} °C/s' ,
                 )
    list_rsq.append(r_sq), list_index.append( [depth, period, variable]), list_slope.append(slope_dt)


    #from nov to may
    period='winter'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion[date_start_regime2:date_final].reset_index(),variable="CT")
    slope_d_day, Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    ax_temp.plot(dt_mooring, Y_pred, linestype_depth[count], label=f'{depth}m, slope: {slope_dt:.2e} °C/s,' ,
                 )
    ax_temp.set_ylabel(
        r'$\theta$ [°C]',
          fontsize=ylabelfontsize)
    list_rsq.append(r_sq), list_index.append( [depth, period, variable]), list_slope.append(slope_dt)




    # divide based on salinity
    variable='Relative Salinity'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion.reset_index(),variable="SA")
    slope_d_day,Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    slope_dt = np.squeeze(slope_d_day)/(24*3600)
    print(f"slope: {slope_dt:.2e} kg/m3/s")
    ax_sal.scatter(dt_mooring, Y, s=2, color=colors[count])

    #until nov
    period='summer'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion[date_start:date_end_regime1].reset_index(),variable="SA")
    slope_d_day, Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    ax_sal.plot(dt_mooring, Y_pred, linestype_depth[count],
                color='black'
                )
    Rsq = dt_confidence_construct_dataframe(df_dt_con, depth, X, Y, slope_dt, period, variable)
    list_rsq.append(r_sq), list_index.append( [depth, period, variable]), list_slope.append(slope_dt)

    #from nov to may
    period='winter'
    X,Y,Z,dt_mooring = make_XYZ(cast_diffusion[date_start_regime2:date_final].reset_index(),variable="SA")
    slope_d_day, Y_pred,slope_dt, r_sq = lin_reg(X, Y, dt_mooring)
    ax_sal.plot(dt_mooring, Y_pred, linestype_depth[count], label=f'{depth}m'      ,color='black'    )
    Rsq = dt_confidence_construct_dataframe(df_dt_con, depth, X, Y, slope_dt, period, variable)
    list_rsq.append(r_sq), list_index.append( [depth, period, variable]), list_slope.append(slope_dt)


    plt.xticks(rotation=45, ha='right')
    ax_sal.set_ylabel(
        r'$S_A$ [$g kg^{-1}$]',
                       fontsize=ylabelfontsize)
    count+=1

# increase font size ylabels
ylabelfontsize=17
for ax in [ax_dens, ax_temp, ax_sal, ax_dens_left]:
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(ylabelfontsize)
# font size legend
ax_sal.legend(fontsize=15)
ax_dens.xaxis.set_minor_locator(mdates.WeekdayLocator())
ax_dens.set_xlim([pd.to_datetime('2018-05-25'), pd.to_datetime('2019-06-15')])
# add a b c to the subplots
for letter, ax in zip(["a)", "b)", "c)"], [ ax_dens, ax_temp, ax_sal]):
    ax.text(-0.1, 1, letter, transform=ax.transAxes,
        fontsize=20, fontweight='bold', va='top', ha='right')
# make all lines black
for ax in [ax_dens, ax_temp, ax_sal, ax_dens_left]:
    for line in ax.get_lines():
        line.set_color('black')


plt.tight_layout()
plt.subplots_adjust(hspace=0.07)

fig.savefig(f"{figpath}/Vertical_diffusion_from_microcat_two_regimes_both_all.jpg", dpi=300)
plt.show()


#apply format_float_scientific to all values in list_slope
list_slope = [np.format_float_scientific(x, precision=3) for x in list_slope]
mi_list_index = pd.MultiIndex.from_tuples(list_index, names = ['depth', 'period', 'variable'])
df_rsq_dt  = pd.DataFrame(index=mi_list_index, data= list_slope, columns=['slope [var/s]'])
df_rsq_dt['R2'] = list_rsq
df_rsq_dt.R2 = df_rsq_dt.R2.round(2)
# reshuffle multi index
df_rsq_dt = df_rsq_dt.reorder_levels(['variable', 'depth', 'period']).sort_index()
# df_rsq_dt.to_csv("intermediate/slope_mooring_r2_dataframe.csv")
# df_rsq_dt.to_latex("intermediate/slope_mooring_r2_dataframe.tex", bold_rows=True, caption="Slope from Mooring GF10", label="tab:slope_mooring_r2_dataframe")

# %%

def lin_reg2(X, Y, label=None):
    ''' Linear regression for Array X and Y
    Creates plot also'''
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X, Y)  # perform linear regression
    Y_pred = linear_regressor.predict(X)  # make predictions
    # plt.plot(X, Y_pred, color='blue')

    r_sq = linear_regressor.score(X,Y)
    # print(f"coefficient of determination: {r_sq}")
    # print(f"intercept: {linear_regressor.intercept_}")
    # print(f"slope: {linear_regressor.coef_}")
    return np.squeeze(linear_regressor.coef_), r_sq, linear_regressor



# %% d rho /d z  calc
df_sorted = pd.read_csv("intermediate/df_sorted.csv")
df_sorted['datetime'] =pd.to_datetime(df_sorted.date) + pd.to_timedelta(df_sorted.time.mod(1), 'days')
df_sorted = to_gsw(df_sorted)


period = "winter"

full_parameter = {'dens':'Density', 'CT':'Temperature', 'SA':'Relative Salinity'}
dates_summer = ['2018-06-06','2018-06-27','2018-07-18', '2018-08-08', '2018-09-06']
dates_winter = ['2018-12-11', '2019-02-05', '2019-03-05']
interval=20

df_dz_comp = pd.DataFrame(columns=['target_depth', 'period', 'min_depth', 'max_depth', 'interval', 'average','min', 'max', 'Rsq'])
df_dz_comp2 = pd.DataFrame(columns=["Depth", "Period", "type","Variable", "value"])
for interval in [40]:
    for relative_density_level in [330, 540]: # Pressure at which to check density for Delta dnesity
        shallow_boundary = relative_density_level-interval/2  #m
        deep_boundary =relative_density_level+interval/2  # m
        for period in [ "summer", "winter"]:
            if period == "summer": 
                date_select = dates_summer # choose period
            elif period == "winter": 
                date_select = dates_winter # choose period
            elif period == "all":
                date_select =  df_sorted.date.unique()
            else: 
                print("no period chosen!!!")
            # df_sorted.date.unique()


            var = 'dens'
            slope_dz = []
            date = []

            df_sorted.loc[(df_sorted.date == '2018-08-09')&(df_sorted["Pressure [dbar]"]==568)] = np.nan
            fig = plt.figure(figsize=[4,6])

            bigX = np.array([])
            bigY = np.array([])
            for i in df_sorted.date.unique() :
                monthly = df_sorted[df_sorted.date == i]
                rel_dens = monthly.loc[monthly["Pressure [dbar]"] ==relative_density_level, var].to_numpy() 
                if len(rel_dens)>0 :
                    monthly['delta_var'] = monthly[var] - rel_dens
                    diff_monthly = monthly.set_index('Pressure [dbar]')[shallow_boundary:deep_boundary].reset_index()
                # if not diff_monthly[var].isnull().values.all():
                    diff_monthly.date = pd.to_datetime(df_sorted['date'])
                    X = diff_monthly["Pressure [dbar]"].values.reshape(-1, 1)  # values converts it into a numpy array
                    Y = diff_monthly["delta_var"].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
                    
                    if i in date_select:
                        print(i, interval, period)
                        slope,r_sq, linregmodel = lin_reg2(X,Y, label=i)
                        slope_dz.append(slope)
                        plt.plot(monthly["delta_var"], monthly["Pressure [dbar]"], "--",label=f"{i}", zorder=3)
                        plt.scatter(diff_monthly["delta_var"], diff_monthly["Pressure [dbar]"], label=f"{slope:.2e}, R2: {r_sq*100:.1f}", zorder=3, s=1)
                        date.append(pd.to_datetime(i))
                        bigX = np.concatenate((bigX, np.squeeze(X)))
                        bigY = np.concatenate((bigY, np.squeeze(Y)))

                    else:
                        # plt.plot(monthly["delta_var"], monthly["Pressure [dbar]"], color='silver' , label = "", zorder=2)
                        print("")
                    # print(f"slope: {slope:.3e} kg/m3/m")
                
            plt.ylabel("Depth [m]")
            plt.xlabel(r"$\Delta$Density [$kg/m^3$ ]")
            plt.ylim([deep_boundary,shallow_boundary])
            # plt.xlim([-.5,0.5])
            plt.title(f"{full_parameter[var]} relative to {relative_density_level:.0f}m, in {period}")
            slope, intercept, slope_ci, intercept_ci , Rsq = lin_reg_confidence(bigX,bigY)
            plt.legend(loc=3)
            df_dz_comp.loc[len(df_dz_comp)] = [relative_density_level,period, shallow_boundary, deep_boundary, interval, np.mean(slope),np.min(slope_ci), np.max(slope_ci), Rsq]
            df_dz_comp2.loc[len(df_dz_comp2)] = [relative_density_level, period, 'middle', full_parameter[var], slope]
            df_dz_comp2.loc[len(df_dz_comp2)] = [relative_density_level, period, 'upper', full_parameter[var], slope_ci[0] ]
            df_dz_comp2.loc[len(df_dz_comp2)] = [relative_density_level, period, 'lower', full_parameter[var], slope_ci[1] ]

            # Plot linear fit per interval
            # fig = plt.figure(figsize=[4,6])
            plt.scatter(bigY,bigX, s=10, label = "for big x")
            plt.plot(np.sort(bigX)*slope + intercept , np.sort(bigX),color='red', label = f"{slope:.2e} +- {(np.max(slope_ci)-slope):.2e}, R2: {Rsq*100:.1f}")
            plt.ylim([deep_boundary, shallow_boundary]  )
            plt.legend()
            # plt.title(f"drho/dz between {shallow_boundary} and {deep_boundary}")
            plt.xlim([min(bigY)*2, max(bigY)*2])
            plt.show()



        
plt.show()

# %% All in one
# All in one
interval = 40

fig, (all_axes) = plt.subplots(2,2, sharey='row', sharex='col', figsize=(7,8))
ax_ori, ax = all_axes[:,0], all_axes[:,1]

# make an empty dictionary
dict_drhodz = dict()
dict_drhodz_BE = dict()

selected_dates = df_sorted.date.unique()[5:-7]
selected_dates = dates_summer
addition = "_summer"

color = iter(cm.plasma(np.linspace(0, 1, len(selected_dates))))
color = dict(zip(selected_dates, color))


for j, relative_density_level in enumerate([330, 530]): # Pressure at which to check density for Delta dnesity
    shallow_boundary = relative_density_level-interval/2  #m
    deep_boundary =relative_density_level+interval/2  # m

    dict_drhodz_date = dict()

    for i in selected_dates :
        SE = []
        monthly = df_sorted[df_sorted.date == i]
        # mean_dens_at_relative_density_level 
        mask_relevant= monthly["Pressure [dbar]"].isin(np.arange(shallow_boundary,deep_boundary))
        if mask_relevant.sum() >0:
            rel_dens = monthly.loc[mask_relevant].sort_values(["Pressure [dbar]"]).iloc[0]['dens']

            # calculate slope
            X = monthly["Pressure [dbar]"].loc[mask_relevant].values.reshape(-1, 1)  # values converts it into a numpy array
            Y = monthly["dens"].loc[mask_relevant].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
                
            slope,r_sq,linregmodel = lin_reg2(X,Y, label=i)
            # get standard error with lin_reg_confidence
            slope, intercept, slope_ci, intercept_ci, rsquared, bse = lin_reg_confidence(X,Y)
            SE.append(bse)            
            signif = ""
            if r_sq>0.99:
                signif = "**"
            elif r_sq>0.95:
                signif = "*"
            ax[j].plot(monthly["dens"].loc[mask_relevant]-rel_dens, monthly["Pressure [dbar]"].loc[mask_relevant] , label = f"{i}   {slope*1000:.2f}{signif} {bse:.0e}",  c=color[i])
            
            ax_ori[j].plot(monthly["dens"].loc[mask_relevant], monthly["Pressure [dbar]"].loc[mask_relevant] ,  c=color[i])
            dict_drhodz_date[i] = slope
    # add slope to dictionary
    dict_drhodz[relative_density_level] = dict_drhodz_date
    # add standard error to dictionary


    ax[j].set_ylim(deep_boundary, shallow_boundary)
    ax[j].legend(title=r"Date  Slope $\times 10^{-3} \ kg \ m^{-4}$     SE", loc='upper left', bbox_to_anchor=(1, 1), fontsize=10)
#  find mean slope of per depth
mean_slope = []
for i in dict_drhodz.keys():
    mean_slope.append(pd.DataFrame(dict_drhodz[i], index=[i]).mean(axis=1))
    print(f"{(mean_slope[-1])}")
#  add in ax as text
for i, axx in enumerate(ax):
    axx.text(0.95, 0.99, "Mean = \n" + f"{(mean_slope[i]*1000).values[0]:.2f}"+r"$\times 10^{-3} \ kg \ m^{-4}$", transform=axx.transAxes, va='top', ha='right')
# set legend on the right outside ax
# set legend title 

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
# set minor x ticks every 0.1
for i in range(2):
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.05))
    ax_ori[i].xaxis.set_minor_locator(MultipleLocator(0.1))

# also add x ticks on top
for i in range(2):
    ax[i].xaxis.set_ticks_position('both')
    ax_ori[i].xaxis.set_ticks_position('both')
    ax[i].tick_params(axis='x', which='both', labelbottom=True)
    ax_ori[i].tick_params(axis='x', which='both', labelbottom=True)
ax_ori[j].set_ylabel("Depth [m]")
plt.tight_layout()
# reduce hspace
plt.subplots_adjust(hspace=0.1)
# add a b c d to the subplots
for axx, letter in zip(fig.axes, ["a)", "b)", "c)", "d)"]): 
    axx.text(0.02, 0.01, letter, transform=axx.transAxes,
        fontsize=15, fontweight='bold', va='bottom', ha='left')

fig.savefig(f"{figpath}/Diffusion/Vertical_diffusion_drhodz_interval_{interval}m{addition}.jpg", dpi=300)


# %% d rho /d z  FIGURE


period = "winter"
shallow_boundary = 500  #m
deep_boundary =520 # m
relative_density_level = np.mean([shallow_boundary, deep_boundary]) # Pressure at which to check density for Delta dnesity

dates_summer = ['2018-06-06','2018-06-27','2018-07-18', '2018-08-08', '2018-09-06']
dates_summer = ['2018-07-18', '2018-08-08', '2018-09-06']
dates_winter = ['2018-12-11', '2019-02-05', '2019-03-05']


fig2, axes = plt.subplots(1,2,sharex=True)

for k, period in enumerate(["summer", "winter"]):
    if period == "summer": 
        date_select = dates_summer # choose period
    elif period == "winter": 
        date_select = dates_winter # choose period
    elif period == "all":
        date_select =  df_sorted.date.unique()
    else: 
        print("no period chosen!!!")
    # df_sorted.date.unique()


    var = 'sigma-dens'
    slope_dz = []
    date = []

    df_sorted.loc[(df_sorted.date == '2018-08-09')&(df_sorted["Pressure [dbar]"]==568)] = np.nan
    fig, ax= plt.subplots(1,1,figsize=[4,6])
    axs = axes[k]


    for i in df_sorted.date.unique() :
        monthly = df_sorted[df_sorted.date == i]
        rel_dens = monthly.loc[monthly["Pressure [dbar]"] ==relative_density_level,var].to_numpy() 
        if len(rel_dens)>0 :
            monthly['delta_var'] = monthly[var] - rel_dens
            diff_monthly = monthly.set_index('Pressure [dbar]')[shallow_boundary:deep_boundary].reset_index()
        # if not diff_monthly[var].isnull().values.all():
            diff_monthly.date = pd.to_datetime(df_sorted['date'])
            X = diff_monthly["Pressure [dbar]"].values.reshape(-1, 1)  # values converts it into a numpy array
            Y = diff_monthly["delta_var"].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column

            if i in date_select:
                slope,r_sq = lin_reg2(X,Y, label=i)

                ax.plot(monthly["delta_var"], monthly["Pressure [dbar]"], "--",label=f"{i}", zorder=3)
                ax.scatter(diff_monthly["delta_var"], diff_monthly["Pressure [dbar]"], label=f"{slope:.1e}, R2: {r_sq:.5f}", zorder=3)
                date.append(pd.to_datetime(i))

                # make a mask where pressure is between shallow_boundary and deep_boundary
                mask = (monthly["Pressure [dbar]"] >= shallow_boundary) & (monthly["Pressure [dbar]"] <= deep_boundary)
                mask = (monthly["Pressure [dbar]"] >= shallow_boundary) & (monthly["Pressure [dbar]"] <= deep_boundary)
                coef1, r1 = lin_reg2(monthly[mask]['Pressure [dbar]'].values.reshape(-1,1),monthly[mask]['sigma-dens'].values.reshape(-1,1) ,label=i)
                axs.plot(monthly[mask]['sigma-dens'], monthly[mask]['Pressure [dbar]'], label=f"{i}, {coef1:.1e}, R2: {r_sq:.1f}", zorder=3)
                axs.set_title(f"{period}")
                                
                slope_dz.append(coef1)
                print(slope, coef1) # the same!


            else:
                ax.plot(monthly["delta_var"], monthly["Pressure [dbar]"], color='silver' , label = "", zorder=2)
                print()
            # print(f"slope: {slope:.3e} kg/m3/m")
        
    ax.set_ylabel("Depth [m]")
    ax.set_xlabel(r"$\Delta$Density [$kg/m^3$ ]")
    ax.set_ylim([shallow_boundary-20,deep_boundary+20])
    ax.set_xlim([-0.1,0.1])

    ax.set_title(f"Density relative to {relative_density_level:.0f}m, in {period}")
    fig.legend(loc=3)
    # fig.savefig(f"{figpath}/Vertical_diffusion_from_CTD_{relative_density_level:.0f}_{period}.jpg", dpi=600 )
    np.mean(slope_dz)
    
    axs.set_ylim([deep_boundary+5, shallow_boundary-5])
    # add slope_dz in axs in scientific notation
    axs.text(0.05, 0.05, f"Mean slope: {np.mean(slope_dz):.1e} kg/m3/m", transform=axs.transAxes)
    axs.legend(loc='upper right')

    # plt.show()
plt.show()

plt.plot( date, slope_dz)
plt.ylabel("Coefficient [kg/m3/m]")
plt.title("Coefficient density change with depth drho/dz 530-550")
plt.show()

middle, _, slope_dz_con, _ , Rsq ,_= lin_reg_confidence(X,Y)
df_dz = pd.DataFrame.from_dict({'type':['middle','lower','upper'], 'value':[middle, slope_dz_con[0], slope_dz_con[1]]})

# %% CALCULATE KAPPA
slope_dt_dens = df_dt.loc[540, 'winter']

kappa = slope_dt_dens * (deep_boundary-shallow_boundary) / np.mean(slope_dz)
print(f"kappa =  {kappa:.4e}")

print(f"drho/dt = {slope_dt_dens:.4e}")
print(f"drho/dz = {np.mean(slope_dz):.4e}")


df_kappa = df_dt_con[df_dt_con.Variable=='Relative Salinity'].set_index(['Depth','Period', 'type', 'Variable'])*(20) / df_dz_comp2.set_index(['Depth','Period', 'type', 'Variable'])
df_kappa.head(20)

df_kappa_uncertainty = (df_kappa.loc[:,:,"lower"]- df_kappa.loc[:,:,"upper"])/2/1e-5

# format df_dz in scientific notation
df_dz_comp2.value = df_dz_comp2.value.apply(lambda x: np.format_float_scientific(x, precision=3))

# %% DOUBLE CHECK WITH T AND S
# =========================




# %% PROFILE KAPPA
# =========================

def make_same_length(profile1, profile2, depth1, depth2, date1, date2):
    if np.max(abs(depth1))>np.max(abs(depth2)):
        max_d2 = np.squeeze(np.where(abs(depth1) == abs(depth2).max()))+1
        profile1= profile1[:max_d2]
        depth1=depth1[:max_d2]
        date1=date1[:max_d2]
    else:
        max_d1 = np.squeeze(np.where(abs(depth2) == abs(depth1).max()))+1
        profile2= profile2[:max_d1]
        depth2=depth2[:max_d1]
        date2=date2[:max_d1]    
    return profile1, profile2, depth1, depth2, date1, date2

def roll_then_resample(profile, ax, step, window):
    ''' takes a rolling mean of window size'''
    profile = pd.Series(index=ax, data=profile)
    df = pd.DataFrame(index=np.arange(min(ax), max(ax)+1))
    df['profile'] = profile
    df = df.rolling(window=window, center=True, min_periods=int(np.round(window/2))).mean().iloc[::step]
    return(df.index.to_numpy(), df.profile.to_numpy())

def calculate_kappa_z(A, drhodz, u,drhodt,dz):
    ''''
    diffusivity coefficient for a certain depth based on Stigebrandt 1989
    A is horizonta area 
    drhodt is an array of values of drho/dt, in dataframe, from shalllow to deep
    dz is an array of dz or a value with the z interval
    '''
    if not (isinstance(drhodt, pd.DataFrame) or isinstance(drhodt, pd.Series)):
        print('your variables are not a data frame')
        return 
    if not (isinstance(A, pd.DataFrame) or isinstance(A,pd.Series)):
        A=pd.Series(index=drhodt.index, data=A)
    sum_drhodt = (drhodt.loc[u:]*A.loc[u:]*dz).sum()
    kappa_z = 1/(A.loc[u]*drhodz.loc[u])*sum_drhodt
    return kappa_z

def give_k_profile(profile1, profile2, depth1, depth2, date1, date2, A, window=1,step=1):
    '''
    Calculates diffusive coefficient profile along depth, according to the budget method
    profile 1 and profile 2 should be two variable profiles (from shallow to deep)
    date1 and date 2 are two dates in pd datetime format
    depth1 and depth2 are depth profiles (positive)
    '''
    if (depth1 !=depth2):
        profile1, profile2, depth1, depth2, date1, date2 = make_same_length(profile1, profile2, depth1, depth2, date1, date2)
        print("Depth profiles are not equal")

    if window>1:
        depth1, profile1 = roll_then_resample(profile1, depth1, window=window, step=step)
        depth2, profile2 = roll_then_resample(profile2, depth2, window=window, step=step)
        date1, date2 = date1[::step], date2[::step]

    dz = depth1[0]-depth1[1]
    depth = depth1[1:-1]
    drhodt = (profile2-profile1) / (date2-date1).dt.total_seconds().values
    drhodz1= (profile1[2:]- profile1[:-2])/(depth1[2:]-depth1[:-2]) # centre difference
    drhodz2= (profile2[2:]- profile2[:-2])/(depth2[2:]-depth2[:-2]) # centre difference

    drhodz = (drhodz1+drhodz2)/2   # not sure if this is correct
    kappa  = pd.Series(index=depth)
    counter =0
    for u in depth:
        counter+=1
        kappa.loc[u] =  calculate_kappa_z(A, drhodz=pd.Series(index=depth,data=drhodz), drhodt=pd.Series(index=depth,data=drhodt[1:-1]), u=u, dz=dz)


    return kappa, depth1[1:-1], (date1+(date2-date1)/2)[0], drhodt[1:-1], drhodz

# %%

# set stylesheet
plt.style.use('seaborn-whitegrid')
from import_crosssection import depth_int, width_int

initial_date ='2018-08-09'
final_date = '2018-10-15'
upper_limit = 250
df_sorted['Sigma_dens_smooth'] = df_sorted.Sigma_dens.copy()

for i in df_sorted.date.unique():
    df_sorted.loc[df_sorted.date == i , 'Sigma_dens_smooth'] = df_sorted.loc[df_sorted.date == i,'Sigma_dens'].rolling(center=True, window=40,).mean()

df1 = df_sorted[df_sorted.date == initial_date ].set_index(['Pressure [dbar]'])[upper_limit:].reset_index()
df2 = df_sorted[df_sorted.date == final_date].set_index(['Pressure [dbar]'])[upper_limit:].reset_index()


fig, axs = plt.subplots(1,5,sharex=False, sharey=True)
axs[0].plot(df1.Sigma_dens_smooth, df1['Pressure [dbar]'], label="Initial profile")
axs[0].plot(df2.Sigma_dens_smooth, df2['Pressure [dbar]'], label="Final profile")
axs[0].plot(df2.Sigma_dens, df2['Pressure [dbar]'], label="Final profile")

axs[0].set_ylim([565,250])

A_area = pd.read_csv("/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/Data/Other/Bathymetry/Area/Area_basin_gf10.csv")
A_area.index.name = 'Pressure [dbar]'


# Construct Area data_frame
A_int = pd.Series(index=depth_int*-1, data=width_int)
A = pd.DataFrame(index=df1["Pressure [dbar]"], data=A_int)[0]
A[A.isnull()] = A.max()
A = A*np.linspace(125,25, len(A.index))

kappa, depth, date, drhodt, drhodz  = give_k_profile(df1.Sigma_dens_smooth.values, df2.Sigma_dens_smooth.values, df1['Pressure [dbar]'].values, df2['Pressure [dbar]'].values, df1.datetime, df2.datetime, window=1, A=A_area['Area_km2'])

axs[1].plot(drhodt, depth)
axs[2].plot(drhodz, depth)
axs[3].plot(kappa, depth)
axs[3].ticklabel_format(axis='x', style='sci', scilimits=(1,2))
axs[4].plot(A_area['Area_km2'], A_area.index)
ax_titles = ["Density", "drho/dt","drho/dz",  r"$\kappa_z$"  ]
axs[4].set_xscale('log')
for i, title in zip(range(4),ax_titles): axs[i].set_xlabel(title)
fig.suptitle(f"{initial_date} -  {final_date}, {date.strftime('%Y-%m-%d')}")


# %%
# use template normal
# set stylesheet
plt.style.use('default')

dates_summer = ['2018-06-27', '2018-08-08','2018-09-06' ,'2018-10-15']
dates_summer_in = dates_summer[:-1]
dates_summer_final = dates_summer[1:]
dates_winter = ['2018-12-11', '2019-02-05','2019-03-05', '2019-04-10',]
dates_winter_in = ['2018-12-11', '2019-02-05']
dates_winter_final = ['2019-02-05', '2019-04-10',]
legend_labels_long = []
legend_labels = ["S1", "S2", "S3", "W1", "W2"]

dfHeatFlux = pd.DataFrame(index=np.arange(600), columns=legend_labels)

fig, axs = plt.subplots(1,5,sharex=False, sharey=True, figsize=[8,5])
z_coordinate = 'Depth'
upper_limit_depth =200
scilim_drhodt = -8
scilim_drhodz = -3
c =0
for period in ['summer', 'winter']:
    if period == 'winter':
        dates_in = dates_winter_in
        dates_final = dates_winter_final
    else:
        dates_in = dates_summer_in
        dates_final = dates_summer_final


    def reindex_and_interpolate(df, new_index):
        return df.reindex(df.index | new_index).interpolate(method='index', limit_direction='both').loc[new_index]

    from scipy.stats import linregress
    def heat_flux(df, kappa, interval=10):
        '''Calclulates heat flux
        df needs to have Depth, CT, and dens
        kappa needs to have almost the same index as df'''
        spec_heat = 3850
        df_heat_flux =  pd.DataFrame(index=df1.index, columns=['Heat flux'])
        for i in range(int(np.floor(len(df)/interval))):
            X = df[i*interval:i*interval+interval].Depth.to_numpy()
            Y = df[i*interval:i*interval+interval].CT.to_numpy()
            slope, intercept, r_value, p_value, std_err = linregress(X, Y)
            dTdz = slope
            dens = df[i*interval:i*interval+interval].dens.mean()
            kappa10m = kappa.reset_index(drop=True)[i*interval:i*interval+interval].mean()
            df_heat_flux.iloc[i*interval+int(interval/2)] = dens*spec_heat*dTdz*kappa10m # from Bendtsen2021
        return df_heat_flux


    # choose two different colors
    # make a list  of 5 discrete colors from 
    import cmcrameri as ccm
    colors = [ccm.cm.batlow_r(i) for i in np.linspace(0, 1, 6)]

    for i in range(len(dates_in)):
        label=  f"{pd.to_datetime(dates_in[i]).strftime('%B %d, %Y')} - {pd.to_datetime(dates_final[i]).strftime('%B %d, %Y')}"
        legend_labels_long.append(label)
    
        df1 = df_sorted[df_sorted.date == dates_in[i] ].set_index([z_coordinate]).reset_index()
        df2 = df_sorted[df_sorted.date == dates_final[i]].set_index([z_coordinate]).reset_index()
        area_depth_level = reindex_and_interpolate(A_area.set_index('Depth'), df1.Depth)['Area_km2']

        kappa, depth, date, drhodt, drhodz  = give_k_profile(df1.Sigma_dens_smooth.values, df2.Sigma_dens_smooth.values, df1[z_coordinate].values, df2[z_coordinate].values, df1.datetime, df2.datetime, window=1, A=area_depth_level )

        upper_limit = int(np.where(np.round(abs(depth) )== upper_limit_depth)[0] )
        axs[0].plot(df1.Sigma_dens.iloc[upper_limit:], df1[z_coordinate].iloc[upper_limit:], ":", label="Initial profile", color=colors[c])
        axs[0].plot(df2.Sigma_dens.iloc[upper_limit:], df2[z_coordinate].iloc[upper_limit:], label="Final profile", color=colors[c])
        axs[1].plot(drhodt[upper_limit:]/10**(scilim_drhodt), depth[upper_limit:], color=colors[c])
        if not ((len(kappa)<(250-100)) or  (np.max(abs(kappa.index))<250) ):
            axs[2].plot(drhodz[upper_limit:]/10**(scilim_drhodz), depth[upper_limit:], color=colors[c])
            axs[3].plot(kappa.iloc[upper_limit:], depth[upper_limit:],label = label, color=colors[c])

        # heat flux
        df_heat_flux = heat_flux(df1, kappa, interval=10)
        df1['Heat flux'] = df_heat_flux
        dfHeatFlux[legend_labels[c]] = df_heat_flux['Heat flux']
        axs[4].plot(df_heat_flux[(~df_heat_flux.isnull()).to_numpy()], df1.loc[(~df_heat_flux.isnull()).to_numpy(),z_coordinate],  color=colors[c], label=label)
        axs[4].set_xlim([-15, 20])

        c+=1

dfHeatFlux = dfHeatFlux.dropna(how='all', axis=0)
dfHeatFlux.to_csv(f"intermediate/VerticalHeatFluxfromdiffusionprofiles.csv")

for i, letter in enumerate(['a)', 'b)', 'c)', 'd)', 'e)']): axs[i].text(0.05, 0.03, letter, transform=axs[i].transAxes,weight='bold')

axs[1].set_ylim([-575,-upper_limit_depth])
axs[1].set_yticklabels([str(abs(int(x))) for x in axs[1].get_yticks()])

axs[3].legend(loc='upper right', bbox_to_anchor=(1, -.1),
            ncol=5, frameon=False)
# save legend labels for later
legend_labels = ["S1", "S2", "S3", "W1", "W2"]
for i, label in enumerate(legend_labels): axs[3].get_legend().get_texts()[i].set_text(label)


for i in range(1,4): 
    axs[i].ticklabel_format(axis='x', style='sci', scilimits=(1,2))
# set xlim for plots
axs[0].set_xlim([27.5, 29.5])
axs[1].set_xlim([-2.5e-8/10**(scilim_drhodt), 3e-9/10**(scilim_drhodt)])
# add text to axs[1], add 10**-8 to the plot in the bottom right in MathText
axs[1].text(0.65, 0.03, rf'${{\times}}10^{{{scilim_drhodt}}}$', transform=axs[1].transAxes)
axs[2].set_xlim([-6.5e-3/10**(scilim_drhodz), -4.6e-3/10**(scilim_drhodz), ])
axs[2].text(0.65, 0.03,rf'${{\times}}10^{{{scilim_drhodz}}}$', transform=axs[2].transAxes)
axs[3].set_xlim([1e-6, 1e-3])
axs[3].set_xscale('log',)
ax_titles = [r"$\sigma$ "+ "\n"+r" [kg m$^{-3}$] ", r"$d\rho$/dt"+ "\n"+r" [kg m$^{-3}$s$^{-1}$]",r"$d\rho$/dz"+ "\n"+r" [kg m$^{-3}$m$^{-1}$]", r"$\kappa_z$"+ "\n"+" [m$^{2}$s$^{-1}]$" , "Q "+ "\n"+r"[Wm$^{-2}$]"]
for i, title in zip(range(len(axs)),ax_titles): axs[i].set_xlabel(title)
axs[0].set_ylabel("Depth [m]")

# format axes
nr_axes = 5
for i in range(nr_axes): axs[i].grid(False)
for i in range(nr_axes): axs[i].spines['right'].set_visible(False)
for i in range(nr_axes): axs[i].spines['top'].set_visible(False)
for i in range(nr_axes): axs[i].xaxis.set_label_position('top')
for i in range(nr_axes): axs[i].tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)

for i in range(3):
    mf = mpl.ticker.ScalarFormatter(useMathText=True)
    mf.set_powerlimits((-8,2))
    axs[i].xaxis.set_major_formatter(mf)
# plt.tight_layout()
# plt.savefig(f"{figpath}/Diffusion/winter_and_summer_profiles_40m_rolling.png", dpi=600, bbox_inches='tight')
# plt.savefig(f"{figpath}/Diffusion/winter_and_summer_profiles_40m_rolling_draft.jpg", dpi=150, bbox_inches='tight')
# plt.savefig(f"{figpath}/Diffusion/winter_and_summer_profiles_40m_rolling.jpg", dpi=300, bbox_inches='tight')




# %%   =====================
# #    TS PLOT


# mpl.rcParams['axes.prop_cycle'] = cycler(linestyle=['-','--','-.',':'])
mpl.rcParams['axes.prop_cycle'] = cycler(linestyle=['-'])


color = iter(cm.rainbow(np.linspace(0, 1, 16)))

all_dates_mooring = ['2018-06-06', '2018-06-27', '2018-07-18', '2018-08-08',
       '2018-08-09', '2018-09-06', '2018-10-15', '2018-10-31',
       '2018-12-11', '2019-02-05', '2019-03-05', '2019-04-10', 
      '2019-05-06' ]
      #'2019-06-04','2019-06-14', '2019-07-01',]

fig, ax = plt.subplots()
TSdiagram(ax, Svals=np.linspace(33.4,33.8,10), Tvals = np.linspace(0.5,3.0,10),  levels= 6)


for i in range(len(all_dates_mooring)):
    c = next(color)
    dfTS = df_sorted[df_sorted.date == all_dates_mooring[i] ].set_index([z_coordinate])[-300:]
    if len(dfTS) >100:
        ax.plot(dfTS.SA, dfTS.CT, label=all_dates_mooring[i], c=c, zorder=4)
        print(all_dates_mooring[i], dfTS.index[-1])



# ax.scatter(both_casts[540].rolling(window=25*6, center=True).mean().SA, both_casts[540].rolling(window=50*6, center=True).mean().CT, s=1,c=both_casts[540].time, cmap=cm.magma)

ax.scatter(
    both_casts[540][:pd.to_datetime('2019-05-06')].SA[0], 
    both_casts[540][:pd.to_datetime('2019-05-06')].CT[0],
    marker="*", color = 'black', s=50, zorder=5, label ='Begin Mooring')
ax.scatter(
    both_casts[540][:pd.to_datetime('2019-04-10')].rolling(window=50*6, center=True).mean().SA, 
    both_casts[540][:pd.to_datetime('2019-04-10')].rolling(window=50*6, center=True).mean().CT, 
    s=1, label="Until April 10", color='grey', zorder=3)
ax.scatter(
    both_casts[540][:pd.to_datetime('2019-05-06')].rolling(window=50*6, center=True).mean().SA, 
    both_casts[540][:pd.to_datetime('2019-05-06')].rolling(window=50*6, center=True).mean().CT, 
    s=1, label="April 10 - May 6", zorder=2)

ax.legend(loc='upper left', bbox_to_anchor=(1, 1.01),
            ncol=1, prop={'size': 10})
ax.set_xlabel('Absolute salinity [g/kg]')
ax.set_xlabel('Absolute salinity [g/kg]')
ax.set_ylabel('Conservative Temperature [C]')
plt.tight_layout()
# fig.savefig(f"{figpath}/Diffusion/diffusion_TS_540m_until_May.png", dpi=1200)

# plt.plot(both_casts[540].rolling(window=1).mean().CT, color='red')
# plt.plot(both_casts[540].rolling(window=50*6, center=True).mean().CT)
    


# %% ============
# TS PLOT
from plot_transect import df_spring,df_all
from plotly.subplots import make_subplots
stylesheet =  'seaborn-ticks',
plt.style.use(stylesheet)
# # plt.rcParams.update({'font.size': 15})
# # make legend, ylabel, yticks fontsize larger
# fsize = 20
# # make legend fontsize larger
# plt.rcParams['ytick.labelsize'] = fsize
# plt.rcParams['xtick.labelsize'] = fsize
# plt.rcParams['axes.labelsize'] = fsize

# make plotly figure 
plotly_figure = make_subplots(rows=1, cols=2)
# add line to plotly figure
tempunit = 'Potential temperature [°C]'
salunit = 'Salinity [PSU]'

sallim = [33.15,33.7]
templim = [1,4]
templim = [-.5, 2.5]

fig, (ax, ax10) = plt.subplots(1, 2, figsize=(6,3), sharex=True, sharey=True)
TSdiagram(ax, Svals=np.linspace(sallim[0], sallim[1],10), Tvals = np.linspace(templim[0], templim[1],11),  levels= 6, alpha=0.2)
TSdiagram(ax10, Svals=np.linspace(sallim[0], sallim[1],10), Tvals = np.linspace(templim[0], templim[1],11),  levels= 6, alpha=0.2)
ax.set_xlim(sallim)
ax.set_ylim(templim)

df_spring_HS = ['HS190123', 'HS190213', 'HS190328', 'HS190423','HS190514']

df_all = df_all.sort_values(["timeJ", "Distance",  "Pressure"])
df_spring = [ 'HS190123','HS190213', 'GF19010', 'GF19011', 'GF19012', 'GF19013', 'GF19017', 'HS190328',
       'GF19024', 'GF19026', 'GF19027', 'GF19029',   'HS190423','GF19033',
       'GF19034', 'GF19035', 'GF19036', 'HS190514']
df_spring_select =[ 'HS190213', 'GF19010', 'GF19011',  'GF19017', 'HS190328',
       'GF19029',   'HS190423','GF19033',
       'GF19034', 'GF19035', 'GF19036', 'HS190514']
df_all['DateTIME'] = df_all.Date.astype('datetime64[D]').copy()
df_other_names = df_all.set_index('DateTIME')['2019-01-01':'2019-07-01'].Name.unique()

# import crameri colormaps
import cmcrameri as ccm


color = iter(cm.rainbow(np.linspace(0, 1, len(df_spring))))

linestyle_st = dict(zip([np.nan, 'GF5', 'GF6', 'GF8', 'GF10', 'Other', 'GF3'], ['--', '-', '-.', ':', '-', '-', '-.']))
df_all['Plot_st'] =df_all['St.']
df_all.loc[df_all['Name'] == "GF19010", 'Plot_st'] = "Close to GF3"

dates = df_all[(df_all.Name.isin(df_spring)&(df_all.Pressure >150))].Date.unique()
cmap = plt.cm.jet  # define the colormap
cmaplist = plt.cm.rainbow(np.linspace(0, 1, len(dates)))
color_date = dict(zip(dates,cmaplist))

for i in df_spring:
# for i in df_other_names:
    df = df_all[(df_all.Name == i)&(df_all.Pressure >150)]
    if len(df) == 0: continue
    Station = df['St.'].iloc[0]
    Date = df.Date.iloc[0]
    # label = f"{df.Date.iloc[0]}, {df['Plot_st'].iloc[0]}"
    label = f"{df.Date.iloc[0]}"
    if Station == 'GF10':
        ax10.plot(df[salunit], df[tempunit], 
            label=label,
             zorder=4, c=color_date[Date], linestyle = linestyle_st[Station]
            )
    elif Station == 'GF3':
        ax.plot(df[salunit], df[tempunit], 
            label=label,
             zorder=4, c=color_date[Date], linestyle = linestyle_st[Station]
            )
    elif Station == 'GF5':
        ax.plot(df[salunit], df[tempunit], 
            label=label,
             zorder=4, c=color_date[Date], linestyle = linestyle_st[Station]
            )
    plotly_figure.add_trace(go.Scatter(x=df.CT, y=df.Depth, mode='lines', name=f"{df.Date.iloc[0]},  {df.Name.iloc[0]}, {df['St.'].iloc[0]}"), row=1, col=1)
    plotly_figure.add_trace(go.Scatter(x=df.SA, y=df.Depth, mode='lines', name=f"{df.Date.iloc[0]},  {df.Name.iloc[0]}, {df['St.'].iloc[0]}"), row=1, col=2)    
    

max_sal = both_casts[540][pd.to_datetime('2019-04-10'):].sal.idxmax()

for single_ax in [ax, ax10]:
    single_ax.scatter(
        both_casts[540].loc[max_sal].sal,both_casts[540].loc[max_sal].temp,zorder=6, s=50, c='black', marker="*", alpha=0.5, label='Cold dense water mass \n April/May 2019')
ax10.plot(
    both_casts[540].rolling(window=50*6, center=True).mean().sal, 
    both_casts[540].rolling(window=50*6, center=True).mean().temp, 
    # s=1, 
    label="Mooring 540 m", zorder=2, c='grey')
# ax10.scatter(
#     both_casts[540][:pd.to_datetime('2019-05-06')].sal[0], 
#     both_casts[540][:pd.to_datetime('2019-05-06')].temp[0],
#     marker="*", color = 'black', s=50, zorder=5, label ='Begin Mooring')


# ax.legend(loc='upper left', bbox_to_anchor=(1, 1.01),
#             ncol=1)
handles, labels = [], []
for single_ax in [ax, ax10]:
    hs, ls = single_ax.get_legend_handles_labels()
    handles.extend(hs)
    labels.extend(ls)
by_label = (dict(zip(labels, handles)))
by_label = dict(sorted(by_label.items()))
fig.legend(by_label.values(), by_label.keys(), loc='lower left', ncol=1, bbox_to_anchor=(1, 0.15), frameon=True)
ax.set_xlabel('Salinity')
ax.set_ylabel('Potential temperature [°C]')
ax10.set_xlabel('Salinity')

# label a b 
ax.text(0.05, 0.95, 'a)', transform=ax.transAxes,fontsize=16, fontweight='bold', va='top')
ax10.text(0.05, 0.95, 'b)', transform=ax10.transAxes,fontsize=16, fontweight='bold', va='top')
# turn off grid
ax.grid(False)
plt.tight_layout()
# MAKE ALL Text larger


fig.savefig(f"{figpath}/Event/TS_CTD_transects_inc_GF3_GF10_selection2.jpg", dpi=600, bbox_inches='tight')
fig.savefig(f"{figpath}/Event/TS_CTD_transects_inc_GF3_GF10_selection2_draft.jpg", dpi=150, bbox_inches='tight')
plotly_figure.update_yaxes(title_text="Depth [m]", row=1, col=1)
plotly_figure.update_xaxes(title_text="CT [°C]", row=1, col=1)
plotly_figure.update_xaxes(title_text="SA [g/kg]", row=1, col=2)
plotly_figure.show()
#%%
def add_timedate_col(df, time_column='timeJ', date_column='Date'):
    '''Creates datetime column from timeJ and Date columns called time_date
    define time_column and date_column if not default, 
    Timezone unaware'''
    df['time_date'] = pd.to_datetime(df_all[date_column])+ pd.to_timedelta(df_all[time_column]%1, 'days')
#%%
names = df_all.Name.unique()
# Select all values in names that start with 'HS'
names_GF3 = [name for name in names if name.startswith('HS')] 
names_GF10 = df_all[df_all['St.'] == 'GF10'].Name.unique()


# %%
pd.options.plotting.backend = "plotly"
var = 'CT'
fig = df_all[(df_all.Name.isin(names_GF10))].plot.line(y='Pressure', x=var,color='Date')
# plt.ylim([350,0])
# plt.show()
fig.update_yaxes(autorange="reversed")
fig.show()

var = 'SA'
fig = df_all[(df_all.Name.isin(df_spring_HS))].plot.line(y='Pressure', x=var,color='Date', hover_data=['Date', 'Name', 'St.'])
fig.update_yaxes(autorange="reversed")
fig.show()
# plt.ylim([350,0])
# plt.show()

var = 'Sigma_dens'
fig = df_all[(df_all.Name.isin(df_spring_HS))].plot.line(y='Pressure', x=var,color='Date')
fig.update_yaxes(autorange="reversed")
fig.show()

#%%
import plotly.graph_objects as go
def add_sigma_contour(fig, Tvals = np.linspace(-1, 7, 100), Svals = np.linspace(26, 35, 100)):
    Tg,Sg = np.meshgrid((Tvals),(Svals))
    sigma = gsw.sigma0(Sg,Tg)
    fig.add_trace(go.Contour(
        x=Sg.flatten(),
        y=Tg.flatten(),
        z=sigma.flatten(),
        contours_coloring='lines',
        ))
    fig.update_yaxes(range=[min(Tvals), max(Tvals)])
    fig.update_xaxes(range=[min(Svals), max(Svals)])
    fig.update_layout(template='plotly_white')
# make contour lines black


fig = df_all[(df_all.Name.isin(names_GF3))].plot.scatter(y='CT', x='SA',color='timeJ', hover_data=['Date', 'Name', 'St.'], color_continuous_scale='Phase')

# add_sigma_contour(fig)


#add traces contour line 
# fig.add_trace(go.Contour(
#     x=df_all[(df_all.Name.isin(names_GF3))].SA,
#     y=df_all[(df_all.Name.isin(names_GF3))].CT,
#     z=df_all[(df_all.Name.isin(names_GF3))].Sigma_dens,
#     contours_coloring='lines',
#     ))
fig.show()

# %%

fig =go.Figure()


# %% ======================================
# Plotting inflow moment
pd.options.plotting.backend = "matplotlib"
plt.style.use('seaborn-ticks')
fig, axs = plt.subplots(2,1, sharex=True, figsize=(10,8))

dates = ['2019-04-01', '2019-06-01']
var = 'CT'
both_casts[540].rolling(window=25*6, center=True).mean()['2019-04':'2019-05'].plot( y=var, use_index= True,ax=axs[0], kind='line')
both_casts[540]['2019-04':'2019-05'].reset_index().plot( y=var, x='date' ,ax=axs[0],kind='scatter', s=1, ylabel='CT [°C]')
# add annotation a) in bottom left
axs[0].annotate('a)', xy=(0.01, 0.05), xycoords='axes fraction', fontsize=20)

var = 'SA'
both_casts[540].rolling(window=25*6, center=True).mean()['2019-04':'2019-05'].plot( y=var , use_index=True, ax=axs[1], kind='line')
both_casts[540]['2019-04':'2019-05'].reset_index().plot( y=var, x='date' ,ax=axs[1],kind='scatter', s=1, ylabel= 'SA [g/kg]')
# add annotation b) in bottom left
axs[1].annotate('b)', xy=(0.01, 0.05), xycoords='axes fraction', fontsize=20)

# use dates on x axis and tilt labels, remove legend
# set m
for ax in axs:
    ax.xaxis.set_major_formatter( mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    # set minor locator every day
    ax.xaxis.set_minor_locator(mdates.DayLocator())

    ax.tick_params(axis='x', labelrotation=45)
    ax.grid(False)
    ax.legend().remove()
    ax.set_xlabel('')
    ax.set_xlim(dates)

# set major ticks every week
# ax.xaxis.set_major_locator(mdates.WedayLocator())
fig.tight_layout()
fig.savefig(f"{figpath}/Event/Mooring_timeseries_spring_zoom_TS.png", dpi=600)


# plot timeseries of 

# %% ======================================
# Plotting inflow moment
pd.options.plotting.backend = "matplotlib"
plt.style.use('seaborn-ticks')
fig, axs = plt.subplots(2,1, sharex=True, figsize=(10,8))

dates = ['2018-10-01', '2019-01-01']
var = 'CT'
both_casts[330].rolling(window=25*6, center=True).mean()[dates[0]:dates[1]].plot( y=var, use_index= True,ax=axs[0], kind='line')
both_casts[330][dates[0]:dates[1]].reset_index().plot( y=var, x='date' ,ax=axs[0],kind='scatter', s=1, ylabel='CT [°C]')
# add annotation a) in bottom left
axs[0].annotate('a)', xy=(0.01, 0.05), xycoords='axes fraction', fontsize=20)

var = 'SA'
both_casts[330].rolling(window=25*6, center=True).mean()[dates[0]:dates[1]].plot( y=var , use_index=True, ax=axs[1], kind='line')
both_casts[330][dates[0]:dates[1]].reset_index().plot( y=var, x='date' ,ax=axs[1],kind='scatter', s=1, ylabel= 'SA [g/kg]')
# add annotation b) in bottom left
axs[1].annotate('b)', xy=(0.01, 0.05), xycoords='axes fraction', fontsize=20)

# use dates on x axis and tilt labels, remove legend
# set m
for ax in axs:
    ax.xaxis.set_major_formatter( mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    # set minor locator every day
    ax.xaxis.set_minor_locator(mdates.DayLocator())

    ax.tick_params(axis='x', labelrotation=45)
    ax.grid(False)
    ax.legend().remove()
    ax.set_xlabel('')
    ax.set_xlim(dates)

# set major ticks every week
# ax.xaxis.set_major_locator(mdates.WedayLocator())
fig.tight_layout()
# fig.savefig(f"{figpath}/Event/Mooring_timeseries_autumn_zoom_TS.png", dpi=600)
# %%
