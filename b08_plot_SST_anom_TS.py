#%% Python file set-up
# 
import os
# NetCDF data handler
import netCDF4 as nc
# numerical computing package
import numpy as np
# basemap toolkit to plot maps
from mpl_toolkits.basemap import Basemap
# command style functions that make matplotlib work like MATLAB
import matplotlib.pyplot as plt
# matplotlib
import matplotlib
# find nearest value
def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()
#
import sys
#
import pickle

os.chdir('/Users/earl/wombat_jra_mom025')
figures_path = '/Users/earl/Dropbox/wombat_jra_mom025/figures/'
data_path = '/Users/earl/Dropbox/Data/wombat_jra_mom025/'

years = range(1958,2015,1)

scriptname = os.path.basename(sys.argv[0])[:-3]

out = {'out':'out'}

print('OK, all set')


def round_to_base(x, base=5):
    return int(base * round(float(x) / base))


#%% Load ocean_month_temp
nc_fid = nc.Dataset(data_path + 
                    'ocean_month_temp_ncrcat_1958_to_2014_ncdiff_ncwa' 
                    + '.nc', 'r')

# same for temperature (1,2,3)
temp_anom = nc_fid.variables['temp'][:,0,:,:]

# get dimensions
lat = nc_fid.variables['yt_ocean'][:]
lon1 = nc_fid.variables['xt_ocean'][:]
lon = lon1 + 360
z = nc_fid.variables['st_ocean'][0]

#
print('ocean_month_temp OK !')


#%%
# Getting back the objects:
with open('areas.pkl') as f:  # Python 3: open(..., 'rb')
    lon_EAC_ext, lat_EAC_ext, lon_EAC_sep, lat_EAC_sep, \
    lon_south_LC, lat_south_LC, lon_LCE, lat_LCE \
                 = pickle.load(f)
    

#%%
lon_EAC_ext_idx = (lon >= lon_EAC_ext[0]) * (lon <= lon_EAC_ext[1])
lon_EAC_sep_idx = (lon >= lon_EAC_sep[0]) * (lon <= lon_EAC_sep[1])
lon_south_LC_idx = (lon >= lon_south_LC[0]) * (lon <= lon_south_LC[1])
lon_LCE_idx = (lon >= lon_LCE[0]) * (lon <= lon_LCE[1])

lat_EAC_ext_idx = (lat <= lat_EAC_ext[0]) * (lat >= lat_EAC_ext[1])
lat_EAC_sep_idx = (lat <= lat_EAC_sep[0]) * (lat >= lat_EAC_sep[1])
lat_south_LC_idx = (lat <= lat_south_LC[0]) * (lat >= lat_south_LC[1])
lat_LCE_idx = (lat <= lat_LCE[0]) * (lat >= lat_LCE[1])


#%%
EAC_ext_temp_anom_mean = [None] * 57
EAC_sep_temp_anom_mean = [None] * 57
south_LC_temp_anom_mean = [None] * 57
LCE_temp_anom_mean = [None] * 57

for y in years:
    idx = years.index(y)
    temp_anom_now = temp_anom[idx,:,:]
    
    EAC_ext_temp_anom = \
    temp_anom_now[:,lon_EAC_ext_idx][lat_EAC_ext_idx,:]
    EAC_ext_temp_anom_mean[idx] = np.ma.mean(EAC_ext_temp_anom)
    
    EAC_sep_temp_anom = \
    temp_anom_now[:,lon_EAC_sep_idx][lat_EAC_sep_idx,:]
    EAC_sep_temp_anom_mean[idx] = np.ma.mean(EAC_sep_temp_anom)
    
    south_LC_temp_anom = \
    temp_anom_now[:,lon_south_LC_idx][lat_south_LC_idx,:]
    south_LC_temp_anom_mean[idx] = np.ma.mean(south_LC_temp_anom)
    
    LCE_temp_anom = \
    temp_anom_now[:,lon_LCE_idx][lat_LCE_idx,:]
    LCE_temp_anom_mean[idx] = np.ma.mean(LCE_temp_anom)
    
    
#%%
EAC_ext_temp_anom_stdev1 = np.std(EAC_ext_temp_anom_mean)
EAC_ext_temp_anom_stdev1p5 = 1.5*np.std(EAC_ext_temp_anom_mean)
EAC_ext_temp_anom_stdev2 = 2*np.std(EAC_ext_temp_anom_mean)

EAC_sep_temp_anom_stdev1 = np.std(EAC_sep_temp_anom_mean)
EAC_sep_temp_anom_stdev1p5 = 1.5*np.std(EAC_sep_temp_anom_mean)
EAC_sep_temp_anom_stdev2 = 2*np.std(EAC_sep_temp_anom_mean)

south_LC_temp_anom_stdev1 = np.std(south_LC_temp_anom_mean)
south_LC_temp_anom_stdev1p5 = 1.5*np.std(south_LC_temp_anom_mean)
south_LC_temp_anom_stdev2 = 2*np.std(south_LC_temp_anom_mean)

LCE_temp_anom_stdev1 = np.std(LCE_temp_anom_mean)
LCE_temp_anom_stdev1p5 = 1.5*np.std(LCE_temp_anom_mean)
LCE_temp_anom_stdev2 = 2*np.std(LCE_temp_anom_mean)

# Saving the objects:
with open('SST_TS_stdev.pkl', 'w') as f:  # Python 3: open(..., 'wb')
    pickle.dump(\
                [EAC_ext_temp_anom_mean, EAC_sep_temp_anom_mean,\
                 south_LC_temp_anom_mean, LCE_temp_anom_mean, \
                 EAC_ext_temp_anom_stdev1, EAC_sep_temp_anom_stdev1,\
                 south_LC_temp_anom_stdev1, LCE_temp_anom_stdev1], f)
    
    
#%% plot surface maps: TEMPERATURE
matplotlib.rcParams.update({'font.size': 14}) 
row = 4
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
# filled contour plot
plot = plt.plot(years, EAC_ext_temp_anom_mean)
plot = plt.plot([1958,2014],\
                [EAC_ext_temp_anom_stdev1,EAC_ext_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [-EAC_ext_temp_anom_stdev1,-EAC_ext_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [EAC_ext_temp_anom_stdev1p5,EAC_ext_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
plot = plt.plot([1958,2014],\
                [-EAC_ext_temp_anom_stdev1p5,-EAC_ext_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
cur_axes = plt.gca()
cur_axes.axes.get_xaxis().set_ticklabels([])
plt.axis([1958, 2014, -1.2, 1.2])
ax.set_title('EAC Extension')
ax.grid(color='black', linestyle='--')
plt.xticks(np.arange(1958,2014,4))

#
ax = fig.add_subplot(row,col,2)
plot = plt.plot(years, EAC_sep_temp_anom_mean)
plot = plt.plot([1958,2014],\
                [EAC_sep_temp_anom_stdev1,EAC_sep_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [-EAC_sep_temp_anom_stdev1,-EAC_sep_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [EAC_sep_temp_anom_stdev1p5,EAC_sep_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
plot = plt.plot([1958,2014],\
                [-EAC_sep_temp_anom_stdev1p5,-EAC_sep_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
cur_axes = plt.gca()
cur_axes.axes.get_xaxis().set_ticklabels([])
plt.axis([1958, 2014, -1.2, 1.2])
ax.set_title('EAC Separation')
ax.grid
ax.grid(color='black', linestyle='--')
plt.xticks(np.arange(1958,2014,4))


#
ax = fig.add_subplot(row,col,3)
plot = plt.plot(years, south_LC_temp_anom_mean)
plot = plt.plot([1958,2014],\
                [south_LC_temp_anom_stdev1,south_LC_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [-south_LC_temp_anom_stdev1,-south_LC_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [south_LC_temp_anom_stdev1p5,south_LC_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
plot = plt.plot([1958,2014],\
                [-south_LC_temp_anom_stdev1p5,-south_LC_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
cur_axes = plt.gca()
cur_axes.axes.get_xaxis().set_ticklabels([])
plt.axis([1958, 2014, -1.2, 1.2])
ax.set_title('South Western Leeuwin')
ax.grid
ax.grid(color='black', linestyle='--')
plt.xticks(np.arange(1958,2014,4))


#
ax = fig.add_subplot(row,col,4)
plot = plt.plot(years, LCE_temp_anom_mean)
plot = plt.plot([1958,2014],\
                [LCE_temp_anom_stdev1,LCE_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [-LCE_temp_anom_stdev1,-LCE_temp_anom_stdev1],\
                linestyle='--',color='green')
plot = plt.plot([1958,2014],\
                [LCE_temp_anom_stdev1p5,LCE_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
plot = plt.plot([1958,2014],\
                [-LCE_temp_anom_stdev1p5,-LCE_temp_anom_stdev1p5],\
                linestyle='--',color='orange')
plt.axis([1958, 2014, -1.2, 1.2])
ax.set_title('Great Australian Bight')
ax.grid
ax.grid(color='black', linestyle='--')
plt.xticks(np.arange(1958,2014,4))


print(str(y) + ' OK !')

    
plt.suptitle(r"Wombat-JRA-MOM025 run. Annual mean SST anomalies (1958 to 2014)", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig1_' + str(y) + '.jpeg', bbox_inches='tight', dpi=200)

