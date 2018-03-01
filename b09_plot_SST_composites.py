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

years = np.arange(1958,2015,1)

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
                 
                 
# Getting back the objects:
with open('SST_TS_stdev.pkl') as f:  # Python 3: open(..., 'rb')
    EAC_ext_temp_anom_mean, EAC_sep_temp_anom_mean,\
                 south_LC_temp_anom_mean, LCE_temp_anom_mean, \
                 EAC_ext_temp_anom_stdev1, EAC_sep_temp_anom_stdev1,\
                 south_LC_temp_anom_stdev1, LCE_temp_anom_stdev1 \
                 = pickle.load(f)
    

#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
Bm = Basemap(projection='mill', llcrnrlat=-50,urcrnrlat=-10,\
llcrnrlon=100,urcrnrlon=170, resolution='c')
Bm.fix_aspect = True
m_aspect = Bm.aspect


#%% 1. EAC Ext Warm
EAC_ext_compos_warm_std1p25_idx = \
EAC_ext_temp_anom_mean >= EAC_ext_temp_anom_stdev1*1.25

years_EAC_ext_compos_warm_std1p25_idx = \
np.where(EAC_ext_compos_warm_std1p25_idx)[0]

years_EAC_ext_compos_warm_std1p25 = \
years[years_EAC_ext_compos_warm_std1p25_idx]

warm_compos_std1p25_idx_EAC_ext = \
np.ma.mean(temp_anom[years_EAC_ext_compos_warm_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, warm_compos_std1p25_idx_EAC_ext[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_EAC_ext, Bm_lat_EAC_ext = Bm(\
[lon_EAC_ext[0],lon_EAC_ext[0],lon_EAC_ext[1],lon_EAC_ext[1],lon_EAC_ext[0]],\
[lat_EAC_ext[0],lat_EAC_ext[1],lat_EAC_ext[1],lat_EAC_ext[0],lat_EAC_ext[0]])

plot = Bm.plot(Bm_lon_EAC_ext,Bm_lat_EAC_ext)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_EAC_ext_compos_warm_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. EAC Extension SST Warm Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig1_.jpeg', bbox_inches='tight', dpi=200)


#%% 2. EAC Ext Cold
EAC_ext_compos_warm_std1p25_idx = \
EAC_ext_temp_anom_mean <= -EAC_ext_temp_anom_stdev1*1.25

years_EAC_ext_compos_warm_std1p25_idx = \
np.where(EAC_ext_compos_warm_std1p25_idx)[0]

years_EAC_ext_compos_warm_std1p25 = \
years[years_EAC_ext_compos_warm_std1p25_idx]

warm_compos_std1p25_idx_EAC_ext = \
np.ma.mean(temp_anom[years_EAC_ext_compos_warm_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, warm_compos_std1p25_idx_EAC_ext[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_EAC_ext, Bm_lat_EAC_ext = Bm(\
[lon_EAC_ext[0],lon_EAC_ext[0],lon_EAC_ext[1],lon_EAC_ext[1],lon_EAC_ext[0]],\
[lat_EAC_ext[0],lat_EAC_ext[1],lat_EAC_ext[1],lat_EAC_ext[0],lat_EAC_ext[0]])

plot = Bm.plot(Bm_lon_EAC_ext,Bm_lat_EAC_ext)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_EAC_ext_compos_warm_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. EAC Extension SST Cold Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig2_.jpeg', bbox_inches='tight', dpi=200)


#%% 3. EAC Sep Warm
EAC_sep_compos_warm_std1p25_idx = \
EAC_sep_temp_anom_mean >= EAC_sep_temp_anom_stdev1*1.25

years_EAC_sep_compos_warm_std1p25_idx = \
np.where(EAC_sep_compos_warm_std1p25_idx)[0]

years_EAC_sep_compos_warm_std1p25 = \
years[years_EAC_sep_compos_warm_std1p25_idx]

warm_compos_std1p25_idx_EAC_sep = \
np.ma.mean(temp_anom[years_EAC_sep_compos_warm_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, warm_compos_std1p25_idx_EAC_sep[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_EAC_sep, Bm_lat_EAC_sep = Bm(\
[lon_EAC_sep[0],lon_EAC_sep[0],lon_EAC_sep[1],lon_EAC_sep[1],lon_EAC_sep[0]],\
[lat_EAC_sep[0],lat_EAC_sep[1],lat_EAC_sep[1],lat_EAC_sep[0],lat_EAC_sep[0]])

plot = Bm.plot(Bm_lon_EAC_sep,Bm_lat_EAC_sep)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_EAC_sep_compos_warm_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. EAC Separation SST Warm Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig3_.jpeg', bbox_inches='tight', dpi=200)


#%% 4. EAC Sep Cold
EAC_sep_compos_cold_std1p25_idx = \
EAC_sep_temp_anom_mean <= -EAC_sep_temp_anom_stdev1*1.25

years_EAC_sep_compos_cold_std1p25_idx = \
np.where(EAC_sep_compos_cold_std1p25_idx)[0]

years_EAC_sep_compos_cold_std1p25 = \
years[years_EAC_sep_compos_cold_std1p25_idx]

cold_compos_std1p25_idx_EAC_sep = \
np.ma.mean(temp_anom[years_EAC_sep_compos_cold_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, cold_compos_std1p25_idx_EAC_sep[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_EAC_sep, Bm_lat_EAC_sep = Bm(\
[lon_EAC_sep[0],lon_EAC_sep[0],lon_EAC_sep[1],lon_EAC_sep[1],lon_EAC_sep[0]],\
[lat_EAC_sep[0],lat_EAC_sep[1],lat_EAC_sep[1],lat_EAC_sep[0],lat_EAC_sep[0]])

plot = Bm.plot(Bm_lon_EAC_sep,Bm_lat_EAC_sep)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_EAC_sep_compos_cold_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. EAC Separation SST Cold Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)


# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig4_.jpeg', bbox_inches='tight', dpi=200)


#%% 5. South LC Warm
south_LC_compos_warm_std1p25_idx = \
south_LC_temp_anom_mean >= south_LC_temp_anom_stdev1*1.25

years_south_LC_compos_warm_std1p25_idx = \
np.where(south_LC_compos_warm_std1p25_idx)[0]

years_south_LC_compos_warm_std1p25 = \
years[years_south_LC_compos_warm_std1p25_idx]

warm_compos_std1p25_idx_south_LC = \
np.ma.mean(temp_anom[years_south_LC_compos_warm_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, warm_compos_std1p25_idx_south_LC[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_south_LC, Bm_lat_south_LC = Bm(\
[lon_south_LC[0],lon_south_LC[0],lon_south_LC[1],lon_south_LC[1],lon_south_LC[0]],\
[lat_south_LC[0],lat_south_LC[1],lat_south_LC[1],lat_south_LC[0],lat_south_LC[0]])

plot = Bm.plot(Bm_lon_south_LC,Bm_lat_south_LC)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_south_LC_compos_warm_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. Leeuwin Current SST Warm Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)

# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig5_.jpeg', bbox_inches='tight', dpi=200)


#%% 6. South LC Cold
south_LC_compos_cold_std1p25_idx = \
south_LC_temp_anom_mean <= -south_LC_temp_anom_stdev1*1.25

years_south_LC_compos_cold_std1p25_idx = \
np.where(south_LC_compos_cold_std1p25_idx)[0]

years_south_LC_compos_cold_std1p25 = \
years[years_south_LC_compos_cold_std1p25_idx]

cold_compos_std1p25_idx_south_LC = \
np.ma.mean(temp_anom[years_south_LC_compos_cold_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, cold_compos_std1p25_idx_south_LC[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_south_LC, Bm_lat_south_LC = Bm(\
[lon_south_LC[0],lon_south_LC[0],lon_south_LC[1],lon_south_LC[1],lon_south_LC[0]],\
[lat_south_LC[0],lat_south_LC[1],lat_south_LC[1],lat_south_LC[0],lat_south_LC[0]])

plot = Bm.plot(Bm_lon_south_LC,Bm_lat_south_LC)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_south_LC_compos_cold_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. Leeuwin Current SST Cold Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)

# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig6_.jpeg', bbox_inches='tight', dpi=200)


#%% 7. LCE Warm
LCE_compos_warm_std1p25_idx = \
LCE_temp_anom_mean >= LCE_temp_anom_stdev1*1.25

years_LCE_compos_warm_std1p25_idx = \
np.where(LCE_compos_warm_std1p25_idx)[0]

years_LCE_compos_warm_std1p25 = \
years[years_LCE_compos_warm_std1p25_idx]

warm_compos_std1p25_idx_LCE = \
np.ma.mean(temp_anom[years_LCE_compos_warm_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, warm_compos_std1p25_idx_LCE[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_LCE, Bm_lat_LCE = Bm(\
[lon_LCE[0],lon_LCE[0],lon_LCE[1],lon_LCE[1],lon_LCE[0]],\
[lat_LCE[0],lat_LCE[1],lat_LCE[1],lat_LCE[0],lat_LCE[0]])

plot = Bm.plot(Bm_lon_LCE,Bm_lat_LCE)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_LCE_compos_warm_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. GAB SST Warm Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)

# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig7_.jpeg', bbox_inches='tight', dpi=200)


#%% 8. LCE Cold
LCE_compos_cold_std1p25_idx = \
LCE_temp_anom_mean <= -LCE_temp_anom_stdev1*1.25

years_LCE_compos_cold_std1p25_idx = \
np.where(LCE_compos_cold_std1p25_idx)[0]

years_LCE_compos_cold_std1p25 = \
years[years_LCE_compos_cold_std1p25_idx]

cold_compos_std1p25_idx_LCE = \
np.ma.mean(temp_anom[years_LCE_compos_cold_std1p25_idx,:,:],axis=0)

matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.close('all') # close all existing figures
fig = plt.figure() # generate figure
fig.set_size_inches(12, 10) # set figure size in inches

# position figure wrt window template
ax = fig.add_subplot(row,col,1)
pos = ax.get_position()
bnd = list(pos.bounds)
magn = 0.04
bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
ax.set_position(bnd)

# colourmap: blue to red because looking at bias.
cmap = plt.get_cmap('seismic')
step = 0.25
# levels to show on colourbar
contf_lvls = np.arange(-2,2+1e-08,step)
    
#cmap = plt.get_cmap('seismic')
#step = 0.02
## levels to show on colourbar
#contf_lvls = np.arange(-0.14,0.14+1e-08,step)            


ax.set_facecolor('grey')

 # meshgrid of lon lats
lons, lats = np.meshgrid(lon, lat)
# use projection template to create x and y axis
Bm_lons, Bm_lats = Bm(lons, lats)
        
        
#
pickle.load(open('map.pickle','rb'))   # load here the above pickle

# draw land outlines
Bm.drawcoastlines(linewidth=0.05)
Bm.fillcontinents(color='white')

# filled contour plot
contf = Bm.contourf(Bm_lons, Bm_lats, cold_compos_std1p25_idx_LCE[:,:], 
                   contf_lvls, cmap=cmap, extend='both')

Bm_lon_LCE, Bm_lat_LCE = Bm(\
[lon_LCE[0],lon_LCE[0],lon_LCE[1],lon_LCE[1],lon_LCE[0]],\
[lat_LCE[0],lat_LCE[1],lat_LCE[1],lat_LCE[0],lat_LCE[0]])

plot = Bm.plot(Bm_lon_LCE,Bm_lat_LCE)

# title ...
ax.set_title('z=' + str(round(z)) + ' m '\
             + str(years_LCE_compos_cold_std1p25))

# meridians. last input is meridians tick label
meridians = map(round_to_base, np.arange(110, 180, 10))
Bm.drawmeridians(meridians, linewidth=0.2, labels=[0,0,0,1])
# parallels. last input is paralles tick label
parallels = map(round_to_base, np.arange(-50, -10, 10))
Bm.drawparallels(parallels, linewidth=0.2, labels=[1,0,0,0])
    
cbar_pos = [bnd[0]-0.08, bnd[1], 0.01, bnd[3]*1] 
cbar_axe = fig.add_axes(cbar_pos)
cbar = plt.colorbar(contf, cax=cbar_axe,
                    orientation='vertical', drawedges=True)
cbar.set_label(r'Temperature anomaly $^{\circ}C$') # units label on colourbar
cbar.dividers.set_linewidth(0.2)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
cbar.ax.tick_params(width=0.2, length= 2)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.yaxis.set_ticks_position('left')
        
plt.suptitle(r"Wombat-JRA-MOM025 run. GAB SST Cold Years composite", 
             y=0.97, fontsize=20)

output_ls = os.listdir(figures_path)

#
if not scriptname:
    scriptname = 'test'

elif scriptname not in output_ls:
    os.mkdir(figures_path + '/' + scriptname)

# save figure
plt.savefig(figures_path + '/' + scriptname + '/' + scriptname[0:3] \
            + '_fig8_.jpeg', bbox_inches='tight', dpi=200)