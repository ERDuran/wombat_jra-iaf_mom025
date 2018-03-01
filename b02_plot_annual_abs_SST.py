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
                    'ocean_month_temp_ncrcat_1958_to_2014' 
                    + '.nc', 'r')

# same for temperature (1,2,3)
abs_temp = nc_fid.variables['temp'][:,0,:,:]

# get dimensions
lat = nc_fid.variables['yt_ocean'][:]
lon1 = nc_fid.variables['xt_ocean'][:]
lon = lon1 + 360
z = nc_fid.variables['st_ocean'][0]

#
print('ocean_month_temp OK !')
            
            
#%% Calculate projection. mill is 'Miller Cylindrical'
# gall is 'Gall Stereographic Equidistant.
# takes some time to run...............
Bm = Basemap(projection='mill', llcrnrlat=-50,urcrnrlat=-10,\
llcrnrlon=100,urcrnrlon=170, resolution='c')

pickle.dump(Bm,open('map.pickle','wb'),-1)  # pickle it 

Bm.fix_aspect = True
m_aspect = Bm.aspect


#%% plot surface maps: TEMPERATURE
matplotlib.rcParams.update({'font.size': 14}) 
row = 1
col = 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

idx = -1
for y in years:
    plt.close('all') # close all existing figures
    fig = plt.figure() # generate figure
    fig.set_size_inches(12, 10) # set figure size in inches
    
    idx = idx + 1

    # position figure wrt window template
    ax = fig.add_subplot(row,col,1)
    pos = ax.get_position()
    bnd = list(pos.bounds)
    magn = 0.04
    bnd = [bnd[0], bnd[1], bnd[2]+magn, bnd[3]+magn*m_aspect]
    ax.set_position(bnd)
    
    # colourmap: blue to red because looking at bias.
    cmap = plt.get_cmap('gist_ncar')
    step = 1
    # levels to show on colourbar
    contf_lvls = np.arange(6,30+1e-08,step)
        
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
    contf = Bm.contourf(Bm_lons, Bm_lats, abs_temp[idx,:,:], 
                       contf_lvls, cmap=cmap, extend='both')
    
    # title ...
    ax.set_title('z=' + str(round(z)) + ' m ' + str(y))
    
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
    cbar.set_label(r'Absolute temperature $^{\circ}C$') # units label on colourbar
    cbar.dividers.set_linewidth(0.2)
    cbar.outline.set_linewidth(0.5)
    cbar.set_ticks(contf_lvls[np.arange(0,np.size(contf_lvls),2)])
    cbar.ax.tick_params(width=0.2, length= 2)
    cbar.ax.yaxis.set_label_position('left')
    cbar.ax.yaxis.set_ticks_position('left')
    
        
    print(str(y) + ' OK !')
            
    plt.suptitle(r"Wombat-JRA-MOM025 run. Absolute annual mean (1958 to 2014)", 
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

