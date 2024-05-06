#from sa_ww import simulated_annealing_opt, data_plants_read

from sa_ww_par import SimulatedAnnealingOptimization
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pgeocode
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib
import time
import numpy as np
'''
This code varies the number of maximum plant in solution (k), it gets the solution
for each of them and save a .csv with the plants in that solution. It also saves
a .csv containing the energy and population served for each solution. It saves a 
scatterplot .png of the energy as a function of k 

'''

#scenario = 'cdph_w_scan'
#dissim_term = True

shapefile_path = '../data/geodata/CA_Counties_TIGER2016.shp'  # polygons for counties
#scenario = 'diss'
#dissim_term = True

scenario = 'all'
dissim_term = False

sa_op = SimulatedAnnealingOptimization(scenario=scenario, dissim_term=dissim_term)  # read optimization class
scn = scenario

if scenario=='all':
    plans1_a = np.load('../output/sol_plants_names_all_ps_0.00_pd_0.00_ds_0.00_svd_1.00_sim_0.00.npy')
    plans1_b = np.load('../output/sol_plants_names_all_ps_1.00_pd_0.00_ds_0.00_svd_0.00_sim_0.00.npy')
    plans1_c = np.load('../output/sol_plants_names_all_ps_0.25_pd_0.25_ds_0.25_svd_0.25_sim_0.00.npy')

else:
    plans1_a = np.load('../output/sol_plants_names_diss_ps_0.00_pd_0.00_ds_0.00_svd_0.00_sim_1.00.npy')
    plans1_b = np.load('../output/sol_plants_names_diss_ps_0.20_pd_0.20_ds_0.20_svd_0.20_sim_0.20.npy')
    plans1_c = np.load('../output/sol_plants_names_diss_ps_0.10_pd_0.10_ds_0.10_svd_0.35_sim_0.35.npy')


def data_SVI():

    # Load data for social social vulnerability index (svi)
    data_svi = pd.read_csv("../data/svi_interactive_map.csv")
    data_svi = data_svi.rename(columns={'COUNTY': 'County'})
    data_svi = data_svi.rename(columns={'RPL_THEMES': 'svi'})
    data_svi = data_svi[['County', 'svi']]
    for i in range(len(data_svi.index)):
        data_svi.at[i, 'County'] = data_svi.at[i, 'County'][0: -7]
    return data_svi

def map_svi(cbar, cbar_svi, plants_name, all, ax, fig):
    mpl.rcParams['font.size'] = 14
    if all:
        city_df= sa_op.data_plants
    else:
        city_df = sa_op.data_plants_read[sa_op.data_plants_read.Plant.isin(plants_name)]

    # read shape file to plot counties polygons
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.rename(columns={'NAME': 'County'})

   # if cbar:
   #     fig, ax = plt.subplots(figsize=(12, 10))
   # else:
   #     fig, ax = plt.subplots(figsize=(8, 9))

    # read and merged SVI data for all California counties
    data_svi = data_SVI()
    gdf = pd.merge(gdf, data_svi, on='County')

    cmap_svi =matplotlib.colormaps.get_cmap('Greens')

    # Set the minimum and maximum values for SVI
    svi_min = sa_op.data_plants_read['svi'].min()
    svi_max = sa_op.data_plants_read['svi'].max()
    divider = make_axes_locatable(ax)
    if cbar:
        #divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.1)
        gdf.plot(column='svi', cmap=cmap_svi, linewidth=0.8, ax=ax, legend=True, cax=cax, vmin=svi_min, vmax=svi_max, alpha=0.85)
    else:
        gdf.plot(column='svi', cmap=cmap_svi, linewidth=0.8, ax=ax, vmin=svi_min, vmax=svi_max, alpha=0.85)

    ax.set_xticks([])
    ax.set_yticks([])
    #ax.set_title("California Social Vulnerability Index (SVI)")

    # Plot city markers with population-based color
    geometry = [Point(lon, lat) for lon, lat in zip(city_df['lng'], city_df['lat'])]
    crs = {'init': 'epsg:4326'}
    city_gdf = gpd.GeoDataFrame(city_df, crs=crs, geometry=geometry)
    city_gdf = city_gdf.to_crs(gdf.crs)

    #cmap_population = plt.cm.get_cmap('hot')
    cmap_population = matplotlib.colormaps.get_cmap('YlOrRd_r')

    norm_population = Normalize(vmin=sa_op.data_plants_read['Population_Served'].min(), vmax=sa_op.data_plants_read['Population_Served'].max())
    sm_population = ScalarMappable(cmap=cmap_population, norm=norm_population)
    sm_population.set_array([])

    city_gdf.plot(ax=ax, marker='o', color=sm_population.to_rgba(city_gdf['Population_Served']), markersize=8)

    if cbar_svi:
    # Create a separate color bar for population served
        #cbar_population = plt.colorbar(sm_population, ax=ax, label='Population Served')
        #cbar_population.ax.set_ylabel('Population Served', fontsize=14)
        cax_population = divider.append_axes("left", size="3%", pad=0.1)
        cbar_population = plt.colorbar(sm_population, ax=ax, cax=cax_population, location='left')
        #cbar_population.yaxis.set_ticks_position('left')

        cbar_population.ax.set_ylabel('Population Served', fontsize=15)
    if cbar:
    # Create a color bar for SVI
        cbar_svi = plt.colorbar(ScalarMappable(cmap=cmap_svi, norm=Normalize(vmin=svi_min, vmax=svi_max)), cax=cax)
        cbar_svi.set_label('Social Vulnerability Index (SVI)', fontsize=15)
    #plt.tight_layout()


#map_svi()

#I want to 2x2 plot that have
#ax[0,0]= map_svi(cbar=False, plants_name=plans1_a, all=True)
#ax[0,1]= map_svi(cbar=True, plants_name=plans1_a, all=False)
#ax[1,0]= map_svi(cbar=False, plants_name=plans1_b, all=False)
#ax[1,1]= map_svi(cbar=True, plants_name=plans1_c, all=False)

# Create figure and axes objects
fig, ax = plt.subplots(2, 2, figsize=(10, 10))
#fig, ax = plt.subplots(2, 2, figsize=(12, 20), constrained_layout=True)

#fig, ax = plt.subplots(2, 2, figsize=(12, 25))

# Plot each subplot
map_svi(cbar=False,  cbar_svi=True, plants_name=plans1_a, all=True, ax=ax[0, 0], fig=fig)
ax[0, 0].text(0.9, 0.95, "A", ha='center', va='center', transform=ax[0, 0].transAxes, fontsize=18)

map_svi(cbar=True, cbar_svi=False, plants_name=plans1_a, all=False, ax=ax[0, 1], fig=fig)
ax[0, 1].text(0.9, 0.95, "B", ha='center', va='center', transform=ax[0, 1].transAxes, fontsize=18)

map_svi(cbar=False, cbar_svi=True, plants_name=plans1_b, all=False, ax=ax[1, 0], fig=fig)
ax[1, 0].text(0.9, 0.95, "C", ha='center', va='center', transform=ax[1, 0].transAxes, fontsize=18)

map_svi(cbar=True, cbar_svi=False, plants_name=plans1_c, all=False, ax=ax[1,1], fig=fig)
ax[1, 1].text(0.9, 0.95, "D", ha='center', va='center', transform=ax[1, 1].transAxes, fontsize=18)

#ax[0, 0].set_rasterized(True)
#ax[0, 1].set_rasterized(True)
#ax[1, 0].set_rasterized(True)
#ax[1, 1].set_rasterized(True)

#plt.savefig('figures/sol_plants_' + scenario+ '.eps', dpi=600)

# Adjust layout
plt.tight_layout()

# Show plot
#plt.show()


#plt.savefig('figures/box_plot_scn_'+scn+  '.eps', dpi=900)

#plt.savefig('figures/sol_plants_' + scenario + '.png', dpi=300)
plt.savefig('figures/sol_plants_' + scenario + '.pdf', dpi=900)