"""
This code generates a geographical visualization of Wastewater Treatment Plant (WWTP) locations in California,
along with the Social Vulnerability Index (SVI) at the county level. It emphasizes the relationship between
WWTPs, their population served, and the social vulnerability of surrounding areas. Data for WWTP locations is
sourced from CDPH WWTPs localities.

"""
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import geopandas as gpd
from shapely.geometry import Point
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pgeocode
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from main_cluster import main_data_cluster

save_path = "../output/"
read_path = "../data/"
shapefile_path = read_path+'geodata/CA_Counties_TIGER2016.shp' # polygons for counties
gdf = gpd.read_file(shapefile_path)
mpl.rcParams['font.size'] = 18

#scenario = 'diss' # cdph plants except those provide by scan
#dissim_term = True

scenario = 'all' # cdph plants except those provide by scan
dissim_term = False

mc=main_data_cluster(scenario, dissim_term)

def data_SVI():
    # Load data for social social vulnerability index (svi)
    data_svi = pd.read_csv(read_path+"svi_interactive_map.csv")
    data_svi = data_svi.rename(columns={'COUNTY': 'County'})
    data_svi = data_svi.rename(columns={'RPL_THEMES': 'svi'})
    data_svi = data_svi[['County', 'svi']]
    for i in range(len(data_svi.index)):
        data_svi.at[i, 'County'] = data_svi.at[i, 'County'][0: -7]
    return data_svi


def map_svi_all1(cbar):
    mpl.rcParams['font.size'] = 18
    city_df = mc.data_sl

    # read shape file to plot counties polygons
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.rename(columns={'NAME': 'County'})

    # read and merged SVI data for all California counties
    data_svi = data_SVI()
    gdf = pd.merge(gdf, data_svi, on='County')

    fig, ax = plt.subplots(figsize=(8, 9))
    if cbar:
        fig, ax = plt.subplots(figsize=(12, 10))


    cmap_svi = plt.cm.get_cmap('Greens')
    #cmap_svi = plt.cm.get_cmap('Blues')

    # Set the minimum and maximum values for SVI
    svi_min = gdf['svi'].min()
    svi_max = gdf['svi'].max()
    if cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        gdf.plot(column='svi', cmap=cmap_svi, linewidth=0.8, ax=ax, legend=True, cax=cax, vmin=svi_min, vmax=svi_max)
    else:
        gdf.plot(column='svi', cmap=cmap_svi, linewidth=0.8, ax=ax, vmin=svi_min, vmax=svi_max)

    ax.set_xticks([])
    ax.set_yticks([])
    #ax.set_title("California Social Vulnerability Index (SVI)")

    # Plot city markers with population-based color
    geometry = [Point(lon, lat) for lon, lat in zip(city_df['lng'], city_df['lat'])]
    crs = {'init': 'epsg:4326'}
    city_gdf = gpd.GeoDataFrame(city_df, crs=crs, geometry=geometry)
    city_gdf = city_gdf.to_crs(gdf.crs)

    cmap_population = plt.cm.get_cmap('YlOrRd_r')
    #cmap_population = plt.cm.get_cmap('PuRd')
    #cmap_population= cmap_population.reversed()
    norm_population = Normalize(vmin=city_gdf['Population_Served'].min(), vmax=city_gdf['Population_Served'].max())
    sm_population = ScalarMappable(cmap=cmap_population, norm=norm_population)
    sm_population.set_array([])

    city_gdf.plot(ax=ax, marker='o', color=sm_population.to_rgba(city_gdf['Population_Served']), markersize=30)

    # Create a separate color bar for population served
    if cbar:
        cbar_population = plt.colorbar(sm_population, ax=ax, label='Population Served')
        cbar_population.ax.set_ylabel('Population Served', fontsize=16)

        # Create a color bar for SVI
        cbar_svi = plt.colorbar(ScalarMappable(cmap=cmap_svi, norm=Normalize(vmin=svi_min, vmax=svi_max)), cax=cax)
        cbar_svi.set_label('Social Vulnerability Index (SVI)', fontsize=16)
    plt.tight_layout()
    if cbar:
        plt.savefig('figures/map_with_svi_and_population_'+mc.scenario+'.png', dpi=300)
    else:
        plt.savefig('figures/map_with_svi_and_population_'+mc.scenario+'_not_bar.png', dpi=300)






def map_svi_all_joint(scenario, dissim_term, cbar, ax, fig):
    mc = main_data_cluster(scenario, dissim_term)

    mpl.rcParams['font.size'] = 14
    city_df = mc.data_sl

    # read shape file to plot counties polygons
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.rename(columns={'NAME': 'County'})

    # read and merged SVI data for all California counties
    data_svi = data_SVI()
    gdf = pd.merge(gdf, data_svi, on='County')

    #fig, ax = plt.subplots(figsize=(8, 9))
    #if cbar:
    #    fig, ax = plt.subplots(figsize=(12, 10))


    cmap_svi = plt.cm.get_cmap('Greens')

    # Set the minimum and maximum values for SVI
    svi_min = gdf['svi'].min()
    svi_max = gdf['svi'].max()
    if cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        gdf.plot(column='svi', cmap=cmap_svi, linewidth=0.8, ax=ax, legend=True, cax=cax, vmin=svi_min, vmax=svi_max)
    else:
        gdf.plot(column='svi', cmap=cmap_svi, linewidth=0.8, ax=ax, vmin=svi_min, vmax=svi_max)

    ax.set_xticks([])
    ax.set_yticks([])
    #ax.set_title("California Social Vulnerability Index (SVI)")

    # Plot city markers with population-based color
    geometry = [Point(lon, lat) for lon, lat in zip(city_df['lng'], city_df['lat'])]
    crs = {'init': 'epsg:4326'}
    city_gdf = gpd.GeoDataFrame(city_df, crs=crs, geometry=geometry)
    city_gdf = city_gdf.to_crs(gdf.crs)

    cmap_population = plt.cm.get_cmap('hot')
    norm_population = Normalize(vmin=city_gdf['Population_Served'].min(), vmax=city_gdf['Population_Served'].max())
    sm_population = ScalarMappable(cmap=cmap_population, norm=norm_population)
    sm_population.set_array([])

    city_gdf.plot(ax=ax, marker='o', color=sm_population.to_rgba(city_gdf['Population_Served']), markersize=8)

    # Create a separate color bar for population served
    if cbar:
        cbar_population = plt.colorbar(sm_population, ax=ax, label='Population Served')
        cbar_population.ax.set_ylabel('Population Served', fontsize=14)

        # Create a color bar for SVI
        cbar_svi = plt.colorbar(ScalarMappable(cmap=cmap_svi, norm=Normalize(vmin=svi_min, vmax=svi_max)), cax=cax)
        cbar_svi.set_label('Social Vulnerability Index (SVI)', fontsize=14)
    plt.tight_layout()
    #if cbar:
    #    plt.savefig('figures/map_with_svi_and_population_'+mc.scenario+'.png', dpi=300)
    #else:
    #    plt.savefig('figures/map_with_svi_and_population_'+mc.scenario+'_not_bar.png', dpi=300)

def plot_all():
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    map_svi_all_joint(scenario='all', dissim_term=False, cbar=False, ax=ax[0], fig=fig)
    ax[0].text(0.9, 0.95, "A", ha='center', va='center', transform=ax[0].transAxes, fontsize=16)

    map_svi_all_joint(scenario='diss', dissim_term=True, cbar=True, ax=ax[1], fig=fig)
    ax[1].text(0.9, 0.95, "B", ha='center', va='center', transform=ax[1].transAxes, fontsize=16)
    plt.savefig('figures/map_with_svi_and_population.png', dpi=300)


#plot_all()
map_svi_all1(cbar=True)