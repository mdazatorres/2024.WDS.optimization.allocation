import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from main_cluster import main_data_cluster
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import pywt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

scenario = 'diss' # cdph plants except those provide by scan
dissim_term = True
# self.scenario = 'cdph_scan' # cdph with scan data
# self.dissim_term = False    # If we can add the dissimilarity term, this when we have wastewater data available

save_path = "../output/"
mc = main_data_cluster(scenario, dissim_term)

data = mc.data_sl
plants = data.Plant.unique()

plans2_a = np.load('../output/sol_plants_names_diss_ps_0.00_pd_0.00_ds_0.00_svd_0.00_sim_1.00.npy')
plans2_b = np.load('../output/sol_plants_names_diss_ps_0.20_pd_0.20_ds_0.20_svd_0.20_sim_0.20.npy')
plans2_c = np.load('../output/sol_plants_names_diss_ps_0.10_pd_0.10_ds_0.10_svd_0.35_sim_0.35.npy')
#plants_sol= plans2_a
def city_plant_dict_create():
    city_counts = {}
    city_plant_dict = data[['Plant', 'City']].drop_duplicates(subset='Plant')
    city_plant_dict.index = city_plant_dict['Plant']
    city_plant_dict = city_plant_dict['City'].to_dict()
    for key, city in city_plant_dict.items():
        if city in city_counts:
            city_counts[city] += 1
            new_city_name = f'{city}_{city_counts[city]}'
            if city=='Tracy':
                new_city_name='Tracy (MH)'
        else:
            city_counts[city] = 1
            new_city_name = city
        city_plant_dict[key] = new_city_name
    return city_plant_dict


def comp_vmin_vmax(metric):
    plants_sol = {'a':plans2_a, 'b': plans2_b, 'c':plans2_c}
    vmax= 0
    for i in ['a', 'b', 'c']:
        sol_plant_series_dict= { k:mc.plant_series_dict[k] for k in plants_sol[i]}
        sol_plant_series_list = list(sol_plant_series_dict.values())
        pairwise_distances_corr = pdist(sol_plant_series_list, metric=metric)
        distance_matrix = squareform(pairwise_distances_corr)
        if distance_matrix.max()>vmax:
            vmax = distance_matrix.max()
    return vmax

def comp_distance_matrix(city_names, metric, case):
    if case=='a':
        plants_sol = plans2_a
    elif case=='b':
        plants_sol = plans2_b
    else:
        plants_sol = plans2_c
    sol_plant_series_dict= { k:mc.plant_series_dict[k] for k in plants_sol}
    sol_plant_series_list = list(sol_plant_series_dict.values())
    pairwise_distances_corr = pdist(sol_plant_series_list, metric=metric)
    distance_matrix = squareform(pairwise_distances_corr)

    plant_names_p = np.vectorize(mc.procces_plant_name)(plants_sol)

    linkage_matrix = linkage(distance_matrix, method='average')
    order = dendrogram(linkage_matrix, no_plot=True)['leaves']
    distance_matrix_sorted = distance_matrix[order, :]
    distance_matrix_sorted = distance_matrix_sorted[:, order]
    #plants_sorted = np.array(plants)[order]
    plants_sorted = np.array(plant_names_p)[order]
    if city_names:
        city_plant_dict = city_plant_dict_create()
        city_names = [city_plant_dict[plant] for plant in plants_sol]
        plants_sorted = np.array(city_names)[order]


    return  distance_matrix_sorted, plants_sorted


def plot_distance_matrix(city_names, metric, case, ax, colorbar):
    vmax = comp_vmin_vmax(metric='correlation')
    distance_matrix_sorted, plants_sorted = comp_distance_matrix(city_names, metric, case)
    #sns.set(font_scale=2)
    #plt.rcParams.update({'font.size': 30})
    #fig, ax = plt.subplots(figsize=(12, 10))

    heatmap=sns.heatmap(distance_matrix_sorted, annot=False, cmap='viridis', xticklabels=plants_sorted,
                yticklabels=plants_sorted, fmt=".2f", cbar_kws={'label': 'Correlation Distance'}, ax=ax, cbar=colorbar,vmin=0,vmax=vmax)
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=30, ha='right', fontsize=14)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation= 30, va='top', fontsize=14)
    plt.tight_layout()
    return heatmap

    #fig.tight_layout()
    #if save:
    #    plt.savefig('figures/dis_matrix_'+scenario+'_' +case +'.png')

colorbar = None


#fig, ax = plt.subplots(2, 2, figsize=(12, 25))






fig, axes = plt.subplots(1, 3, figsize=(18, 6))

heatmap = plot_distance_matrix(city_names=True, metric='correlation', case='a', ax=axes[0], colorbar=colorbar)
colorbar = heatmap.collections[0].colorbar
plot_distance_matrix(city_names=True, metric='correlation', case='b', ax=axes[1], colorbar=colorbar)
plot_distance_matrix(city_names=True, metric='correlation', case='c', ax=axes[2], colorbar=colorbar)

axes[0].text(-0.2, 0.97, "A", ha='center', va='center', transform=axes[0].transAxes, fontsize=20)
axes[1].text(-0.2, 0.97, "B", ha='center', va='center', transform=axes[1].transAxes, fontsize=20)
axes[2].text(-0.2, 0.97, "C", ha='center', va='center', transform=axes[2].transAxes, fontsize=20)

plt.subplots_adjust(right=0.940, top=0.985, bottom=0.20,left=0.08,wspace=0.480, hspace=0.30)

cbar_ax = fig.add_axes([0.943, 0.205, 0.015, 0.781])  # [left, bottom, width, height]
cbar = plt.colorbar(heatmap.collections[0], cax=cbar_ax, orientation='vertical', label='Correlation Distance')
cbar.ax.tick_params(labelsize=14)
cbar.ax.yaxis.label.set_size(14)

plt.savefig('figures/dis_matrix_'+scenario+'_' +'a_b_c'+'.png')



