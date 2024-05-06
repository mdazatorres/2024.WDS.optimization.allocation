import pandas as pd
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
scenario =mc.scenario
data = mc.data_sl
plants = data.Plant.unique()


def city_plant_dict_create():
    city_counts = {}
    city_plant_dict = data[['Plant', 'City']].drop_duplicates(subset='Plant')
    city_plant_dict.index = city_plant_dict['Plant']
    city_plant_dict = city_plant_dict['City'].to_dict()
    for key, city in city_plant_dict.items():
        if city in city_counts:
            city_counts[city] += 1
            new_city_name = f'{city}_{city_counts[city]}'
            if city == 'Tracy':
                new_city_name = 'Tracy (MH)'
        else:
            city_counts[city] = 1
            new_city_name = city

        city_plant_dict[key] = new_city_name
    return city_plant_dict



def plot_distance_matrix(save, city_names, metric):
    pairwise_distances_corr = pdist(mc.plant_series_list, metric=metric)
    distance_matrix = squareform(pairwise_distances_corr)

    plant_names_p = np.vectorize(mc.procces_plant_name)(plants)

    linkage_matrix = linkage(distance_matrix, method='average')
    order = dendrogram(linkage_matrix, no_plot=True)['leaves']
    distance_matrix_sorted = distance_matrix[order, :]
    distance_matrix_sorted = distance_matrix_sorted[:, order]
    #plants_sorted = np.array(plants)[order]
    plants_sorted = np.array(plant_names_p)[order]
    if city_names:
        city_plant_dict = city_plant_dict_create()
        city_names = [city_plant_dict[plant] for plant in plants]
        plants_sorted = np.array(city_names)[order]

    sns.set(font_scale=1.5)
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(figsize=(12, 10))

    sns.heatmap(distance_matrix_sorted, annot=False, cmap='viridis', xticklabels=plants_sorted,
                yticklabels=plants_sorted, fmt=".2f", cbar_kws={'label': 'Correlation Distance'})
    fig.tight_layout()
    if save:
        dis_matrix_df = pd.DataFrame(distance_matrix, index=plants, columns=plants)
        #dis_matrix_df.to_csv(save_path + "dis_matrix.csv", index=True)
        dis_matrix_df.to_csv(save_path + "dis_matrix_"+scenario+'.csv', index=True)
        plt.savefig('figures/dis_matrix_'+scenario+'.png')
    return dis_matrix_df

#for i in range(21):
#    print(len(mc.plant_series_list[i]), plants[i])


ds_matrix=plot_distance_matrix(save=True, city_names=True, metric ='correlation')

def plot_distance_matrix_den(save, city_names, metric ):
    pairwise_distances_corr = pdist(mc.plant_series_list, metric=metric)
    distance_matrix = squareform(pairwise_distances_corr)

    plant_names_p = np.vectorize(mc.procces_plant_name)(plants)

    linkage_matrix = linkage(distance_matrix, method='average')
    order = dendrogram(linkage_matrix, no_plot=True)['leaves']
    distance_matrix_sorted = distance_matrix[order, :]
    distance_matrix_sorted = distance_matrix_sorted[:, order]
    #plants_sorted = np.array(plants)[order]
    plants_sorted = np.array(plant_names_p)[order]

    # If you have city names, create a dictionary to map plants to cities
    if city_names:
        city_plant_dict = city_plant_dict_create()
        city_names = [city_plant_dict[plant] for plant in plants]
        plants_sorted = np.array(city_names)[order]


    # Create a DataFrame for the heatmap
    heatmap_df = pd.DataFrame(distance_matrix_sorted, index=plants_sorted, columns=plants_sorted)

    # Set up the plot
    sns.set(font_scale=1.5)
    plt.rcParams.update({'font.size': 10})

    # Use sns.clustermap to create the heatmap with dendrogram
    cl=sns.clustermap(heatmap_df, method='average', metric=metric, cmap='viridis',
                   annot=False, figsize=(12, 12), cbar_kws={'label': 'Distance'},
                   row_cluster=True, col_cluster=True,cbar_pos=(0.06, 0.82, 0.04, 0.15))
    #plt.setp(cl.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    #plt.tight_layout()
    if save:
        plt.savefig('figures/dis_matrix_den_' +scenario+'.eps', dpi=900)

    plt.show()

#plot_distance_matrix_den(save=True, city_names=True, metric ='cosine')
plot_distance_matrix_den(save=True, city_names=True, metric ='correlation')
#plot_distance_matrix(save=True, city_names=True, metric ='correlation')
####














