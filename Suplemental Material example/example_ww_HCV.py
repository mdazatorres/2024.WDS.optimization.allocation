import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import numpy as np
import pywt
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import pywt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns
from scipy.spatial.distance import correlation


data_ww = pd.read_csv('data_ww_cases_full.csv')
data_county= pd.read_csv('data_cases_county.csv')
data_county= data_county[data_county.County=='Yolo'].reset_index()
data_county['cases_av'] = data_county['cases'].rolling(window=7, center=False).mean()
data_county['pos_rate_av'] = data_county['pos_rate'].rolling(window=7, center=False).mean()
data_county= data_county[(data_county.date>='2022-05-07') & (data_county.date<='2022-09-28')].reset_index()

cols = ['SampleDate', 'positives', 'pos_rate','N1_Concentration_Merged', 'N2_Concentration_Merged',
       'N/PMMoV', 'PV_Concentration_Merged', 'City', 'NormalizedConc_crude', 'NormalizedConc', 'Type', 'N_gene',
       'S_gene']
data_ww = data_ww[cols]
data_ww['SampleDate'] = pd.to_datetime(data_ww['SampleDate'])
#data_ww = data_ww[~data_ww.isin(['Modesto (sludge, Eurofins)', 'Merced (sludge, Eurofins)', 'Emerald Lift Station','UCDavis', 'UCDavis (sludge)', 'Merced', 'Modesto'])]
data_ww = data_ww[~data_ww.isin(['Modesto (sludge, Eurofins)', 'Merced (sludge, Eurofins)', 'Emerald Lift Station', 'UCDavis', 'UCDavis (sludge)', 'Merced', 'Modesto'])]

city_influent=['Davis',  'Winters', 'Woodland', 'Esparto', 'Los Banos', 'Turlock']
city_sludge= ['Davis (sludge)', 'Winters (solids)', 'Woodland (solids)', 'Esparto (solids)', 'Los Banos (solids)', 'Turlock (solids)']


city_eur = ['City of Davis (sludge, Eurofins)']
city1 = 'Davis (sludge)'
city2 = 'City of Davis (sludge, Eurofins)'
cities= ['Davis (sludge)', 'City of Davis (sludge, Eurofins)']
d_Davis1 = data_ww[data_ww.City=='Davis (sludge)']
d_Davis2 = data_ww[data_ww.City=='City of Davis (sludge, Eurofins)']
d_Davis3 = data_ww[data_ww.City=='Davis']


#sns.lineplot(data=d_Davis1 , x='SampleDate', y='NormalizedConc_crude')
#sns.lineplot(data=d_Davis2,  x='SampleDate', y='NormalizedConc_crude')
#sns.lineplot(data=d_Davis3,  x='SampleDate', y='NormalizedConc_crude')

def trim_fun(x):
       x = x.dropna()
       if len(x) <= 2:
              return np.mean(x)
       else:
              x1 = x.sort_values().ravel()
              return np.mean(x1[1:-1])

def plot(city1, city2):
       df1= data_ww[data_ww.City == city1]
       df2 = data_ww[data_ww.City == city2]
       sns.lineplot(data=df1, x='SampleDate', y='NormalizedConc_crude')
       sns.lineplot(data=df2, x='SampleDate', y='NormalizedConc_crude')

       plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
       plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))  # Format to show year and month

# Rotate the x-axis labels for better readability
       plt.xticks(rotation=45)
'''
i=1
#plot(city1=city_influent[i], city2=city_sludge[i])
#sns.lineplot(data=data_ww, x='SampleDate', y='NormalizedConc_crude', hue='City')

'''
#

def wavelet_transform(data, all=True):
       wavelet = 'db4' #'db2'#'rbio3.7' #'db1'
       level = 4
       coeffs = pywt.wavedec(data, wavelet, level=level)
       if all:
              return np.concatenate(coeffs[:-4])  # Flatten and concatenate the coefficients
       else:
              return coeffs[0]


def wavelet_transform_data(data):
       transformed_data = data.groupby('City')['SC2_N_norm_PMMoV_av10'].apply(wavelet_transform).reset_index()
       plant_series_dict = {}

       # Iterate through unique county names
       for plant in data['City'].unique():
              transf_plant_data = transformed_data[(transformed_data['City'] == plant)]
              concentration_series = transf_plant_data['SC2_N_norm_PMMoV_av10'].values[0]
              plant_series_dict[plant] = concentration_series
       plant_series_list = list(plant_series_dict.values())

       return plant_series_dict, plant_series_list

def data_process(sl_inf, Davis, data_ww):
       if sl_inf:
              data_ww= data_ww[data_ww.City!='City of Davis (sludge, Eurofins)']
       if Davis:
              data_ww = data_ww[data_ww.City.isin(['City of Davis (sludge, Eurofins)', 'Davis (sludge)'])]

       min_dates_by_plant = data_ww.groupby('City')['SampleDate'].min()
       max_dates_by_plant = data_ww.groupby('City')['SampleDate'].max()

       min_date = min_dates_by_plant.max()
       max_date = max_dates_by_plant.min()

       df=data_ww[(data_ww['SampleDate']>min_date) & (data_ww['SampleDate']<max_date)]

       agg_funcs = {'N1_Concentration_Merged': 'mean', 'N2_Concentration_Merged':'mean',
              'N/PMMoV':'mean', 'PV_Concentration_Merged':'mean', 'NormalizedConc_crude':'mean',
              'NormalizedConc':'mean', 'Type':'first', 'N_gene':'mean', 'S_gene':'mean','positives':'mean',
                    'pos_rate':'mean'}

       df = df.groupby(['City', 'SampleDate']).agg(agg_funcs).reset_index()
       df['SC2_N_norm_PMMoV_av10'] = df.groupby('City')['NormalizedConc_crude'].rolling(window=10, min_periods=1,  center=False).apply(lambda x: trim_fun(x)).interpolate().reset_index(level=0, drop=True)
       df = df[['City', 'SampleDate','positives', 'NormalizedConc_crude', 'pos_rate','SC2_N_norm_PMMoV_av10']]
       return df

d1 = data_process(sl_inf=True, Davis=False, data_ww=data_ww)
d2 = data_process(sl_inf=False, Davis= True, data_ww=data_ww)
d3 = data_process(sl_inf=False, Davis= False, data_ww=data_ww)

#sns.lineplot(data=d2, x='SampleDate', y='SC2_N_norm_PMMoV_av10', hue='City')

plant_series_dict1, plant_series_list1 = wavelet_transform_data(d3)
plant_series_list = plant_series_list1
plant_series_dict = plant_series_dict1
city_names = d3.City.unique()

def plot_distance_matrix_den(metric='correlation'):
    pairwise_distances_corr = pdist(plant_series_list, metric=metric)
    distance_matrix = squareform(pairwise_distances_corr)

    #plant_names_p = np.vectorize(procces_plant_name)(plants)

    linkage_matrix = linkage(distance_matrix, method='average')
    order = dendrogram(linkage_matrix, no_plot=True)['leaves']
    distance_matrix_sorted = distance_matrix[order, :]
    distance_matrix_sorted = distance_matrix_sorted[:, order]
    #plants_sorted = np.array(plants)[order]
    plants_sorted = np.array(city_names)[order]

    # Create a DataFrame for the heatmap
    heatmap_df = pd.DataFrame(distance_matrix_sorted, index=plants_sorted, columns=plants_sorted)

    # Set up the plot
    sns.set(font_scale=1.5)
    plt.rcParams.update({'font.size': 10})

    # Use sns.clustermap to create the heatmap with dendrogram
    cl=sns.clustermap(heatmap_df, method='average', metric=metric, cmap='viridis',
                   annot=False, figsize=(12, 12), cbar_kws={'label': 'Distance'},
                   row_cluster=True, col_cluster=True, cbar_pos=(0.06, 0.82, 0.04, 0.15))
    #plt.setp(cl.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    #plt.tight_layout()
    #if save:
    #    plt.savefig('figures/dis_matrix_den_' +scenario+'.png')

    plt.show()

#plot_distance_matrix_den()
#df_1= 'Davis (sludge)'



def plot_all_coeff(data, cities):
    #'coif3' # haar 4 #rbio3.7 4
    #i=105
    wavelet = 'db4' #'db2'
    coeff_labels = ['cA', r'$cD_4$', r'$cD_3$', r'$cD_2$', r'$cD_1$']
    #print(wavelet)
    level = 4
    fig, ax = plt.subplots(3,2, figsize=(10, 8))
    ax_flat = ax.flatten()
    Coeff=[]; ts=[]
    Coeff_l1=[]; Coeff_l2=[]; Coeff_l3=[]
    labels={'City of Davis (sludge, Eurofins)':'Lab 2', 'Davis (sludge)':'Lab 1'}
    #ax_flat[0].plot(coeffs1[i], label=city)
    for city in cities:
        df = data[data.City == city][5:]
        coeffs1 = pywt.wavedec(df['SC2_N_norm_PMMoV_av10'], wavelet, level=level)
        #coeffs_all= coeffs1.copy()

        Coeff.append(np.concatenate(coeffs1))
        Coeff_l1.append(np.concatenate(coeffs1[:-1]))
        Coeff_l2.append(np.concatenate(coeffs1[:-2]))
        Coeff_l3.append(np.concatenate(coeffs1[:-3]))
        ts.append(df['SC2_N_norm_PMMoV_av10'])
        ax_flat[0].plot(df['SC2_N_norm_PMMoV_av10'].values, label=labels[city])
        # Plot coefficients
        for i in range(len(coeffs1)):
            ax_flat[i+1].plot(coeffs1[i], label=city)
    #ax_flat[0].plot(ts[0], label=city)
    #ax_flat[0].plot(ts[1], label=city)
    ax_flat[0].legend()
    for i, label in enumerate(coeff_labels):
            i = i+1
            ax = ax_flat[i]
            ax.set_title(f'{label} Coefficients')
    plt.tight_layout()
    plt.savefig('figures/Davis_SCAN_eurofins.png')

    print("wa_1", correlation(Coeff_l1[0], Coeff_l1[1]))
    print("wa_2", correlation(Coeff_l2[0], Coeff_l2[1]))
    #print("wa_3", correlation(Coeff_l3[0], Coeff_l3[1]))
    print("wa", correlation(Coeff[0], Coeff[1]))
    print("ts", correlation(ts[0], ts[1]))


#plot_all_coeff(data=d2, cities=['City of Davis (sludge, Eurofins)', 'Davis (sludge)'])

