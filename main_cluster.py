from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import KMeans
from geopy.geocoders import Nominatim
import pandas as pd
import numpy as np
from shapely.geometry import Point
import geopandas as gpd
import pywt
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import ward, fcluster
from scipy.cluster.hierarchy import dendrogram

from scipy.ndimage import gaussian_filter
import scipy.cluster.hierarchy as sch



class main_data_cluster:
    def __init__(self, scenario, dissim_term):
        #self.scenario = 'cdph_w_scan' # cdph plants except those provide by scan
        #self.dissim_term = True
        #self.scenario = 'cdph_scan' # cdph with scan data
        #self.dissim_term = False    # If we can add the dissimilarity term, this when we have wastewater data available

        self.scenario = scenario
        self.dissim_term = dissim_term

        self.save_path = "../output/"
        self.Data = {'all': pd.read_csv(self.save_path + "data_ww_geodata_CA_cdph_scan_eur.csv"),
                     'diss': pd.read_csv(self.save_path + "data_ww_geodata_CA_cdph_scan_eur.csv")}

        self.data = self.Data[self.scenario]

        #"data_ww_geodata_CA_cdph_w_scan.csv"
        self.shapefile_path = '../data/geodata/CA_Counties_TIGER2016.shp'
        self.gdf = gpd.read_file(self.shapefile_path)
        self.data['Date'] = pd.to_datetime(self.data['Date'])
        #self.data = self.data[self.data.Date >= '2022-05-01']

        self.data = self.smothed_data(trimmed=True)
        self.Min_date = {'all':'2023-02-01', 'diss': '2022-03-21'}
        self.Max_date = {'all': '2024-02-01', 'diss': '2023-05-21'}
        self.min_date = self.Min_date[self.scenario]
        self.max_date = self.Max_date[self.scenario]
        if self.dissim_term:
            self.data_sl = self.data_select(self.data)  # data same length
        else:
            self.data_sl = self.data

        self.plant_series_dict, self.plant_series_list = self.wavelet_transform_data(self.data_sl)
       # self.data_full = self.data[~self.data.Plant.isin(self.Plants_finished_before)]  # Full data just eliminated plants to stop to collect ww

    def trim_fun(self, x):
        x = x.dropna()
        if len(x)<=2:
            return np.mean(x)
        else:
            x1 = x.sort_values().ravel()
            return np.mean(x1[1:-1])

    def smothed_data(self, trimmed=True):
        #'SC2_Norm'
        col= 'SC2_Norm'
        if trimmed:
            #self.data['SC2_N_norm_PMMoV_av10'] = self.data.groupby('Plant')['SC2_N_gc_g_dry_weight'].rolling(window=10, min_periods=1, center=False).apply(lambda x: self.trim_fun(x)).interpolate().reset_index(level=0, drop=True)
            self.data['SC2_N_norm_PMMoV_av10'] = self.data.groupby('Plant')['SC2_Norm'].rolling(window=10, min_periods=1, center=False).apply(lambda x: self.trim_fun(x)).interpolate().reset_index(level=0, drop=True)

            self.data['SC2_N_norm_PMMoV_av10'] = self.data.groupby('Plant')['SC2_N_norm_PMMoV_av10'].fillna(0)
        else:
            #self.data['SC2_N_norm_PMMoV_av10'] = self.data.groupby('Plant')['SC2_N_gc_g_dry_weight'].rolling(window=10, center=False).mean().reset_index(level=0, drop=True)
            self.data['SC2_N_norm_PMMoV_av10'] = self.data.groupby('Plant')['SC2_Norm'].rolling(window=10, center=False).mean().reset_index(level=0, drop=True)

            self.data['SC2_N_norm_PMMoV_av10'] = self.data.groupby('Plant')['SC2_N_norm_PMMoV_av10'].fillna(0)
        return self.data


    def wavelet_transform(self, data, all=True):
        wavelet = 'db4'
        level = 4
        coeffs = pywt.wavedec(data, wavelet, level=level)
        if all:
            return np.concatenate(coeffs[:-2])  # Flatten and concatenate the coefficients
        else:
            return coeffs[0]
    def wavelet_transform_data(self, data):
        transformed_data = data.groupby('Plant')['SC2_N_norm_PMMoV_av10'].apply(self.wavelet_transform).reset_index()
        plant_series_dict = {}

        # Iterate through unique county names
        for plant in data['Plant'].unique():
            transf_plant_data = transformed_data[(transformed_data['Plant'] == plant)]
            concentration_series = transf_plant_data['SC2_N_norm_PMMoV_av10'].values[0]
            plant_series_dict[plant] = concentration_series
        plant_series_list = list(plant_series_dict.values())

        return plant_series_dict, plant_series_list

    def exclude_plants(self, data):
        data = data[data.Plant != 'CALEXICO_CTY']  # This plant does no have data for PMMOV sars cov2 concentration
        data = data[data.Plant != 'Oxnard_WWTP']  # This plant does no have data for PMMOV sars cov2 concentration
        min_dates_by_plant = data.groupby('Plant')['Date'].min()
        max_dates_by_plant = data.groupby('Plant')['Date'].max()
        data_exc = data[['Date', 'Plant', 'Population_Served', 'County', 'City']]
        self.Plants_started_after = min_dates_by_plant[min_dates_by_plant > self.min_date].index  # Plants_started_after.append(pd.Index(['Oxnard_WWTP']))
        end_dates = data_exc.groupby('Plant')['Date'].max().reset_index()
        end_dates.rename(columns={'Date': 'End_Date'}, inplace=True)
        data_exc = data_exc[data_exc['Plant'].isin(self.Plants_started_after)].drop_duplicates(subset='Plant')
        data_exc = data_exc.merge(end_dates, on='Plant')
        data_exc.rename(columns={'Date': 'Init_Date'}, inplace=True)
        self.data_exc = data_exc[['Init_Date', 'End_Date', 'Plant', 'County', 'City', 'Population_Served']]
        self.latex_table = self.data_exc.to_latex(index=False)
        return min_dates_by_plant, max_dates_by_plant

    def data_select(self, data):
        # We kept data sets with data between min and max date in order to have a same number of size in our analysis
        min_dates_by_plant, max_dates_by_plant = self.exclude_plants(data)

        df = data[(data.Date >= self.min_date) & (data.Date <= self.max_date)]
        self.selected_plants = min_dates_by_plant[(min_dates_by_plant <= self.min_date) & (max_dates_by_plant >= self.max_date)].index
        if self.dissim_term:
            pd.DataFrame({'Plant': self.selected_plants}).to_csv(self.save_path + "wwtps_CA_sel_"+self.scenario+'.csv')
        filt_df = df[df['Plant'].isin(self.selected_plants)]
        #unique_dates_per_plant = filt_df.groupby('Plant')['Date'].nunique()
        return filt_df


    def plot_data(self, data):

        sns.set_style("darkgrid")

        # Get the unique 'Plant' values in the DataFrame
        unique_plants = data['Plant'].unique()
        num_subplots = len(unique_plants)

        # Create a 3x10 grid of subplots
        ncols = 3
        nrows = int(np.ceil(data['Plant'].unique().shape[0]/ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(18, 10))

        for idx, plant in enumerate(unique_plants):
            # Filter the DataFrame for the specific 'Plant'
            plant_data = data[data['Plant'] == plant]

            # Calculate the subplot row and column for the current 'Plant'
            col = idx // nrows
            row = idx % nrows

            # Plot the data for the current 'Plant' in the corresponding subplot
            ax = axes[row, col]
            sns.lineplot(data=plant_data, x='Date', y='SC2_N_norm_PMMoV_av10', ax=ax)
            ax.set_title(plant)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

        # Remove any empty subplots if there are fewer than 30 unique plants
        for i in range(num_subplots, 30):
            col = i // nrows
            row = i % nrows
            fig.delaxes(axes[row, col])

        # Add a common legend for all subplots
        #fig.legend(loc='upper right')

        # Set spacing between subplots
        plt.tight_layout()

        sns.set_style("darkgrid")

    def procces_plant_name(self, plant_name):
        words = plant_name.split()
        plants_ = ['West County Wastewater District',
                   'Central Marin Sanitation Agency', 'Central Marin Sanitation Agency - West Railroad'
                   ]
        if len(words) == 1:
            # If it's a single word, keep the full name
            return words[0]
        elif words[0] == "City" and words[1] == "of":
            # If it starts with "City of," separate and add the next word
            if words[2] == 'San' or words[2] == 'Santa' or words[2] == 'Paso':
                if plant_name=='City of Santa Cruz WTF - County Influent':
                    return words[2] + ' ' + words[3] + ' CI'
                else:
                    return words[2] + ' ' + words[3]
            else:
                return words[2]
        elif words[0] == "Los" or words[0] == "Las":
            return words[0] + ' ' + words[1]
        elif plant_name == 'Sewerage Agency of Southern Marin Wastewater Treatment Plant':
            return 'Southern Marin WTP'
        elif np.isin(plant_name, plants_):
            return words[0] + ' ' + words[1] + ' ' + words[-1]
        elif words[0][0] == '[':
            return words[0] + ' ' + words[1]
        elif plant_name == 'Sewer Authority Mid-Coastside':
            return 'SA Mid-Coastside'
        else:
            abbreviated = [word[0] for word in words[1:]]
            return words[0] + ' ' + "".join(abbreviated)


scenario = 'diss' # cdph plants except those provide by scan
dissim_term = True

mc=main_data_cluster(scenario, dissim_term)
mc.data_sl

#plt.plot(data1.Date, data1.SC2_Norm)
#plt.plot(data1.Date, data1.SC2_N_gc_g_dry_weight)
#plt.plot(data1.Date, data1.SC2_N_norm_PMMoV_av10)

#plt.show()
#mc.plot_data(mc.data_sl)