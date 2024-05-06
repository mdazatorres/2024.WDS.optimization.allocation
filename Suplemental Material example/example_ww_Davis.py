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

cities=['Davis (sludge)', 'City of Davis (sludge, Eurofins)']
#d_Davis1 = data_ww[data_ww.City=='Davis (sludge)']
#d_Davis2 = data_ww[data_ww.City=='City of Davis (sludge, Eurofins)']
#d_Davis3 = data_ww[data_ww.City=='Davis']


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

#plot(city1='Davis (sludge)', city2='City of Davis (sludge, Eurofins)')



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

def data_process(data_ww):
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

data_p = data_process(data_ww=data_ww)


#sns.lineplot(data=data_p, x='SampleDate', y='SC2_N_norm_PMMoV_av10', hue='City')




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
    labels={ 'Davis (sludge)':'Lab 1','City of Davis (sludge, Eurofins)':'Lab 2'}
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
        ax_flat[0].set_title('Original Scale')
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
    print("wa_3", correlation(Coeff_l3[0], Coeff_l3[1]))
    print("wa", correlation(Coeff[0], Coeff[1]))
    print("ts", correlation(ts[0], ts[1]))


plot_all_coeff(data=data_p, cities=cities)



def reconstructed_signal(df, wavelet='db4', level=4):

    coeffs = pywt.wavedec(df['SC2_N_norm_PMMoV_av10'], wavelet, level=level)
    coeffs_l1 = coeffs.copy()
    coeffs_l1[-1] = np.zeros(len(coeffs[-1]))
    coeffs_l2 =  coeffs_l1.copy()
    coeffs_l2[-2] = np.zeros(len(coeffs[-2]))
    coeffs_l3 = coeffs_l2.copy()
    coeffs_l3[-3] = np.zeros(len(coeffs[-3]))
    reconstructed_signal1 = pywt.waverec(coeffs_l1, 'db4')
    reconstructed_signal2 = pywt.waverec(coeffs_l2, 'db4')
    reconstructed_signal3 = pywt.waverec(coeffs_l3, 'db4')

    data_county_ = data_county[5:]

    correlation0 = np.corrcoef(df['SC2_N_norm_PMMoV_av10'], data_county_['cases_av'])[0, 1]
    correlation1 = np.corrcoef(reconstructed_signal1, data_county_['cases_av'])[0, 1]
    correlation2 = np.corrcoef(reconstructed_signal2, data_county_['cases_av'])[0, 1]
    correlation3 = np.corrcoef(reconstructed_signal3, data_county_['cases_av'])[0, 1]

    print(correlation0, correlation1, correlation2, correlation3)

    return reconstructed_signal1, reconstructed_signal2, reconstructed_signal3


def plot_reconstructed_signal(data_p):
    labels = {'Davis (sludge)': 'Lab 1', 'City of Davis (sludge, Eurofins)': 'Lab 2'}
    matplotlib.rcParams.update({'font.size': 15})
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20, 5))

    for i, city in enumerate(cities):
        df = data_p[data_p.City == city][5:]
        r1, r2, r3 = reconstructed_signal(df)

        #plt.plot(df.SampleDate, df['SC2_N_norm_PMMoV_av10'], label='Original Signal', color='k', marker='o', markersize=3, linestyle='')
        axes[i].plot(df.SampleDate, df['SC2_N_norm_PMMoV_av10'], label='Original signal', color='grey', lw=2)
        axes[i].plot(df.SampleDate, r1, label='Rec. signal without $cD_1$', color='red', lw=1)
        axes[i].plot(df.SampleDate, r2, label='Rec. signal without  $cD_{1,2}$', color='blue', lw=1 )
        axes[i].plot(df.SampleDate, r3, label='Rec. signal without $cD_{1,2,3}$', color='green', lw=1 )
        axes[i].xaxis.set_major_locator(mdates.MonthLocator())
        axes[i].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))  # Format to show year and month
        axes[i].set_xlabel('Date', fontsize=16)
        axes[i].set_title(labels[df.City.iloc[0]], fontsize=20)
    axes[0].legend(loc='best', fontsize=15)
    axes[0].set_ylabel('N/PMMoV', fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.21)
    plt.savefig('figures/Comparing_signals.png')

#plot_reconstructed_signal(data_p)
def plot_cases_reconstructed_signal(data=data_p):
    labels = {'Davis (sludge)':'Lab 1', 'City of Davis (sludge, Eurofins)':'Lab 2'}
    matplotlib.rcParams.update({'font.size': 18})
    df1 = data[data.City == cities[0]][5:]
    df2 = data[data.City == cities[1]][5:]
    data_county_ = data_county[5:]
    r1, r2, r3 = reconstructed_signal(df1)
    r1_, r2_, r3_ = reconstructed_signal(df2)

    fig, ax = plt.subplots(num=1, figsize=(12, 5))
    ax_p = ax.twinx()
    p1, = ax.plot(df1.SampleDate, data_county_['cases_av'],  linewidth=1.5, color='black', label='Cases')
    p2, = ax_p.plot(df1.SampleDate, r2, linewidth=1.5, color='blue', label='Rec. Signal '+ labels[cities[0]])
    ax_p2 = ax.twinx()

    p3, = ax_p2.plot(df1.SampleDate, r2_, linewidth=1.5, color='red', label='Rec. Signal '+labels[cities[1]])
    ax_p2.spines['right'].set_position(('outward', 80))  # Adjust the distance from the first y-axis

    ax.legend(handles=[p1, p2, p3], fontsize=20)  # loc=5 woodland
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
    #ax_p.set_ylabel('N/PMMoV')
    ax_p2.set_ylabel('N/PMMoV', fontsize=20)
    ax.set_ylabel('Cases', fontsize=20)
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonthday=1, interval=1))
    ax.set_xlabel('Date', fontsize=20)
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right", rotation_mode="anchor")
    fig.tight_layout()
    plt.savefig('figures/Comparing_signals_cases.png')
plot_cases_reconstructed_signal(data=data_p)
# plt.plot(data_county_['cases_av'],reconstructed_signal2,'o')

#plt.plot(df1['SC2_N_norm_PMMoV_av10'].values)
#plt.plot(df2['SC2_N_norm_PMMoV_av10'].values)
    #plt.plot(reconstructed_signal3)

#City of Merced wastewater concentrations for 10-day average for N/PMMoV compared to weekly average of Merced County cases per 100 k population, weekly average county hospitalizations with 14-day lag, and weekly average county ICU patients with 9-day lag.

