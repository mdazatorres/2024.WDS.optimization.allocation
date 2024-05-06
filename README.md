## Optimizing Spatial Distribution of Wastewater-Based Disease Surveillance to Advance Health Equity




Files:
A) Create_data
We create all the data set to use in the project. In the folder output will be saved all the data to use. We create data sets for our problem  considering the three available data sets: SCAN, CDPH and CDC


create_data_ww_cpdh.py
In this code, we perform a preprocessing of wastewater datasets sourced from the CDPH. We incorporate geographical
 information such as latitude, longitude, city, and county, which is derived from the zip codes associated with
  each wastewater treatment plant. This data set is ready to use in the cluster problem and also in the optimization take just the needed information.

This code creates a dataset for use in the allocation optimization problem by merging two data sets and using a dictionary:

This date set is based on CPDH wastewater location.
1) A dataset containing information about cities monitoring wastewater, their corresponding population served, and other variables (from CPDH).
2) A dataset with information about cities in California, including longitude, latitude, area, density, and other variables needed for the optimization problem.

3) A dictionary mapping city names to counties
4) A data set with  the Social Vulnerability Index to county level
Input:
 CA_Counties_TIGER2016.shp (shape file for county polygons)
us-cities-table-for-california.csv Load data related to cities in California, including population, longitude, latitude, density, etc.

create_data_ww_cpdh_not_scan.py
This data set is based on create_data_ww_cpdh_scan.py the difference is here we selected the plants that are scan to deleted in the data set. We can obtained all the previous data sets here.

create_data_ww_cdc.py
In this code, we preprocess the WW dataset provided by CDC. 
Note: this dataset does not contain precise information about the location of the WWTPs.
Read.
county_fips.csv
cdc-nwss-restricted-data-set-wastewater.csv (data from CDC)


B) Optimization_spatial

There is two folders with similar files, Optimization_spatial_all_CA is just a version in class of all of the following files.

example.py
This code solves an optimization allocation problem using synthetic data for cities. The process involves generating data for cities, including population density and a distance matrix. Then, it formulates the problem as an LP model to find the optimal locations for wastewater treatment plants. The LP model aims to maximize the population density coverage while penalizing the proximity of selected cities.
1) It generates synthetic data for cities, including population density and distance matrix, and then
2) formulates the problem as a Linear Programming (LP) model. The objective is to maximize the total population
density covered while penalizing the proximity of the selected cities. The LP problem is subject to constraints
on the minimum population coverage and a maximum number of cities to be included in the solution.

sa_ww.py
Optimization allocation problem using Simulated Annealing

testing.py
This code is used to run sa_ww.py for different parameters. It can be used to vary the weights of the minimizing function and find the solution of k plants for each k=1, ….50, and it saves .csv files containing the solution of cities given by the algorithm. They are saved in Results/csv_files/Solution_plants/ with the format solution_{k}_{w1}_{w2}_{w3}_{w4}.csv, where k is the number of plants, w1, w2, w3 are the pop served weight, pop density weight,prox weight, svi weight, respectively. 
It also saves a .csv file containing the population served and energy of each k in Results/csv_files/Population_and_energy/ as {w1}_{w2}_{w3}.csv, where w1, w2, w3 are the same as before.

map_all_cities.py
This code plots a map of all cities in California.It reads city data from a CSV file containing the latitude,
longitude, and name of each city. The code  create an interactive map with markers representing each city.



C. Cluster and other approach
clusters_euc_dtw.py
In this code we create the clusters for ww signal across california to wwtps level to be used in the optimization problem. We use the locations provide by cdph. Here we use the euclidean distance to compute clusters or the fast DTW algorithm. The solutions are plotted in the California map together with ski and the population served.

clusters_wavlets.py
Compute clusters using wavlets, we save a similarity matrix to be used in the optimization problem

map_ww_cdph.py
Plot the map by California with the cdph data.

map_ww_scan.py
Plot the map by California with the SCAN data.

D. Dissmilarity matrix
Here, we create the dissimilarity matrix that is used in the optimization problem.

dissimilarity_matrix.py
We create the dissimilarity matrix using the correlation distance from the wavelet coefficients of the wastewater data.

main_cluster.py
We set all the variables for the problem and the computation of the wavelets coefficients in this main file

map_ww_cdph_select_plants
Map of the California plants used in the solution.




D. Data 

California WWTP List.xlsx
A list with all wwtps across California, Colleen shared this data set with us

svi_interactive_map.csv 
This contains sci to county level. This dat set was download from 

Geodata
CA_Counties_TIGER2016.*
These are the polygons for the counties in california for doing the map.

cal_cities_lat_long.csv
Contains latitude and longitude from the cities across california

California_Counties_Incorporated_Cities.csv
This data sets contains a column with the cities in california and the name of their respective county

Econda.xlsx and demographic data are the data that Steven sent us about some demographics variables in California a county level

us-cities-table-for-california.csv
Population and density for california cities

