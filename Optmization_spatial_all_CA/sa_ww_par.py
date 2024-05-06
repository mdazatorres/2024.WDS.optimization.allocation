import random
#import pulp
import pandas as pd
from geopy.distance import geodesic
#import random
import math


class SimulatedAnnealingOptimization:
    def __init__(self, scenario, dissim_term):
        #self.scenario = 'cdph_w_scan'
        #self.dissim_term = True
        #self.scenario = 'cdph_scan' # scenario1
        #self.dissim_term = False

        self.scenario = scenario # scenario1
        self.dissim_term = dissim_term


        #self.Data_path = {'cdph_w_scan':"../output/data_optimization_cdph_w_scan.csv",'cdph_scan': "../output/data_optimization_cdph_scan.csv"}
        #self.data_path = self.Data_path[self.scenario]
        self.data_path = "../output/data_optimization_cdph_scan_eur.csv"
        self.data_plants_read = pd.read_csv(self.data_path)

        self.data_plants = pd.DataFrame(data=self.data_plants_read)
        #self.data_plants = self.data_plants[self.sel_wwtps['Plant']]
        if self.dissim_term:
            self.sel_wwtps = pd.read_csv('../output/wwtps_CA_sel_'+self.scenario+'.csv')
            self.data_plants = self.data_plants[self.data_plants.Plant.isin(self.sel_wwtps['Plant'])].reset_index(drop=True)

        self.plants = self.data_plants['Plant'].tolist()
        self.n = len(self.plants)
        self.cities = dict(zip(self.plants, self.data_plants['City'].tolist()))
        self.counties = dict(zip(self.plants, self.data_plants['County'].tolist()))
        self.population_density = dict(zip(self.plants, self.data_plants['density'].tolist()))
        self.sum_population_density = sum([self.population_density[plant] for plant in self.plants])
        self.svi = dict(zip(self.plants, self.data_plants['svi'].tolist()))
        self.lat = dict(zip(self.plants, self.data_plants['lat'].tolist()))
        self.lng = dict(zip(self.plants, self.data_plants['lng'].tolist()))

        self.population_served = dict(zip(self.plants, self.data_plants['Population_Served'].tolist()))
        self.sum_population_served = sum([self.population_served[plant] for plant in self.plants])
        self.population_city = dict(zip(self.plants, self.data_plants['pop2023'].tolist()))
        self.total_distance,  self.distance_matrix = self.dist_matrix()
        self.sum_svi = sum([self.svi[plant] for plant in self.plants])
        if self.dissim_term:
            self.sim_matrix = pd.read_csv('../output/dis_matrix_'+self.scenario+'.csv',index_col=0)
            self.sum_sim_matrix = 0.5 * self.sim_matrix.values.sum()
            self.prod_list_sim_matrix = [(index, column) for index in self.sim_matrix.index for column in self.sim_matrix.columns]
        self.cache = {}

    def dist_matrix(self):
        distance_matrix = {(plant1, plant2): geodesic((self.lat[plant1], self.lng[plant1]), (self.lat[plant2], self.lng[plant2])).miles
                           for plant1 in self.plants for plant2 in self.plants if plant1 != plant2}
        total_distance = 0.5 * sum([distance_matrix[i, j] for i, j in distance_matrix.keys() if i != j])
        return total_distance, distance_matrix

    def simulated_annealing_opt(self, minimum_population_served, minimum_cities, minimum_counties,
                                 max_plants_in_solution, weight_population_served, weight_population_density,
                                 weight_proximity, weight_svd, weight_sim):
        # All the weights grouped in one list
        weights = [weight_population_served, weight_population_density, weight_proximity, weight_svd,weight_sim]

        # Number of total minimum_population_served
        min_population_served = minimum_population_served * self.sum_population_served

        # Create a binary variable for each plant
        c_list = [0] * self.n
        X = dict(zip(self.plants, c_list))
        X_new = dict(zip(self.plants, c_list))

        # Create binary variables for the presence of edges (connections) between each pair of plants
        edges = {(plant1, plant2): 0 for plant1 in self.plants for plant2 in self.plants if plant1 != plant2}

        # Initialize indices of plants in the solution
        solution_indices = self.random_list(max_plants_in_solution, self.n)

        # Set up the variable X
        for i in range(max_plants_in_solution):
            X[self.plants[solution_indices[i]]] = 1
            X_new[self.plants[solution_indices[i]]] = 1

        # Set up the edges
        for i in range(max_plants_in_solution):
            for j in range(i):
                edges[(self.plants[solution_indices[i]], self.plants[solution_indices[j]])] = 1
                edges[(self.plants[solution_indices[j]], self.plants[solution_indices[i]])] = 1

        # Initial temperature and energy and minimum temperature
        temperature = 100000
        minimum_temperature = 0.000001

        # Objective function: maximize the total population density while penalizing proximity
        Energy = []
        energy = self.energy_function(weights, X, edges)
        Energy.append(energy)
        # Simulated annealing
        new_edges = edges.copy()

        while temperature > minimum_temperature:
            accepted_solutions = 0
            itr = 0
            while accepted_solutions < 10 * self.n and itr < 100 * self.n:
                new_solution_indices = self.new_permutation(solution_indices, self.n)
                a, b = self.find_first_difference(solution_indices, new_solution_indices)
                X_new[self.plants[a]] = 0
                X_new[self.plants[b]] = 1
                for i in range(max_plants_in_solution):
                    c = new_solution_indices[i]
                    if c != a:
                        new_edges[(self.plants[c], self.plants[a])] = 0
                        new_edges[(self.plants[a], self.plants[c])] = 0
                    if c != b:
                        new_edges[(self.plants[c], self.plants[b])] = 1
                        new_edges[(self.plants[b], self.plants[c])] = 1
                new_energy = self.energy_function(weights, X_new, new_edges)
                delta_E = energy - new_energy
                if self.acceptance_function(delta_E, temperature):
                    [new_population_served, new_cities_count, new_counties_count] = self.count(X_new)
                    if min_population_served <= new_population_served and minimum_cities <= new_cities_count \
                            and minimum_counties <= new_counties_count:
                        X = X_new.copy()
                        energy = new_energy
                        accepted_solutions += 1
                        solution_indices = new_solution_indices.copy()
                        edges = new_edges.copy()
                        Energy.append(energy)
                X_new = X.copy()
                new_edges = edges.copy()
                itr += 1
            temperature *= 0.1

        [total_population_served, final_cities_count, final_counties_count] = self.count(X)


        # In case a constraint is not satisfied 
        if total_population_served < min_population_served:
            print(f"It was not possible to have more than {minimum_population_served*100}% of population_served")
        if  final_cities_count < minimum_cities:
            print(f"It was not possible to have more than {minimum_cities} cities")
        if final_counties_count < minimum_counties:
            print(f"It was not possible to have more than {minimum_counties} counties")


        selected_plants = [plant for plant in self.plants if X[plant] == 1] 
        selected_lat = [self.lat[plant] for plant in selected_plants]
        selected_lng = [self.lng[plant] for plant in selected_plants]
        selected_pop_served = [self.population_served[plant] for plant in selected_plants]
        selected_pop_city = [self.population_city[plant] for plant in selected_plants]
        selected_cities = [self.cities[plant] for plant in selected_plants]
        selected_svi = [self.svi[plant] for plant in selected_plants]
        dict_plants = {'Name': selected_plants, 'City': selected_cities, 'Latitude': selected_lat,
                       'Longitude': selected_lng, 'Population Served': selected_pop_served,
                       'Population City': selected_pop_city, 'SVI': selected_svi}
        return [dict_plants, total_population_served, energy, Energy]
        #return [dict_plants, total_population_served, energy]

    def energy_function(self, weights, X, edges):
        '''
        Minimizing function
        '''
        # Check if the result is already cached
        cache_key = tuple(sorted(X.items()))
        if cache_key in self.cache:
            return self.cache[cache_key]

        #0.5 * (weights[2] / self.total_distance)
        energy = ((weights[0] / self.sum_population_served) * sum([self.population_served[plant] * X[plant] for plant in self.plants])
                + ((weights[1]) / self.sum_population_density) * sum(
                    [self.population_density[plant] * X[plant] for plant in self.plants])
                + 0.5 * (weights[2] / self.total_distance) * sum(
                    [self.distance_matrix[i, j] * edges[i, j] for i, j in self.distance_matrix.keys()])
                + (weights[3] / self.sum_svi) * sum([self.svi[plant] * X[plant] for plant in self.plants]))
        if self.dissim_term:
            energy = energy + 0.5 * (weights[4] / self.sum_sim_matrix) * sum([self.sim_matrix.loc[i, j] * edges[i, j] for i, j in self.distance_matrix.keys()])

        self.cache[cache_key] = energy
        return energy

    def count(self, X):
        '''
        This function returns a list cotaining the population served, the number of cities
        covered, and the number of counties covered by a given solution X
        '''
        new_population_served = sum([self.population_served[plant] * X[plant] for plant in self.plants])
        cities_count = len(set([self.cities[plant] for plant in self.plants if X[plant] == 1]))
        counties_count = len(set([self.cities[plant] for plant in self.plants if X[plant] == 1]))
        return [new_population_served, cities_count, counties_count]


    def new_permutation(self, lst, n):
        '''
        Given a list lst of k < n different numbers from 0 to n-1, it replaces a random number from this list
        by a new random number (between 0 and n-1) that was not in the list previously
        '''
        k = len(lst)
        index_to_change = random.randint(0, k - 1)
        new_number = lst[index_to_change]

        while new_number in lst:
            new_number = random.randint(0, n - 1)

        new_list = lst.copy()
        new_list[index_to_change] = new_number
        return new_list


    def random_list(self, k, n):
        '''
        It returns a list of k random different numbers from 0 to n-1
        '''
        if k >= n:
            raise ValueError("k must be less than n.")

        lst = random.sample(range(n), k)
        return lst


    def find_first_difference(self,lst1, lst2):
        '''
        It finds the first integers that are in a different lst1 and lst2
        '''
        set_lst1 = set(lst1)
        set_lst2 = set(lst2)

        for a in set_lst1:
            if a not in set_lst2:
                for b in set_lst2:
                    if b not in set_lst1:
                        return a, b
        return None, None

    '''
    Acceptance function for Simulated Annealing
    '''

    def acceptance_function(self, delta_E, temperature):
        if delta_E < 0:
            return True
        else:
            r = random.uniform(0, 1)
            if r < math.exp((-delta_E) / temperature):
                return True
            else:
                return False
#sa=SimulatedAnnealingOptimization()


"""
array(['AmerValley_Quincy',
       'City of Paso Robles Wastewater Treatment Plant',
       'Esparto Wastewater Treatment Facility', 'LASAN_Hyp',
       'SFPUC_Ocean', 'SFPUC_SE'], dtype=object)
"""
