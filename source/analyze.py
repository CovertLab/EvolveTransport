from __future__ import absolute_import, division, print_function

import os
import csv

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from source.kinetic_flux_model import KineticFluxModel



class Analyze(object):

	def __init__(self, config, fitness_function):

		self.out_dir = config.get('out_dir', None)
		self.parameter_out_dir = config.get('parameter_out_dir', None)
		self.saved_param_file = config.get('saved_param_file', None)
		self.parameter_analytics = config.get('parameter_analytics', False)
		self.mutation_variance = config.get('mutation_variance', 1.0)
		self.seed = config.get('seed', None)
		self.replicate_id = config.get('replicate_id', None)
		self.save_fitness_threshold = config.get('fitness_threshold')
		self.wcm_sim_data = config.get('wcm_sim_data', None)
		self.exchange_molecules = config.get('exchange_molecules', None)

		self.fitness_function = fitness_function
		self.kinetic_model = fitness_function.kinetic_model
		self.reactions = fitness_function.reactions
		self.parameter_indices = fitness_function.parameter_indices
		self.km_range = fitness_function.km_range
		self.kcat_range = fitness_function.kcat_range

	# Analyses -- TODO -- this should not be part of plots
	def parameters(self, final_population, final_fitness):

		# get indices of individualus with fitness higher than 0.95
		top_fit_indices = [index for index, value in enumerate(final_fitness.values()) if value >= self.save_fitness_threshold]
		top_fit_parameters = [final_population[index] for index in top_fit_indices]

		# TODO -- convert genotype to phenotype

		# save top parameters to 'best_parameters' file
		if not os.path.exists(self.parameter_out_dir):
			os.mkdir(self.parameter_out_dir)

		with open(os.path.join(self.parameter_out_dir, self.saved_param_file), 'a') as tsv_file:
			writer = csv.writer(tsv_file)
			for parameter in top_fit_parameters:
				writer.writerow(parameter)
		tsv_file.close()

		# Plot parameter space
		if self.parameter_analytics:
			# gather all saved parameter values
			with open(os.path.join(self.parameter_out_dir, self.saved_param_file), 'r') as tsv_file:
				read = csv.reader(tsv_file)
				best_parameters = list(read)
			tsv_file.close()
			self.parameters(best_parameters)

		print('parameters analyzed')