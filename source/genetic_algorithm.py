from __future__ import absolute_import, division, print_function

import numpy as np


class GeneticAlgorithm(object):

	def __init__(self, config, fitness_function):

		# configuration
		self.population_size = config.get('population_size', None)
		self.rank_based = config.get('rank_based', False)
		self.number_elitist = config.get('number_elitist', False)
		self.enforce_bounds = config.get('enforce_bounds', True)
		self.diagnose_error = config.get('diagnose_error', False)
		self.mutation_variance = config.get('mutation_variance', 1.0)
		self.max_generations = config.get('max_generations', 1000)
		self.max_fitness = config.get('max_fitness', 0.99)
		self.initial_parameters = config.get('initial_parameters', None)
		self.temperature = config.get('temperature', 1.0)
		self.stochastic_acceptance = config.get('stochastic_acceptance', False)

		# fitness function
		self.fitness_function = fitness_function
		self.genome_size = self.fitness_function.genome_size
		self.reactions = self.fitness_function.reactions
		self.parameter_indices = self.fitness_function.parameter_indices

		# initialize the population
		self.population = self.initialize_population(self.population_size)

		# initialize dictionaries for the output
		self.fitness = {individual: 0 for individual, genome in self.population.iteritems()}
		self.total_error = {individual: 0 for individual, genome in self.population.iteritems()}
		self.diagnosis = {individual: {} for individual, genome in self.population.iteritems()}

		# if self.diagnose_error:
		self.saved_mutation_diagnosis = []

	def initialize_population(self, population_size):
		'''
		fill the population dicitonary with {individual: [parameters]}
		for each individual in the population

		'''

		population = {}
		for ind in xrange(population_size):
			population[ind] = self.initialize_genome()

		return population

	def initialize_genome(self):

		genome = np.random.uniform(0, 1, self.genome_size)

		# initialize defined parameters to target value
		for rxn, target_parameters in self.initial_parameters.iteritems():
			for parameter, phenotypic_target in target_parameters.iteritems():
				if phenotypic_target:

					# get index of this phenotypic target by looking up reaction's transporters,
					# and the parameters associated with them.
					transporters = self.reactions[rxn]['transporters']
					for transporter in transporters:
						param_indices = self.parameter_indices[rxn][transporter]

						if 'kcat' in parameter:
							param_idx = param_indices[parameter]
							gene_value = self.fitness_function.phenotype_transform[param_idx]['pheno_to_geno'](phenotypic_target)
							genome[param_idx] = gene_value

						else:  # km, uses molecule id as parameter name
							param_idx = param_indices['kms'][parameter]
							gene_value = self.fitness_function.phenotype_transform[param_idx]['pheno_to_geno'](phenotypic_target)
							genome[param_idx] = gene_value
		return genome

	def evolve(self):
		generations = 0
		top_fit = 0
		saved_error = []
		saved_fitness = []
		saved_penality_diagnosis = []

		# genetic algorithm loop
		while generations < self.max_generations and top_fit < self.max_fitness:

			# evaluate fitness of each individual
			for individual, genome in self.population.iteritems():
				# TODO -- clean up diagnose_error option
				if self.diagnose_error:
					self.total_error[individual], self.diagnosis[individual] = self.fitness_function.evaluate(genome, self.diagnose_error)
				else:
					self.total_error[individual] = self.fitness_function.evaluate(genome)

			# get fitness from error
			self.fitness = {individual: self.fitness_from_error(total_error) for individual, total_error in self.total_error.iteritems()}

			# save values for analysis
			saved_fitness.append(self.fitness.values())
			saved_error.append(self.total_error.values())
			top_fit = max(self.fitness.values())

			# save error contributions to top individual
			if self.diagnose_error:
				top_index = self.fitness.values().index(top_fit)
				top_diagnosis = self.diagnosis[top_index]
				saved_penality_diagnosis.append(top_diagnosis)

			print('gen ' + str(generations) + ' fittest: ' + str(top_fit)) # TODO -- print error for error terms

			# repopulate based on fitness -- higher fitness gets higher selection
			self.population = self.repopulate(self.population, self.fitness)

			generations += 1

		if top_fit >= self.max_fitness:
			print('Success!')
		else:
			print('Did not find solution')

		saved_diagnosis = {
			'penality_diagnosis': saved_penality_diagnosis,
			'mutation_diagnosis': self.saved_mutation_diagnosis,
		}

		return self.population, self.fitness, saved_error, saved_fitness, saved_diagnosis

	def repopulate(self, population, fitness):
		'''
		population is a dictionary with {id: genome}
		fitness is a dictionary with {id: fitness value}
		'''

		new_population = {}

		if self.rank_based:
			# order individuals by rank
			# TODO -- make sure this is still working
			rank_ordered = sorted(fitness, key=fitness.get, reverse=True)
			total_rank_fitness = sum([1.0 / n for n in xrange(1, self.population_size + 1)])

			# selection fitness is 1/rank normalized to 1
			selection_fitness = {individual: 1.0 / rank / total_rank_fitness for rank, individual in enumerate(rank_ordered, 1)}
		else:
			# normalized fitness
			total_fitness = np.sum(fitness.values())
			selection_fitness = {individual: value / total_fitness for individual, value in fitness.items()}

		repopulation_index = 0

		# elitist selection: re-seed the top individuals
		top_indices = np.argsort(fitness.values())[-self.number_elitist:]
		for idx in top_indices:
			new_population[repopulation_index] = population[idx]
			repopulation_index += 1

		while len(new_population) < self.population_size:

			enforce_bounds = 0

			# Selection
			# TODO -- use numpy multimodal instead, for fitness-proportionate selection
			selection_index = 0
			total = selection_fitness[selection_index]
			rand_selection = np.random.uniform(0, 1)
			while total < rand_selection:
				selection_index += 1
				total += selection_fitness[selection_index]

			# Mutation
			# TODO -- vectorize these steps, apply mutations to the entire population at once.
			genome = population[selection_index]
			genome_fitness = fitness[selection_index]

			# gaussian distance
			magnitude = np.random.normal(0, self.mutation_variance)

			# random unit vector
			direction = [np.random.normal(0, 1) for i in range(len(genome))]
			direction_mag = np.sum(x ** 2 for x in direction) ** 0.5
			mutation = [magnitude * x / direction_mag for x in direction]

			# apply mutation
			new_genome = np.array([x + y for x, y in zip(genome, mutation)])

			# enforce bounds on genome
			if self.enforce_bounds:
				# clip at bounds
				# new_genome[new_genome >= 1.0] = 1.0
				# new_genome[new_genome <= 0.0] = 0.0

				# if parameter is not in range, initialize it randomly within range
				out_of_range = np.where(np.logical_or(new_genome <= 0.0, new_genome >= 1.0))
				new_genome[out_of_range] = np.random.uniform(0.0, 1.0)

				enforce_bounds = float(len(out_of_range[0])) / len(new_genome)

			# TODO -- check if mutation improves objective. Accept improvements. Take reductions with p(d_obj)

			# genome_fitness
			new_genome_error = self.fitness_function.evaluate(new_genome)
			new_genome_fitness = self.fitness_from_error(new_genome_error)

			# compare fitness
			fitness_change = new_genome_fitness - genome_fitness

			if self.stochastic_acceptance:
				# save all fitness improvements
				if fitness_change > 0:
					new_population[repopulation_index] = new_genome
					repopulation_index += 1

				# save fitness reductions with probability based on exponential function.
				elif np.random.rand() < np.exp(fitness_change / self.temperature):
					new_population[repopulation_index] = new_genome
					repopulation_index += 1
			else:
				# repopulate without boltzmann acceptance function
				new_population[repopulation_index] = new_genome
				repopulation_index += 1

			if self.diagnose_error and len(self.saved_mutation_diagnosis) < 5000:

				# TODO -- if self.stochastic_acceptance is False, don't save acceptance_p
				diagnosis = {
					'mutation': mutation,
					'fitness_change': fitness_change,
					'acceptance_p': np.exp(fitness_change / self.temperature),
					'enforce_bounds': enforce_bounds,
				}

				self.saved_mutation_diagnosis.append(diagnosis)

		return new_population

	def fitness_from_error(self, error):
		return 1 / (1 + error)
