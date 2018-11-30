from __future__ import absolute_import, division, print_function

from source.kinetic_flux_model import KineticFluxModel
from source.fitness_function import FitnessFunction
from source.genetic_algorithm import GeneticAlgorithm


class ConfigureEvolution(object):
	'''
	Sets up and runs the genetic algorithm for a given condition
	'''

	def __init__(self, config):


		self.all_reactions = config.get('all_reactions', None)
		include_reactions = config.get('initial_reactions', None)

		self.kinetic_model_config = {
			'km_range': config['km_range'],
			'kcat_range': config['kcat_range'],
			'wcm_sim_data': config['wcm_sim_data'],
			'set_baseline': config['set_baseline'],
			}

		self.evaluator_config = {
			'conditions': config['conditions'],
			}

		self.ga_config = {
			'population_size': config['population_size'],
			'rank_based': config['rank_based'],
			'number_elitist': config['number_elitist'],
			'mutation_variance': config['mutation_variance'],
			'max_fitness': config['max_fitness'],
			'diagnose_error': config['diagnose_error'],
			'seed_parameters': config['seed_parameters'], # TODO -- this can be passed to GA from fitness function.
			'temperature': config['temperature'],
			'stochastic_acceptance': config['stochastic_acceptance'],
			}

		# initialize reactions
		self.reactions = {}
		self.add_reactions(include_reactions)


	# def run_evolution(self, condition, n_generations):
	def run_evolution(self, n_generations):

		# make the kinetic transport model with baseline concentrations
		self.kinetic_model = KineticFluxModel(self.kinetic_model_config, self.reactions)

		# configure the fitness function, passing in the kinetic model
		self.fitness_function = FitnessFunction(self.evaluator_config, self.kinetic_model)

		# configure the genetic algorithm, passing in a fitness function
		self.genetic_algorithm = GeneticAlgorithm(self.ga_config, self.fitness_function, n_generations)

		# run the genetic algorithm
		results = self.genetic_algorithm.evolve()

		return results


	def add_reactions(self, add_reactions):

		new_reactions = {reaction: self.all_reactions[reaction] for reaction in add_reactions}
		self.reactions.update(new_reactions)
