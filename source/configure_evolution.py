from __future__ import absolute_import, division, print_function

from source.kinetic_flux_model import KineticFluxModel
from source.fitness_function import FitnessFunction
from source.genetic_algorithm import GeneticAlgorithm


class ConfigureEvolution(object):
	'''
	Sets up and runs the genetic algorithm for a given condition
	'''

	def __init__(self, config, initial_reactions, condition):

		self.kinetic_model_config = config.get('kinetic_model_config', None)
		self.evaluator_config = config.get('evaluator_config', None)
		self.ga_config = config.get('ga_config', None)
		self.all_reactions = config.get('all_reactions', None)

		# initialize reactions
		self.reactions = {}
		self.add_reactions(initial_reactions)

		self.condition = condition


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

	def reconfigure(self, config):

		pass


	def map_parameters(self):

		pass
