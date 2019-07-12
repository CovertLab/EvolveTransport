from __future__ import absolute_import, division, print_function

import os
import datetime

import numpy as np

from source.configure_evolution import ConfigureEvolution
from source.analyze import Analyze
from source.visualize import Visualize
from source import data


# options
PARAMETER_ANALYTICS = False
RANK_BASED_SELECTION = False
DIAGNOSE_ERROR = True

# threshold for saving successful parameters
SAVE_FITNESS_THRESHOLD = 0.95

# simulation parameters
TIME_TOTAL = 1.0  # seconds
TIME_STEP = 0.1  # seconds

# genetic algorithm parameters
POPULATION_SIZE = 100
DEFAULT_N_GENERATIONS = 101
FITNESS_MAX = 0.999
NUMBER_ELITIST = 2
STOCHASTIC_ACCEPTANCE = True

# for staging
ACCEPTANCE_TEMPERATURE = 0.3
MUTATION_VARIANCE = 0.005  # 0.001 # 0.1 # TODO -- make this default, and adjust mutation variance in conditions


def reactions_from_exchange(include_exchanges):
	include_reactions = []
	for reaction_id, specs in data.ALL_REACTIONS.iteritems():
		reaction_molecules = specs['stoichiometry'].keys()

		for exchange in include_exchanges:
			if exchange in reaction_molecules:
				# add the reaction
				include_reactions.append(reaction_id)

	return include_reactions


# set allowable parameter ranges
# A concentration of one molecule per E. coli cell is roughly 1 nM (1e-9 M),
# while water, the most abundant species, has a concentration of about 50 M.
PARAM_RANGES = {
	'km': [
		1e-9,  # 1e-9,  # units in M
		1e-1    # 1e1 units in M
	],
	'kcat': [
		1e-2,  # 1e-2 gives average kcat of about 100 w/ upper kcat of 1e6
		1e5    # 1e5 catalase is around 1e5/s
	],
	}

# Set up conditions
'''
A condition is a dictionary that includes 'initial_concentrations', 'targets', 'penalties' as keys. 
Each of these is itself a dictionary. 
 
'initial_concentrations' has {molecule: concentration}, and sets this conditions concentrations to those listed. 
If unlisted, it uses WCM concentrations.
 
'targets' has these four sub dictionaries: 'reaction_fluxes', 'exchange_fluxes', 'concentrations', 'parameters', 
each which has a reaction, molcule, or parameter id as entries. 

'penalities' has the same four entries as 'targets', but with each having a single penalty term.
  
'''

# Saved conditions
TEST_SHARED_TRANSPORTER = True
TEST_LEUCINE = False
TEST_PIPERNO = False

BASELINE_CONCS = {}
INCLUDE_EXCHANGE = []

if TEST_LEUCINE:

	INCLUDE_EXCHANGE = ['LEU[p]']  # Piperno data: ['GLY[p]', 'ILE[p]', 'MET[p]', 'PHE[p]']
	initial_reactions = reactions_from_exchange(INCLUDE_EXCHANGE)

	C1 = {
		'initial_concentrations': {
			'LEU[p]': 1e-4,
		},
		'targets': {
			'exchange_fluxes': {
				'LEU[p]': 0.5e-5,
			},
		},
		'penalties': {
			'exchange_fluxes': 10.0
		}
	}

	C2 = {
		'initial_concentrations': {
			'LEU[p]': 1e-3,
		},
		'targets': {
			'exchange_fluxes': {
				'LEU[p]': 0.5e-4,
			},
		},
		'penalties': {
			'exchange_fluxes': 10.0
		}
	}

	CONDITIONS = [C1, C2]

if TEST_SHARED_TRANSPORTER:

	initial_reactions = ['RXN0-5202', 'TRANS-RXN-62B']

	BASELINE_CONCS = {
		'PROTON[p]': 1e-2,
		}

	C1 = {
		'initial_concentrations': {
			'GLY[p]': 1e-4,
			'L-ALPHA-ALANINE[p]': 1e-5,  # lower
		},
		'targets': {
			'reaction_fluxes': {
				'TRANS-RXN-62B': 0.5e-5,
				'RXN0-5202': 0.5e-5,
			},
		},
		'penalties': {
		  'reaction_fluxes': 10.0
		}
	}
	C2 = {
		'initial_concentrations': {
			'GLY[p]': 1e-1,
			'L-ALPHA-ALANINE[p]': 1e-3,  # higher
		},
		'targets': {
			'reaction_fluxes': {
				'TRANS-RXN-62B': 1e-4,
				'RXN0-5202': 1e-3,
			},
		},
		'penalties': {
			'reaction_fluxes': 10.0
		}
	}

	CONDITIONS = [C1, C2]

if TEST_PIPERNO:

	INCLUDE_EXCHANGE = ['GLY[p]']  # ['GLY[p]'] #['GLY[p]', 'ILE[p]', 'MET[p]', 'PHE[p]']
	initial_reactions = reactions_from_exchange(INCLUDE_EXCHANGE)

	BASELINE_CONCS = {
		'PROTON[p]': 1e-2,
		}

	# make conditions from data
	CONDITIONS = []
	for flux_id in INCLUDE_EXCHANGE:
		target_data = data.target_definition[flux_id]

		for target in target_data:

			initial_concentration = float(target['substrate_concentration'])
			target_flux = float(target['flux'])

			condition = {
				'initial_concentrations': {
					flux_id: initial_concentration,
					# 'CYCA-MONOMER' : 0.0, # turn off 'TRANS-RXN-62B' by setting transporter concentration to 0
					# 'CPLX0 - 7654' : 0.0, # turn off 'TRANS-RXN0-537' by setting transporter concentration to 0
				},
				'targets': {
					'exchange_fluxes': {flux_id: - target_flux},  # need negative flux, because uptake removes from [p]
					'parameters': {
						'TRANS-RXN0-537': {
							'CPLX0-7654': {'kcat_f': 1e-2, 'GLY[p]': 1e-2}, # low kcat to turn off rxn
						},
					}
				},
				'penalties': {
					'exchange_fluxes': 1,#1e6,
				},
			}

			CONDITIONS.append(condition)


# # get parameters initialization values from targets, if not included here they are set to random.
# INITIAL_PARAMETERS = {}
# for condition in CONDITIONS:
# 	if 'parameters' in condition['targets']:
# 		params = condition['targets']['parameters']
# 		INITIAL_PARAMETERS.update(params)



class Main(object):

	def __init__(self):

		# set random seed.
		# seed at a constant to repeat searches
		self.seed = np.random.randint(2 ** 32 - 1)
		print('seed = ' + str(self.seed))
		np.random.seed(self.seed)

		# replicate id for naming output figures
		self.replicate_id = self.get_replicate_id()

		self.evo_config = {
			'all_reactions': data.ALL_REACTIONS,
			'include_reactions': [],

			# for kinetic model config
			'km_range': PARAM_RANGES['km'],
			'kcat_range': PARAM_RANGES['kcat'],
			'wcm_sim_data': data.wcm_sim_out,
			'set_baseline': BASELINE_CONCS,

			# for fitness function config
			'conditions': None, # CONDITIONS, # TODO -- this should be set within the stage.

			# for genetic algorithm config
			'n_generations': DEFAULT_N_GENERATIONS,
			'population_size': POPULATION_SIZE,
			'rank_based': RANK_BASED_SELECTION,
			'number_elitist': NUMBER_ELITIST,
			'mutation_variance': MUTATION_VARIANCE,
			'max_fitness': FITNESS_MAX,
			'diagnose_error': DIAGNOSE_ERROR,
			'seed_parameters': None,
			'temperature': ACCEPTANCE_TEMPERATURE,
			'stochastic_acceptance': STOCHASTIC_ACCEPTANCE,
		}

		self.visualize_config = {
			'out_dir': data.PLOTOUTDIR,
			'parameter_out_dir': data.PARAMOUTDIR,
			'saved_param_file': data.PARAM_FILE,
			'parameter_analytics': PARAMETER_ANALYTICS,
			'mutation_variance': MUTATION_VARIANCE, # TODO -- this should get mutation variance from the evo_config. changes through the stages.
			'seed': self.seed,
			'replicate_id': self.replicate_id,
			'fitness_threshold': SAVE_FITNESS_THRESHOLD,
			'exchange_molecules': INCLUDE_EXCHANGE,
			'wcm_sim_data': data.wcm_sim_out, # TODO -- initial concentrations are already available in the kinetic model
			}

	def main(self):

		stages = {
			1: {
			'n_generations': 50,
			'include_reactions': initial_reactions,
			'seed_results_from': [],
			'add_reactions': [],
			'conditions': CONDITIONS,
			'mutation_variance': 0.05,
			},
			2: {
			'n_generations': 50,
			# 'include_reactions': initial_reactions,
			'seed_results_from': [1],
			'add_reactions': ['RXN0-5202'],
			'conditions': CONDITIONS,
			'mutation_variance': 0.001,
			},
		}

		phenotype_summaries = {}
		all_results = {}

		for stage_id, stage in stages.iteritems():

			self.conditions = stage['conditions']
			add_reactions = stage['add_reactions']
			seed_results_from = stage['seed_results_from']

			# update reactions
			include_reactions = self.evo_config['include_reactions']
			include_reactions.extend(add_reactions)
			self.evo_config['include_reactions'] = include_reactions

			# update evo_config
			self.evo_config.update(stage)

			# initialize seed_parameters
			seed_parameters = {}

			# add target parameters
			target_parameters = {}
			for condition in stage['conditions']:
				if 'parameters' in condition['targets']:
					params = condition['targets']['parameters']
					target_parameters.update(params)
			seed_parameters.update(target_parameters)

			# add parameters from previous stages.
			stages_parameters = {}
			for stage in seed_results_from:
				params = phenotype_summaries[stage]
				stages_parameters.update(params)

			# TODO -- this could overwrite target parameters. Is this what you want?
			seed_parameters.update(stages_parameters)

			# seed parameters
			self.evo_config['seed_parameters'] = seed_parameters

			# configure evolution and run
			self.configuration = ConfigureEvolution(self.evo_config)
			results = self.configuration.run_evolution()

			# save top phenotype's parameters
			# TODO -- should be able to save the top N phenotypes, reseed all of them.
			final_population = results['final_population']
			final_fitness = results['final_fitness']
			top_phenotype = self.get_top_phenotype(final_population, final_fitness)
			top_phenotype_summary = self.configuration.kinetic_model.get_phenotype_summary(top_phenotype)
			phenotype_summaries[stage_id] = top_phenotype_summary

			# save results
			all_results[stage_id] = results


		# configure plotting
		## TODO -- should this use kinetic_model instead of fitness_function?
		self.analyze = Analyze(self.visualize_config, self.configuration.fitness_function)
		self.plot = Visualize(self.visualize_config, self.configuration.fitness_function)

		# Visualization and Analysis
		self.visualize(all_results)


	def seed_parameters_from_targets(self):

		pass


	def get_top_phenotype(self, population, fitness):
		top_index = fitness.values().index(max(fitness.values()))
		top_genotype = population[top_index]
		top_phenotype = self.configuration.fitness_function.get_phenotype(top_genotype)

		return top_phenotype


	# TODO -- this should be in visualize. set with self.visualize_config
	def visualize(self, all_results):

		# TODO -- use all results to run visualize. for evolution in particular
		results = all_results[1]

		final_population = results['final_population']
		final_fitness = results['final_fitness']
		saved_error = results['saved_error']
		saved_fitness = results['saved_fitness']
		saved_diagnosis = results['saved_diagnosis']

		self.analyze.parameters(final_population, final_fitness)

		self.plot.evolution(saved_error, saved_fitness, saved_diagnosis)

		# simulate the best individual
		top_phenotype = self.get_top_phenotype(final_population, final_fitness)

		run_for = 1
		sim_output = self.configuration.kinetic_model.run_simulation(top_phenotype, run_for)

		# plot simulation of best individual
		self.plot.sim_out(sim_output, top_phenotype, self.conditions) # TODO -- why does sim_out care about conditions?

		# plot best individual across all conditions
		self.plot.conditions(top_phenotype, self.conditions)

		# plot michaelis menten curves for best individual across all reactions
		self.plot.michaelis_menten(top_phenotype)


	def get_replicate_id(self):

		now = datetime.datetime.now()
		time_stamp = now.strftime('%m-%d_%H:%M:%S')

		replicate_nums = os.listdir(data.PLOTOUTDIR)
		# remove files that don't match pattern
		replicate_nums = [name for name in replicate_nums if '__' in name]
		if not replicate_nums:
			replicate_num = 1
		else:
			replicate_nums = [name.replace(name[0:name.find('__')+2], '') for name in replicate_nums]
			replicate_nums = [name.replace(name[name.find('.'):], '') for name in replicate_nums]
			replicate_nums = [int(name) for name in replicate_nums]
			replicate_num = max(replicate_nums) + 1

		replicate_id = (time_stamp + '__' + str(replicate_num))

		return replicate_id


if __name__ == '__main__':
	Main().main()
