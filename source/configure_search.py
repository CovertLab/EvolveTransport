from __future__ import absolute_import, division, print_function

import os
import json
import argparse
import datetime

import numpy as np

from source.kinetic_flux_model import KineticFluxModel
from source.evaluation import FitnessEvaluator
from source.genetic_algorithm import GeneticAlgorithm
from source.visualize import Visualize

from source import data


## options
ENFORCE_BOUNDS = True
PARAMETER_ANALYTICS = False
RANK_BASED_SELECTION = False
DIAGNOSE_ERROR = True

# threshold for saving successful parameters
SAVE_FITNESS_THRESHOLD = 0.95


# simulation parameters
TIME_TOTAL = 1.0 # seconds
TIME_STEP = 0.1 # seconds

# genetic algorithm parameters
POPULATION_SIZE = 100
MAX_GENERATIONS = 101
FITNESS_MAX = 0.999
NUMBER_ELITIST = 2
STOCHASTIC_ACCEPTANCE = False

# for staging
ACCEPTANCE_TEMPERATURE = 0.3
MUTATION_VARIANCE = 0.1 #0.001 #0.1 # TODO -- make this default, and allow adjustable mutation variance in conditions


# set allowable parameter ranges
# A concentration of one molecule per E. coli cell is roughly 1 nM (1e-9 M),
# while water, the most abundant species, has a concentration of about 50 M.
PARAM_RANGES = {
	'km': [
		1e-9, #1e-9,  # units in M
		1e0  #1e1 units in M
	],
	'kcat': [
		1e-2, # gives average kcat of about 100 w/ upper kcat of 1e6
		1e5   # catalase is around 1e5 /s
	],
	}



# Set up conditions
'''
A condition is a dictionary that includes "initial_concentrations", "targets", "penalties" as keys. 
Each of these is itself a dictionary. 
 
"initial_concentrations" has {molecule : concentration}, and sets this conditions concentrations to those listed. 
If unlisted, it uses WCM concentrations.
 
"targets" has these four sub dictionaries: "reaction_fluxes", "exchange_fluxes", "concentrations", "parameters", 
each which has a reaction, molcule, or parameter id as entries. 

"penalities" has the same four entries as "targets", but with each having a single penalty term.
  
'''

# Saved conditions
TEST_SHARED_TRANSPORTER = False
TEST_LEUCINE = True


with open(data.CONDITIONS_FILE, "r") as f:
	conditions_dict = json.loads(f.read())


	# INCLUDE_REACTIONS = [
	# 	'TRANS-RXN0-265-HIS//HIS.9.',
	# 	'TRANS-RXN0-265-TRP//TRP.9.',
	# 	'TRANS-RXN0-265-PHE//PHE.9.'
	# ]

if TEST_LEUCINE:

	INCLUDE_REACTIONS = [
		'TRANS-RXN-126B',
		'TRANS-RXN0-569-LEU//LEU.9.',
		'TRANS-RXN0-270',
		'ABC-35-RXN',
	]

	BASELINE_CONCS = {
		# 'PROTON[p]' : 1e-2,
		}

	C1 = {
		"initial_concentrations": {
			'LEU[p]' : 1e-4,
		},
		"targets": {
		  "exchange_fluxes": {
			'LEU[p]': 0.5e-5,
		  },
		},
		"penalties": {
		  "exchange_fluxes": 10.0
		}
	}

	C2 = {
		"initial_concentrations": {
			'LEU[p]' : 1e-3,
		},
		"targets": {
		  "exchange_fluxes": {
			'LEU[p]': 0.5e-4,
		  },
		},
		"penalties": {
		  "exchange_fluxes": 10.0
		}
	}

	CONDITIONS = [C1,C2]


if TEST_SHARED_TRANSPORTER:

	INCLUDE_REACTIONS = ['RXN0-5202', 'TRANS-RXN-62B']

	BASELINE_CONCS = {
		'PROTON[p]' : 1e-2,
		}



	C1 = {
		"initial_concentrations": {
			'GLY[p]' : 1e-4,
			'L-ALPHA-ALANINE[p]' : 1e-5, #lower
		},
		"targets": {
		  "reaction_fluxes": {
			'TRANS-RXN-62B': 0.5e-5,
			'RXN0-5202': 0.5e-5,
		  },
		},
		"penalties": {
		  "reaction_fluxes": 10.0
		}
	}

	C2 = {
		"initial_concentrations": {
			'GLY[p]' : 1e-1,
			'L-ALPHA-ALANINE[p]' : 1e-3, # higher
		},
		"targets": {
		  "reaction_fluxes": {
			'TRANS-RXN-62B': 1e-4,
			'RXN0-5202' : 1e-3,
		  },
		},
		"penalties": {
		  "reaction_fluxes": 10.0
		}
	}

	CONDITIONS = [C1, C2]




# TODO -- add check that all reactions listed in conditions actually exist in INCLUDE_REACTIONS
# CONDITIONS = [conditions_dict['C3']]
# CONDITIONS = [conditions_dict['C1'], conditions_dict['C2']]



# get parameters initialization values from targets, if not included here they are set to random.
INITIAL_PARAMETERS = {}
for condition in CONDITIONS:
	if 'parameters' in condition['targets']:
		params = condition['targets']['parameters']
		INITIAL_PARAMETERS.update(params)



class TransportEstimation(object):

	def __init__(self):

		# set random seed to make search deterministic
		self.seed = np.random.randint(2 ** 32 - 1)
		print('seed = ' + str(self.seed))
		np.random.seed(1)

		# replicate id is used to name outputs of this replicate
		self.replicate_id = self.get_replicate_id()

		self.kinetic_model_config = {
			'km_range': PARAM_RANGES['km'],
			'kcat_range': PARAM_RANGES['kcat'],
			'wcm_sim_data' : data.wcm_sim_out,
			'set_baseline' : BASELINE_CONCS,
			}

		self.evaluator_config = {
			'conditions' : CONDITIONS,
			}

		self.ga_config = {
			'population_size' : POPULATION_SIZE,
			'rank_based' : RANK_BASED_SELECTION,
			'number_elitist' : NUMBER_ELITIST,
			'enforce_bounds' : ENFORCE_BOUNDS,
			'mutation_variance' : MUTATION_VARIANCE,
			'max_generations' : MAX_GENERATIONS,
			'max_fitness' : FITNESS_MAX,
			'diagnose_error' : DIAGNOSE_ERROR,
			'initial_parameters' : INITIAL_PARAMETERS,
			'temperature' : ACCEPTANCE_TEMPERATURE,
			'stochastic_acceptance' : STOCHASTIC_ACCEPTANCE,
			}

		self.plot_config = {
			'out_dir' : data.OUTDIR,
			'parameter_out_dir' : data.PARAMOUTDIR,
			'saved_param_file' : data.PARAM_FILE,
			'parameter_analytics' : PARAMETER_ANALYTICS,
			'mutation_variance': MUTATION_VARIANCE,
			'seed' : self.seed,
			'replicate_id' : self.replicate_id,
			'fitness_threshold' : SAVE_FITNESS_THRESHOLD,
			'wcm_sim_data' : data.wcm_sim_out, # is this needed? initial concentrations are available in the kinetic model
			}



	def main(self):

		# TODO -- staging should be done here.
		# allow passing state between stages, seed population in new GA
		# new reaction definitions based on conditions. new parameter indices.
		# for stage in STAGES:





		# initialize reactions
		self.reactions = {reaction: data.ALL_REACTIONS[reaction] for reaction in INCLUDE_REACTIONS}

		self.conditions = CONDITIONS

		# make the kinetic transport model with baseline concentrations
		self.kinetic_model = KineticFluxModel(self.kinetic_model_config, self.reactions)

		# configure the fitness function, passing in the kinetic model
		self.fitness_function = FitnessEvaluator(self.evaluator_config, self.kinetic_model)

		# configure the genetic algorithm, passing in a fitness function
		self.genetic_algorithm = GeneticAlgorithm(self.ga_config, self.fitness_function)

		# configure plotting
		self.plot = Visualize(self.plot_config, self.fitness_function) ## TODO -- should this use kinetic_model instead of fitness_function?

		# run the genetic algorithm
		final_population, final_fitness, saved_error, saved_fitness, saved_diagnosis = self.genetic_algorithm.evolve()







		## Visualization and Analysis
		# TODO -- make a separate analysis class
		self.plot.parameter_analysis(final_population, final_fitness)

		self.plot.evolution(saved_error, saved_fitness, saved_diagnosis)

		# get best individual's parameters and simulate it
		top_index = final_fitness.values().index(max(final_fitness.values()))
		top_genotype = final_population[top_index]
		top_phenotype = self.fitness_function.get_phenotype(top_genotype)

		run_for = 1
		sim_output = self.kinetic_model.run_simulation(top_phenotype, run_for)

		# plot simulation of best individual
		self.plot.sim_out(sim_output, top_phenotype, self.conditions) # TODO -- why does sim_out care about conditions?

		# # TODO -- make parameter plot of conditions, rather parameters in sim_out
		# plot best individual across all conditions
		self.plot.conditions(top_phenotype, self.conditions)

		# plot michaelis menten curves for best individual across all reactions
		self.plot.michaelis_menten(top_phenotype)


	def get_replicate_id(self):

		now = datetime.datetime.now()
		time_stamp = now.strftime('%m-%d_%H:%M:%S')

		replicate_nums = os.listdir(data.OUTDIR)
		# remove files that don't match pattern
		replicate_nums = [name for name in replicate_nums if '__' in name]
		if not replicate_nums:
			replicate_num = 1
		else:
			replicate_nums = [name.replace(name[0:name.find("__")+2], '') for name in replicate_nums]
			replicate_nums = [name.replace(name[name.find("."):], '') for name in replicate_nums]
			replicate_nums = [int(name) for name in replicate_nums]
			replicate_num = max(replicate_nums) + 1

		replicate_id = (time_stamp + '__' + str(replicate_num))

		return replicate_id


if __name__ == "__main__":
	# parser = argparse.ArgumentParser(description='evolve parameters for transport')
	# parser.add_argument('--simout', help='directory of sim out data', default='out/manual/condition_000002/000000/generation_000000/000000/simOut')
	# args = parser.parse_args()
	TransportEstimation().main()
