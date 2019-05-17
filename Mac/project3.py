from __future__ import absolute_import, division, print_function

import datetime
import os

import numpy as np

from source import data
#from source.analyze import Analyze
from source.configure_evolution import ConfigureEvolution
#from source.visualize import Visualize

# options
PARAMETER_ANALYTICS = False
RANK_BASED_SELECTION = False
DIAGNOSE_ERROR = True

# threshold for saving successful parameters
SAVE_FITNESS_THRESHOLD = 0.95 # TODO -- track threshold



# genetic algorithm parameters
POPULATION_SIZE = 100
DEFAULT_N_GENERATIONS = 101
FITNESS_MAX = 0.999
NUMBER_ELITIST = 2
STOCHASTIC_ACCEPTANCE = True

# for staging
ACCEPTANCE_TEMPERATURE = 0.3
MUTATION_VARIANCE = 0.005  # 0.001 # 0.1 # TODO -- make this default, and adjust mutation variance in conditions


"process reaction from exchange file"



# set allowable parameter ranges
# A concentration of one molecule per E. coli cell is roughly 1 nM (1e-9 M),
# while water, the most abundant species, has a concentration of about 50 M.
PARAM_RANGES = {
    'km': [
        1e-9,  # 1e-9,  # units in M
        1e-1  # 1e1 units in M
    ],
    'kcat': [
        1e-2,  # 1e-2 gives average kcat of about 100 w/ upper kcat of 1e6
        1e5  # 1e5 catalase is around 1e5/s
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
AA_Exchange= True


BASELINE_CONCS = {}
INCLUDE_EXCHANGE = []



if AA_Exchange:
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
            'conditions': None,

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

    def main(self):

        stage = {
                'n_generations': 50,
                'include_reactions': initial_reactions, #TODO -- check where we use initial_reaction
                'seed_results_from': [],
                'add_reactions': [],
                'conditions': CONDITIONS,  # will be 72 reaction
                'mutation_variance': 0.05,
        }
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

        # stages_parameters = {}
        # for stage in seed_results_from:
        #     params = phenotype_summaries[stage]
        #     stages_parameters.update(params)

        # TODO -- this could overwrite target parameters. Is this what you want?
        #seed_parameters.update(stages_parameters)

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


    def get_top_phenotype(self, population, fitness):
        top_index = fitness.values().index(max(fitness.values()))
        top_genotype = population[top_index]
        top_phenotype = self.configuration.fitness_function.get_phenotype(top_genotype)

        return top_phenotype

    def get_replicate_id(self):

        now = datetime.datetime.now()
        time_stamp = now.strftime('%m-%d_%H:%M:%S')

        replicate_nums = os.listdir(data.PLOTOUTDIR)
        # remove files that don't match pattern
        replicate_nums = [name for name in replicate_nums if '__' in name]
        if not replicate_nums:
            replicate_num = 1
        else:
            replicate_nums = [name.replace(name[0:name.find('__') + 2], '') for name in replicate_nums]
            replicate_nums = [name.replace(name[name.find('.'):], '') for name in replicate_nums]
            replicate_nums = [int(name) for name in replicate_nums]
            replicate_num = max(replicate_nums) + 1

        replicate_id = (time_stamp + '__' + str(replicate_num))

        return replicate_id



















if __name__ == '__main__':
    Main().main()
