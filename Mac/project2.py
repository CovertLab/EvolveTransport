#Building fitness function

# find fitness/error of each transport reaction, molecules sum across all molecules


from project1 import get_flux
from source.fitness_function import *
from source.kinetic_flux_model import *

"""
total_error = 0
reaction_flux = 5*(10**-6)
flux = get_flux()[0]
for rxn,target_flux in flux.iteritems():
    error = 10 * abs( reaction_flux- target_flux)
    total_error += error

print(total_error)
"""

import json
import os

import numpy as np

from source.kinetic_flux_model import KineticFluxModel



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
ARAM_FILE = 'best_parameters.tsv'
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(BASE_DIR, 'out')
DATADIR = os.path.join(BASE_DIR, 'data')
WCM_SIMDATA_FILE = os.path.join(DATADIR, 'wcm_sim_data.json')

BASELINE_CONCS = {
        'PROTON[p]': 1e-2,
    }



with open(WCM_SIMDATA_FILE, 'r') as f:
    wcm_sim_out = json.loads(f.read())

kinetic_model_config = {
            # for kinetic model config
            'km_range': PARAM_RANGES['km'],
            'kcat_range': PARAM_RANGES['kcat'],
            'wcm_sim_data': wcm_sim_out,
            'set_baseline': BASELINE_CONCS,
    # for fitness function config
        }


# get a flux at a specific instant in time



# from data.py


REACTIONS_FILE = os.path.join(
    DATADIR, 'aa_transport_reactions.json'
)

with open(REACTIONS_FILE, 'r') as f:
    reactions = json.loads(f.read())


# make kinetic flux model

kinetic_model = KineticFluxModel(kinetic_model_config, reactions)


# make fitness function config

config = {'conditions': [{'penalties': {'reaction_fluxes': 10.0}, 'targets': {'reaction_fluxes': {'RXN0-5202': 5e-06, 'TRANS-RXN-62B': 5e-06}}, 'initial_concentrations': {'L-ALPHA-ALANINE[p]': 1e-05, 'GLY[p]': 0.0001}}, {'penalties': {'reaction_fluxes': 10.0}, 'targets': {'reaction_fluxes': {'RXN0-5202': 0.001, 'TRANS-RXN-62B': 0.0001}}, 'initial_concentrations': {'L-ALPHA-ALANINE[p]': 0.001, 'GLY[p]': 0.1}}]}

# pass kinetic model and config to fitness function
fitness_function= FitnessFunction(config, kinetic_model)

genome_size = kinetic_model.n_parameters

genome = np.random.uniform(0, 1, genome_size)
# evaluate fitness function for a single random genotype
error = fitness_function.evaluate(genome)

fitness = 1 / (1 + error)

print(fitness)