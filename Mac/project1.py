import json
import os

import numpy as np

from source.kinetic_flux_model import KineticFluxModel

def get_flux():

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


    kinetic_model = KineticFluxModel(kinetic_model_config, reactions)


    parameters = np.random.rand(220,1)
    concentrations = kinetic_model.initialize_state()

    fluxes = kinetic_model.get_fluxes(parameters, concentrations)
    return fluxes