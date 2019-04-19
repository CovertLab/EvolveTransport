import numpy as np
from source.kinetic_flux_model import KineticFluxModel

kinetic_model_config = {}
reactions = {}

kinetic_model = KineticFluxModel(kinetic_model_config, reactions)


parameters = np.array([])
concentrations = {}

fluxes = kinetic_model.get_fluxes(parameters, concentrations)
