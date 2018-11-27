from __future__ import absolute_import, division, print_function

import scipy.constants as constants
import numpy as np

class KineticFluxModel(object):

	def __init__(self, config, reactions):
		'''
		Args:
			reactions: a dictionary with reaction ids assigned to subdictionary with a 'transporters' list,
				stoichiometry, and reversibility

		parameter_indices: a dictionary, organized by transporter -- each with 'reaction_cofactors' for
			each reaction that it is a part of, 'partition' for all of its bounded forms, and 'parameter_indices'
			for the indices at which each of its parameters can be assigned in an array.

		rate_laws: a dict, with a key for each reaction id, and then subdictionaries with each reaction's transporters
			and their rate law function. These rate laws are used directly from within this dictionary
		'''

		self.time_step = 0.1
		self.avogadro = constants.Avogadro

		self.km_range = config.get('km_range', None)
		self.kcat_range = config.get('kcat_range', None)
		self.wcm_sim_data = config.get('wcm_sim_data', None)
		self.set_baseline = config.get('set_baseline', None)

		self.reactions = reactions
		self.parameter_indices, self.transport_configuration, self.n_parameters = self.make_parameter_indices(reactions)

		self.rate_laws = self.make_rate_laws(self.reactions, self.transport_configuration, self.parameter_indices)

		self.baseline_concentrations = self.initialize_state(self.set_baseline)

	def make_parameter_indices(self, reactions):

		# initialize parameter dictionary for all reactions and transporters
		transport_configuration = {}
		parameter_indices = {}
		parameter_index = 0

		for reaction, specs in reactions.iteritems():
			transporters = specs['transporters']
			parameter_indices[reaction] = {}

			# initialize transporters' entries
			for transporter in transporters:
				if transporter not in parameter_indices[reaction]:
					parameter_indices[reaction][transporter] = {
						'kms': {}
					}
				if transporter not in transport_configuration:
					transport_configuration[transporter] = {
						'partition': [],
						'reaction_cofactors': {},
					}

		# get parameters for all reactions
		for reaction, specs in reactions.iteritems():
			stoich = specs['stoichiometry']
			transporters = specs['transporters']
			reversibility = specs['reversibility']

			# get sets of cofactors driving this reaction
			forward_cofactors = [mol for mol, coeff in stoich.iteritems() if coeff < 0]
			cofactors = [forward_cofactors]

			if reversibility:
				reverse_cofactors = [mol for mol, coeff in stoich.iteritems() if coeff > 0]
				cofactors.append(reverse_cofactors)

			# get partition, reactions, and parameter indices for each transporter, and save to transport_configuration dictionary
			for transporter in transporters:

				# get competition for this transporter from other reactions
				competing_reactions = [rxn for rxn, specs2 in reactions.iteritems() if
						(rxn is not reaction) and (transporter in specs2['transporters'])]

				competitors = []
				for reaction2 in competing_reactions:
					stoich2 = reactions[reaction2]['stoichiometry']
					reactants2 = [mol for mol, coeff in stoich2.iteritems() if coeff < 0]
					competitors.append(reactants2)

				# partition includes competitors and cofactors.
				partition = competitors + cofactors

				transport_configuration[transporter]['partition'] = partition
				transport_configuration[transporter]['reaction_cofactors'][reaction] = cofactors

				# add kms to parameter indices
				partitioned_molecules = set([molecule for part in partition for molecule in part])

				# add each new parameter to parameter_indices
				sharing_reactions = competing_reactions + [reaction]
				for molecule in partitioned_molecules:
					new_param = False

					# assign same km indices to all reactions that use this same transporter
					for rx in sharing_reactions:
						if molecule not in parameter_indices[rx][transporter]['kms']:
							parameter_indices[rx][transporter]['kms'][molecule] = parameter_index
							new_param = True

					# only increase index if a new param has been assigned
					if new_param:
						parameter_index += 1

				# each reaction gets a kcat_f for a given combination of cofactors:
				parameter_indices[reaction][transporter]['kcat_f'] = parameter_index
				parameter_index += 1

				if reversibility:
					parameter_indices[reaction][transporter]['kcat_r'] = parameter_index
					parameter_index += 1

		return parameter_indices, transport_configuration, parameter_index


	# Make rate laws
	def make_rate_laws(self, reactions, transport_configuration, parameter_indices):

		rate_laws = {reaction: {} for reaction, specs in reactions.iteritems()}

		# make rate law for each reaction
		for reaction, specs in reactions.iteritems():
			stoichiometry = specs['stoichiometry']
			transporters = specs['transporters']

			# rate law for each transporter
			for transporter in transporters:

				cofactors_sets = transport_configuration[transporter]['reaction_cofactors'][reaction]
				partition = transport_configuration[transporter]['partition']

				param_idxs = parameter_indices[reaction][transporter]
				kcat_indices = {param:idx for param, idx in param_idxs.iteritems() if param != 'kms'}
				km_indices = param_idxs['kms']

				rate_law = self.generate_rate_law(
					stoichiometry,
					transporter,
					cofactors_sets,
					partition,
					kcat_indices,
					km_indices
				)

				# save the rate law for each transporter in this reaction
				# reactions[reaction]['rate_laws'][transporter] = rate_law
				rate_laws[reaction][transporter] = rate_law

		return rate_laws

	def cofactor_numerator(self, concentration, km):
		def term():
			return concentration / km if km else 0

		return term

	def cofactor_denominator(self, concentration, km):
		def term():
			return 1 + concentration / km if km else 1

		return term

	def generate_rate_law(self, stoichiometry, transporter, cofactors_sets, partition, kcat_indices, km_indices):
		'''

		Args:
			cofactors: a list with the required cofactors , each pair needs a kcat.
			partition: a list of lists. each sublist is the set of cofactors for a given partition.
				[[C1, C2],[C3, C4], [C5]

		Returns: a kinetic rate law for the reaction, with arguments for concentrations and parameters,
			and returns flux.

		'''

		def rate_law(concentrations, parameters):

			# construct numerator
			transporter_concentration = concentrations[transporter]
			kcat_f = parameters[kcat_indices['kcat_f']]

			kcat_r = None
			if 'kcat_r' in kcat_indices:
				kcat_r = parameters[kcat_indices['kcat_r']]

			numerator = 0
			for cofactors in cofactors_sets:
				# if reversible, determine direction by looking at stoichiometry
				if kcat_r:
					coeffs = [stoichiometry[mol] for mol in cofactors]  # coeffs should be of the same sign
					if coeffs[0] > 0:
						kcat = -kcat_r  # use reverse rate
					else:
						kcat = kcat_f
				else:
					kcat = kcat_f

				# multiply the affinities of all cofactors
				term = np.prod([
					self.cofactor_numerator(
						concentrations[molecule],
						parameters[km_indices[molecule]]
					)() for molecule in cofactors
				])
				numerator += kcat * term

			numerator *= transporter_concentration

			# construct denominator, with all competing terms in the partition
			# denominator starts at +1 for the unbound state
			denominator = 1
			for cofactors in partition:
				# multiply the affinities of all cofactors in this partition
				term = np.prod([
					self.cofactor_denominator(
						concentrations[cofactor],
						parameters[km_indices[cofactor]]
					)() for cofactor in cofactors
				])

				denominator += term - 1

			flux = numerator / denominator

			return flux

		return rate_law

	# use rate laws to calculate flux
	def get_fluxes(self, parameters, concentrations):

		reaction_fluxes = {rxn: 0.0 for rxn, specs in self.reactions.iteritems()}
		exchange_fluxes = {mol: 0.0 for mol in concentrations.keys()}

		# loop through all reactions, save reaction flux and molecule flux.
		for rxn, specs in self.reactions.iteritems():
			transporters = specs['transporters']
			stoich = specs['stoichiometry']

			for transporter in transporters:
				flux = self.rate_laws[rxn][transporter](concentrations, parameters)
				reaction_fluxes[rxn] += flux

				for molecule, coefficient in stoich.iteritems():
					if coefficient < 0:
						exchange_fluxes[molecule] -= flux
					elif coefficient > 0:
						exchange_fluxes[molecule] += flux

		return reaction_fluxes, exchange_fluxes

	def run_simulation(self, parameters, run_for):

		timeline = np.arange(0, run_for + self.time_step, self.time_step)

		# initialize concentrations
		# concentrations = self.initialize_state()
		concentrations = self.baseline_concentrations.copy()
		reaction_fluxes, exchange_fluxes = self.get_fluxes(parameters, concentrations)

		# create dicts for saved timeseries, starting with initial concentrations and fluxes at t = 0
		saved_concentrations = {mol: [conc] for mol, conc in concentrations.iteritems()}
		saved_reaction_fluxes = {rxn: [flux] for rxn, flux in reaction_fluxes.iteritems()}
		saved_exchange_fluxes = {molecule: [flux] for molecule, flux in exchange_fluxes.iteritems()}

		for time in timeline[1:]:
			reaction_fluxes, exchange_fluxes = self.get_fluxes(parameters, concentrations)

			for rxn, flux in reaction_fluxes.iteritems():
				concentration_change = flux * self.time_step  # / self.coefficient
				substrates = self.reactions[rxn]['stoichiometry']

				for substrate, coefficient in substrates.iteritems():
					if coefficient < 0:
						concentrations[substrate] -= concentration_change
					elif coefficient > 0:
						concentrations[substrate] += concentration_change

			# save state
			for molecule, timeseries in saved_concentrations.iteritems():
				timeseries.append(concentrations[molecule])

			for rxn, timeseries in saved_reaction_fluxes.iteritems():
				timeseries.append(reaction_fluxes[rxn])

			for molecule, timeseries in saved_exchange_fluxes.iteritems():
				timeseries.append(exchange_fluxes[molecule])

		sim_output = {
			'concentrations': saved_concentrations,
			'reaction_fluxes': saved_reaction_fluxes,
			'exchange_fluxes': saved_exchange_fluxes,
			'time': timeline,
		}

		return sim_output

	def initialize_state(self, set_concentrations={}):
		''' set all initial undefined molecular concentrations to their initial concentrations in the WCM'''

		concentrations = set_concentrations.copy()  # copy.deepcopy(CONDITIONS[0]['initial_concentrations'])

		time_index = int(len(self.wcm_sim_data['time']) / 2)  # get midpoint of timeseries
		cell_volume_fL = self.wcm_sim_data['volume'][time_index]  # [fL]
		cell_volume = cell_volume_fL / 1e15  # convert to L

		initial_concentrations = {}
		for molecule, series in self.wcm_sim_data.iteritems():
			if molecule not in ['time', 'cell_mass', 'volume']:
				# convert counts to molar concentrations
				initial_concentrations[molecule] = series[time_index] / self.avogadro / cell_volume  # [fM]

		# get all substrates in REACTIONS that are not yet set
		for rxn, specs in self.reactions.iteritems():
			# substrates = specs['substrates'].values()
			substrates = specs['stoichiometry'].keys()
			transporters = specs['transporters']

			# loop through substrates
			for substrate in substrates:
				# if substrate is not in concentrations dict
				if substrate not in concentrations.keys() or concentrations[substrate] is None:
					concentrations[substrate] = initial_concentrations[substrate]

			# loop through transporters
			# this is different from the substrate loop in that it looks for the transporter name within the string
			for transporter in transporters:
				# if substrate is not in concentrations dict
				if transporter not in concentrations.keys() or concentrations[transporter] is None:
					transporter_id = [mol_id for mol_id in initial_concentrations.keys() if transporter in mol_id]
					concentrations[transporter] = initial_concentrations[transporter_id[0]]

		return concentrations
