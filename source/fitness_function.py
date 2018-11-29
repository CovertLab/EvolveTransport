from __future__ import absolute_import, division, print_function

import numpy as np



class FitnessFunction(object):
	'''
	FitnessFunction takes a kinetic model, and builds a fitness function around it, which can be used by
	the genetic algorithm class.

	'''

	def __init__(self, config, kinetic_model):

		self.conditions = config.get('conditions')

		self.kinetic_model = kinetic_model
		self.reactions = kinetic_model.reactions
		self.parameter_indices = kinetic_model.parameter_indices
		self.baseline_concentrations = kinetic_model.baseline_concentrations
		self.km_range = kinetic_model.km_range
		self.kcat_range = kinetic_model.kcat_range
		self.genome_size = kinetic_model.n_parameters

		# make mapping from genotype to phenotype
		self.phenotype_transform = self.make_phenotype_transform()

	def evaluate(self, genotype, diagnose_error=False):

		phenotype = self.get_phenotype(genotype)
		diagnosis = {}

		total_error = 0.0

		# loop through all conditions
		for number, condition in enumerate(self.conditions):
			targets = condition['targets']
			penalties = condition['penalties']
			initial_concentrations = condition['initial_concentrations']

			# replace baseline concentrations with initial concentrations for this condition
			concentrations = self.baseline_concentrations.copy()
			for mol_id, init_conc in initial_concentrations.iteritems():
				concentrations[mol_id] = init_conc

			# get condition - dependent fluxes
			reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(phenotype, concentrations)

			error_terms = {term: 0.0 for term, value in penalties.iteritems()}
			# concentration error
			if 'concentrations' in penalties:
				for molecule, target_conc in targets['concentrations'].iteritems():

					# error = penalties['concentrations'] * np.log( abs(concentrations[molecule] - target_conc) + 1)
					# error = penalties['concentrations'] * (concentrations[molecule] - target_conc)**2
					error = penalties['concentrations'] * abs(concentrations[molecule] - target_conc)

					error_terms['concentrations'] += error
					total_error += error

			# reaction flux error
			if 'reaction_fluxes' in penalties:
				for rxn, target_flux in targets['reaction_fluxes'].iteritems():

					# log abs of flux, b/c flux can be negative and also potentially spans many scales
					# error = penalties['reaction_fluxes'] * np.log( abs(reaction_fluxes[rxn] - target_flux) + 1)
					# error = penalties['reaction_fluxes'] * (reaction_fluxes[rxn] - target_flux)**2
					error = penalties['reaction_fluxes'] * abs(reaction_fluxes[rxn] - target_flux)

					error_terms['reaction_fluxes'] += error
					total_error += error

			# molecular flux error
			if 'exchange_fluxes' in penalties:
				for molecule, target_flux in targets['exchange_fluxes'].iteritems():

					# log abs of flux, b/c flux can be negative and also potentially spans many scales
					# error = penalties['exchange_fluxes'] * np.log( abs(exchange_fluxes[molecule] - target_flux) + 1)
					# error = penalties['exchange_fluxes'] * (exchange_fluxes[molecule] - target_flux)**2
					# error = penalties['exchange_fluxes'] * abs(exchange_fluxes[molecule] - target_flux)
					error = penalties['exchange_fluxes'] * (1 - exchange_fluxes[molecule]/target_flux)**2

					error_terms['exchange_fluxes'] += error
					total_error += error

			# parameter error -- use gene values for comparison
			if 'parameters' in penalties:
				for rxn, parameters in targets['parameters'].iteritems():
					for param_id, parameter_target in parameters.iteritems():
						if parameter_target:

							import ipdb
							ipdb.set_trace()
							# TODO -- make sure param_indices points correctly
							index = self.parameter_indices[rxn][param_id]
							# parameter_value = phenotype[index]

							gene_value_target = self.phenotype_transform[index]['pheno_to_geno'](parameter_target)
							gene_value = genotype[index]

							# error = penalties['parameters'] * (gene_value - gene_value_target) ** 2
							error = penalties['parameters'] * abs(gene_value - gene_value_target)

							error_terms['parameters'] += error
							total_error += error

			# steady state error --
			if 'steady_state' in penalties:
				for molecule, steady_state in targets['steady_state'].iteritems():
					# for param_id, steady_state in steady_states.iteritems():

					# TODO -- get derivative at steady_state
					# set concentration to steady state
					steady_state_conc = concentrations.copy()
					steady_state_conc[molecule] = steady_state

					# get fluxes
					ss_reaction_fluxes, ss_exchange_fluxes = self.kinetic_model.get_fluxes(phenotype, steady_state_conc)

					# molecule's current flux
					molecule_flux = ss_exchange_fluxes[molecule]

					# penalize for distance from steady state (flux = 0)
					error_terms['steady_state'] = penalties['steady_state'] * molecule_flux ** 2
					total_error += error_terms['steady_state']

			if diagnose_error:
				diagnosis[number] = {term: error_terms[term] for term in penalties.keys()}

		if diagnose_error:
			diagnosis['total'] = total_error
			return total_error, diagnosis
		else:
			return total_error

	def get_phenotype(self, genotype):
		# convert all genes to param values using geno_to_pheno function

		phenotype = np.empty(self.genome_size)
		for index, gene in enumerate(genotype):
			phenotype[index] = self.phenotype_transform[index]['geno_to_pheno'](gene)

		return phenotype

	def make_phenotype_transform(self):
		# create a list that maps each parameter to an index in the parameter array

		phenotype_transform = [None] * self.genome_size

		for reaction, transporters in self.parameter_indices.iteritems():
			for transporter, params in transporters.iteritems():
				for param_type, indices in params.iteritems():

					# kms
					if param_type is 'kms':
						bounds = self.km_range

						for param, idx in indices.iteritems():
							g_to_p = self.make_genotype_to_phenotype(bounds[0], bounds[1])
							p_to_g = self.make_phenotype_to_genotype(bounds[0], bounds[1])

							mapping = {
								'bounds': bounds,
								'geno_to_pheno': g_to_p,
								'pheno_to_geno': p_to_g,
							}

							phenotype_transform[idx] = mapping

					# kcats (one at a time)
					else:
						bounds = self.kcat_range

						# for param, idx in indices.iteritems():
						g_to_p = self.make_genotype_to_phenotype(bounds[0], bounds[1])
						p_to_g = self.make_phenotype_to_genotype(bounds[0], bounds[1])

						mapping = {
							'bounds': bounds,
							'geno_to_pheno': g_to_p,
							'pheno_to_geno': p_to_g,
						}

						phenotype_transform[indices] = mapping

		return phenotype_transform

	def make_phenotype_to_genotype(self, min_phenotype, max_phenotype):
		a = min_phenotype
		b = max_phenotype / min_phenotype

		def p_to_g(x):
			return np.log(x / a) / np.log(b)

		return p_to_g

	def make_genotype_to_phenotype(self, min_phenotype, max_phenotype):
		a = min_phenotype
		b = max_phenotype / min_phenotype

		def g_to_p(x):
			return a * b ** x

		return g_to_p
