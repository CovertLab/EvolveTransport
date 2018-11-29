from __future__ import absolute_import, division, print_function

import os
import csv

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from source.kinetic_flux_model import KineticFluxModel



class Visualize(object):
	'''
	Visualize class for plotting the evolution of parameters by a genetic algorithm.

	To initialize this object, pass a config and a fitness function constructed in evaluation.py
	'''

	def __init__(self, config, fitness_function):

		self.out_dir = config.get('out_dir', None)
		self.parameter_out_dir = config.get('parameter_out_dir', None)
		self.saved_param_file = config.get('saved_param_file', None)
		self.parameter_analytics = config.get('parameter_analytics', False)
		self.mutation_variance = config.get('mutation_variance', 1.0)
		self.seed = config.get('seed', None)
		self.replicate_id = config.get('replicate_id', None)
		self.save_fitness_threshold = config.get('fitness_threshold')
		self.wcm_sim_data = config.get('wcm_sim_data', None)
		self.exchange_molecules = config.get('exchange_molecules', None)

		self.fitness_function = fitness_function
		self.kinetic_model = fitness_function.kinetic_model
		self.reactions = fitness_function.reactions
		self.parameter_indices = fitness_function.parameter_indices
		self.km_range = fitness_function.km_range
		self.kcat_range = fitness_function.kcat_range

		# for michaelis-menten plot
		self.transport_configuration = fitness_function.kinetic_model.transport_configuration


	def parameters(self, best_parameters):

		best_param_array = np.array(best_parameters, np.float64)

		pca = PCA(n_components=2)
		pca.fit(best_param_array)
		best_params_reduced = pca.transform(best_param_array)

		# todo -- cluster analysis
		plt.figure(figsize=(8.5, 8.5))
		plt.scatter(best_params_reduced[:, 0], best_params_reduced[:, 1])
		plt.xlabel('PC1')
		plt.ylabel('PC2')

		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)
		fig_name = ('param_space_' + self.replicate_id)
		plt.savefig(os.path.join(self.out_dir, fig_name))

		print('parameters plot saved')

	# Plots
	def conditions(self, phenotype, conditions):

		run_for = 1.0
		rows = 4
		n_conditions = len(conditions)

		plt.figure(figsize=(5 * n_conditions + 2, 2 * rows + 1))

		for condition_index, cond in enumerate(conditions):

			# run simulation
			set_concentrations = cond['initial_concentrations']
			initial_concentrations = self.kinetic_model.initialize_state(set_concentrations)
			sim_output = self.fitness_function.kinetic_model.run_simulation(
				phenotype,
				run_for,
			)

			concentrations_timeseries = sim_output['concentrations']
			reaction_fluxes_timeseries = sim_output['reaction_fluxes']
			exchange_fluxes_timeseries = sim_output['exchange_fluxes']
			time_timeseries = sim_output['time']

			# look at all targets for this condition
			condition_targets = cond['targets']

			row_index = 0
			for target_type, targets in condition_targets.iteritems():

				if target_type == 'exchange_fluxes':
					for molecule, target_flux in targets.iteritems():
						exchange_flux = exchange_fluxes_timeseries[molecule]
						concentration = concentrations_timeseries[molecule]

						plt.subplot(rows, n_conditions, n_conditions * row_index + condition_index + 1)  # (n_conditions, columns, columns * condition_index + index)
						plt.plot(time_timeseries, exchange_flux, 'g', label=('ex flux, avg = %.2e' % np.mean(exchange_flux)))
						plt.axhline(y=target_flux, linestyle='--', color='r', label=('target = %.2e' % target_flux))
						plt.title('cond_' + str(condition_index) + ' ex_flux target: ' + molecule, y=1.15)
						plt.ylabel('flux (M/s)')
						plt.xlabel('t (s)')
						plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

						row_index += 1

						plt.subplot(rows, n_conditions, n_conditions * row_index + condition_index + 1)
						plt.plot(time_timeseries, concentration, 'c')
						plt.title(molecule + ' concentration')
						plt.ylabel('conc (M)')
						plt.xlabel('t (s)')

						row_index += 1

				if target_type == 'reaction_fluxes':

					for reaction, target_flux in targets.iteritems():
						reaction_flux = reaction_fluxes_timeseries[reaction]

						plt.subplot(rows, n_conditions, n_conditions * row_index + condition_index + 1)
						plt.plot(time_timeseries, reaction_flux, 'b', label='reaction flux')
						plt.axhline(y=target_flux, linestyle='--', color='r', label=('target flux = %.2e' % target_flux))
						plt.title('cond_' + str(condition_index) + ' rxn_flux target: ' + reaction, y=1.15)
						plt.ylabel('flux (M/s)')
						plt.xlabel('t (s)')
						# plt.yscale('log')
						plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

						row_index += 1

					for molecule, init_conc in set_concentrations.iteritems():
						concentration = concentrations_timeseries[molecule]

						plt.subplot(rows, n_conditions, n_conditions * row_index + condition_index + 1)
						plt.plot(time_timeseries, concentration, 'c')
						plt.title('cond_' + str(condition_index) + ': ' + molecule + ' (M)')
						plt.ylabel('conc')
						plt.xlabel('t (s)')

						row_index += 1


				# if target_type == 'parameters':
				#
				# 	# conc_range = PARAM_RANGES['km']
				# 	n_sample_concentrations = 1000000
				# 	conc_range = np.linspace(self.km_range[0], self.km_range[1], n_sample_concentrations)
				# 	flux_range = np.empty_like(conc_range)
				#
				# 	# show michaelis-menten for each reaction
				# 	for reaction_id, parameters in targets.iteritems():
				#
				# 		# create new reaction and parameter dicts
				# 		reaction = {}
				#
				# 		# configure transport model to this reaction alone
				# 		reaction[reaction_id] = self.reactions[reaction_id]
				#
				# 		kinetic_model = KineticFluxModel(self.kinetic_model.config, reaction)
				#
				# 		reactants = [mol for mol, coeff in reaction[reaction_id]['stoichiometry'].iteritems() if coeff < 0]
				# 		a1 = reactants[0]
				# 		concentrations = self.kinetic_model.initialize_state()
				#
				# 		for idx, conc in enumerate(conc_range):
				# 			concentrations[a1] = conc
				# 			reaction_fluxes, exchange_fluxes = kinetic_model.get_fluxes(phenotype, concentrations)
				# 			flux_range[idx] = reaction_fluxes[reaction_id]
				#
				# 		# plot M-M curve for this reaction
				# 		plt.subplot(rows, n_conditions, n_conditions * row_index + condition_index + 1)
				# 		plt.plot(conc_range, flux_range)
				#
				# 		# put target parameters on top of M-M
				# 		for param, value in parameters.iteritems():
				# 			if 'kcat' in param:
				# 				transporter_id = reaction[reaction_id]['transporters']
				# 				t_conc = concentrations[transporter_id[0]]  # uses concentration of first transporter
				# 				vmax = value * t_conc
				# 				plt.axhline(
				# 					y=vmax, linestyle='--', color='r',
				# 					label='target vmax, w/ kcat = ' + str(value)
				# 				)
				# 			else:
				# 				plt.axvline(x=value, linestyle='--', color='k', label=('target km = ' + str(value)))
				#
				# 		plt.xscale('log')
				# 		plt.xlabel(a1 + ' concentration (M)')
				# 		plt.ylabel('flux (M/s)')
				# 		plt.title('cond_' + str(condition_index) + ' param targets, ' + reaction_id, y=1.15)
				# 		# plt.title(reaction_id)
				# 		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
				#
				# 		row_index += 1


		# if PIPERNO:
		#
		# 	piperno = {
		# 		'ILE[p] v': [1.2824824035087214e-07, 1.0923345329533425e-07, 1.0707756748696281e-07,
		# 					 9.040275407599245e-08, 6.92480888500132e-08, 6.01151597095376e-08,
		# 					 4.585349771652952e-08, 3.180622373151381e-08],
		# 		'PHE[p] c': [5.494395897471078e-05, 2.561190625241187e-05, 1.585982454485426e-05,
		# 					 8.165867714023554e-06, 3.4301223726749635e-06, 1.981888501084805e-06,
		# 					 6.187237910660243e-07, 7.075661815780681e-07],
		# 		'GLY[p] c': [6.070602258792052e-05, 2.992264880671636e-05, 1.4947046162025882e-05,
		# 					 7.328896930821278e-06, 3.4201747691041166e-06, 2.1578581780449557e-06],
		# 		'MET[p] c': [5.430010890910807e-05, 2.6998396377700222e-05, 1.3544091037723712e-05,
		# 					 6.867533937621667e-06, 2.8559055320681584e-06, 7.950364888389618e-07,
		# 					 3.7758015967618875e-07],
		# 		'ILE[p] c': [5.8812890722144566e-05, 3.006568848564886e-05, 1.6831602507482144e-05,
		# 					 8.293128435567832e-06, 3.4400699762458105e-06, 1.8875577775681566e-06,
		# 					 4.297193232199383e-07, 5.106722350378625e-07],
		# 		'GLY[p] v': [9.211728643857689e-08, 8.292432865586342e-08, 8.211594294039347e-08,
		# 					 6.365227599073463e-08, 5.631620420791254e-08, 3.581042001945691e-08],
		# 		'PHE[p] v': [8.010927118558473e-08, 8.094166835794985e-08, 7.542189177883924e-08,
		# 					 7.434995173887732e-08, 6.278214652896286e-08, 5.476346332869303e-08,
		# 					 3.5373068483152454e-08, 2.645395564931458e-08],
		# 		'MET[p] v': [6.160672854380995e-08, 4.8826916462045214e-08, 3.686177881065156e-08,
		# 					 3.043299707765325e-08, 2.2874533467873714e-08, 1.6642988702228408e-08,
		# 					 1.1963041413093294e-08]
		# 	}
		#
		# 	plt.subplot(rows, n_conditions, n_conditions * (rows -1) + 1)
		# 	plt.plot(piperno['ILE[p] c'], piperno['ILE[p] v'], 'rs', linewidth=1, markersize=3, label='isoleucine target')
		# 	plt.plot(piperno['GLY[p] c'], piperno['GLY[p] v'], 'r^', linewidth=1, markersize=6  )  # , label='glycine target')
		# 	plt.plot(piperno['PHE[p] c'], piperno['PHE[p] v'], 'rP', linewidth=1, markersize=3, label='phenylalanine target')
		# 	plt.plot(piperno['MET[p] c'], piperno['MET[p] v'], 'ro', linewidth=1, markersize=3, label='methionine target')
		#
		# 	plt.title('uptake')
		# 	plt.ylabel('exchange flux (M/s/gDCW)')
		# 	plt.xlabel('concentration (M)')
		# 	plt.ticklabel_format(style='sci', axis='x', scilimits=(0 ,0))

		plt.subplots_adjust(hspace=2.0, wspace=2.0)
		# plt.tight_layout()

		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)
		fig_name = ('conditions_' + self.replicate_id)
		plt.savefig(os.path.join(self.out_dir, fig_name))  # format='pdf'

		print('condition plots saved')

	def sim_out(self, sim_output, parameters, conditions):

		concentrations_timeseries = sim_output['concentrations']
		reaction_fluxes_timeseries = sim_output['reaction_fluxes']
		exchange_fluxes_timeseries = sim_output['exchange_fluxes']
		time_timeseries = sim_output['time']

		# all reaction targets, for adding to parameter plot
		# TODO -- this should be in parameters plot, not sim out
		target_params = {}
		for specs in conditions:
			if 'targets' in specs and 'parameters' in specs['targets']:
				rxns = specs['targets']['parameters']
				for rxn, params in rxns.iteritems():
					target_params[rxn] = params

		columns = 4
		rows = max([len(concentrations_timeseries), len(reaction_fluxes_timeseries), len(parameters)])

		plt.figure(figsize=(4*columns, 1*rows))
		# plot concentrations over time
		index = 1
		for molecule, timeseries in concentrations_timeseries.iteritems():
			plt.subplot(rows, columns, columns * index - (columns - 1))
			plt.plot(time_timeseries, timeseries, 'c')
			plt.title(molecule)
			if index < len(concentrations_timeseries):
				plt.tick_params(
					axis='x',  # changes apply to the x-axis
					which='both',  # both major and minor ticks are affected
					bottom=False,  # ticks along the bottom edge are off
					top=False,  # ticks along the top edge are off
					labelbottom=False)  # labels along the bottom edge are off
			else:
				plt.ylabel('Conc (M)')
				plt.xlabel('Time (s)')
			index += 1

		# plot reaction fluxes over time
		index = 1
		for reaction, timeseries in reaction_fluxes_timeseries.iteritems():
			plt.subplot(rows, columns, columns * index - (columns - 2))
			plt.plot(time_timeseries, timeseries, 'b')
			plt.title(reaction)
			if index < len(reaction_fluxes_timeseries):
				plt.tick_params(
					axis='x',  # changes apply to the x-axis
					which='both',  # both major and minor ticks are affected
					bottom=False,  # ticks along the bottom edge are off
					top=False,  # ticks along the top edge are off
					labelbottom=False)  # labels along the bottom edge are off
			else:
				plt.ylabel('Flux (M/s)')
				plt.xlabel('Time (s)')
			index += 1

		if 'TARGET_FLUXES' in globals():
			# plot time series of exchange fluxes in TARGET_FLUXES
			index = 1
			for molecule, timeseries in exchange_fluxes_timeseries.iteritems():
				if molecule in TARGET_FLUXES:
					plt.subplot(rows, columns, columns * index -(columns -3))
					plt.plot(time_timeseries, timeseries, 'g')
					plt.title(molecule)
					if index < len(TARGET_FLUXES):
						plt.tick_params(
							axis='x',  # changes apply to the x-axis
							which='both',  # both major and minor ticks are affected
							bottom=False,  # ticks along the bottom edge are off
							top=False,  # ticks along the top edge are off
							labelbottom=False)  # labels along the bottom edge are off
					else:
						plt.ylabel('Flux (M/s)')
						plt.xlabel('Time (s)')
					index += 1

		# plot parameters
		for rxn, params in self.parameter_indices.iteritems():
			for transporter, param_indices in params.iteritems():
				for param_type, params in param_indices.iteritems():

					# plot kms
					if 'km' in param_type:
						for param, param_idx in params.iteritems():
							bounds = self.km_range

							param_value = parameters[param_idx]

							plt.subplot(rows, columns, columns * (param_idx + 1) - (columns - 3))
							plt.axvline(x=bounds[0])
							plt.axvline(x=bounds[1])
							plt.axhline(y=0.5)

							# plot target parameters if defined in target_params
							if rxn in target_params and param in target_params[rxn]:
								target_value = target_params[rxn][param]
								plt.axvline(x=target_value, linewidth=4, color='r')

							plt.plot(param_value, 0.5, 'bo', markersize=10)

							info = (transporter + ': km_' + param)
							plt.title(info)
							plt.ylim(0, 1)
							plt.xscale('log')
							plt.tick_params(
								left=False,
								right=False,
								bottom=False,
								top=False,
								labelleft=False)

					# plot kcats
					else:
						bounds = self.kcat_range
						param_value = parameters[params]

						plt.subplot(rows, columns, columns * (params + 1) - (columns - 3))
						plt.axvline(x=bounds[0])
						plt.axvline(x=bounds[1])
						plt.axhline(y=0.5)

						# import ipdb; ipdb.set_trace()
						# plot target parameters if defined in target_params
						if rxn in target_params and param_type in target_params[rxn]:
							target_value = target_params[rxn][param_type]
							plt.axvline(x=target_value, linewidth=4, color='r')

						plt.plot(param_value, 0.5, 'bo', markersize=10)

						info = (rxn + ': ' + param_type)
						plt.title(info)
						plt.ylim(0, 1)
						plt.xscale('log')
						plt.tick_params(
							left=False,
							right=False,
							bottom=False,
							top=False,
							labelleft=False)

		plt.subplots_adjust(hspace=1.5, wspace=1.0)
		plt.tight_layout()

		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)
		fig_name = ('sim_' + self.replicate_id)
		plt.savefig(os.path.join(self.out_dir, fig_name))

		print('top simulation plot saved')

	def michaelis_menten(self, all_parameters):

		set_parameters = False
		set_concentrations = False

		test_transporter = True
		test_cofactor = True
		test_competitor = True

		columns = 1 + sum([test_transporter, test_cofactor, test_competitor])

		n_vary = 10
		n_samples = 100
		n_rxns = len(self.parameter_indices)
		rows = 2*n_rxns + 2  # extra row for each reaction header

		cmap = plt.cm.get_cmap('Spectral')
		colors = [cmap(float(idx) / n_vary) for idx in range(n_vary)]

		baseline_concentrations = self.kinetic_model.baseline_concentrations

		if set_concentrations:
			new_concs = {
				# 'CYCA-MONOMER': 1e0,
				# 'L-ALPHA-ALANINE[p]': 1e-2,
				# 'L-ALPHA-ALANINE[c]': 1e-5,
				# 'GLY[p]': 1e-2,
				# 'GLY[c]': 1e-2,
				# 'PROTON[p]': 1e-2,
				# 'PROTON[c]': 1e-2,
			}
			# set new baseline_concentrations
			baseline_concentrations.update(new_concs)

		if set_parameters:
			new_params = {
				'RXN0-5202': {
					'CYCA-MONOMER': {
						# 'kcat_f': 1e0,
						'kms': {
							# 'L-ALPHA-ALANINE[p]': 1e-3,
							# 'GLY[p]': 1e-3,
							# 'PROTON[p]': 1,
						}
					}
				},
				'TRANS-RXN-62B': {
					'CYCA-MONOMER': {
						# 'kcat_f': 1e2,
					}
				}
			}

			for rxn, trps in new_params.iteritems():
				for trp, pars in trps.iteritems():
					for par, par_type in pars.iteritems():
						if 'km' in par:
							for mol, val in par_type.iteritems():
								idx = self.parameter_indices[rxn][trp][par][mol]
								all_parameters[idx] = val
						else:
							idx = self.parameter_indices[rxn][trp][par]
							all_parameters[idx] = par_type

		plt.figure(figsize=(6*columns, 3*rows))
		plot_number = 1
		row_number = 0
		for reaction_id, specs in self.reactions.iteritems():
			transporters = specs['transporters']
			stoich = specs['stoichiometry']
			parameters = self.parameter_indices[reaction_id]

			# TODO -- set a1 to amino acid... or show all?
			reactants = [mol for mol, coeff in stoich.iteritems() if coeff < 0]
			products = [mol for mol, coeff in stoich.iteritems() if coeff > 0]

			a1_set = False
			for mol in self.exchange_molecules:
				if mol in reactants:
					a1 = mol
					a1_set = True

			if not a1_set:
				a1 = reactants[0]

			# get cofactor
			b1 = None
			if len(reactants) > 1:
				cofactors = [x for x in reactants if x != a1]
				b1 = cofactors[0]

			# plot info in whole row
			param_values = {}
			for trans, params in self.parameter_indices[reaction_id].iteritems():
				param_values[trans] = {}
				param_values[trans]['kms'] = {}
				for param_type, params in params.iteritems():
					if 'km' in param_type:
						for param, idx in params.iteritems():
							param_values[trans]['kms'][param] = all_parameters[idx]
					else:
						param_values[trans][param_type] = all_parameters[params]

			plt.subplot(rows, columns, plot_number)
			plt.text(0.02, 0.6, 'reaction: ' + reaction_id, weight='bold')
			plt.text(0.02, 0.45, 'reactants: %s' % reactants)
			plt.text(0.02, 0.3, 'products: %s' % products)
			plt.text(0.02, 0.15, 'transporters: %s' % transporters)
			plt.text(0.02, 0.0, 'parameters: %s' % param_values[transporters[0]])
			plt.axis('off')
			plot_number += columns
			row_number += 1

			# michaelis menten by sampling substrate concentrations
			for transporter in transporters:

				concentrations = baseline_concentrations.copy()
				conc_values = np.logspace(-9, 0, num=n_samples, endpoint=True, base=10)

				flux_values = np.empty_like(conc_values)
				for idx, conc in enumerate(conc_values):
					concentrations[a1] = conc
					reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
					flux_values[idx] = reaction_fluxes[reaction_id]

				# plot M-M curve for this reaction
				plt.subplot(rows, columns, plot_number)
				plt.plot(conc_values, flux_values)

				# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
				plt.xscale('log')
				plt.xlabel(a1 + ' concentration (M)')
				plt.ylabel('flux (M/s)')
				plt.title('transporter: %s' % transporter)

				plot_number += 1

				if test_transporter:
					concentrations = baseline_concentrations.copy()
					conc_values = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
					transporter_concs = np.logspace(-4, 1, num=n_vary, endpoint=True, base=10)

					plt.subplot(rows, columns, plot_number)
					for index, transporter_conc in enumerate(transporter_concs):
						concentrations[transporter] = transporter_conc

						flux_values = np.empty_like(conc_values)
						for idx, conc in enumerate(conc_values):

							concentrations[a1] = conc
							reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
							flux_values[idx] = reaction_fluxes[reaction_id]

						# plot M-M curve for this reaction
						plt.plot(conc_values, flux_values,
												color = colors[index],
												label = ('conc = %.2e' % (transporter_conc)),
												)

					plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
					plt.xscale('log')
					plt.xlabel(a1 + ' concentration (M)')
					plt.ylabel('flux (M/s)')
					plt.title('transporter: %s' % transporter)

					plot_number += 1

				if test_cofactor:

					concentrations = baseline_concentrations.copy()
					conc_values = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
					cofactor_concs = np.logspace(-8, 1, num=n_vary, endpoint=True, base=10)

					if b1 is not None:
						plt.subplot(rows, columns, plot_number)
						for index, cofactor_conc in enumerate(cofactor_concs):
							concentrations[b1] = cofactor_conc

							flux_values = np.empty_like(conc_values)
							for idx, conc in enumerate(conc_values):

								concentrations[a1] = conc
								reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
								flux_values[idx] = reaction_fluxes[reaction_id]

							# plot M-M curve for this reaction
							plt.plot(conc_values, flux_values,
											color = colors[index],
											label = ('conc = %.2e' % (cofactor_conc)),
											)

						plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
						plt.xscale('log')
						plt.xlabel(a1 + ' concentration (M)')
						plt.ylabel('flux (M/s)')
						plt.title('cofactor: %s' % b1)

					plot_number += 1

				if test_competitor:
					# get competitor
					rxns_transporter = self.transport_configuration[transporter]['reaction_cofactors'].keys()
					competing_rxns = [trpr for trpr in rxns_transporter if trpr not in reaction_id]
					competitor = None
					for rx in competing_rxns:
						competitors = self.transport_configuration[transporter]['reaction_cofactors'][rx]
						competitor = competitors[0][0]

					if competitor is not None:
						concentrations = baseline_concentrations.copy()
						conc_values = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
						competitor_concs = np.logspace(-8, 1, num=n_vary, endpoint=True, base=10)

						plt.subplot(rows, columns, plot_number)
						for index, competitor_conc in enumerate(competitor_concs):
							concentrations[competitor] = competitor_conc

							flux_values = np.empty_like(conc_values)
							for idx, conc in enumerate(conc_values):

								concentrations[a1] = conc
								reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
								flux_values[idx] = reaction_fluxes[reaction_id]

							# plot M-M curve for this reaction
							plt.plot(conc_values, flux_values,
													color = colors[index],
													label = ('conc = %.2e' % (competitor_conc)),
													)

						plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
						plt.xscale('log')
						plt.xlabel(a1 + ' concentration (M)')
						plt.ylabel('flux (M/s)')
						plt.title('competitor: %s' % (competitor))

					plot_number += 1

				row_number += 1

			plot_number = row_number * columns + 1

		plt.subplots_adjust(hspace=0.5, wspace=1.5)

		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)
		fig_name = ('MM_' + self.replicate_id)
		plt.savefig(os.path.join(self.out_dir, fig_name), bbox_inches='tight')

		print('michaelis-menten plot saved')

	def evolution(self, saved_error, saved_fitness, saved_diagnosis):

		penality_diagnosis = saved_diagnosis['penality_diagnosis']
		mutation_diagnosis = saved_diagnosis['mutation_diagnosis']
		n_bins = 10
		max_shown_gens = 20
		max_gens_labels = 10

		n_saved_gens, population_size = np.shape(saved_error)

		# define shown generations for fitness distribution plots
		if n_saved_gens >= max_shown_gens:
			nth_gen_shown = int(n_saved_gens / max_shown_gens)
		else:
			nth_gen_shown = 1
		shown_gens = [gen for gen in xrange(0, n_saved_gens, nth_gen_shown)]

		# make labels
		if n_saved_gens >= max_gens_labels:
			nth_gen_label = int(n_saved_gens / max_gens_labels)
		else:
			nth_gen_label = 1
		labeled_gens = [gen for gen in xrange(0, n_saved_gens, nth_gen_label)]
		gen_label = ['gen ' + str(gen) for gen in labeled_gens]

		hist_fit_by_gen = [saved_fitness[gen] for gen in labeled_gens]
		plot_fit_by_gen = [saved_fitness[gen] for gen in shown_gens]
		plot_error_by_gen = [saved_error[gen] for gen in shown_gens]

		# fitness_ordered_gens: a 2d array of [gens, population], with population ordered by rank of fitness
		error_ordered_gens = np.array([sorted(gen, reverse=False) for gen in plot_error_by_gen])
		fitness_ordered_gens = np.array([sorted(gen, reverse=True) for gen in plot_fit_by_gen])
		selection_fitness_ordered_gens = np.array([fitness / np.sum(fitness) for fitness in fitness_ordered_gens])

		top_fitness = [max(fit) for fit in saved_fitness]
		avg_fitness = [sum(fit) / len(fit) for fit in saved_fitness]

		plt.figure(figsize=(20, 20))

		# plot fitness over gens
		plt.subplot(3, 2, 1)
		plt.plot(top_fitness, linewidth=3, label='top fitness')
		plt.plot(avg_fitness, 'r', linewidth=3, label='mean fitness')
		plt.title('seed = ' + str(self.seed) + ' mutation variance = ' + str(self.mutation_variance))
		plt.ylabel('Fitness (1/1+error)')
		plt.xlabel('Generation')
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

		# plot histograms of fitness over several generations
		plt.subplot(3, 2, 3)
		cmap = plt.cm.get_cmap('Spectral')
		colors = [cmap(float(idx) / len(hist_fit_by_gen)) for idx in range(len(hist_fit_by_gen))]

		plt.hist(hist_fit_by_gen, bins=n_bins, color=colors, label=gen_label)
		plt.title('Fitness distribution for population with n=' + str(population_size))
		plt.ylabel('Counts')
		plt.xlabel('Fitness: 1/(1+error)')
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

		# plot relative fitness in a stacked bar graph for all gens in shown_gens
		plt.subplot(3, 2, 2)
		cmap = plt.cm.get_cmap('Spectral')
		colors = [cmap(float(idx) / len(shown_gens)) for idx in range(len(shown_gens))]

		for idx, selection_fitness in enumerate(selection_fitness_ordered_gens):
			cum_fitness = np.cumsum(selection_fitness)

			if shown_gens[idx] in labeled_gens:
				plt.plot(cum_fitness, linewidth=3, color=colors[idx], label=('gen ' + str(shown_gens[idx])))
			else:
				plt.plot(cum_fitness, linewidth=3, color=colors[idx])

		plt.plot([0, population_size], [1./population_size, 1.], 'k--', label='equal selection')
		plt.title('Selection fitness')
		plt.ylabel('Cumulative selection fitness')
		plt.xlabel('Individuals, ordered by rank (most to least fitness)')
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

		# plot error for a error-ordered population at each gen in shown_gens
		plt.subplot(3, 2, 4)
		cmap = plt.cm.get_cmap('Spectral')
		colors = [cmap(float(idx) / len(shown_gens)) for idx in range(len(shown_gens))]

		for idx, gen in enumerate(error_ordered_gens):
			if shown_gens[idx] in labeled_gens:
				plt.plot(gen, linewidth=3, color=colors[idx], label=('gen ' + str(shown_gens[idx])))
			else:
				plt.plot(gen, linewidth=3, color=colors[idx])

		plt.title('Total error for population')
		plt.ylabel('Total error')
		plt.xlabel('Individuals, ordered by rank (least to most error)')
		# plt.yscale('log')
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

		if penality_diagnosis:
			total_error = [diagnosis['total'] for generation, diagnosis in enumerate(penality_diagnosis)]

			error_contributions = {
				condition: {} for condition, terms in penality_diagnosis[0].iteritems() if condition != 'total'}
			for condition in error_contributions:
				error_contributions[condition] = {term: [] for term, value in penality_diagnosis[0][condition].iteritems()}

			# go through each generation's diagnosis, and concatenate errors from each condition and term
			for time, diagnosis in enumerate(penality_diagnosis):
				for condition, terms in diagnosis.iteritems():
					if condition != 'total':
						for term, value in terms.iteritems():
							error_contributions[condition][term].append(value)

			plt.subplot(3, 2, 6)
			plt.plot(total_error, linewidth=3, label='total error')
			for condition, terms in error_contributions.iteritems():
				for term, series in terms.iteritems():
					term_label = ('condition ' + str(condition) + ': ' + term)
					plt.plot(series, linewidth=3, label=term_label)

			plt.title('Error contributions to top individual')
			plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			plt.ylabel('Error contribution')
			plt.xlabel('Generation')
		# plt.yscale('log')

		plt.subplots_adjust(hspace=0.5, wspace=1.0)

		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)
		fig_name = ('GA_' + self.replicate_id)
		plt.savefig(os.path.join(self.out_dir, fig_name), bbox_inches='tight')  # format='pdf'

		print('evolution plot saved')

		# if mutation_diagnosis:
		#
		# 	plt.figure(figsize=(10, 20))
		#
		# 	p_enforce_bounds = np.mean \
		# 		([mutation['enforce_bounds'] for mutation in mutation_diagnosis])  # / float(len(mutation_diagnosis))
		# 	acceptance_p = np.array \
		# 		([[mutation['fitness_change'], mutation['acceptance_p']] for mutation in mutation_diagnosis if (mutation['fitness_change'] < 0)])
		# 	fitness_change = [mutation['fitness_change'] for mutation in mutation_diagnosis]
		#
		# 	# analyze mutation directions
		# 	# mutation_directions = np.array([mutation['mutation'] for mutation in mutation_diagnosis])
		# 	good_mutations = np.array \
		# 		([mutation['mutation'] for mutation in mutation_diagnosis if mutation['fitness_change' ] >0])
		# 	bad_mutations = np.array \
		# 		([mutation['mutation'] for mutation in mutation_diagnosis if mutation['fitness_change' ] <0])
		# 	# TODO -- learn labeled data (bad vs good directions)

			# TODO -- fitness change needs some range for histogram to work

			#
			# plt.subplot(2, 1, 1)
			# plt.hist(fitness_change, bins=50)
			# plt.title('mutation variance = ' + str(self.mutation_variance) + '. probability of enforce bounds = ' + str
			# 	(p_enforce_bounds))
			# plt.xlabel('fitness changes resulting from mutation')
			# plt.ylabel('counts')
			#
			# plt.subplot(2, 1, 2)
			# plt.scatter(acceptance_p[: ,0], acceptance_p[: ,1])
			# # plt.hist(acceptance_p, bins=50) #, bins=np.logspace(np.log10(0.1),np.log10(1.0), 50))
			# plt.xlabel('fitness')
			# plt.ylabel('acceptance probability')
			#
			#
			# # pca = PCA(n_components=2)
			# # pca.fit(mutation_directions)
			# #
			# # good_mutations_reduced = pca.transform(good_mutations)
			# # bad_mutations_reduced = pca.transform(bad_mutations)
			# #
			# # plt.subplot(2, 1, 2)
			# # plt.scatter(good_mutations_reduced[:, 0], good_mutations_reduced[:, 1])
			# # plt.scatter(bad_mutations_reduced[:, 0], bad_mutations_reduced[:, 1])
			# # plt.xlabel('PC1')
			# # plt.ylabel('PC2')
			# #
			# # plt.subplots_adjust(hspace=0.5)
			#
			# if not os.path.exists(self.out_dir):
			# 	os.mkdir(self.out_dir)
			# fig_name = ('mutation_' + self.replicate_id)
			# plt.savefig(os.path.join(self.out_dir ,fig_name), bbox_inches='tight')
			#
			# print('mutation plot saved')
