import os
import numpy as np
import matplotlib.pyplot as plt

out_dir = os.path.join(
	os.path.split(__file__)[0],
	'out'
	)
# 	v = exchange flux (M/s/gDCW)
# 	c = concentration (M)
piperno = {
	'ILE[p] v': [1.2824824035087214e-07, 1.0923345329533425e-07, 1.0707756748696281e-07,
				 9.040275407599245e-08, 6.92480888500132e-08, 6.01151597095376e-08,
				 4.585349771652952e-08, 3.180622373151381e-08],
	'PHE[p] c': [5.494395897471078e-05, 2.561190625241187e-05, 1.585982454485426e-05,
				 8.165867714023554e-06, 3.4301223726749635e-06, 1.981888501084805e-06,
				 6.187237910660243e-07, 7.075661815780681e-07],
	'GLY[p] c': [6.070602258792052e-05, 2.992264880671636e-05, 1.4947046162025882e-05,
				 7.328896930821278e-06, 3.4201747691041166e-06, 2.1578581780449557e-06],
	'MET[p] c': [5.430010890910807e-05, 2.6998396377700222e-05, 1.3544091037723712e-05,
				 6.867533937621667e-06, 2.8559055320681584e-06, 7.950364888389618e-07,
				 3.7758015967618875e-07],
	'ILE[p] c': [5.8812890722144566e-05, 3.006568848564886e-05, 1.6831602507482144e-05,
				 8.293128435567832e-06, 3.4400699762458105e-06, 1.8875577775681566e-06,
				 4.297193232199383e-07, 5.106722350378625e-07],
	'GLY[p] v': [9.211728643857689e-08, 8.292432865586342e-08, 8.211594294039347e-08,
				 6.365227599073463e-08, 5.631620420791254e-08, 3.581042001945691e-08],
	'PHE[p] v': [8.010927118558473e-08, 8.094166835794985e-08, 7.542189177883924e-08,
				 7.434995173887732e-08, 6.278214652896286e-08, 5.476346332869303e-08,
				 3.5373068483152454e-08, 2.645395564931458e-08],
	'MET[p] v': [6.160672854380995e-08, 4.8826916462045214e-08, 3.686177881065156e-08,
				 3.043299707765325e-08, 2.2874533467873714e-08, 1.6642988702228408e-08,
				 1.1963041413093294e-08]
}


molecule_ids = [
	'L-Alanine',
	'Glycine',
	'L-Isoleucine',
	'L-Leucine',
	'L-Valine',
	'L-Methionine',
	'L-Phenylalanine',
	'L-Tryptophan',
	'D-Alanine',
]

molecule_concs = {
	'L-Alanine': 2.04, #mmoles/kg, wet weight
	'Glycine': 0.82,
	'Serine': 0.04,
	'Proline': 0.12,
	'L-Isoleucine': 0.30,
	'L-Leucine': 0.22,
	'L-Valine': 0.22,
	'L-Methionine': 0.17,
	'L-Phenylalanine': 0.83,
	'Tyrosine': 0.75,
	'Aspartate': 0.10,
	'Threonine': 0.09,
	'Aspartate plus Glutamate': 0.20,
	'Glutamate': 3.06,
	# 'L-Tryptophan',
	# 'D-Alanine',
}


kinetic_params = {
	'L-Alanine': {'km': 9.2, 'vmax': 1.17},  # km in uM, vmax in mmol/30s/kgWCW
	'Glycine': {'km': 3.8, 'vmax': 1.08},  # km in uM, vmax in mmol/30s/kgWCW
	'L-Isoleucine': {'km': 1.22, 'vmax': 0.96},  # km in uM, vmax in mmol/30s/kgWCW
	'L-Leucine': {'km': 1.07, 'vmax': 1.58},  # km in uM, vmax in mmol/30s/kgWCW
	'L-Valine': {'km': 8.0, 'vmax': 1.41},  # also {'km': 0.7, 'vmax': 0.6} km in uM, vmax in mmol/30s/kgWCW
	'L-Methionine': {'km': 2.27, 'vmax': 0.39},  # km in uM, vmax in mmol/30s/kgWCW
	'L-Phenylalanine': {'km': 0.72, 'vmax': 0.75},  # km in uM, vmax in mmol/30s/kgWCW
	'L-Tryptophan': {'km': 0.9, 'vmax': 0.59},  # km in uM, vmax in mmol/30s/kgWCW
	'D-Alanine': {'km': 8.3, 'vmax': 0.36},  # km in uM, vmax in mmol/30s/kgWCW
}

inhibitor_params = {
	'L-Leucine': {'L-Valine': 5, 'L-Isoleucine': 1},
	'L-Isoleucine': {'L-Leucine': 1, 'L-Valine': 4},
	'L-Alanine': {'Glycine': 9, 'D-Alanine': 14},
	'Glycine': {'L-Alanine': 23},
	'L-Phenylalanine': {'L-Tryptophan': 2},
	'L-Methionine': {},
}


def uM_to_M(conc):
	return conc * 1e-6

def flux(exchange_molecule, concentrations, kinetic_params, inhibitor_params):

	A1 = concentrations[exchange_molecule]
	kinetics = kinetic_params[exchange_molecule]
	km = uM_to_M(kinetics['km']) # convert from (uM) to (M)
	vmax = kinetics['vmax'] * 1e-3/30/1000/0.3 # convert from (mmol/30s/kgWCW) to (M/s/gDCW)
	inhibitors = inhibitor_params.get(exchange_molecule, None)

	inhibition = inhibition_term(concentrations, inhibitors)

	numerator = vmax * A1
	denomenator = A1 + km * inhibition

	flux = numerator / denomenator

	return flux

def inhibition_term(concentrations, inhibitors):

	inhibition = 1
	if inhibitors:
		for inhibitor, k_i in inhibitors.iteritems():
			conc_i = concentrations[inhibitor]
			term = 1 + conc_i/(uM_to_M(k_i)) # convert from (uM) to (M)
			inhibition *= term

	return inhibition

# Initialize concentrations
concentrations = {molecule: 1e-6 for molecule in molecule_ids}
fluxes = {molecule: flux(molecule, concentrations, kinetic_params, inhibitor_params) for molecule in molecule_ids}

# run sweep
plt.figure(figsize=(10, 10))
concentration_sweep = np.linspace(0.01e-6, 60e-6, num=50, endpoint=True)
for mol_idx, molecule in enumerate(molecule_ids):
	conc_test = concentrations.copy()
	flux_sweep = [0] * len(concentration_sweep)

	for conc_idx, conc in enumerate(concentration_sweep):
		conc_test[molecule] = conc
		flux_sweep[conc_idx] = flux(molecule, conc_test, kinetic_params, inhibitor_params)

	plt.plot(concentration_sweep, flux_sweep, linewidth=2, label=molecule)



plt.plot(piperno['ILE[p] c'], piperno['ILE[p] v'], 'rs', linewidth=1, markersize=6, label='isoleucine target')
plt.plot(piperno['GLY[p] c'], piperno['GLY[p] v'], 'r^', linewidth=1, markersize=6, label='glycine target')
plt.plot(piperno['PHE[p] c'], piperno['PHE[p] v'], 'rP', linewidth=1, markersize=6, label='phenylalanine target')
plt.plot(piperno['MET[p] c'], piperno['MET[p] v'], 'ro', linewidth=1, markersize=6, label='methionine target')


plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Concentration (M)')
plt.ylabel('Flux (M/s/gDCW)')
# 	plt.title('uptake')
# 	plt.ylabel('exchange flux (M/s/gDCW)')
# 	plt.xlabel('concentration (M)')
# 	plt.ticklabel_format(style='sci', axis='x', scilimits=(0 ,0))


fig_name = ('conc_sweep_inhibition')
plt.savefig(os.path.join(out_dir, fig_name), bbox_inches='tight')



# make lineweaver-burk plot
n_plots = len(kinetic_params)
plt.figure(figsize=(5, 2*n_plots))

plot_number = 1
for molecule_id, params in kinetic_params.iteritems():

	# km = uM_to_M(params['km']) # convert from (uM) to (M)
	# vmax = params['vmax'] * 1e-3/30/1000/0.3 # convert from (mmol/30s/kgWCW) to (M/s/gDCW)
	km = params['km']
	vmax = params['vmax']

	x_intercept = -1/km
	y_intercept = 1/vmax
	slope = km/vmax

	x_end = 50
	y_end = y_intercept + slope * x_end

	x = [x_intercept, 0, x_end]
	y = [0, y_intercept, y_end]

	plt.subplot(n_plots, 1, plot_number)
	plt.plot(x, y, label=molecule_id)
	plt.axvline(x=0.0, color='k', linestyle='--')
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.xlabel('(1 / [S])')
	plt.ylabel('(1 / v)')

	plot_number +=1

plt.subplots_adjust(hspace=1.0, wspace=1.0)

fig_name = ('lineweaver_burk')
plt.savefig(os.path.join(out_dir, fig_name), bbox_inches='tight')

#
# import ipdb;
# ipdb.set_trace()