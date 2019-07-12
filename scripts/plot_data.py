import os
import csv
import matplotlib.pyplot as plt
import numpy as np

# OUTDIR = os.path.join(
# 	os.path.split(__file__)[0],
# 	'out/plot_out'
# 	)
OUTDIR = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'out/plot_out'))

GLYCINE_DATA = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'data/piperno_oxender_1968/glycine.csv'))
ISOLEUCINE_DATA = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'data/piperno_oxender_1968/isoleucine.csv'))
METHIONINE_DATA = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'data/piperno_oxender_1968/methionine.csv'))
PHENYLALANINE_DATA = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'data/piperno_oxender_1968/phenylalanine.csv'))

EXCHANGE_FLUX_PARAM_DATA = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'data/piperno_oxender_1968/exchange_flux_params.csv'))

with open(GLYCINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)
	glycine_concs = []
	glycine_uptake = []
	for row in reader:
		glycine_concs.append(float(row[0]))
		glycine_uptake.append(float(row[1]))

with open(ISOLEUCINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)
	isoleucine_concs = []
	isoleucine_uptake = []
	for row in reader:
		isoleucine_concs.append(float(row[0]))
		isoleucine_uptake.append(float(row[1]))

with open(METHIONINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)
	methionine_concs = []
	methionine_uptake = []
	for row in reader:
		methionine_concs.append(float(row[0]))
		methionine_uptake.append(float(row[1]))

with open(PHENYLALANINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)
	phenylalanine_concs = []
	phenylalanine_uptake = []
	for row in reader:
		phenylalanine_concs.append(float(row[0]))
		phenylalanine_uptake.append(float(row[1]))

with open(EXCHANGE_FLUX_PARAM_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)
	exchange_param = {}
	for row in reader:
		exchange_param[row[0]] = {'km' : float(row[1]), 'vmax' : float(row[2])}


def get_flux(conc, km, vmax):
	return vmax * conc / (conc + km)



plt.figure(figsize=(8, 8))
#
# # convert units
# isoleucine_concs = [conc * 1e-6 for conc in isoleucine_concs]  # convert to M
# glycine_concs = [conc * 1e-6 for conc in glycine_concs]  # convert to M
# phenylalanine_concs = [conc * 1e-6 for conc in phenylalanine_concs]  # convert to M
# methionine_concs = [conc * 1e-6 for conc in methionine_concs]  # convert to M
#
# isoleucine_uptake = [flux * 1e-6 * (1.0 / 30) * (1.0 / 0.3) for flux in isoleucine_uptake]  # to mol AA / s / g dry cell mass
# glycine_uptake = [flux * 1e-6 * (1.0 / 30) * (1.0 / 0.3) for flux in glycine_uptake]  # to mol AA / s / g dry cell mass
# phenylalanine_uptake = [flux * 1e-6 * (1.0 / 30) * (1.0 / 0.3) for flux in phenylalanine_uptake]  # to mol AA / s / g dry cell mass
# methionine_uptake = [flux * 1e-6 * (1.0 / 30) * (1.0 / 0.3) for flux in methionine_uptake]  # to mol AA / s / g dry cell mass
#
# # plot data only
# # plt.subplot(3, 1, 1)
# plt.plot(isoleucine_concs, isoleucine_uptake, 's--', linewidth=2, markersize=12, label='isoleucine')
# plt.plot(glycine_concs, glycine_uptake, '^--', linewidth=2, markersize=12, label='glycine')
# plt.plot(phenylalanine_concs, phenylalanine_uptake, 'o--', linewidth=2, markersize=12, label='phenylalanine')
# plt.plot(methionine_concs, methionine_uptake, '.--', linewidth=2, markersize=12, label='methionine')
#
# plt.title('concentration dependance of amino acid uptake')
# plt.ylabel('uptake (mol AAs / 1s / g DCW)')
# plt.xlabel('concentration (M)')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))


# save_piperno = {
# 	'ILE[p] c' : isoleucine_concs,
# 	'GLY[p] c' : glycine_concs,
# 	'PHE[p] c': phenylalanine_concs,
# 	'MET[p] c' : methionine_concs,
# 	'ILE[p] v': isoleucine_uptake,
# 	'GLY[p] v' : glycine_uptake,
# 	'PHE[p] v' : phenylalanine_uptake,
# 	'MET[p] v' : methionine_uptake,
# 	}

# import ipdb; ipdb.set_trace()

# # plot michaelis menten curves
# # plt.subplot(3, 1, 2)
# concentrations = np.arange(0.0, 60.0, 1)
# for molecule, params in exchange_param.iteritems():
# 	km = params['km']
# 	vmax = params['vmax']
#
# 	uptake = []
# 	for conc in concentrations:
# 		uptake.append(get_flux(conc, km, vmax))
#
# 	plt.plot(concentrations, uptake, label=molecule)
#
# plt.title('concentration dependance of amino acid uptake')
# plt.ylabel('uptake (umol AAs / 30s / g wet cells)')
# plt.xlabel('concentration (uM) (1e-6M)')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))



#
# plot data with michaelis menten curves
# plt.subplot(3, 1, 3)
concentrations = np.arange(0.0, 60.0, 1)
for molecule, params in exchange_param.iteritems():
	km = params['km']
	vmax = params['vmax']

	uptake = []
	for conc in concentrations:
		uptake.append(get_flux(conc, km, vmax))

	plt.plot(concentrations, uptake, linewidth=2, label=molecule)

plt.plot(isoleucine_concs, isoleucine_uptake, 's', linewidth=2, markersize=8, label='isoleucine data')
plt.plot(glycine_concs, glycine_uptake, '^', linewidth=2, markersize=8, label='glycine data')
plt.plot(phenylalanine_concs, phenylalanine_uptake, 'o', linewidth=2, markersize=8, label='phenylalanine data')
plt.plot(methionine_concs, methionine_uptake, '.', linewidth=2, markersize=8, label='methionine data')

plt.title('concentration dependance of amino acid uptake')
plt.ylabel('uptake (umol AAs / 30s / g wet cells)')
plt.xlabel('concentration (uM) (1e-6M)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))




if not os.path.exists(OUTDIR):
	os.mkdir(OUTDIR)
fig_name = ('update_data')
plt.savefig(os.path.join(OUTDIR, fig_name), bbox_inches='tight')
