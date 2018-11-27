from __future__ import absolute_import, division, print_function

import os
import csv
import json


PARAM_FILE = 'best_parameters.tsv'

PARAMOUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out/saved_parameters'
	)

OUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out', 'plot_out'
	)

WCM_SIMDATA_FILE = os.path.join(
	os.path.split(__file__)[0],
	'data', 'wcm_sim_data.json'
	)

with open(WCM_SIMDATA_FILE, "r") as f:
	wcm_sim_out = json.loads(f.read())

CONDITIONS_FILE = os.path.join(
	os.path.split(__file__)[0],
	'data', 'conditions.json'
	)


# Piperno and Oxender 1968 target data TODO -- put in conditions.
GLYCINE_DATA = os.path.join(
	os.path.split(__file__)[0],
	'data', 'piperno_oxender_1968', 'glycine.csv'
	)

ISOLEUCINE_DATA = os.path.join(
	os.path.split(__file__)[0],
	'data', 'piperno_oxender_1968', 'isoleucine.csv'
	)

METHIONINE_DATA = os.path.join(
	os.path.split(__file__)[0],
	'data', 'piperno_oxender_1968', 'methionine.csv'
	)

PHENYLALANINE_DATA = os.path.join(
	os.path.split(__file__)[0],
	'data', 'piperno_oxender_1968', 'phenylalanine.csv'
	)

target_definition = {'GLY[p]': [], 'ILE[p]': [], 'MET[p]': [], 'PHE[p]': []}


# TODO -- write function for units conversion
with open(GLYCINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0]) # in uM
		flux = float(row[1]) # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6 # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3) # to mol AA / s / g dry cell mass

		target = {'substrate_concentration' : conc, 'flux' : flux}
		target_definition['GLY[p]'].append(target)

with open(ISOLEUCINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0]) # in uM
		flux = float(row[1]) # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6 # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3) # to mol AA / s / g dry cell mass

		target = {'substrate_concentration': conc, 'flux': flux}
		target_definition['ILE[p]'].append(target)

with open(METHIONINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0]) # in uM
		flux = float(row[1]) # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6 # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3) # to mol AA / s / g dry cell mass

		target = {'substrate_concentration': conc, 'flux': flux}
		target_definition['MET[p]'].append(target)

with open(PHENYLALANINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0]) # in uM
		flux = float(row[1]) # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6 # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3) # to mol AA / s / g dry cell mass

		target = {'substrate_concentration': conc, 'flux': flux}
		target_definition['PHE[p]'].append(target)

## Initialize the reactions
# load all reactions from file
REACTIONS_FILE = os.path.join(
	os.path.split(__file__)[0],
	'data/aa_transport_reactions.json'
	# 'data/aa_transport_reactions_curated.json'
	)

with open(REACTIONS_FILE, "r") as f:
	ALL_REACTIONS = json.loads(f.read())

# # TODO -- this should replace REACTIONS_FILE entirely
# REACTIONS_FILE_RAW = os.path.join(
# 	os.path.split(__file__)[0],
# 	'data/aa_transport_reactions.json'
# 	)
#
# with open(REACTIONS_FILE_RAW, "r") as f:
# 	ALL_REACTIONS_RAW = json.loads(f.read())

# define parameter for each kinetic rate law
REACTION_PARAMS = {
	'uniport' : ['kcat', 'km'],
	'uniport_reversible' : ['kcat_f', 'kcat_r', 'km_A1', 'km_A2'],
	'symport' : ['kcat', 'km_A', 'km_B'],
	'symport_reversible' : ['kcat_f', 'kcat_r', 'km_A1', 'km_A2', 'km_B1', 'km_B2'],
	}

