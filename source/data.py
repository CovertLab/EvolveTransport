from __future__ import absolute_import, division, print_function

import os
import csv
import json


PARAM_FILE = 'best_parameters.tsv'

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(BASE_DIR, 'out')
DATADIR = os.path.join(BASE_DIR, 'data')

# output locations
PARAMOUTDIR = os.path.join(OUTDIR, 'saved_parameters')
PLOTOUTDIR = os.path.join(OUTDIR, 'plot_out')

# Initialize the reactions
REACTIONS_FILE = os.path.join(
	DATADIR, 'aa_transport_reactions.json'
	)

with open(REACTIONS_FILE, 'r') as f:
	ALL_REACTIONS = json.loads(f.read())

# data
WCM_SIMDATA_FILE = os.path.join(DATADIR, 'wcm_sim_data.json')
CONDITIONS_FILE = os.path.join(DATADIR, 'conditions.json')

with open(WCM_SIMDATA_FILE, 'r') as f:
	wcm_sim_out = json.loads(f.read())


amino_acids = [
	'L-ALPHA-ALANINE',
	'ARG',
	'ASN',
	'L-ASPARTATE',
	'CYS',
	'GLT',
	'GLN',
	'GLY',
	'HIS',
	'ILE',
	'LEU',
	'LYS',
	'MET',
	'PHE',
	'PRO',
	'SER',
	'THR',
	'TRP',
	'TYR',
	'L-SELENOCYSTEINE',
	'VAL'
]


# Piperno and Oxender 1968 target data
GLYCINE_DATA = os.path.join(
	DATADIR, 'piperno_oxender_1968', 'glycine.csv'
	)

ISOLEUCINE_DATA = os.path.join(
	DATADIR, 'piperno_oxender_1968', 'isoleucine.csv'
	)

METHIONINE_DATA = os.path.join(
	DATADIR, 'piperno_oxender_1968', 'methionine.csv'
	)

PHENYLALANINE_DATA = os.path.join(
	DATADIR, 'piperno_oxender_1968', 'phenylalanine.csv'
	)

target_definition = {'GLY[p]': [], 'ILE[p]': [], 'MET[p]': [], 'PHE[p]': []}


# TODO -- write function for units conversion
with open(GLYCINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0])  # in uM
		flux = float(row[1])  # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6  # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3)  # to mol AA / s / g dry cell mass

		target = {'substrate_concentration': conc, 'flux': flux}
		target_definition['GLY[p]'].append(target)

with open(ISOLEUCINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0])  # in uM
		flux = float(row[1])  # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6    # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3)  # to mol AA / s / g dry cell mass

		target = {'substrate_concentration': conc, 'flux': flux}
		target_definition['ILE[p]'].append(target)

with open(METHIONINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0])  # in uM
		flux = float(row[1])  # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6    # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3)  # to mol AA / s / g dry cell mass

		target = {'substrate_concentration': conc, 'flux': flux}
		target_definition['MET[p]'].append(target)

with open(PHENYLALANINE_DATA) as target_file:
	reader = csv.reader(target_file)
	next_line = next(reader)

	for row in reader:
		conc = float(row[0])  # in uM
		flux = float(row[1])  # in umol AA/ 30 s/ g wet cell
		conc = conc * 1e-6    # convert to M
		flux = flux * 1e-6 * (1.0/30) * (1.0/0.3)  # to mol AA / s / g dry cell mass

		target = {'substrate_concentration': conc, 'flux': flux}
		target_definition['PHE[p]'].append(target)

with open(CONDITIONS_FILE, 'r') as f:
	SET_CONDITIONS = json.loads(f.read())