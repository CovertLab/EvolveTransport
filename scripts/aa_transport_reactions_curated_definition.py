"""
make_aa_transport_reactions
- makes a tsv file for aa transport reactions
- makes a json file with a dict that has {aa: [rxns]}, with each amino acid given a list of transport reactions
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import os
import cPickle
import csv
import re
import json
import numpy as np
import ast

directory = os.path.join('reconstruction', 'ecoli', 'flat')

# transport reactions
with open(os.path.join(directory, 'aa_transport_reactions.tsv'), 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	headers = reader.next()
	wcTransport = np.array(list(reader))
	rxn_id_header = headers.index('reaction id')
	stoic_header = headers.index('stoichiometry')
	reversible_header = headers.index('is reversible')
	transporter_header = headers.index('catalyzed by')
	types_header = headers.index('type')

	rxn_id = (wcTransport[:, rxn_id_header]).tolist()
	stoichiometries = (wcTransport[:, stoic_header]).tolist()
	reversibility = (wcTransport[:, reversible_header]).tolist()
	catalysts = (wcTransport[:, transporter_header]).tolist()
	types = (wcTransport[:, types_header]).tolist()

	# make catalysts a list of lists
	catalysts_lists = []
	for catalyst in catalysts:
		catalysts_lists.append(ast.literal_eval(catalyst))
	# make stoic a list of dicts
	stoic_dicts = []
	for st in stoichiometries:
		stoic_dicts.append(ast.literal_eval(st))

	rxn_dict = dict(zip(rxn_id, stoic_dicts))
	reversibility_dict = dict(zip(rxn_id, reversibility))
	transporter_dict = dict(zip(rxn_id, catalysts_lists))
	types_dict = dict(zip(rxn_id, types))


# get amino acids, remove compartment from id to make amino_acids, initialize amino acids summary dictionary
sim_data = cPickle.load(open("out/manual/kb/simData_Fit_1.cPickle", "rb"))
amino_acids_compartments = sim_data.moleculeGroups.aaIDs
amino_acids = [re.sub("[[@*&?].*[]@*&?]", "", aa) for aa in amino_acids_compartments]


reactions = {}
for rxn, stoich in rxn_dict.iteritems():
	type = types_dict[rxn]

	transporter = transporter_dict[rxn]

	if type == 'ABC':
		type = 'symport_reversible'

		for mol, coefficient in stoich.items():
			mol_no_comp = re.sub("[[@*&?].*[]@*&?]", "", mol)

			# don't include ATP hydrolysis substrates
			if mol_no_comp  not in ['PI', 'WATER', 'PROTON']:
				is_aa = False
				for aa in amino_acids:
					if aa == mol_no_comp:
						is_aa = True

				if coefficient == -1 and is_aa:
					A1 = mol
				elif coefficient == 1 and is_aa:
					A2 = mol
				elif coefficient == -1:
					B1 = mol
				elif coefficient == 1:
					B2 = mol
		specs = {'A1': A1, 'A2': A2, 'B1': B1, 'B2': B2}

	elif type == 'symporter':
		type = 'symport'

		for mol, coefficient in stoich.items():
			mol_no_comp = re.sub("[[@*&?].*[]@*&?]", "", mol)

			is_aa = False
			for aa in amino_acids:
				if aa == mol_no_comp:
					is_aa = True

			if coefficient == -1 and is_aa:
				A1 = mol
			elif coefficient == 1 and is_aa:
				A2 = mol
			elif coefficient == -1:
				B1 = mol
			elif coefficient == 1:
				B2 = mol
		specs = {'A1': A1, 'A2': A2, 'B1': B1, 'B2': B2}

	elif type == 'Membrane_Electrochemical_Gradient':
		type = 'uniport'
		for mol, coefficient in stoich.items():
			if coefficient == -1:
				A1 = mol
			elif coefficient == 1:
				A2 = mol
		specs = {'A1' : A1, 'A2' : A2}

	elif type == 'Unknown_Mechanism':
		type = 'uniport'
		for mol, coefficient in stoich.items():
			if coefficient == -1:
				A1 = mol
			if coefficient == 1:
				A2 = mol
		specs = {'A1' : A1, 'A2' : A2}

	elif type == 'Channels_and_Pores':
		type = 'uniport'
		for mol, coefficient in stoich.items():
			if coefficient == -1:
				A1 = mol
			if coefficient == 1:
				A2 = mol
		specs = {'A1' : A1, 'A2' : A2}

	else:
		raise Exception('Unknown reaction type: {}'.format(type))

	reactions[rxn] = {'type' : type, 'substrates' : specs, 'transporter': transporter}

# save aa_to_rxns dict to json
with open('user/aa_transport_reactions_curated.json', 'w') as outfile:
	json.dump(reactions, outfile, indent=2)
