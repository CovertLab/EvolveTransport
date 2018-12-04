'''

run this with by passing path to simout:
python prototypes/aa_transport_estimation/scripts/initial_state_from_WCM.py --simout out/manual/condition_000002/000000/generation_000000/000000/simOut

'''

import os, cPickle
import argparse
import json
import copy

from wholecell.io.tablereader import TableReader

OUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out'
	)


TIME_STEP = 0.1 # seconds

## Initialize the reactions
# load all reactions from file
REACTIONS_FILE = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'data/aa_transport_reactions.json'))
# REACTIONS_FILE = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'data/aa_transport_reactions_reduced.json'))

with open(REACTIONS_FILE, "r") as f:
	ALL_REACTIONS = json.loads(f.read())




class getInitialFlux(object):

	def main(self, args):

		self.nth_timestep = 100

		self.time, self.cell_mass, self.cell_volume = self.initialize_from_WCM(args.simout)
		self.counts = self.initialize_counts(args.simout)

		self.wcm_sim_data = self.counts.copy()
		self.wcm_sim_data['time'] = self.time
		self.wcm_sim_data['cell_mass'] = self.cell_mass
		self.wcm_sim_data['volume'] = self.cell_volume

		# save aa_to_rxns dict to json
		with open('user/wcm_sim_data.json', 'w') as outfile:
			json.dump(self.wcm_sim_data, outfile, indent=2)


	## Initialization
	def initialize_from_WCM(self, simOutDir):

		# get time
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		# get mass, density, and compute volume
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cell_mass = mass.readColumn("cellMass")
		# dry_cell_mass = mass.readColumn("dryMass")
		mass.close()
		density = 1100
		cell_volume = cell_mass / density

		# Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		# coefficient = dry_cell_mass[0] / cell_mass[0] * density * (TIME_STEP)

		# convert to list to jsonify data
		time = time.tolist()
		cell_mass = cell_mass.tolist()
		cell_volume = cell_volume.tolist()

		return time, cell_mass, cell_volume

	def initialize_counts(self, simOutDir):
		''' set all initial undefined molecular concentrations to their initial concentrations in the WCM'''

		counts = {}
		indices = {}

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		molecule_ids = bulkMolecules.readAttribute("objectNames")

		# get indices of all molecules in self.reactions
		for rxn, specs in ALL_REACTIONS.iteritems():

			substrates = specs['stoichiometry'].keys()
			transporters = specs['transporters']

			# loop through substrates
			for substrate in substrates:
				substrate_index = molecule_ids.index(substrate)
				indices[substrate] = substrate_index

			# loop through transporters
			for transporter in transporters:
				transporter_id = [id for id in molecule_ids if transporter in id]
				transporter_index = molecule_ids.index(transporter_id[0])
				indices[transporter] = transporter_index

				# # if transporter is not in counts dict
				# if transporter not in counts.keys() or counts[transporter] is None:
				#
				# 	#get transporter id with added compartment tag
				# 	transporter_id = [id for id in molecule_ids if transporter in id]
				# 	transporter_index = molecule_ids.index(transporter_id[0])
				# 	indices[transporter] = transporter_index

		for molecule, molecule_index in indices.iteritems():

			print('reading ' + molecule)

			time_series = bulkMolecules.readColumn("counts")[:, molecule_index]

			counts[molecule] = time_series.tolist()


		return counts



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='evolve parameters for transport')
	parser.add_argument('--simout', help='directory of sim out data', default='out/manual/wildtype_000000/000000/generation_000000/000000/simOut')
	args = parser.parse_args()
	getInitialFlux().main(args)