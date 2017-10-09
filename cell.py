"""
create cells and activate subsets of synapses
"""
# imports
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import stims

# create cell
class Cell_Migliore_2005:
	"""
	pyramidal neuron based on Migliore et al. 2005
	"""
	def __init__(self,p):
		self.p = p
		self._init_geometry(p)
		self.insert_mech(p)

	def _init_geometry(self,p):
		"""
		create cell geometry at hoc top level
		"""
		h.load_file('geo5038804.hoc')   # load cell class geometry from Migliore 2005
		h.load_file('fixnseg.hoc')  	# set discretization based on dlambda rule (set dlambda in hoc file)
		
		self.geo = {}
		self.syns = {}
		self.geo['soma'] = h.soma
		self.geo['axon'] =  h.axon
		self.geo['basal'] = h.dendrite
		self.geo['apical_trunk'] = h.user5
		self.geo['apical_tuft'] = h.apical_dendrite
		

		h.celsius = p['celsius']
		# set soma as origin for distance measurements
		h.distance(sec = self.geo['soma'][0])

	def insert_mech(self,p):
		# loop over trees
		for tree_key,tree in self.geo.iteritems():
			
			# create sub-dictionary for different types
			self.syns[tree_key] = {
			'ampa' : [],
			'nmda' : [],
			'clopath' : []
			}

			# loop over sections
			for sec_i,sec in enumerate(tree):
				
				# add list to store synapses for each section
				for syn_key,syn in self.syns[tree_key].iteritems():
					syn.append([])

				# common passive biophysics for all sections
				sec.insert('pas')
				sec.g_pas = 1/p['RmAll']				# passive conductance (S/cm2)
				sec.e_pas = p['Vrest']				# leak reversal potential (mV)
				sec.cm = p['Cm'] 						# specific capacitance (uf/cm2)
				sec.Ra = p['RaAx'] 					# axial resistance (ohm cm) 

				# axon
				if tree_key == 'axon':
					# active biophysics (parmameters: Migliore 2005)
					sec.insert('nax')						# voltage gated sodium
					sec.gbar_nax = p['gna']*p['AXONM']
					sec.insert('kdr')						# delayed rectifier potassium
					sec.gkdrbar_kdr = p['gkdr']
					sec.insert('kap')						# a-type potassium
					sec.gkabar_kap = p['KMULTP']
					sec.ena = p['ena']					# sodium reversal potential (see _init_parameters)
					sec.ek = p['ek']						# potassium reversal potential (see _init_parameters)

				# soma
				elif tree_key == 'soma':
					# soma biophysics (parmameters: Migliore 2005)
					sec.insert('na3')
					sec.gbar_na3 = p['gna']				# voltage gated sodium
					sec.insert('hd')
					sec.ghdbar_hd = p['ghd']				# h-current
					sec.vhalfl_hd = p['vhalfl_prox']		# h-current activation threshold
					sec.insert('kdr')
					sec.gkdrbar_kdr = p['gkdr']			# delayed rectifier potassium
					sec.insert('kap')
					sec.gkabar_kap = p['KMULTP']			# scaling factor for a-type potassium current
					sec.ena = p['ena']					# sodium reversal potential 
					sec.ek = p['ek']						# potassium reversal potential 

				# dendrites
				elif ((tree_key == 'basal') or 
				(tree_key == 'apical_trunk') or 
				(tree_key == 'apical_tuft')):
					# active biophysics (parameters: Migliore 2005)
				    sec.insert('hd')
				    sec.ghdbar_hd = p['ghd']		# h-current
				    sec.insert('na3')
				    sec.gbar_na3 = p['gna']		# voltage gated sodium
				    sec.insert('kdr')
				    sec.gkdrbar_kdr = p['gkdr']	# delayed rectifier potassium
				    sec.insert('kap')
				    sec.gkabar_kap = 0			# a-type potassium proximal
				    sec.insert('kad')
				    sec.gkabar_kad = 0			# a-type potassium distal
				    sec.ena = p['ena']			# sodium reversal potential 
				    sec.ek = p['ek']				# potassium reversal potential 

				    # mechanisms that vary with distance from soma
				    # loop over segments
				    for seg_i,seg in enumerate(sec):
				    	
				    	# distance from soma
				    	seg_dist = h.distance(seg.x,sec=sec)
				    	
				    	# h current
				    	seg.ghdbar_hd = p['ghd']*(1+3*seg_dist/100)
				    	
				    	# A-type potassium
				        if seg_dist > 100:	# distal
				            seg.vhalfl_hd = p['vhalfl_dist']
				            seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*seg_dist/100)
				        
				        else:	# proximal
				            seg.vhalfl_hd = p['vhalfl_prox']
				            seg.gkabar_kap = p['KMULTP']*(1+p['ka_grad']*seg_dist/100)

				        # loop over synapse types
				        for syn_key,syn in self.syns[tree_key].iteritems():
				        	if syn_key == 'ampa':
				        		syn[sec_i].append(h.Exp2Syn(sec(seg.x)))
				        		syn[sec_i][seg_i].tau1 = p['tau1_ampa']
				        		syn[sec_i][seg_i].tau2 = p['tau2_nmda']
				        		syn[sec_i][seg_i].i = p['i_ampa']
				        	elif syn_key == 'nmda':
				        		syn[sec_i].append(h.Exp2SynNMDA(sec(seg.x)))
				        		syn[sec_i][seg_i].tau1 = p['tau1_nmda']
				        		syn[sec_i][seg_i].tau2 = p['tau2_nmda']
				        	elif syn_key == 'clopath':
				        		syn[sec_i].append(h.STDPSynCCNon(sec(seg.x)))

# activate a specified subset of synapses with NetStim/NetCon
class Syn_act:
	"""
	Activate a specific subset of synpases with NetCon objects
	"""
	def __init__(self,syns,p,stim):
		# store netcon objects [dic of synapse type][section][segment][list of netstim objects]
		self.nc = {}
		
		# loop over synapse types
		for tree_key,tree in syns.iteritems():
			if tree_key == p['tree']:
				self.nc[tree_key] = {}

				for syntype_key,syn_type in syns[tree_key].iteritems():
					self.nc[tree_key][syntype_key] = []

					# loop over active sections
					for sec_i,sec in enumerate(p['sec_idx']):
						self.nc[tree_key][syntype_key].append([])
						
						# loop over active segments
						for seg_i,seg in enumerate(p['seg_idx'][sec_i]):
							self.nc[tree_key][syntype_key][sec_i].append([])

							print tree_key,syntype_key
							print sec,seg,len(syns[tree_key][syntype_key])
							# loop over stimulation bursts
							for syn_stim_i,syn_stim in enumerate(stim):
								self.nc[tree_key][syntype_key][sec_i][seg_i].append(
									h.NetCon(syn_stim,syns[tree_key][syntype_key][sec][seg],0,0,p['w_list'][sec_i][seg_i]))

# set procedure if called as a script
if __name__ == "__main__":
	cell_1 = Cell_Migliore_2005()
