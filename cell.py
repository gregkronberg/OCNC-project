"""
create cells and activate subsets of synapses
"""
# imports
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import stims

# create cell
class CellMigliore2005:
	""" pyramidal neuron based on Migliore et al. 2005

	An instance of this object will creates a cell (hoc objects) at the top level of the hoc interpreter using the hoc files in _init_geometry.  The .geo attribute contains a python mapping to these hoc objects.  The geo object is organized as geo['section tree'][section](segment location)

	the syns attribute creates a container for synapse objects that are added to each segment in the hoc cell.  syns is organized as syns['section tree']['synapse type'][section][segment number]
	
	"""
	def __init__(self,p):
		# initialize geometry
		self.geometry(p)
		# insert membrane mechanisms
		self.mechanisms(p)

	def geometry(self,p):
		""" create cell geometry at hoc top level
		
		"""
		# load cell geometry into hoc interpreter
		h.load_file('geo5038804.hoc')  
		# set discretization based on dlambda rule (set dlambda in hoc file) 
		h.load_file('fixnseg.hoc')  	
		# dictionary for storing geometry ['tree'][sec](seg location)
		self.geo = {}
		# dictionary for storing synapse objects ['tree']['type'][sec][seg]
		self.syns = {}
		# add section trees to geometry dictionary
		self.geo['soma'] = h.soma
		self.geo['axon'] =  h.axon
		self.geo['basal'] = h.dendrite
		self.geo['apical_trunk'] = h.user5
		self.geo['apical_tuft'] = h.apical_dendrite
		
		# set temperature in hoc
		h.celsius = p['celsius']
		# set soma as origin for distance measurements
		h.distance(sec = self.geo['soma'][0])

	def mechanisms(self,p):
		""" insert membrane mechanisms

		self.syns is updated to store an object for each synapse.  It is organized as ['tree']['synapse type'][section][segment].  Note that the last index will depend on how the cell is discretized as the number segments changes in each sections 

		the parameters for each membrane mechanism  are store in a dictionary called p.  See the param module for details.
		"""

		# loop over trees
		for tree_key,tree in self.geo.iteritems():
			
			# create sub-dictionary for different types of synapses
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
				# passive conductance (S/cm2)
				sec.g_pas = 1/p['RmAll']			
				# leak reversal potential (mV)	
				sec.e_pas = p['Vrest']				
				# specific capacitance (uf/cm2)
				sec.cm = p['Cm'] 			
				# axial resistance (ohm cm) 		
				sec.Ra = p['RaAx'] 
										


				# axon active bipophysics
				if tree_key == 'axon':
					# voltage gated sodium
					sec.insert('nax')						
					sec.gbar_nax = p['gna']*p['AXONM']
					# delayed rectifier potassium
					sec.insert('kdr')						
					sec.gkdrbar_kdr = p['gkdr']
					# a-type potassium
					sec.insert('kap')						
					sec.gkabar_kap = p['KMULTP']
					# sodium reversal potential 
					sec.ena = p['ena']		
					# potassium reversal potential 
					sec.ek = p['ek']
					
				# soma active biophysics
				elif tree_key == 'soma':
					# voltage gated sodium
					sec.insert('na3')
					sec.gbar_na3 = p['gna']	
					# h-current			
					sec.insert('hd')
					sec.ghdbar_hd = p['ghd']				
					sec.vhalfl_hd = p['vhalfl_prox']
					# delayed rectifier potassium		
					sec.insert('kdr')
					sec.gkdrbar_kdr = p['gkdr']	
					# a-type potassium		
					sec.insert('kap')
					sec.gkabar_kap = p['KMULTP']
					# sodium reversal potential 
					sec.ena = p['ena']		
					# potassium reversal potential 
					sec.ek = p['ek']			

				# dendrites active biophysics
				elif ((tree_key == 'basal') or 
				(tree_key == 'apical_trunk') or 
				(tree_key == 'apical_tuft')):
					# h-current
				    sec.insert('hd')
				    sec.ghdbar_hd = p['ghd']
				    # voltage gated sodium		
				    sec.insert('na3')
				    sec.gbar_na3 = p['gna']	
				    # delayed rectifier potassium	
				    sec.insert('kdr')
				    sec.gkdrbar_kdr = p['gkdr']	
				    # a-type potassium proximal
				    sec.insert('kap')
				    sec.gkabar_kap = 0	
				    # a-type potassium distal		
				    sec.insert('kad')
				    sec.gkabar_kad = 0	
				    # sodium reversal potential 
				    sec.ena = p['ena']
				    # potassium reversal potential
				    sec.ek = p['ek']		

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

class CellKim2015:
	pass

class Syn_act:
	"""Activate a specific subset of synpases with NetCon objects
	
	arguments: 
	syns = nested dictionary containing synapse objects organized as ['section tree']['synapse type'][section number][segment number]

	p =  dictionary of parameters including sec_idx and seg_idx lists, including section and segment numbers to be activated.  p also contains 'w_list', which sets the weight for each of the activated synapses. w_list has the same dimensions as seg_idx

	stim = list of NetStim objects to be connected (via NetCon) to the synapses in syns that designated by sec_idx and seg_idx. The stim list will be iterated through and a distinct NetCon object will be created for each NetStim for each activated synapse.  The need for multiple NetStim objects arises from complicated stimulation patterns, like theta bursts, which are not easily programmed with a single NetStim

	The created NetCon objects are referenced by the nc object, which is organized the same way as syns, namely ['section tree']['synapse type'][section number][segment number][NetStim object]
	"""
	def __init__(self, syns, p, stim):
		# store netcon objects ['tree']['syn type'][section][segment][list of netstim objects]
		self.nc = {}
		
		# loop over synapse types
		for tree_key,tree in syns.iteritems():
			if tree_key in p['tree']:
				self.nc[tree_key] = {}

				for syntype_key,syn_type in syns[tree_key].iteritems():
					self.nc[tree_key][syntype_key] = []

					# loop over active sections
					for sec_i,sec in enumerate(p['sec_idx']):
						self.nc[tree_key][syntype_key].append([])
						
						# loop over active segments
						for seg_i,seg in enumerate(p['seg_idx'][sec_i]):
							self.nc[tree_key][syntype_key][sec_i].append([])

							# loop over stimulation bursts
							for syn_stim_i,syn_stim in enumerate(stim):
								print syns[tree_key][syntype_key][sec][seg]
								print p['w_list'][sec_i][seg_i]
								self.nc[tree_key][syntype_key][sec_i][seg_i].append(
									h.NetCon(syn_stim,syns[tree_key][syntype_key][sec][seg],0,0,p['w_list'][sec_i][seg_i]))

# set procedure if called as a script
if __name__ == "__main__":
	cell_1 = Cell_Migliore_2005()
