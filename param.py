# parameters 

from neuron import h
import numpy as np
import cell 
import stims

class exp_3:
	"""
	Replicate TBS results from Rahman 2018

	Activate a random set of synapses on the apical dendritic tuft with TBS and record the resulting plasticity
	"""
	def __init__(self, syn_frac=0.4, w_mean=.0018, w_std=.0002, w_rand=True, exp='exp_4'):
		"""set parameters in dictionary p
		
		p cannot contain any hoc objects, as this will be pickled and stored with each experiment so that the parameters can be retrieved
		
		arguments:
		syn_frac = fraction of synapses to activated

		w_mean = mean synaptic weight

		w_std = standard deviation of synaptic weights drawn from gaussian distribution

		w_rand = if True synaptic weights are drawn from gaussian(w_mean,w_std). If False, weights are uniformly w_mean for all synapses  

		exp = name of the experiment for storing and retrieving data
		"""
	
		# set default parameters  
		self.p = {
			'experiment':exp,
			'data_folder':'Data/'+exp+'/',
			'fig_folder':'png figures/'+exp+'/',
			'syn_frac':syn_frac,
			'trial':0,
			'trial_id':0,
			'w_rand':w_rand,
			'w_std':w_std,
			'w_mean':w_mean, # mean synaptic weight (microsiemens or micro-ohms)
			'tree':'apical_trunk',
			'w_list':[],
			'sec_list':[],
			'seg_list':[],
			'sec_idx': [],
			'seg_idx':[],
			'field_angle':0,
			'field':[-20,0,20],
			'field_color':['r','k','b'],
			'dt' : .025,
			'warmup': 30,
			'tstop' : 70,#5*1000/5 + 30 + 5*1000/100 +30
			'bursts':1,
			'pulses':4,
			'pulse_freq':100,
			'burst_freq':5,

			# clopath synapse parameters
			'clopath_delay_steps': 1,
			'clopath_tau_0':12, # time constant (ms) for low passed membrane potential for depression
			'clopath_tau_r' : 30, # time constant (ms) for low pass filter presynaptic variable
			'clopath_tau_y': 10, # time constant (ms) for low pass filter post membrane potential for potentiation
			'clopath_A_m':.0001, # depression magnitude parameter (mV^-1)
			'clopath_A_p': .001*1.01, # amplitude for potentiation (mV^-2)
			'clopath_tetam':-41, # depression threshold (mV)
			'clopath_tetap':-38, # potentiation threshold (mV)

			# ampa synapse parameters
			'tau1_ampa' : 0.2,
			'tau2_ampa' : 2,
			'i_ampa' : 0.18,

			# nmda synapse parameters
			'tau1_nmda' : 1,
			'tau2_nmda' : 50,

			# Parameters from Migliore 2005 (signal propogation in oblique dendrites)
			
			# conductances reported as (nS/um2) in paper, but need to be in (mho/cm2)
			# conversion 10*(nS/um2) = (mho/cm2)
			# *** units in paper are a typo, values are already reported in (mho/cm2) ***
			'Vrest' : -65,				# resting potential (mV)
			'gna' : .25,				# peak sodium conductance (mho/cm2)
			'ena' : 55,					# sodium reversal potential (mV)
			'AXONM' : 5,				# multiplicative factor for axonal conductance
			'gkdr' : 0.1,				# delayed rectifier potassium peak conductance (mho/cm2)
			'ek' : -90,					# potassium reversal potential
			'celsius' : 35.0,  				# temperature (degrees C)
			'KMULT' : 0.3,			# multiplicative factor for distal A-type potassium conductances
			'KMULTP' : 0.3,				# multiplicative factor for proximal A-type potassium conductances
			'ghd' : 0.0005,			# peak h-current conductance (mho/cm2)
			'ehd' : -30,					# h-current reversal potential (mV)
			'vhalfl_prox' : -73,			# activation threshold for proximal a-type potassium (mV)
			'vhalfl_dist' : -81,			# activation threshold for distal a-type potassium (mV)
			'RaAll' : 150,				# axial resistance, all compartments (ohm*cm)
			'RaAx' : 50,					# axial resistance, axon (ohm*cm)					
			'RmAll' : 28000,			# specific membrane resistance (ohm/cm2)
			'Cm' : 1,					# specific membrane capacitance (uf/cm2)
			'ka_grad' : 1,				# slope of a-type potassium channel gradient with distance from soma 

			}

			# fraction of synapses that are randomly activated
		
		# load cell
		self.cell = cell.CellMigliore2005(self.p)

		# randomly choose active segments 
		self.choose_seg_rand(syn_list = self.cell.syns[self.p['tree']]['ampa'], syn_frac=syn_frac)
		
		# set weights for active segments
		self.set_weights(seg_idx=self.p['seg_idx'], w_mean=w_mean, w_std = w_std, w_rand=w_rand)

	def choose_seg_rand(self, syn_list, syn_frac):
		""" choose random segments to activate

		arguments:
		syn_list = list of all synapses to be chosen from ogranized as [section number][segment number]

		syn_frac = fraction of synapses to be chosen, with equal probability for all synapses

		updates the parameter dictionary according to the chosen synapses
		"""

		# list all segments as [[section,segment]] 
		segs_all = [[sec_i,seg_i] for sec_i,sec in enumerate(syn_list) for seg_i,seg in enumerate(syn_list[sec_i])]

		# choose segments to activate
		segs_choose = np.random.choice(len(segs_all),int(syn_frac*len(segs_all)),replace=False)

		# list of active sections (contains duplicates)
		sec_list = [segs_all[a][0] for a in segs_choose]
		
		# list of active segments
		seg_list = [segs_all[a][1] for a in segs_choose]

		# uniqure list of active sections
		sec_idx  = list(set(sec_list))
		
		# list of active segments as [unique section][segments]
		seg_idx = []
		for sec in sec_idx:
			seg_idx.append([seg_list[sec_i] for sec_i,sec_num in enumerate(sec_list) if sec_num==sec])

		# update parameter dictionary
		self.p['sec_list'] = sec_list
		self.p['seg_list'] = seg_list
		self.p['sec_idx'] = sec_idx
		self.p['seg_idx'] = seg_idx

	def set_weights(self, seg_idx, w_mean=.0018, w_std=.0002, w_rand=True):
		w_list = []
		# loop over sections
		for sec_i,sec in enumerate(seg_idx):
			w_list.append([])
			# loop over segments
			for seg_i,seg in enumerate(seg_idx[sec_i]):
				# if weights are randomized
				if w_rand:
					# choose from normal distribution
					w_list[sec_i].append(np.random.normal(w_mean,w_std))
				# otherwise set all weights to the same
				else:
					w_list[sec_i].append(w_mean)

		# update parameter dictionary
		self.p['w_list']=w_list

# set procedure if called as a script
if __name__ == "__main__":
	pass