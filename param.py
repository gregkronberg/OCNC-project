# parameters 

from neuron import h
import numpy as np
import cell 
import stims

sec_idx = [5,6,7]
seg_idx = []
for a in sec_idx:
	seg_idx.append([-1])
params = {
	'w_ampa':5*.00018, # ampa weight (microsiemens or micro-ohms)
	'w_nmda':5*.00018, # nmda weight (microsiemens or micro-ohms)
	'sec_idx': sec_idx,
	'seg_idx':seg_idx,
	'plot_idx':sec_idx,
	'field_angle':0,
	'field':[-20,0,20],
	'field_color':['b','k','r']}

# param = {
# 	'subtree_ampa':subtree_ampa # neuronal subtree_ampa (soma, apical, basal, etc.)
# 		self.subtree_nmda = subtree_nmda
# 		self.sec_idx = sec_idx	# list of section indeces with active synapses
# 		self.seg_idx = seg_idx	# nested list of active segment indices within each section, e.g [section][segment]
# 		self.w_ampa = w_ampa	# maximum ampa conductance (microsiemens or micro-omhos)
# 		self.w_nmda = w_nmda	# maximum ampa conductance (microsiemens or micro-omhos)
# 		self.choose_syn_ampa = [self.subtree_ampa[a] for a in self.sec_idx]
# 		self.choose_syn_nmda = [self.subtree_nmda[a] for a in self.sec_idx]
# 		self.active = []

# 		# default stim objects
# 		self.warm_up = 50   # warm up time (ms)
# 		self.pulse_freq = 100
# 		self.burst_freq = 5
# 		self.bursts = 1
# 		self.pulses = 1
# 		self.stim  = []
# 		self.stim.append(h.NetStim())
# 		self.stim[0].start = self.warm_up
# 		self.stim[0].interval = 1000/self.pulse_freq
# 		self.stim[0].noise  = 0 
# 		self.stim[0].number = self.pulses