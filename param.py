# parameters 

from neuron import h
import numpy as np
import cell 
import stims


class exp1:
	def __init__(self):

		self.sec_idx = [5,6,7]
		self.seg_idx = []
		for a in self.sec_idx:
			self.seg_idx.append([-1])
		self.params = {
			'experiment':'exp1',
			'w_ampa':5*.00018, # ampa weight (microsiemens or micro-ohms)
			'w_nmda':5*.00018, # nmda weight (microsiemens or micro-ohms)
			'sec_idx': self.sec_idx,
			'seg_idx':self.seg_idx,
			'plot_idx':self.sec_idx,
			'field_angle':0,
			'field':[-20,0,20],
			'field_color':['b','k','r']}
