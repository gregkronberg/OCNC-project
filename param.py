# parameters 

from neuron import h
import numpy as np
import cell 
import stims


class exp_default:
	def __init__(self):

		self.sec_idx = [0]
		self.seg_idx = []
		for a in self.sec_idx:
			self.seg_idx.append([-1])
		self.params = {
			'experiment':'default',
			'trial':0,
			'w_ampa':3*.00018, # ampa weight (microsiemens or micro-ohms)
			'w_nmda':3*.00018, # nmda weight (microsiemens or micro-ohms)
			'sec_idx': self.sec_idx,
			'seg_idx':self.seg_idx,
			'plot_idx':self.sec_idx,
			'field_angle':0,
			'field':[-20,0,20],
			'field_color':['b','k','r'],
			}

class exp_1:
	def __init__(self):

		self.sec_idx = [5,6,7]
		self.seg_idx = []
		for a in self.sec_idx:
			self.seg_idx.append([-1])
		self.params = {
			'experiment':'exp_1',
			'trial':0,
			'w_ampa':4*.00018, # ampa weight (microsiemens or micro-ohms)
			'w_nmda':4*.00018, # nmda weight (microsiemens or micro-ohms)
			'sec_idx': self.sec_idx,
			'seg_idx':self.seg_idx,
			'plot_idx':self.sec_idx,
			'field_angle':0,
			'field':[-20,0,20],
			'field_color':['b','k','r'],
			}

class exp_2:
	"""
	choose random sections on apical tuft, specifying number of sections
	"""
	def __init__(self,num_sec):
		# self.cell_temp = cell.Cell()
		self.sec_idx = np.random.choice(70,num_sec)
		self.seg_idx = []
		for a in self.sec_idx:
			self.seg_idx.append([0])
		self.params = {
			'experiment':'exp_2',
			'trial':0,
			'w_ampa':.00018, # ampa weight (microsiemens or micro-ohms)
			'w_nmda':.00018, # nmda weight (microsiemens or micro-ohms)
			'sec_idx': self.sec_idx,
			'seg_idx':self.seg_idx,
			'plot_idx':self.sec_idx,
			'field_angle':0,
			'field':[-20,0,20],
			'field_color':['b','k','r'],
			}
