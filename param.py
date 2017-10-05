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
			'plot_sec_idx':self.sec_idx,
			'plot_seg_idx':self.seg_idx,
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
			'plot_sec_idx':self.sec_idx,
			'plot_seg_idx':self.seg_idx,
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
			'plot_sec_idx':self.sec_idx,
			'plot_seg_idx':self.seg_idx,
			'field_angle':0,
			'field':[-20,0,20],
			'field_color':['b','k','r'],
			'dt' : .025,
			'tstop' : 100
			}

class exp_3:
	"""
	Replicate TBS results from Rahman 2018

	Activate a random set of synapses on the apical dendritic tuft with TBS and record the resulting plasticity
	"""
	def __init__(self,syn_frac=0.4):
		
		# fraction of synapses that are randomly activated
		self.syn_frac = syn_frac
		# load cell
		self.cell = cell.Cell_Migliore_2005()
		# list all segments as [[section,segment]] 
		self.segs_all = ([[sec,seg] for sec in range(len(self.cell.syn_a_tuft_ampa)) for seg in range(len(self.cell.syn_a_tuft_ampa[sec]))])
		# choose percentage of segments to activate
		self.segs_choose = np.random.choice(len(self.segs_all),int(self.syn_frac*len(self.segs_all)),replace=False)

		# list of sections
		self.sec_list = [self.segs_all[a][0] for a in self.segs_choose]
		# list of segments
		self.seg_list = [self.segs_all[a][1] for a in self.segs_choose]
		# uniqure list of sections
		self.sec_idx  = list(set(self.sec_list))
		# list of segments as [unique section][segments]
		self.seg_idx = []
		for sec in self.sec_idx:
			self.seg_idx.append([self.seg_list[sec_i] for sec_i,sec_num in enumerate(self.sec_list) if sec_num==sec])
		
		# set parameters.  Cannot contain any hoc objects, as this will be pickled for data storage
		self.params = {
			'experiment':'exp_3',
			'data_folder':'Data/exp_3/',
			'fig_folder':'png figures/exp_3/',
			'syn_frac':self.syn_frac,
			'trial':0,
			'w_rand':False,
			'w_std':.0005,
			'w_ampa':.0018, # ampa weight (microsiemens or micro-ohms)
			'w_nmda':.0018, # nmda weight (microsiemens or micro-ohms)
			'sec_list':self.sec_list,
			'seg_list':self.seg_list,
			'sec_idx': self.sec_idx,
			'seg_idx':self.seg_idx,
			'plot_sec_list':self.sec_list,
			'plot_seg_list':self.seg_list,
			'plot_sec_idx':self.sec_idx,
			'plot_seg_idx':self.seg_idx,
			'field_angle':0,
			'field':[-20,0,20],
			'field_color':['r','k','b'],
			'dt' : .025,
			'warmup': 60,
			'bursts':1,
			'pulses':4,
			'pulse_freq':100,
			'burst_freq':5,
			'tstop' : 100,#5*1000/5 + 30 + 5*1000/100 +30
			
			# clopath parameters
			'clopath_delay_steps': 1,
			'clopath_tau_0':12, # time constant (ms) for low passed membrane potential for depression
			'clopath_tau_r' : 30, # time constant (ms) for low pass filter presynaptic variable
			'clopath_tau_y': 20, # time constant (ms) for low pass filter post membrane potential for potentiation
			'clopath_A_m':.0001, # depression magnitude parameter (mV^-1)
			'clopath_A_p': .004, # amplitude for potentiation (mV^-2)
			'clopath_tetam':-41, # depression threshold (mV)
			'clopath_tetap':-38, # potentiation threshold (mV)

			}

# set procedure if called as a script
if __name__ == "__main__":
	pass