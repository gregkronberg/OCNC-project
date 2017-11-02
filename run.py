"""
docstring
"""
# imports
from neuron import h
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cell 
import itertools as it
import stims
import pickle
import param
import os

# load standard runtime settings.  this is critical.
h.load_file("stdrun.hoc")

# run control
class Run():
	"""
	Arguments list:
	p - dictionary of parameters

	each experiment appends a list to the appropriate key in the data dictionary
	data is organized as data['type of data'][experiments][sections][time series data vector]
	details of each experiment are tracked via the data['detail'][experiment number], e.g. data['field'][2]
	"""
	def __init__(self, p):

		# create cell
		# self.cell1 = cell.CellMigliore2005(p)
		self.cell1 = cell.PyramidalCell(p) #p['cell']
		self.update_clopath( p, syns=self.cell1.syns[p['tree']]['clopath'])
		self.activate_synapses(p)
		self.recording_vectors(p)
		self.run_sims(p)

	# update clopath parameters
	def update_clopath(self, p, syns):
		# iterate over parameters
		for parameter_key,parameter in p.iteritems():
			# if it is a clopath learning rule parameter
			if 'clopath_' in parameter_key:
				# get parameter name
				p_clopath = parameter_key[len('clopath_'):]
				# iterate over sections
				for sec_i,sec in enumerate(syns):
					# iterate over segments
					for seg_i,seg in enumerate(syns[sec_i]):
						# set synapse parameter value
						setattr(syns[sec_i][seg_i], p_clopath, p['clopath_'+p_clopath])

	# activate synapses
	def activate_synapses(self,p):
		bipolar = stims.Bipolar()
		bipolar.tbs(bursts=p['bursts'], warmup=p['warmup'], pulses=p['pulses'], pulse_freq=p['pulse_freq'])
		self.stim = bipolar.stim
		self.nc = cell.Syn_act(p=p, syns=self.cell1.syns, stim=self.stim)

	def shape_plot(self,p):
		# highlight active sections
		self.shapeplot = h.PlotShape()
		
		# create section list of active sections
		self.sl = h.SectionList()    # secetion list of included sections
		for sec_i,sec in enumerate(p['sec_idx']):
			self.sl.append(sec=self.cell1.geo[p['tree']][sec])
			self.shapeplot.color(2, sec=self.cell1.geo[p['tree']][sec])

	def recording_vectors(self,p):
		# set up recording vectors
		self.rec =  {}
		self.data = {}
		
		# loop over trees
		for tree_key, tree in self.cell1.geo.iteritems():
			self.rec[tree_key+'_v'] = []
			self.rec[tree_key+'_w'] = []
			self.data[tree_key + '_v'] = []
			self.data[tree_key + '_w'] = []
			
			# loop over sections
			for sec_i,sec in enumerate(tree):
				self.rec[tree_key+'_v'].append([])
				self.rec[tree_key+'_w'].append([])
				
				# loop over segments
				for seg_i,seg in enumerate(tree[sec_i]):
					
					# determine relative segment location in (0-1) 
					seg_loc = (seg_i+1)/(self.cell1.geo[tree_key][sec_i].nseg+1)
					# record voltage
					self.rec[tree_key+'_v'][sec_i].append(h.Vector())
					self.rec[tree_key+'_v'][sec_i][seg_i].record(
						tree[sec_i](seg_loc)._ref_v)
					
					# if segment is in dendrite
					if ((tree_key == 'basal') or 
					(tree_key == 'apical_trunk') or 
					(tree_key == 'apical_tuft')):
						
						# record clopath weight
						self.rec[tree_key+'_w'][sec_i].append(h.Vector())
						self.rec[tree_key+'_w'][sec_i][seg_i].record(
							self.cell1.syns[tree_key]['clopath'][sec_i][seg_i]._ref_gbar)
		# time
		self.data['t'] = []
		self.rec['t'] = h.Vector()
		self.rec['t'].record(h._ref_t)

	def run_sims(self,p):
		# data organized as ['tree']['polarity'][section][segment]
		# loop over dcs fields
		for f_i,f in enumerate(p['field']):

			# insert extracellular field
			stims.DCS(cell=0, field_angle=p['field_angle'], intensity=f)
			
			# run time
			h.dt = p['dt']
			h.tstop = p['tstop']

			# run simulation
			h.run()

			# store recording vectors as arrays
			# loop over trees
			for tree_key,tree in self.rec.iteritems():
				# add list for each field polarity
				self.data[tree_key].append([])

				if tree_key != 't':
					# loop over sections
					for sec_i,sec in enumerate(self.rec[tree_key]):
						self.data[tree_key][f_i].append([])
						
						# loop over segments
						for seg_i,seg in enumerate(sec):
							if '_v' in tree_key:
								self.data[tree_key][f_i][sec_i].append(np.array(self.rec[tree_key][sec_i][seg_i]))
							
							if ((tree_key == 'basal_w') or 
							(tree_key == 'apical_trunk_w') or 
							(tree_key == 'apical_tuft_w')):
								self.data[tree_key][f_i][sec_i].append(np.array(self.rec[tree_key][sec_i][seg_i]))

			self.data['t'].append([])
			self.data['t'][f_i] = np.array(self.rec['t'])
		self.data['p'] = p

def save_data(data):	# save data
	p = data['p']
	# delete cell hoc object (can't be pickled)
	p['cell']=[]
	# check if folder exists with experiment name
	if os.path.isdir(p['data_folder']) is False:
		os.mkdir(p['data_folder'])

	with open(p['data_folder']+'data_'+
		p['experiment']+
		'_trial_'+str(p['trial'])+
		'_weight_'+str(p['w_mean'])+
		'_synfrac_'+str(p['syn_frac'])+
		'_'+p['trial_id']+
		'.pkl', 'wb') as output:

		pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

# procedures to be initialized if called as a script
if __name__ =="__main__":
	plot_sections(None,None)

