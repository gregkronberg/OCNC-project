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
	def __init__(self,p):

		# create cell
		self.cell1 = cell.Cell_Migliore_2005(p)
		self.update_clopath(p,syns=self.cell1.syns['apical_tuft']['clopath'])
		self.activate_synapses(p)
		self.recording_vectors(p)
		self.run_sims(p)

	# update clopath parameters
	def update_clopath(self,p,syns):
		for sec_i,sec in enumerate(syns):
			for seg_i,seg in enumerate(syns[sec_i]):
				syns[sec_i][seg_i].delay_steps = p['clopath_delay_steps']
				syns[sec_i][seg_i].tau_0 = p['clopath_tau_0']
				syns[sec_i][seg_i].tau_r = p['clopath_tau_r']
				syns[sec_i][seg_i].tau_y = p['clopath_tau_y']
				syns[sec_i][seg_i].A_m = p['clopath_A_m']
				syns[sec_i][seg_i].A_p = p['clopath_A_p']
				syns[sec_i][seg_i].tetam = p['clopath_tetam']
				syns[sec_i][seg_i].tetap = p['clopath_tetap']

	# activate synapses
	def activate_synapses(self,p):
		self.stim = stims.tbs(bursts=p['bursts'],warmup=p['warmup'],pulses=p['pulses']).stim
		self.nc = cell.Syn_act(p=p, syns=self.cell1.syns, stim=self.stim)

	def shape_plot(self,p):
		# highlight active sections
		self.shapeplot = h.PlotShape()
		
		# create section list of active sections
		self.sl = h.SectionList()    # secetion list of included sections
		for seec_i,sec in enumerate(p['seg_idx']):
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
		
		# dendrite voltage (sections chosen with 'plot_sec_idx' in parameter module)
		# plot_sec_idx is a list organized as [sections]
		# plot_seg_idx is [sections][segments]

	def run_sims(self,p):
		# loop over dcs fields
		for f_i,f in enumerate(p['field']):

			
			# insert extracellular field
			stims.dcs(cell=0,field_angle=p['field_angle'],intensity=f)
			
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
	if os.path.isdir(p['data_folder']) is False:
		os.mkdir(p['data_folder'])

	with open(p['data_folder']+
		'data_'+p['experiment']+
		'_trial_'+str(p['trial'])+
		'_weight_'+str(p['w_mean'])
		+'_synfrac_'+str(p['syn_frac'])+
		'.pkl', 'wb') as output:

		pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

# plot control
def plot_sections(data_file,fig_folder):
	"""
	creates a single figure with subplots for each section to be recorded from
	"""
	# open data
	data_file = p['data_folder']+'data_'+p['experiment']+'_trial_'+str(p['trial'])+'_weight_'+str(p['w_ampa'])

	pkl_file = open(data_file, 'rb')
	
	data = pickle.load(pkl_file)
	
	# number of segments to plot
	n_seg = len(data['params']['plot_seg_idx'])+1
	# number of elements n in each dimension of n x n plot grid 
	n = int(np.ceil(np.sqrt(n_seg)))
	# create figure
	fig = plt.figure(figsize=(10, 10))
	# setup figure grid
	gs = gridspec.GridSpec(n, n, wspace=0.10, hspace=0.05, left=0.1, right=0.95, bottom=0.1, top=0.95)
	# dictionary for storing individual plots
	ax = {}
	# rows and columns
	rows = np.arange(0, n, 1, dtype=int)
	cols = np.arange(0, n, 1, dtype=int)
	# loop over grid elements
	for k, (i, j) in enumerate(it.product(rows, cols)):
		if k < n_seg-2:
			axh = "section-{:03d}".format(data['params']['sec_idx'][k])
			ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
			ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
			for exp in range(len(data['t'])):
				plot_color = data['field_color'][exp]
				# ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
				# plot dendritic voltage
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['dend'][exp][k]),plot_color)
		
		elif k==n_seg-2:
			axh = "section-{}".format('soma')
			ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
			ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
			for exp in range(len(data['t'])):
				plot_color = data['field_color'][exp]
				# plot soma voltage
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
				# plot weight changes
				# ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['weight'][exp][0]),plot_color)


		elif k==n_seg-1:
			axh = "section-{}".format('soma')
			ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
			ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
			for exp in range(len(data['t'])):
				plot_color = data['field_color'][exp]
				# plot soma voltage
				# ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
				# plot weight changes
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['weight'][exp][0]),plot_color)
				# plot dendritic voltage
				# ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['dend'][exp][k]),plot_color)

	fig.savefig(p['fig_folder']+data['params']['experiment']+'_syn_'+str(len(data['params']['sec_idx']))+
		'_trial_'+str(data['params']['trial'])+
		'_weight_'+str(data['params']['w_ampa'])+'.png', dpi=250)
	plt.close(fig)

# create and save shape plot
def shapeplot():
		# shape plots
		shapeplot.append(h.PlotShape())
		shapeplot[cnt].color_list(sl,1)

		shapeplot[cnt].variable('v')
		h.fast_flush_list.append(shapeplot[cnt])
		shapeplot[cnt].exec_menu('View = plot')
		shapeplot[cnt].exec_menu('Shape Plot')
		shapeplot[cnt].scale(-65, 10)
		shapeplot[cnt].colormap(1,255,255,0)
		
		shapeplot[cnt].fastflush()

		# pause simulation to view shape plot at specific time
		h.load_file('interrupts_shapeflush.hoc')

		# save shape plot
		pwm = h.PWManager()
		pwm.printfile('shapeplot_'+p['experiment']+'_trial_'+str(p['trial'])+'.eps', 0, 0)

# procedures to be initialized if called as a script
if __name__ =="__main__":
	plot_sections(None,None)

