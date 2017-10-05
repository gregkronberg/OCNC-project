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

		# initialize data
		self.data = {'t':[],
			'soma':[],
			'dend':[],
			'weight':[],
			'field':[],
			'field_color':[],
			'params':p,
				}


		# create cell
		self.cell1 = cell.Cell_Migliore_2005()

		# update clopath parameters
		for a in range(len(self.cell1.syn_a_tuft_clopath)):
			for b in range(len(self.cell1.syn_a_tuft_clopath[a])):
				self.cell1.syn_a_tuft_clopath[a][b].delay_steps = p['clopath_delay_steps']
				self.cell1.syn_a_tuft_clopath[a][b].tau_0 = p['clopath_tau_0']
				self.cell1.syn_a_tuft_clopath[a][b].tau_r = p['clopath_tau_r']
				self.cell1.syn_a_tuft_clopath[a][b].tau_y = p['clopath_tau_y']
				self.cell1.syn_a_tuft_clopath[a][b].A_m = p['clopath_A_m']
				self.cell1.syn_a_tuft_clopath[a][b].A_p = p['clopath_A_p']
				self.cell1.syn_a_tuft_clopath[a][b].tetam = p['clopath_tetam']
				self.cell1.syn_a_tuft_clopath[a][b].tetap = p['clopath_tetap']

		# activate synapses
		self.tuft_act = []	# list of activations (each burst x subtree combination gets an entry)
		for syn_stim_i,syn_stim in enumerate(stims.tbs(bursts=p['bursts']).stim):
			# randomize weights
			if p['w_rand']:
			# activate synapses
				self.tuft_act.append(cell.Syn_act(syn_stim,self.cell1.syn_a_tuft_ampa,self.cell1.syn_a_tuft_nmda,self.cell1.syn_a_tuft_clopath,p['sec_idx'],p['seg_idx'],np.random.normal(loc=p['w_ampa'],scale=p['w_std']),np.random.normal(loc=p['w_nmda'],scale=p['w_std'])))
			# dont randomize weights
			else:
				# activate synapses
				self.tuft_act.append(cell.Syn_act(syn_stim,self.cell1.syn_a_tuft_ampa,self.cell1.syn_a_tuft_nmda,self.cell1.syn_a_tuft_clopath,p['sec_idx'],p['seg_idx'],p['w_ampa'],p['w_nmda']))


		# highlight active sections
		self.shapeplot = h.PlotShape()
		
		# create section list of active sections
		self.sl = h.SectionList()    # secetion list of included sections
		for section in p['sec_idx']:
			self.sl.append(sec = self.cell1.dend_a_tuft[section])
			self.shapeplot.color(2,sec = self.cell1.dend_a_tuft[section])
		
		# run time
		h.dt = p['dt']
		h.tstop = p['tstop']

		# set up recording vectors
		# time
		self.t_rec = h.Vector()
		self.t_rec.record(h._ref_t)
		
		# soma voltage
		self.soma_rec = h.Vector()
		self.soma_rec.record(self.cell1.soma(0.5)._ref_v)
		
		# dendrite voltage (sections chosen with 'plot_sec_idx' in parameter module)
		# plot_sec_idx is a list organized as [sections]
		# plot_seg_idx is [sections][segments]
		
		
		# loop over dcs fields
		cnt=-1
		shapeplot = []
		self.dend_rec=[]
		self.weight_rec=[]
		self.dend_arr=[]
		self.weight_arr=[]
		for f_i,f in enumerate(p['field']):
			cnt +=1
			self.dend_rec.append([])
			self.dend_arr.append([])
			self.weight_rec.append([])
			self.weight_arr.append([])
			for sec_i,sec in enumerate(p['plot_sec_idx']):
				self.dend_rec[cnt].append([])
				self.dend_arr[cnt].append([])
				self.weight_rec[cnt].append([])
				self.weight_arr[cnt].append([])
				for seg_i,seg in enumerate(p['plot_seg_idx'][sec_i]):
					# check if segment exists
					if seg <= self.cell1.dend_a_tuft[sec].nseg:

						# determine relative segment location in (0-1) 
						seg_loc = (seg+1)/(self.cell1.dend_a_tuft[sec].nseg+1)
						# print cnt
						self.dend_rec[cnt][sec_i].append(h.Vector())
						self.dend_arr[cnt][sec_i].append([])
						self.dend_rec[cnt][sec_i][seg_i].record(self.cell1.dend_a_tuft[sec](seg_loc)._ref_v)

					if seg < len(self.cell1.syn_a_tuft_clopath[sec]):
						self.weight_rec[cnt][sec_i].append(h.Vector())
						self.weight_arr[cnt][sec_i].append([])
						self.weight_rec[cnt][sec_i][seg_i].record(self.cell1.syn_a_tuft_clopath[sec][seg]._ref_gbar)

			# insert extracellular field
			stims.dcs(cell=0,field_angle=p['field_angle'],intensity=f)
			
			# run simulation
			h.run()

			# store recording vectors as arrays
			for a in range(len(self.dend_rec[cnt])):
				for b in range(len(self.dend_rec[cnt][a])):
					self.dend_arr[cnt][a][b] = np.array(self.dend_rec[cnt][a][b])
					self.weight_arr[cnt][a][b] = np.array(self.weight_rec[cnt][a][b])
			
			# store data
			self.data['t'].append(np.array(self.t_rec))
			self.data['soma'].append(np.array(self.soma_rec))
			self.data['field'].append(f)
			self.data['field_color'].append(p['field_color'][f_i])	
		
		self.data['dend'] = self.dend_arr
		self.data['weight'] = self.weight_arr

def save_data(data,p):	# save data
	if os.path.isdir(p['data_folder']) is False:
		os.mkdir(p['data_folder'])

	with open(p['data_folder']+'data_'+p['experiment']+'_trial_'+str(p['trial'])+'_weight_'+str(p['w_ampa'])
		+'_synfrac_'+str(p['syn_frac'])+'.pkl', 'wb') as output:

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

