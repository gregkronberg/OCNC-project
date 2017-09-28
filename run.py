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

# load standard runtime settings.  this is critical.
h.load_file("stdrun.hoc")

# run control
def run(p):
	"""
	Arguments list:
	p - dictionary of parameters

	each experiment appends a list to the appropriate key in the data dictionary
	data is organized as data['type of data'][experiments][sections][time series data vector]
	details of each experiment are tracked via the data['detail'][experiment number], e.g. data['field'][2]
	"""

	# initialize data
	data = {'t':[],
		'soma':[],
		'dend':[],
		'weight':[],
		'field':[],
		'field_color':[],
		'params':p,
			}

	# create cell
	cell1 = cell.Cell_Migliore_2005()

	# activate synapses
	tuft_act = []	# list of activations (each burst x subtree combination gets an entry)
	for syn_stim_i,syn_stim in enumerate(stims.tbs(bursts=p['bursts']).stim):

	# activate synapses
		tuft_act.append(cell.Syn_act(syn_stim,cell1.syn_a_tuft_ampa,cell1.syn_a_tuft_nmda,cell1.syn_a_tuft_clopath,p['sec_idx'],p['seg_idx'],p['w_ampa'],p['w_nmda']))

	# highlight active sections
	shapeplot = h.PlotShape()
	
	# create section list of active sections
	sl = h.SectionList()    # secetion list of included sections
	for section in p['sec_idx']:
		sl.append(sec = cell1.dend_a_tuft[section])
		shapeplot.color(2,sec=cell1.dend_a_tuft[section])
	
	# run time
	h.dt = p['dt']
	h.tstop = p['tstop']

	# set up recording vectors
	# time
	t_rec = h.Vector()
	t_rec.record(h._ref_t)
	
	# soma voltage
	soma_rec = h.Vector()
	soma_rec.record(cell1.soma(0.5)._ref_v)
	
	# dendrite voltage (sections chosen with 'plot_sec_idx' in parameter module)
	# plot_sec_idx is a list organized as [sections]
	# plot_seg_idx is [sections][segments]
	dend_rec=[]
	for sec_i,sec in enumerate(p['plot_sec_idx']):
		dend_rec.append([])
		for seg_i,seg in enumerate(p['plot_seg_idx'][sec_i]):
			# check if segment exists
			if seg <= cell1.dend_a_tuft[sec].nseg:

				# determine relative segment location in (0-1) 
				seg_loc = (seg+1)/(cell1.dend_a_tuft[sec].nseg+1)
				dend_rec[sec_i].append(h.Vector())
				dend_rec[sec_i][seg_i].record(cell1.dend_a_tuft[sec](seg_loc)._ref_v)

	# clopath weight update
	weight_rec=[]
	for sec_i,sec in enumerate(p['plot_sec_idx']):
		weight_rec.append([])
		for seg_i,seg in enumerate(p['plot_seg_idx'][sec_i]):
			if seg < len(cell1.syn_a_tuft_clopath[sec]):
				weight_rec[sec_i].append(h.Vector())
				weight_rec[sec_i][seg_i].record(cell1.syn_a_tuft_clopath[sec][seg]._ref_gbar)
	
	# loop over dcs fields
	cnt=-1
	shapeplot = []
	for f_i,f in enumerate(p['field']):
		cnt +=1
		
		# insert extracellular field
		stims.dcs(cell=0,field_angle=p['field_angle'],intensity=f)
		
		# run simulation
		h.run()

		# store data
		data['t'].append(np.array(t_rec))
		data['soma'].append(np.array(soma_rec))
		data['dend'].append(np.array(dend_rec))
		data['weight'].append(np.array(weight_rec))
		data['field'].append(f)
		data['field_color'].append(p['field_color'][f_i])	

	# save data
	data_folder = 'Data/'
	with open(data_folder+'data_'+p['experiment']+'_trial_'+str(p['trial'])+'_weight_'+str(p['w_ampa'])
		+'.pkl', 'wb') as output:
	# 
		pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

# plot control
def plot_sections(data_file):
	"""
	creates a single figure with subplots for each section to be recorded from
	"""
	# open data
	pkl_file = open(data_file, 'rb')
	
	data = pickle.load(pkl_file)

	plot_folder = 'png figures/'
	
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
		if k < n_seg-1:
			axh = "section-{:03d}".format(data['params']['sec_idx'][k])
			ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
			ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
			for exp in range(len(data['t'])):
				plot_color = data['field_color'][exp]
				# ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
				# plot dendritic voltage
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['dend'][exp][k]),plot_color)
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

	fig.savefig(plot_folder+data['params']['experiment']+'_syn_'+str(len(data['params']['sec_idx']))+
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

