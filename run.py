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

	# synaptic stimulation
	syn_stim = stims.tbs().stim[0]

	# activate synapses
	tuft_act  = cell.Syn_act(syn_stim,cell1.syn_a_tuft_ampa,cell1.syn_a_tuft_nmda,cell1.syn_a_tuft_clopath,p['sec_idx'],p['seg_idx'],p['w_ampa'],p['w_nmda'])

	# highlight active sections
	shapeplot = h.PlotShape()
	
	# create section list of active sections
	sl = h.SectionList()    # secetion list of included sections
	for section in p['sec_idx']:
		sl.append(sec = cell1.dend_a_tuft[section])
		shapeplot.color(2,sec=cell1.dend_a_tuft[section])

	# sl.printnames()
	
	# run time
	h.dt = .025
	h.tstop = 100

	#======================================================================================
	# setup recording vectors
	#======================================================================================
	# time
	t_rec = h.Vector()
	t_rec.record(h._ref_t)
	
	# soma voltage
	soma_rec = h.Vector()
	soma_rec.record(cell1.soma(0.5)._ref_v)
	
	# dendrite voltage (sections chosen with 'plot_idx' in parameter module)
	dend_rec=[]
	for sec_i,sec in enumerate(p['plot_idx']):
		dend_rec.append(h.Vector())
		dend_rec[sec_i].record(cell1.dend_a_tuft[sec](0)._ref_v)

	# clopath weight update
	weight_rec=[]
	for sec_i,sec in enumerate(p['plot_idx']):
		weight_rec.append(h.Vector())
		weight_rec[sec_i].record(cell1.syn_a_tuft_clopath[sec][0].gbar)
	
	# loop over dcs fields
	cnt=-1
	shapeplot = []
	for f_i,f in enumerate(p['field']):
		cnt +=1
		# insert extracellular field
		stims.dcs(cell=0,field_angle=p['field_angle'],intensity=f)
		
		# shape plots
		# shapeplot.append(h.PlotShape())
		# shapeplot[cnt].color_list(sl,1)

		# shapeplot[cnt].variable('v')
		# h.fast_flush_list.append(shapeplot[cnt])
		# shapeplot[cnt].exec_menu('View = plot')
		# shapeplot[cnt].exec_menu('Shape Plot')
		# shapeplot[cnt].scale(-65, 10)
		# shapeplot[cnt].colormap(1,255,255,0)
		
		# shapeplot[cnt].fastflush()

		# pause simulation to view shape plot at specific time
		# h.load_file('interrupts_shapeflush.hoc')
		
		# run simulation
		h.run()

		# store data
		data['t'].append(np.array(t_rec))
		data['soma'].append(np.array(soma_rec))
		data['dend'].append(np.array(dend_rec))
		data['weight'].append(np.array(weight_rec))
		data['field'].append(f)
		data['field_color'].append(p['field_color'][f_i])	
	
		# for sec in h.allsec():
		# 	print sec(0.5).e_extracellular

	# save data
	data_folder = 'Data/'
	with open(data_folder+'data_'+p['experiment']+'_syn_'+str(len(p['sec_idx']))+'_trial_'+str(p['trial'])+'_weight_'+str(p['w_ampa'])
		+'.pkl', 'wb') as output:
	# 
		pickle.dump(data, output,protocol=pickle.HIGHEST_PROTOCOL)

	# save shape plot
	# pwm = h.PWManager()
	# pwm.printfile('shapeplot_'+p['experiment']+'_trial_'+str(p['trial'])+'.eps', 0, 0)

def plot_sections(data_file):
	"""
	creates a single figure with subplots for each section to be recorded from
	"""
	# open data
	pkl_file = open(data_file, 'rb')
	
	data = pickle.load(pkl_file)

	plot_folder = 'png figures/'
	
	n_sec = len(data['params']['plot_idx'])+1
	n = int(np.ceil(np.sqrt(n_sec)))
	fig = plt.figure(figsize=(10, 10))
	gs = gridspec.GridSpec(n, n, wspace=0.10, hspace=0.05, left=0.1, right=0.95, bottom=0.1, top=0.95)
	ax = {}
	rows = np.arange(0, n, 1, dtype=int)
	cols = np.arange(0, n, 1, dtype=int)
	for k, (i, j) in enumerate(it.product(rows, cols)):
		if k < n_sec-1:
			axh = "section-{:03d}".format(data['params']['sec_idx'][k])
			ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
			ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
			for exp in range(len(data['t'])):
				plot_color = data['field_color'][exp]
				# ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['dend'][exp][k]),plot_color)
		elif k==n_sec-1:
			axh = "section-{}".format('soma')
			ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
			ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
			for exp in range(len(data['t'])):
				plot_color = data['field_color'][exp]
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
				# ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['dend'][exp][k]),plot_color)

	fig.savefig(plot_folder+data['params']['experiment']+'_syn_'+str(len(data['params']['sec_idx']))+
		'_trial_'+str(data['params']['trial'])+
		'_weight_'+str(data['params']['w_ampa'])+'.png', dpi=250)
	plt.close(fig)

# procedures to be initialized if called as a script
if __name__ =="__main__":
	plot_sections(None,None)

