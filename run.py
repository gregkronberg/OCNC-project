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
from param import params as p
h.load_file("stdrun.hoc")

# initialize data
data = {'t':[],
	'soma':[],
	'dend':[],
	'field':[],
	'field_color':[]
			}

# run control
def run():
	"""
	docstring

	each experiment appends each field in the data dictionary
	data is organized as data['type of data'][experiments][sections][data vector]
	details of each experiment can be accessed from data['params'][experiment number]
	"""
	global cell1, tuft_act, data
	# create cell
	cell1 = cell.Cell()
	# activate synapses
	tuft_act  = cell.Syn_act(cell1.syn_a_tuft_ampa,cell1.syn_a_tuft_nmda,p['sec_idx'],p['seg_idx'],p['w_ampa'],p['w_nmda'])
	
	# run time
	h.dt = .025
	h.tstop = 60

	# setup recording vectors
	t_rec = h.Vector()
	t_rec.record(h._ref_t)
	soma_rec = h.Vector()
	soma_rec.record(cell1.soma(0.5)._ref_v)
	dend_rec=[]
	for sec_i,sec in enumerate(p['plot_idx']):
		dend_rec.append(h.Vector())
		dend_rec[sec_i].record(cell1.dend_a_tuft[sec](0)._ref_v)	
	
	shapeplot = []
	# loop over dcs fields
	cnt=-1
	for f_i,f in enumerate(p['field']):
		cnt +=1
		# insert extracellular field
		stims.dcs(cell=0,field_angle=p['field_angle'],intensity=f)
		
		shapeplot.append(h.PlotShape())
		shapeplot[cnt].variable('v')
		shapeplot[cnt].exec_menu('Shape Plot')
		shapeplot[cnt].scale(-65, -50)
		shapeplot[cnt].flush()
		# run simulation
		h.run()

		# store data
		data['t'].append(np.array(t_rec))
		data['soma'].append(np.array(soma_rec))
		data['dend'].append(np.array(dend_rec))
		data['field'].append(f)
		data['field_color'].append(p['field_color'][f_i])	
	
		for sec in h.allsec():
			print sec(0.5).e_extracellular

	# pickle data
	output = open('data_'+p['experiment']+'_.pkl', 'wb')
	# 
	pickle.dump(data, output)

	output.close()
		

	# # convert to numpy arrays
	# for vec in data:
	# 	data[vec] = np.array(data[vec])

def create_plot(data_file):
	pkl_file = open(data_file, 'rb')
	
	data = pickle.load(pkl_file)
	
	n_sec = len(p['plot_idx'])
	n = int(np.ceil(np.sqrt(n_sec)))
	fig = plt.figure(figsize=(10, 10))
	gs = gridspec.GridSpec(n, n, wspace=0.10, hspace=0.05, left=0.1, right=0.95, bottom=0.1, top=0.95)
	ax = {}
	rows = np.arange(0, n, 1, dtype=int)
	cols = np.arange(0, n, 1, dtype=int)
	for k, (i, j) in enumerate(it.product(rows, cols)):
		if k < n_sec:
			print(i, j)
			axh = "section-{:03d}".format(p['sec_idx'][k])
			ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
			ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
			for exp in range(len(data['t'])):
				plot_color = data['field_color'][exp]
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
				ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['dend'][exp][k]),plot_color)

	fig.savefig(p['experiment']+'.png', dpi=200)
	plt.close(fig)


    	# axes are listed here
    # for sec_i,sec in enumerate(p['plot_idx']):
    # 	# create the figure object
    # 	fig.append(plt.figure()) 
    # 	fpng.append('name'+'_section_'+str(sec)+'.png')
    	
    # 	fig[sec_i].savefig(fpng[sec_i], dpi=250)
    # 	print("Figure {} saved.".format(fpng))
    # 	# and finally do some clean up
    # 	plt.show(fig[sec_i])




    

if __name__ =="__main__":
	# run()
	create_plot(None,None)
# plt.plot(t_arr,soma_arr,'k')
# plt.plot(t_arr,dend_arr)
# plt.show(fig)
