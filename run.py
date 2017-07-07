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
import itertools as iter
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
	
	# run time
	h.dt = .025
	h.tstop = 60
	# loop over dcs fields
	cnt=-1
	for f_i,f in enumerate(p['field']):
		cnt +=1
		cell1 = cell.Cell()
		# activate synapses
		tuft_act  = cell.Syn_act(cell1.syn_a_tuft_ampa,cell1.syn_a_tuft_nmda,p['sec_idx'],p['seg_idx'],p['w_ampa'],p['w_nmda'])
		# insert extracellular field
		stims.dcs(cell=0,field_angle=p['field_angle'],intensity=f)
		# create data vectors
		data['t'].append(h.Vector())
		data['t'][cnt].record(h._ref_t)
		data['soma'].append(h.Vector())
		data['soma'][cnt].record(cell1.soma(0.5)._ref_v)
		data['dend'].append([])
		data['field'].append(f)
		data['field_color'].append(p['field_color'][f_i])
		# loop over dendritic sections that were stimulated
		for sec_i,sec in enumerate(p['plot_idx']):
			# add recording vector for each dendritic section
			data['dend'][cnt].append(h.Vector())
			data['dend'][cnt][sec_i].record(cell1.dend_a_tuft[sec](0)._ref_v)
		for sec in h.allsec():
			print sec(0.5).e_extracellular
		h.run()

	# convert to numpy arrays
	for vec in data:
		data[vec] = np.array(data[vec])

def create_plot(data,name):
    n_sec = len(p['plot_idx'])

    n = int(np.ceil(np.sqrt(n_sec)))
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(n, n, wspace=0.10, hspace=0.05, left=0.1, right=0.95, bottom=0.1, top=0.95)

    ax = {}
    rows = np.arange(0, n, 1, dtype=int)
    cols = np.arange(0, n, 1, dtype=int)

    for k, (i, j) in enumerate(it.product(rows, cols)):
    	print(i, j)
    	axh = "test-{:03d}".format(k)
    	ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
    	ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
    	if k < n_sec:
    		for exp in range(len(data['t'])):
	    		plot_color = data['field_color'][exp]
		    	ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['soma'][exp]),plot_color)
		    	ax[axh].plot(np.transpose(data['t'][exp]), np.transpose(data['dend'][exp][k]),plot_color)

    fig.savefig('testing.png', dpi=200)
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
