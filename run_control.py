
"""
run control
"""
# imports
from neuron import h
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import cell
# import itertools as it
# import stims
# import pickle
import param
import run
import time
h.load_file("stdrun.hoc")


def exp_2(num_secs,trials,weights):
	for tri in range(trials):
		for secs in num_secs:
			for w in weights:
				p = param.exp_2(num_sec=secs).params
				print p['sec_idx']
				p['trial']=tri
				p['w_ampa']=w*p['w_ampa']
				p['w_nmda']=w*p['w_nmda']
				run.run(p)
				data_folder = 'Data/'
				run.plot_sections(data_folder+'data_'+p['experiment']+'_syn_'+str(len(p['sec_idx']))+
					'_trial_'+str(p['trial'])+'_weight_'+str(p['w_ampa'])+'.pkl')

def exp_3(trials,weights):
	# loop over trials
	for tri in range(trials):
		# loop over weights
		for w in weights:
			# choose fraction of synapses to be activated
			syn_frac = np.random.normal(loc=.4, scale=.1) # chosen from gaussian
			# load rest of parameters from parameter module
			p = param.exp_3(syn_frac=syn_frac).params
			p['trial']=tri
			p['w_ampa']=w*p['w_ampa']
			p['w_nmda']=w*p['w_nmda']
			p['w_rand']=True
			start = time.time() # start timer
			sim = run.Run(p)	# run simulation
			end = time.time() # end timer
			print 'trial'+ str(tri) + ' duration:' + str(end -start) # print simulation time
			run.save_data(sim.data,p)

if __name__ =="__main__":
	exp_3(trials=200,weights = [1])
