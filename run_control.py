
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
import uuid
import analysis
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

def exp_3(trials=1,weights = [[.0018,.0002]]):
	# loop over trials
	for tri in range(trials):
		# loop over weights
		for w in weights:
			# choose fraction of synapses to be activated
			syn_frac = np.random.normal(loc=.1, scale=.1) # chosen from gaussian
			
			# load rest of parameters from parameter module
			p = param.exp_3(syn_frac=syn_frac, w_mean=w[0], w_std=w[1], w_rand=True, exp='exp_5').p
			
			print p['seg_idx']
			# store trial number
			p['trial']=tri
			
			# create unique identifier for each trial
			p['trial_id'] = str(uuid.uuid4())
			
			# start timer
			start = time.time() 
			
			# run simulation
			sim = run.Run(p)	

			# end timer
			end = time.time() 

			# print trial and simulation time
			print 'trial'+ str(tri) + ' duration:' + str(end -start) 
			
			# save data for eahc trial
			run.save_data(sim.data)

if __name__ =="__main__":
	exp_3(trials=1,weights = [[.002,.002]])
	dw = analysis.Weights(param.exp_3(exp='exp_5').p)
	print np.mean(dw.w_end_all,axis=1)/.05 
	# analysis.Spikes(param.exp_3().p)
	# analysis.Voltage(param.exp_3().p)
