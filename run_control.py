
"""
run control
"""
# imports
from neuron import h
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import cell
# import itertools as it
# import stims
# import pickle
import param
import run
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
	for tri in range(trials):
		for w in weights:
			p = param.exp_3(syn_frac=0.4).params
			p['trial']=tri
			p['w_ampa']=w*p['w_ampa']
			p['w_nmda']=w*p['w_nmda']
			run.run(p)
			data_folder = 'Data/'
			run.plot_sections(data_folder+'data_'+p['experiment']+
			'_trial_'+str(p['trial'])+'_weight_'+str(p['w_ampa'])+'.pkl')
		



if __name__ =="__main__":
	exp_3(trials=1,weights = [12])
