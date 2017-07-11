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
				p = param.exp_2(secs).params
				print secs
				# print p['sec_idx']
				p['trial']=tri
				p['w_ampa']=w*p['w_ampa']
				p['w_nmda']=w*p['w_nmda']
				print p['w_nmda']
				run.run(p)
				run.plot_sections('data_'+p['experiment']+'_trial_'+str(p['trial'])+'.pkl')

if __name__ =="__main__":
	exp_2(num_secs=[10],trials=1,weights = [5])
