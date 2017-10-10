"""
parameter optimization
"""
import param
import run
import analysis
import numpy as np
from scipy.optimize import curve_fit
import uuid
import shutil
import os

def func(fields, syn_frac, weights, A_p):
	directory = 'Data/exp_8/'
	if os.path.isdir(directory) is True:
		for file in os.listdir(directory):
			os.remove(os.path.join(directory,file))
		# os.rmdir(directory)
	p = param.exp_3(syn_frac=.1*syn_frac,w_mean =.001*weights,w_rand=False,exp='exp_8').p
	# field_i = [field_i for field_i,field in enumerate(p['field']) if fields==field]  
	for trial in range(3):
		p['clopath_A_p'] = .01*A_p
		p['trial'] = trial
		p['trial_id'] = str(uuid.uuid4())
		# run simulation
		sim = run.Run(p)
		run.save_data(sim.data)
		# print trial and simulation time
		print 'trial'+ str(trial) 

	W = analysis.Weights(p)
	dw_mean = np.divide(np.mean(W.w_end_all,axis=1),np.mean(W.w_start_all,axis=1))
	print dw_mean
	return dw_mean

p0 = [1,1,1]
upper = (4,10,10)
lower = (.5,.1,.01)
bounds = [lower,upper]
xdata = [-20,0,20]
ydata = np.array([1.46,1.25,1.22])
p_opt,p_cov = curve_fit(func,xdata,ydata,p0=p0,bounds=bounds)

if __name__=='__main__':
	directory = 'Data/exp_8/'
	p0 = [1,1,1]
	upper = (4,10,10)
	lower = (.5,.1,.01)
	bounds = [lower,upper]
	xdata = [-20,0,20]
	ydata = np.array([1.46,1.25,1.22])
	p_opt,p_cov = curve_fit(func,xdata,ydata,p0=p0,bounds=bounds)
	with open(directory+'parameters_'+'.pkl', 'wb') as output:
			pickle.dump(p_opt, output,protocol=pickle.HIGHEST_PROTOCOL)

