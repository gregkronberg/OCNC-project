"""
parameter optimization
"""
import param
import run
import analysis
import numpy as np
import scipy.optimize
import uuid
import shutil
import os
import pickle

class Optimize:
	def __init__(self):
		self.x0 = [1,1,1]
		self.upper = (4,10,10)
		self.lower = (.5,.1,.01)
		self.trials = 3
		self.exp = 'exp_opt'
		self.target=np.array([1.46,1.25,1.22])
		self.x_history = []
		self.xtol=0.01
		self.ftol=0.01
		self.maxiter = 100
		self.opt = scipy.optimize.fmin(self.error_func,self.x0,args=(self.trials, self.exp, self.target), xtol=self.xtol, ftol = self.ftol, maxiter=self.maxiter, callback=self.store_params)
		self.save_params()

	def error_func(self,x,*args):
		syn_frac = 0.1*x[0]
		weights = .001*x[1]
		A_p = .01*x[2]
		trials = args[0]
		exp = args[1]
		target = args[2]

		p = param.exp_3(syn_frac=syn_frac, w_mean=weights, w_rand=False, exp=exp).p

		for trial in range(trials):
			p['clopath_A_p'] = A_p
			p['trial'] = trial
			p['trial_id'] = str(uuid.uuid4())
			# run simulation
			sim = run.Run(p)
			run.save_data(sim.data)
			# print trial and simulation time
			print 'trial'+ str(trial) 

		W = analysis.Weights(p)
		dw_mean = np.divide(np.mean(W.w_end_all,axis=1),np.mean(W.w_start_all,axis=1))

		error = 0
		for dw_i,dw in enumerate(dw_mean):
			error += (dw - target[dw_i])**2	

		# clear directory
		directory = p['data_folder']
		if os.path.isdir(directory) is True:
			for file in os.listdir(directory):
				try:
					os.remove(os.path.join(directory,file))
				except OSError, e:  ## if failed, report it back to the user ##
	        			print ("Error: %s - %s." % (e.filename,e.strerror))
		

		print x
		print dw_mean
		return error

	def store_params(self, xk):
		self.x_history.append(xk)

	def save_params(self):
		directory = 'Data/'+self.exp+'/'
		with open(directory+'parameters_'+'.pkl', 'wb') as output:
			pickle.dump(self.opt, output,protocol=pickle.HIGHEST_PROTOCOL)

if __name__=='__main__':
	Optimize()

