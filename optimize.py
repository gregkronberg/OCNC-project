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
	"""
	"""
	def __init__(self):
		# initial parameter guess
		self.x0 = [1,1,1]
		# upper parameter bound 
		self.upper = (4,10,10)
		# lower bound
		self.lower = (.5,.1,.01)
		# number of trials to average over
		self.trials = 5
		# folder to store data
		self.exp = 'exp_opt'
		# experimental data to fit
		self.target=np.array([1.46,1.25,1.22])
		# list to store parameters after each iteration
		self.x_history = []
		# tolerance for parameter changes
		self.xtol=0.01
		# tolerance for error function
		self.ftol=0.01
		# maximum number of iterations
		self.maxiter = 1000
		# run optimization
		self.opt = scipy.optimize.fmin(self.error_func,self.x0,
			args=(self.trials, self.exp, self.target), 
			xtol=self.xtol, 
			ftol = self.ftol, 
			maxiter=self.maxiter, 
			callback=self.store_params)
		# save optimization result
		self.save_params()

	def error_func(self,x,*args):
		""" error function 
		"""
		# scaled parameters 
		syn_frac = 0.1*x[0]
		weights = .001*x[1]
		A_p = .01*x[2]
		
		# number of trials
		trials = args[0]
		# name of experiment for storing data
		exp = args[1]
		# experimental data to be fit
		target = args[2]

		# load all simulation parameters
		p = param.exp_3(syn_frac=syn_frac, w_mean=weights, w_rand=False, exp=exp).p

		# clear existing data from directory
		# GOOGLE DRIVE MUST BE DISABLED, it interferes with file i/o
		directory = p['data_folder']
		# check if directory exists
		if os.path.isdir(directory) is True:
			# loop over files in directory
			for file in os.listdir(directory):
				# only delete raw data files
				if 'data' in file:
					try:
						# remove data
						os.remove(os.path.join(directory,file))
					# if failed, report it back to the user 
					except OSError, e:  
		        			print ("Error: %s - %s." % (e.filename,e.strerror))
		# loop over trials
		for trial in range(trials):
			# update parameters
			# clopath potentiation amplitude
			p['clopath_A_p'] = A_p
			#trial number
			p['trial'] = trial
			# unique identifier for each trial
			p['trial_id'] = str(uuid.uuid4())
			
			# run simulation
			sim = run.Run(p)
			
			# save data 
			run.save_data(sim.data)
			
			# print trial number
			print 'trial'+ str(trial) 

		# get weight changes
		W = analysis.Weights(p)
		# average overage overall all synapses, all trials
		dw_mean = np.divide(np.mean(W.w_end_all,axis=1),np.mean(W.w_start_all,axis=1))

		# compute least square error
		error = 0
		# loop over polarities
		for dw_i,dw in enumerate(dw_mean):
			# sum of square difference
			error += (dw - target[dw_i])**2	

		print x
		print dw_mean
		return error

	def store_params(self, xk):
		self.x_history.append(xk)
		directory = 'Data/'+self.exp+'/'
		with open(directory+'parameter_history'+'.pkl', 'wb') as output:
			pickle.dump(self.x_history, output,protocol=pickle.HIGHEST_PROTOCOL)

	def save_params(self):
		directory = 'Data/'+self.exp+'/'
		with open(directory+'parameter_optimum'+'.pkl', 'wb') as output:
			pickle.dump(self.opt, output,protocol=pickle.HIGHEST_PROTOCOL)

if __name__=='__main__':
	Optimize()

