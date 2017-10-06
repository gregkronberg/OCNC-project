"""
analysis
"""
from __future__ import division
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import itertools as it
import os
import cPickle as pickle
import param


class Weights():
	"""
	measure weight change at group of synapses

	saves initial and final weights at each active synapse across all simulated neurons (dcs polarity x synapses)
	"""
	def __init__(self,p):
		self.n_pol = len(p['field'])
		self.group_dw(p)
		self.save_dw(p)
		self.plot_dw_all(p,self.w_end)
		self.plot_dw_mean(p,self.w_end)

	def group_dw(self,p):
		# arrays for storing all weight changes across neurons
		self.w_end_all = np.empty([self.n_pol,0])
		self.w_start_all = np.empty([self.n_pol,0])
		# loop over experiments (neurons)
		for data_file in os.listdir(p['data_folder']):
			# check for proper data file format
			if 'data' in data_file:
				# load data file
				pkl_file = open(p['data_folder']+data_file, 'rb')
				data = pickle.load(pkl_file)

				self.p = data['params']
				
				self.n_act_seg = len(self.p['seg_list'])
				# measure weight changes for individual neuron
				self.measure_dw(data)
				# add individual neuron to the group
				self.w_end_all = np.append(self.w_end_all,self.w_end,axis=1)
				self.w_start_all = np.append(self.w_start_all,self.w_start,axis=1)

	def measure_dw(self,data):
		# set up arrays to record final and initial weights at each active synapse
		self.w_end = np.empty([self.n_pol,self.n_act_seg]) # polarity x segments
		self.w_start = np.empty([self.n_pol,self.n_act_seg]) # polarity x segments

		# find active synapses (all recorded segments were active)
		# measure weight change at each active synapse
		for a in range(len(data['weight'])): # loop fields
			cnt = -1
			for b in range(len(data['weight'][a])): # loop over sections
				for c in range(len(data['weight'][a][b])): # loop over segemnts
					cnt+=1
					self.w_end[a,cnt] = data['weight'][a][b][c][-1]
					self.w_start[a,cnt] = data['weight'][a][b][c][0]

	def save_dw(self,p):
		with open(p['data_folder']+'dw_all_'+p['experiment']+'.pkl', 'wb') as output:
			pickle.dump(self.w_end, output,protocol=pickle.HIGHEST_PROTOCOL)

	def plot_dw_all(self,p,dw):
		# create figure
		self.fig_dw_all = plt.figure()
		# loop over n_pol
		for field_i,field in enumerate(p['field']):
			# plot
			plt.plot(field_i*np.ones(len(dw[field_i,:])),dw[field_i,:],p['field_color'][field_i]+'.')
		# save figure
		self.fig_dw_all.savefig(p['data_folder']+'fig_dw_all'+'.png', dpi=250)
		plt.close(self.fig_dw_all)
	
	def plot_dw_mean(self,p,dw):
		# determine stats
		dw_mean = np.mean(dw,axis=1)
		dw_std = np.std(dw,axis=1)
		dw_sem = stats.sem(dw,axis=1)
		# create figure
		self.fig_dw_mean = plt.figure()
		# loop over n_pol
		for field_i,field in enumerate(p['field']):
			# plot
			plt.errorbar(field_i,dw_mean[field_i],yerr=dw_sem[field_i],color = p['field_color'][field_i],fmt='.')
		# save figure
		self.fig_dw_mean.savefig(p['data_folder']+'fig_dw_mean'+'.png', dpi=250)
		plt.close(self.fig_dw_mean)

class Spikes():
	"""
	detect spikes and determine where they originated
	"""
	def __init__(self,p):
		self.initialize_vectors(p)
		self.group_spikes(p)
		self.spike_start(p)
		self.save_spikes(p)
		self.plot_spike_hist_soma(self.spiket_soma,p)
		self.plot_spike_hist_dend(self.spike_dend_init,p)

	def initialize_vectors(self,p):
		self.n_pol = len(p['field'])
		# initialize lists
		self.spiket_soma = [] # soma spike times [polarity list][spikes array]
		self.spiket_dend = [] # dendrite spike times [polarity list][spikes array]
		self.sec_list = [] # keep track of sections (same dimensions as spiket_dend)
		self.seg_list = [] # keep track of segments (same dimensions as spiket_dend)
		self.cell_list_soma = []
		self.cell_list_dend = []
		self.win_list_soma = []
		# loop over polarity
		for pol in range(self.n_pol):
			self.spiket_soma.append(np.empty([1,0]))
			self.spiket_dend.append(np.empty([1,0]))
			self.sec_list.append([])
			self.seg_list.append([])
			self.cell_list_soma.append([])
			self.cell_list_dend.append([])
			self.win_list_soma.append([])

	def group_spikes(self,p):
		cell_num = -1 	# track which cell number
		for data_file in os.listdir(p['data_folder']):
			# check for proper data file format
			if 'data' in data_file:
				cell_num+=1
				# load data file
				pkl_file = open(p['data_folder']+data_file, 'rb')
				data = pickle.load(pkl_file)
				# get parameters from specific experiment
				self.p = data['params']
				self.n_act_seg = len(self.p['seg_list'])
				self.measure_spikes(data,self.p,cell_num)
				
	def measure_spikes(self,data,p,cell_num=0):
		# determine windows for spike time detection [window number][min max]
		bursts = range(p['bursts'])
		pulses = range(p['pulses'])
		burst_freq = p['burst_freq']
		pulse_freq = p['pulse_freq']
		fs = 1./p['dt']
		warmup = p['warmup']
		# for each input pulse their is a list containing the window start and stop time [window number][start,stop]
		window =  [[warmup*fs+(burst)*1000*fs/burst_freq+(pulse)*1000*fs/pulse_freq,warmup*fs+(burst)*1000*fs/burst_freq+(pulse+1)*1000*fs/pulse_freq] for burst in bursts for pulse in pulses]
		
		# detect spikes for individual neurons
		for pol in range(self.n_pol):
			
			# detect soma spikes
			self.soma_spikes = self.detect_spikes(np.array(data['soma'][pol]))['times']
			if self.soma_spikes.size!=0:
				
				# add spike times to array
				self.spiket_soma[pol] = np.append(self.spiket_soma[pol],self.soma_spikes/fs,axis=1)
				# track cell number
				for spike_i,spike in enumerate(self.soma_spikes[0]):
					# detect the window that the spike occurred in, indexed by the onset time of the window
					spike_soma_win = [win[0] for win in window if spike >= win[0] and spike < win[1] ]
					self.cell_list_soma[pol].append(cell_num)
					self.win_list_soma[pol].append(spike_soma_win)
			
			# detect dendritic spikes and track location
			cnt=-1
			sec_list_track = []
			seg_list_track = []
			for sec in range(len(data['dend'][pol])): # loop over sections
				for seg in range(len(data['dend'][pol][sec])): # loop over segemnts
					cnt+=1
					# detect spikes
					dend_spikes = self.detect_spikes(np.array(data['dend'][pol][sec][seg]))['times']
					if dend_spikes.size!=0:
						# add spike times to array
						self.spiket_dend[pol] = np.append(self.spiket_dend[pol],dend_spikes,axis=1)
						# spiket_dend_track = np.append(spiket_dend_track,dend_spikes,axis=1)
						# for each spike store the section, segment, cell number in the appropriate list
						for spike in dend_spikes[0,:]:
							self.sec_list[pol].append(p['seg_idx'][sec])
							self.seg_list[pol].append(p['seg_idx'][sec][seg])
							self.cell_list_dend[pol].append(cell_num)
							sec_list_track.append(p['seg_idx'][sec])
							seg_list_track.append(p['seg_idx'][sec][seg])

	def spike_start(self,p):
		# determine windows for spike time detection [window number][min max]
		bursts = range(p['bursts'])
		pulses = range(p['pulses'])
		burst_freq = p['burst_freq']
		pulse_freq = p['pulse_freq']
		fs = 1./p['dt']
		warmup = p['warmup']
		# for each input pulse their is a list containing the window start and stop time [window number][start,stop]
		window =  [[warmup*fs+(burst)*1000*fs/burst_freq+(pulse)*1000*fs/pulse_freq,warmup*fs+(burst)*1000*fs/burst_freq+(pulse+1)*1000*fs/pulse_freq] for burst in bursts for pulse in pulses]

		# numpy array for storing minumum spike time for each cell 
		self.spike_dend_init = []	# minimum spike time for each cell
		self.spike_dend_init_sec = []	# section where first spike occured
		self.spike_dend_init_seg = []	# segment where first spike occured
		self.spike_dend_init_cell = [] # keep track of cell number
		self.spike_dend_init_win = [] # timing of presynaptic input 
		# loop over polarities
		for pol in range(self.n_pol):
			# list all cells with a dendritic spike
			cells = list(set(self.cell_list_dend[pol]))
			# numpy array for storing minumum spike time for each cell 
			self.spike_dend_init.append([])	# minimum spike time for each cell
			self.spike_dend_init_sec.append([])	# section where first spike occured
			self.spike_dend_init_seg.append([])	# segment where first spike occured
			self.spike_dend_init_cell.append([]) # keep track of cell number
			self.spike_dend_init_win.append([]) # keep track of cell number
			# loop over cells
			for cell_i,cell in enumerate(cells):
				# print len(self.spiket_dend[pol][0,:])
				# print len(self.cell_list_dend[pol])
				# for each cell list all dendritic spike times
				spiket_dend = [spikes for spike_i,spikes in enumerate(self.spiket_dend[pol][0,:]) if self.cell_list_dend[pol][spike_i]==cell]
				# keep track of the index for in the full list of spike times
				spikei_dend = [spike_i for spike_i,spikes in enumerate(self.spiket_dend[pol][0,:]) if self.cell_list_dend[pol][spike_i]==cell]
				# loop over spike windows
				for win in window:
					# return spikes for this cell that fit the window
					spiket_dend_win = [spike for spike in spiket_dend if spike >= win[0] and spike < win[1]]
					# keep track of indeces
					spikei_dend_win = [spike for spike_i,spike in enumerate(spikei_dend) if spiket_dend[spike_i] >= win[0] and spiket_dend[spike_i] < win[1]]
					# check if spike was found
					if spiket_dend_win:
						# print min(spiket_dend_win)/fs
						# index in full list for first spike in current time window in current cell
						spike_idx = spikei_dend_win[spiket_dend_win.index(min(spiket_dend_win))]
						# store minimum spike time and keep track of section, segment, and cell number
						self.spike_dend_init[pol].append(min(spiket_dend_win)/fs)
						self.spike_dend_init_sec[pol].append(self.sec_list[pol][spike_idx])
						self.spike_dend_init_seg[pol].append(self.seg_list[pol][spike_idx])
						self.spike_dend_init_cell[pol].append(cell)
						self.spike_dend_init_win[pol].append(win[0])
			
			self.spike_dend_init[pol] = np.array([self.spike_dend_init[pol]])
			print self.spike_dend_init[pol].shape

	def save_spikes(self,p):
		with open(p['data_folder']+'spiket_soma_all_'+p['experiment']+'.pkl', 'wb') as output:
			pickle.dump(self.spiket_soma, output,protocol=pickle.HIGHEST_PROTOCOL)

		with open(p['data_folder']+'spiket_dend_all_'+p['experiment']+'.pkl', 'wb') as output:
			pickle.dump(self.spiket_dend, output,protocol=pickle.HIGHEST_PROTOCOL)

	def plot_spike_hist_soma(self,data,p):
		warmup = p['warmup']
		finish = p['tstop']
		bins = np.linspace(warmup,finish ,(finish-warmup)/p['dt'])
		self.fig_spike_hist_soma = plt.figure()
		for pol in range(self.n_pol):
			plt.hist(data[pol][0,:],bins=bins,color = p['field_color'][pol])
		plt.title('Somatic spike time histogram')
		plt.xlabel('spike onset (ms)')
		plt.ylabel('count')
		# save figure
		self.fig_spike_hist_soma.savefig(p['data_folder']+'fig_spike_hist_soma'+'.png', dpi=250)
		plt.close(self.fig_spike_hist_soma)

	def plot_spike_hist_dend(self,data,p):
		warmup = p['warmup']
		finish = p['tstop']
		bins = np.linspace(warmup,finish ,(finish-warmup)/p['dt'])
		self.fig_spike_hist_dend = plt.figure()
		for pol in range(self.n_pol):
			plt.hist(data[pol][0,:],bins=bins,color = p['field_color'][pol])
		plt.title('Dendritic spike time histogram')
		plt.xlabel('spike onset (ms)')
		plt.ylabel('count')
		self.fig_spike_hist_dend.savefig(p['data_folder']+'fig_spike_hist_dend'+'.png', dpi=250)
		plt.close(self.fig_spike_hist_dend)

	def detect_spikes(self,data,threshold=-20):
		spike_times = np.asarray(np.where(np.diff(np.sign(data-threshold))>0))
		spike_train = np.zeros([1,len(data)])
		for time in spike_times:
			spike_train[0,time] = 1 
		# detect indeces where vector crosses threshold in the positive direction
		return {'times':spike_times,'train':spike_train}
		

class Voltage():
	def __init__(self,p):
		self.n_pol = len(p['field'])
		cell_num=-1
		for data_file in os.listdir(p['data_folder']):
			# check for proper data file format
			if 'data' in data_file:
				cell_num+=1
				# load data file
				pkl_file = open(p['data_folder']+data_file, 'rb')
				data = pickle.load(pkl_file)
				self.fig_dend_trace = plt.figure()
				for pol in range(self.n_pol):
					print data['dend'][pol][0][0][200]
					plt.plot(data['t'][pol],data['dend'][pol][10][0],color = p['field_color'][pol])
					# plt.plot(data['t'][pol],np.array(data['soma'][pol]),color = p['field_color'][pol])
				self.fig_dend_trace.savefig(p['data_folder']+'fig_dend_trace'+'cellnum'+str(cell_num)+'.png', dpi=250)
				plt.close(self.fig_dend_trace)


if __name__ =="__main__":
	# Weights(param.exp_3().params)
	Spikes(param.exp_3().params)
	# Voltage(param.exp_3().params)




