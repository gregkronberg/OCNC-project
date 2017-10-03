"""
analysis
"""
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
		self.p = p
		self.n_pol = len(self.p['field'])
		self.n_act_seg = len(self.p['seg_list'])
		self.group_dw(p)
		self.save_dw(p)
		self.plot_dw_all(p)
		self.plot_dw_mean(p)

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
				# measure weight changes for individual neuron
				self.measure_dw(data)
				# add individual neuron to the group
				self.w_end_all = np.append(self.w_end_all,self.w_end,axis=1)
				self.w_start_all = np.append(self.w_start_all,self.w_start,axis=1)

	def measure_dw(self,data):
		# set up arrays to record final and initial weights at each active synapse
		self.w_end = np.empty([n_pol,n_act_seg]) # polarity x segments
		self.w_start = np.empty([n_pol,n_act_seg]) # polarity x segments

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
		plt.close(fig_dw_all)
	
	def plot_dw_mean(self,p,dw):
		# determine stats
		dw_mean = np.mean(dw,axis=1)
		dw_std = np.std(dw,axis=1)
		dw_sem = np.sem(dw,axis=1)
		# create figure
		self.fig_dw_mean = plt.figure()
		# loop over n_pol
		for field_i,field in enumerate(p['field']):
			# plot
			plt.errorbar(field_i,dw_mean[field_i],yerr=dw_sem[field_i],color = p['field_color'][field_i],fmt='.')
		# save figure
		self.fig_dw_all.savefig(p['data_folder']+'fig_dw_mean'+'.png', dpi=250)
		plt.close(fig_dw_mean)

class Spikes():
	"""
	detect spikes and determine where they originated
	"""
	def __init__(self,p):
		self.p = p
		self.n_pol = len(self.p['field'])
		self.n_act_seg = len(self.p['seg_list'])
		self.group_spikes(p)


	def group_spikes(self,p):
		cell_num = -1 
		for data_file in os.listdir(p['data_folder']):
			# check for proper data file format
			if 'data' in data_file:
				cell_num+=1
				# load data file
				pkl_file = open(p['data_folder']+data_file, 'rb')
				data = pickle.load(pkl_file)

				self.measure_spikes(data, cell_num)
				

	def measure_spikes(self,data,cell_num=0):
		# detect spikes for individual neurons
		# initialize lists
		self.spiket_soma = [] # soma spike times [polarity list][spikes array]
		self.spiket_dend = [] # dendrite spike times [polarity list][spikes array]
		self.sec_list = [] # keep track of sections (same dimensions as spiket_dend)
		self.seg_list = [] # keep track of segments (same dimensions as spiket_dend)
		self.cell_list_soma = []
		self.cell_list_dend = []
		# loop over polarity
		for pol in range(self.n_pol):
			self.spiket_soma.append(np.empty([]))
			self.spiket_dend.append(np.empty([]))
			self.sec_list.append([])
			self.seg_list.append([])
			self.cell_list_dend.append([])
			self.cell_list_dend.append([])

			# detect soma spikes
			soma_spikes = self.detect_spikes(np.array(data['soma'][pol]))
			# add spike times to array
			self.spiket_soma[pol] = np.append(self.spiket_soma[pol],soma_spikes,axis=0)
			# track cell number
			for spike in soma_spikes:
				self.cell_list_soma[pol].append(cell_num)
			# loop over dendritic segments
			cnt=-1
			for sec in range(len(data['dend'][pol])): # loop over sections
				for seg in range(len(data['dend'][pol][sec])): # loop over segemnts
				cnt+=1
				# detect spikes
				dend_spikes = self.spiket_dend[cnt][pol,:] = self.detect_spikes(np.array(data['dend'][pol][sec][seg])) 
				# add spike times to array
				self.spiket_dend[pol] = np.append(self.spiket_dend[pol],dend_spikes,axis=0)
				# for each spike store the section, segment, cell number in the appropriate list
				for spike in dend_spikes:
					self.sec_list[pol].append(p['seg_idx'][sec])
					self.seg_list[pol].append(p['seg_idx'][sec][seg])
					self.cell_list_dend[pol].append(cell_num)

	def plot_spike_hist_soma(data,p):
		warmup = p['warmup']./dt
		finish = p['tstop']./dt
		bins = np.linspace(warmup,finish ,(finish-warmup)/2)
		fig_spike_hist_soma = plt.figure()
		for pol in range(self.n_pol):
			plt.hist(self.spiket_soma[pol],bins=bins,color = p['field_color'][pol])
		# save figure
		self.fig_spike_hist_soma.savefig(p['data_folder']+'fig_spike_hist_soma'+'.png', dpi=250)
		plt.close(fig_spike_hist_soma)

	def plot_spike_hist_soma(data,p):
		warmup = p['warmup']./dt
		finish = p['tstop']./dt
		bins = np.linspace(warmup,finish ,(finish-warmup)/2)
		fig_spike_hist_dend = plt.figure()
		for pol in range(self.n_pol):
			plt.hist(self.spiket_dend[pol],bins=bins,color = p['field_color'][pol])
		self.fig_spike_hist_dend.savefig(p['data_folder']+'fig_spike_hist_dend'+'.png', dpi=250)
		plt.close(fig_spike_hist_dend)

	def detect_spikes(self,data,threshold=-20):
		# detect indeces where vector crosses threshold in the positive direction
		return np.where(np.diff(np.sign(data-threshold))>0)
		# loop over section


def exp_2_spike_analysis():
	data_all = {
		'spike_t_soma':[],
		'spike_t_dend':[],
		'n_input':[],
		'field':[],
		'weight':[],
	}
	dpath = 'Data/'
	data_list  = os.listdir(dpath)
	for data_name in data_list:
		data_pkl = open(dpath+data_name, 'rb')
		data = pickle.load(data_pkl)
		n_exp = len(data['soma'])
		for a in range(n_exp):
			data_soma = data['soma'][a]
			soma_cross = np.where(np.diff(np.sign(data_soma)))[0]
			if len(soma_cross)>0:
				soma_cross_min = soma_cross[0]
			else:
				soma_cross_min=0
			print soma_cross_min
			n_sec = len(data['dend'][a])
			dend_cross = []
			dend_cross_sec=[]
			dend_cross_sec_min=[]
			for b in range(n_sec):
				data_dend = data['dend'][a][b]
				dend_cross_sec.append(np.where(np.diff(np.sign(data_dend)))[0])
				if len(dend_cross_sec[b])>0:
					dend_cross_sec_min.append(dend_cross_sec[b][0])
				else:
					dend_cross_sec_min.append(0)
			dend_cross_sec_min_np = np.array(dend_cross_sec_min)
			if np.count_nonzero(dend_cross_sec_min_np)>0:
				dend_cross_min = np.min(dend_cross_sec_min_np[np.nonzero(dend_cross_sec_min_np)])
			else:
				dend_cross_min = 0
			data_all['spike_t_soma'].append(soma_cross_min)
			data_all['spike_t_dend'].append(dend_cross_min)
			data_all['n_input'].append(n_sec)
			data_all['field'].append(data['field'][a])
			data_all['weight'].append(data['params']['w_ampa'])

	data_all['spike_t_soma'] = np.array(data_all['spike_t_soma'])
	data_all['spike_t_dend']=np.array(data_all['spike_t_dend'])
	data_all['n_input']=np.array(data_all['n_input'])
	data_all['field']=np.array(data_all['field'])
	data_all['weight']=np.array(data_all['weight'])
	data_all['spike_source']=[]
	for a in range(len(data_all['spike_t_soma'])):
		if data_all['spike_t_soma'][a]!=0:
			if data_all['spike_t_soma'][a]<=data_all['spike_t_dend'][a]:
				data_all['spike_source'].append(1)
			else:
				data_all['spike_source'].append(0)
		else:
			data_all['spike_source'].append(-1)
	data_all['spike_source']=np.array(data_all['spike_source'])
	
	with open('data_all_exp_2'+'.pkl', 'wb') as output:
	# 
		pickle.dump(data_all, output,protocol=pickle.HIGHEST_PROTOCOL)

def plot_spikes(data_file):
	"""
	creates a single figure with subplots for each section to be recorded from
	"""
	pkl_file = open(data_file, 'rb')
	
	data_all = pickle.load(pkl_file)
	
	cathodal_dend_timing = data_all['spike_t_dend'][np.where(data_all['field']== -20)] 
	# print np.where(data_all['field']== -20)
	# print data_all['spike_t_dend']
	# print cathodal_dend_timing
	control_dend_timing = data_all['spike_t_dend'][np.where(data_all['field']== 0)] 
	anodal_dend_timing = data_all['spike_t_dend'][np.where(data_all['field']== 20)] 
	

	cathodal_soma_timing = data_all['spike_t_soma'][np.where(data_all['field']== -20)]
	cathodal_soma_timing = .025*(cathodal_soma_timing[np.where(cathodal_soma_timing!= 0)] -1200.)
	cathodal_soma_timing_mean = np.mean(cathodal_soma_timing)
	cathodal_soma_timing_std = np.std(cathodal_soma_timing)/np.sqrt(cathodal_soma_timing.size)  
	control_soma_timing = data_all['spike_t_soma'][np.where(data_all['field']== 0)]
	control_soma_timing = .025*(control_soma_timing[np.where(control_soma_timing!= 0)] -1200.)
	control_soma_timing_mean = np.mean(control_soma_timing)  
	control_soma_timing_std = np.std(control_soma_timing)/np.sqrt(control_soma_timing.size)  
	anodal_soma_timing = data_all['spike_t_soma'][np.where(data_all['field']== 20)]
	anodal_soma_timing = .025*(anodal_soma_timing[np.where(anodal_soma_timing!= 0)]-1200.)
	print anodal_soma_timing
	anodal_soma_timing_mean = np.mean(anodal_soma_timing)
	anodal_soma_timing_std = np.std(anodal_soma_timing)/np.sqrt(anodal_soma_timing.size)     

	cathodal_source = data_all['spike_source'][np.where(data_all['field']== -20)] 
	control_source = data_all['spike_source'][np.where(data_all['field']== 0)] 
	anodal_source = data_all['spike_source'][np.where(data_all['field']== 20)]

	total = control_source.size
	cathodal_count = cathodal_source[np.where(cathodal_source!= -1)] 
	control_count = control_source[np.where(control_source!= -1)]
	anodal_count = anodal_source[np.where(anodal_source!= -1)] 
	counts = np.array([cathodal_count.size*1./total,control_count.size*1./total,anodal_count.size*1./total])
	ind = np.arange(3)
	width = .5
	# fig.hist(cathodal_dend_timing)
	fig = plt.figure(figsize=(10, 10))
	n=1
	gs = gridspec.GridSpec(n, n, wspace=0.10, hspace=0.05, left=0.1, right=0.95, bottom=0.1, top=0.95)
	ax = {}
	rows = np.arange(0, n, 1, dtype=int)
	cols = np.arange(0, n, 1, dtype=int)
	for k, (i, j) in enumerate(it.product(rows, cols)):
		# axh = "spike_prob"#"section-{:03d}".format(data['params']['sec_idx'][k])
		axh = "spike_time"
		ax[axh] = fig.add_subplot(gs[i:i+1, j:j+1])
		# ax[axh].text(0.05, 0.90, axh, transform=ax[axh].transAxes)
		plot_color = 'b'
		ax[axh].bar(0,cathodal_soma_timing_mean,width,color=plot_color,yerr=cathodal_soma_timing_std)
		# ax[axh].bar(0,counts[0],width,color=plot_color)
		# ax[axh].hist(cathodal_source,bins=4)
		# print cathodal_source
		# plt.gca().set_xlim([1200,1600])
		# plt.gca().set_ylim([0,50])
		plot_color = 'k'
		ax[axh].bar(1,control_soma_timing_mean,width,color=plot_color,yerr=control_soma_timing_std)
		# ax[axh].bar(1,counts[1],width,color=plot_color)
		# ax[axh].hist(anodal_source,bins=4)
		# plt.gca().set_xlim([1200,1600])
		plt.gca().set_ylim([4.6,4.9])
		plot_color = 'r'
		ax[axh].bar(2,anodal_soma_timing_mean,width,color=plot_color,yerr=anodal_soma_timing_std)
		# ax[axh].bar(2,counts[2],width,color=plot_color)
		# ax[axh].hist(control_source,bins=4)
		ax[axh].set_xticks([0,1,2])
		ax[axh].set_xticklabels(['+DC','control','-DC'],size=35)
		# ax[axh].set_yticks([0,.1,.2,.3,.4,.5])
		# ax[axh].set_yticklabels([0,.1,.2,.3,.4,.5],size=20)
		ax[axh].set_ylabel('spike timing (ms)',size=35)
	fig.savefig('data_spike_exp_2'+'.png', dpi=250)
	plt.close(fig)

if __name__ =="__main__":
	Weights(param.exp_3().params)
	# exp_2_spike_analysis()
	# plot_spikes('data_all_exp_2'+'.pkl')



