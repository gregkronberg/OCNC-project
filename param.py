# parameters 

from neuron import h
import numpy as np
import cell 
import stims
# active 1 synapse on random branches (choose number of branches)

class Active_syn:
	def __init__(self):
		self.warm_up = 50   # warm up time (ms)
		self.pulse_freq = 100
		self.burst_freq = 5
		self.bursts = 1
		self.pulses = 1
		self.total_time = warm_up + a*1000/burst_freq + pulse_freq*pulses+ 20
		self.len_t = total_time/dt+2
	def active(self,subtree,sec_idx,seg_idx,w_ampa,w_nmda)
		#======================================
		# synaptic stimulation
		#======================================
		self.choose_syn_ampa = [subtree[a] for a in sec_idx]
		self.choose_syn_nmda = [subtree[a] for a in sec_idx]
		choose_seg = range(0,segs)
		dend_vec = [] # setup recording vector
		active = []
		for a in range(len(self.choose_syn_ampa)):
		    dend_vec.append(h.Vector())
		    dend_vec[a].record(choose_sec[a](0)._ref_v)
		    active.append([])
		    for b in range(len(choose_syn_ampa[a])):
		        if b in choose_seg: # choose segments within each section
		            active[a].append(1)
		        else:
		            active[a].append(0)

					             # #####################################
			## Theta Burst
			######################################
			# set up NetStim, NetCon objects
			# each burst is a NetStim object
			stim = []
			nc_ampa = []
			nc_nmda = []
			for a in range(0,bursts):
			    stim.append(h.NetStim())
			    stim[a].start = warm_up + a*1000/burst_freq
			    stim[a].interval = 1000/pulse_freq
			    stim[a].noise  = 0 
			    stim[a].number = pulses
			    
			    nc_ampa.append([])
			    nc_nmda.append([])
			    for b in range(len(choose_syn_ampa)):
			        nc_ampa[a].append([])
			        nc_nmda[a].append([])
			        cnt=-1
			        for c in range(len(choose_syn_ampa[b])):
			            cnt+=1
			            delay = (len(choose_syn_ampa[b]) - c)*5
			            if active[b][c]==1:
			                nc_ampa[a][b].append( h.NetCon(stim[a],choose_syn_ampa[b][c],0,delay,w_ampa)) # nested netcon objects ([stim number][section][segment])
			                nc_nmda[a][b].append( h.NetCon(stim[a],choose_syn_nmda[b][c],0,delay,w_nmda)) # nested netcon objects ([stim number][section][segment])
			            else:
			                nc_ampa[a][b].append( h.NetCon(stim[a],choose_syn_nmda[b][c],0,delay,0)) # nested netcon objects ([stim number][section][segment])
			                nc_nmda[a][b].append( h.NetCon(stim[a],choose_syn_nmda[b][c],0,delay,0)) # nested netcon objects ([stim number][section][segment])
			#            nc[a][b][cnt].delay = (len(choose_syn[b]) - c)*5
class Choose_sec:
	def __init__(self,sec_tot,):
		self.sec_n = 1
		self.sec_i = np.random.choice(secs,nsec) 