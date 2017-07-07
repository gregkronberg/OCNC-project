"""
docstring
"""
# imports
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import stims
# load neuron's runtime routine - this is critical
h.load_file("stdrun.hoc")

class Cell():
	"""
	docstring
	"""
	def __init__(self):
		#=============================================================================
		# Parameters from Migliore 2005 (signal propogation in oblique dendrites)
		#=============================================================================
		Vrest   = -65
		dt      = 0.1
		gna     =  .025
		AXONM   = 5
		gkdr    = 0.01
		celsius = 35.0  
		KMULT   =  0.1
		KMULTP  = 0.03
		ghd     = 0.00005
		vhalfl_prox  = -73
		vhalfl_dist = -81
		RaAll   = 150
		RmAll   = 28000
		ka_grad = 1
		w_ampa = 4*.00018 # ampa weight (microsiemens or micro-omhos)
		w_nmda = 4*.00018 # nmda weight (microsiemens or micro-omhos)

		#==============================================================================
		# load morphology
		#==============================================================================
		h.load_file('geo5038804.hoc')   # load cell class geometry
		h.load_file('fixnseg.hoc')  # load functions required for d_lambda rule used in cell class hoc file

		#h("geom_nseg()")
		self.soma = h.soma[0]
		self.axon = h.axon[0]
		self.dend_b = h.dendrite
		self.dend_a_trunk = h.user5
		self.dend_a_tuft = h.apical_dendrite

		#==============================================================================
		#%% soma
		#==============================================================================
		# soma biophysics Hodgkin Huxley
		self.soma.insert('hh')
		self.soma.gnabar_hh = 0.12		# Sodium conductance in S/cm2
		self.soma.gkbar_hh = 0.036      # Potassium conductance in S/cm2
		self.soma.gl_hh = 0.0003        # Leak conductance in S/cm2
		self.soma.el_hh = -54.3         # Reversal potential in mV

		# passive biophysics
		self.soma.insert('pas')
		self.soma.g_pas = 5.0/RmAll		# passive conductance (S/cm^2)
		self.soma.e_pas = -65			# leak reversal potential (mV) 
		self.soma.Ra = RaAll                               

		# soma biophysics (parmameters: Migliore 2005)
		self.soma.insert('na3')
		self.soma.gbar_na3 = gna
		self.soma.insert('hd')
		self.soma.ghdbar_hd = ghd
		self.soma.vhalfl_hd = vhalfl_prox
		self.soma.insert('kdr')
		self.soma.gkdrbar_kdr = gkdr
		self.soma.insert('kap')
		self.soma.gkabar_kap = KMULTP

		# set soma as origin for distance measurements
		h.distance(sec = self.soma)

		#==============================================================================
		#%% Basal Dendrites
		#==============================================================================    
		self.syn_b_ampa = []
		self.syn_b_nmda = []
		# loop over dendrite sections
		for a in range(0,len(self.dend_b)):      
		    # passive biophysics
		    self.dend_b[a].Ra = RaAll          # axial resistance (units?)
		    self.dend_b[a].insert('pas')       # insert passive mechanism
		    self.dend_b[a].g_pas = 5.0/RmAll   # passive conductance (S/cm^2)
		    self.dend_b[a].e_pas = -65         # leak reversal potential (mV)

		    # active biophysics (parameters: Migliore 2005)
		    # h current
		    self.dend_b[a].insert('hd')
		    self.dend_b[a].ghdbar_hd = ghd
		    # sodium
		    self.dend_b[a].insert('na3')
		    self.dend_b[a].gbar_na3 = gna
		    # delayed rectifier potassium
		    self.dend_b[a].insert('kdr')
		    self.dend_b[a].gkdrbar_kdr = gkdr
		    # A-type potassium (proximal)
		    self.dend_b[a].insert('kap')
		    self.dend_b[a].gkabar_kap = 0
		    # A-type potassium (distal)
		    self.dend_b[a].insert('kad')
		    self.dend_b[a].gkabar_kad = 0
		    
		    # channels that vary with distance from soma (parameters: Migliore 2005)
		    # loop over segments     
		    cnt=-1
		    self.syn_b_ampa.append([])
		    self.syn_b_nmda.append([])
		    for seg in self.dend_b[a]:
		        cnt+=1
		        # distance from soma
		        seg_dist = h.distance(seg.x,sec=self.dend_b[a])
		        # h current
		        seg.ghdbar_hd = ghd*(1+3*seg_dist/100)
		        # A-type potassium
		        if seg_dist > 100:	# distal
		            seg.vhalfl_hd = vhalfl_dist
		            seg.gkabar_kad = KMULT*(1+ka_grad*seg_dist/100)
		        else:	# proximal
		            seg.vhalfl_hd = vhalfl_prox
		            seg.gkabar_kap = KMULT*(1+ka_grad*seg_dist/100)

			    # AMPA synapses (params: Kim 2015)
			    self.syn_b_ampa[a].append(h.Exp2Syn(self.dend_b[a](0.5)))
			    self.syn_b_ampa[a][cnt].tau1 = 0.2
			    self.syn_b_ampa[a][cnt].tau2 = 2
			    self.syn_b_ampa[a][cnt].i = 0.18

			    # NMDA synapses (params: Kim 2015)
			    self.syn_b_nmda[a].append(h.Exp2SynNMDA(self.dend_b[a](0.5)))
			    self.syn_b_nmda[a][cnt].tau1 = 1
			    self.syn_b_nmda[a][cnt].tau2 = 50

		#==============================================================================
		#%% Apical Dendrites
		#==============================================================================  
		self.syn_a_trunk_ampa = []
		self.syn_a_trunk_nmda = []
		# loop over dendrite sections
		for a in range(0,len(self.dend_a_trunk)):      
		    # passive biophysics (params: Migliore 2005)
		    self.dend_a_trunk[a].Ra = RaAll          # axial resistance (units?)
		    self.dend_a_trunk[a].insert('pas')       # insert passive mechanism
		    self.dend_a_trunk[a].g_pas = 5.0/RmAll   # passive conductance (S/cm^2)
		    self.dend_a_trunk[a].e_pas = -65         # leak reversal potential (mV)

		    # active biophysics (params Migliore 2005)
		    # h current
		    self.dend_a_trunk[a].insert('hd')
		    self.dend_a_trunk[a].ghdbar_hd = ghd
		    #sodium current
		    self.dend_a_trunk[a].insert('na3')
		    self.dend_a_trunk[a].gbar_na3 = gna
		    # delayed rectifier potassium
		    self.dend_a_trunk[a].insert('kdr')
		    self.dend_a_trunk[a].gkdrbar_kdr = gkdr
		    # A-type potassium proximal
		    self.dend_a_trunk[a].insert('kap')
		    self.dend_a_trunk[a].gkabar_kap = 0    # set to zero, then change if the synapse meets the "proximal" criteria
		    # A-type potassium distal
		    self.dend_a_trunk[a].insert('kad')
		    self.dend_a_trunk[a].gkabar_kad = 0    # set to zero, then change if the synapse meets the "distal" criteria
		        
		    # channels that vary with distance from soma (params: Migliore 2005)
		    # loop over segments    
		    cnt=-1
		    self.syn_a_trunk_ampa.append([])
		    self.syn_a_trunk_nmda.append([])
		    for seg in self.dend_a_trunk[a]:
		        cnt+=1
		        # distance from soma
		        seg_dist = h.distance(seg.x,sec=self.dend_a_trunk[a])
				# h current
		        seg.ghdbar_hd = ghd*(1+3*seg_dist/100)
		        # A-type potassium 
		        if seg_dist > 100:	# distal
		            seg.vhalfl_hd = vhalfl_dist
		            seg.gkabar_kad = KMULT*(1+ka_grad*seg_dist/100)
		        else: # proximal
		            seg.vhalfl_hd = vhalfl_dist
		            seg.gkabar_kap = KMULT*(1+ka_grad*seg_dist/100)

				# AMPA synapses (params: Kim 2015)
		        self.syn_a_trunk_ampa[a].append(h.Exp2Syn(self.dend_a_trunk[a](0.5)))
		        self.syn_a_trunk_ampa[a][cnt].tau1 = 0.2
		        self.syn_a_trunk_ampa[a][cnt].tau2 = 2
		        self.syn_a_trunk_ampa[a][cnt].i = 0.18
		    
		        # NMDA synapses (params: Kim 2015)
		        self.syn_a_trunk_nmda[a].append(h.Exp2SynNMDA(self.dend_a_trunk[a](0.5)))
		        self.syn_a_trunk_nmda[a][cnt].tau1 = 1
		        self.syn_a_trunk_nmda[a][cnt].tau2 = 50

		#==============================================================================
		#%% Apical tuft dendrites
		#==============================================================================  
		self.syn_a_tuft_ampa = []
		self.syn_a_tuft_nmda = []
		# loop over dendrite sections
		for a in range(0,len(self.dend_a_tuft)):      
		    # passive biophysics
		    self.dend_a_tuft[a].Ra = RaAll          # axial resistance (units?)
		    self.dend_a_tuft[a].insert('pas')       # insert passive mechanism
		    self.dend_a_tuft[a].g_pas = 5.0/RmAll   # passive conductance (S/cm^2)
		    self.dend_a_tuft[a].e_pas = -65         # leak reversal potential (mV)

		    # active biophysics (params: Migliore 2005)
		    # h current
		    self.dend_a_tuft[a].insert('hd')
		    self.dend_a_tuft[a].ghdbar_hd = ghd
		    #sodium current
		    self.dend_a_tuft[a].insert('na3')
		    self.dend_a_tuft[a].gbar_na3 = gna
		    # delayed rectifier potassium
		    self.dend_a_tuft[a].insert('kdr')
		    self.dend_a_tuft[a].gkdrbar_kdr = gkdr
		    # A-type potassium proximal
		    self.dend_a_tuft[a].insert('kap')
		    self.dend_a_tuft[a].gkabar_kap = 0    # set to zero, then change if the synapse meets the "proximal" criteria
		    # A-type potassium distal
		    self.dend_a_tuft[a].insert('kad')
		    self.dend_a_tuft[a].gkabar_kad = 0    # set to zero, then change if the synapse meets the "distal" criteria
		        
		    # channels that vary with distance from soma
		    # loop over segments    
		    cnt=-1
		    self.syn_a_tuft_ampa.append([])
		    self.syn_a_tuft_nmda.append([])
		    for seg in self.dend_a_tuft[a]:
		        cnt+=1
		        # distance from soma
		        seg_dist = h.distance(seg.x,sec=self.dend_a_tuft[a])
		        # h current
		        seg.ghdbar_hd = ghd*(1+3*seg_dist/100)
		        # A-type potassium current
		        if seg_dist > 100:	# distal
		            seg.vhalfl_hd = vhalfl_dist
		            seg.gkabar_kad = KMULT*(1+ka_grad*seg_dist/100)
		        else:	# proximal
		            seg.vhalfl_hd = vhalfl_dist
		            seg.gkabar_kap = KMULT*(1+ka_grad*seg_dist/100)

		        #==============================================================================
		        #%%  synapses
		        #==============================================================================
		        # AMPA synapses (Kim et al. 2015 Elife)
		        self.syn_a_tuft_ampa[a].append(h.Exp2Syn(self.dend_a_tuft[a](seg.x)))
		        self.syn_a_tuft_ampa[a][cnt].tau1 = 0.2
		        self.syn_a_tuft_ampa[a][cnt].tau2 = 2
		        self.syn_a_tuft_ampa[a][cnt].i = 0.18
		    
		        # NMDA synapses (Kim et al. 2015 Elife)
		        self.syn_a_tuft_nmda[a].append(h.Exp2SynNMDA(self.dend_a_tuft[a](seg.x)))
		        self.syn_a_tuft_nmda[a][cnt].tau1 = 1
		        self.syn_a_tuft_nmda[a][cnt].tau2 = 50

# activate a specified subset of synapses with NetStim/NetCon
class Syn_act:
	"""
	docstring
	"""
	def __init__(self,subtree_ampa,subtree_nmda,sec_idx,seg_idx,w_ampa,w_nmda):
		
		# choose active segments
		self.subtree_ampa  = subtree_ampa # neuronal subtree_ampa (soma, apical, basal, etc.)
		self.subtree_nmda = subtree_nmda
		self.sec_idx = sec_idx	# list of section indeces with active synapses
		self.seg_idx = seg_idx	# nested list of active segment indices within each section, e.g [section][segment]
		self.w_ampa = w_ampa	# maximum ampa conductance (microsiemens or micro-omhos)
		self.w_nmda = w_nmda	# maximum ampa conductance (microsiemens or micro-omhos)
		self.choose_syn_ampa = [self.subtree_ampa[a] for a in self.sec_idx]
		self.choose_syn_nmda = [self.subtree_nmda[a] for a in self.sec_idx]
		self.active = []

		# default stim objects
		self.warm_up = 30   # warm up time (ms)
		self.pulse_freq = 100
		self.burst_freq = 5
		self.bursts = 1
		self.pulses = 1
		self.stim  = []
		self.stim.append(h.NetStim())
		self.stim[0].start = self.warm_up
		self.stim[0].interval = 1000/self.pulse_freq
		self.stim[0].noise  = 0 
		self.stim[0].number = self.pulses
		
		self.nc_ampa = []
		self.nc_nmda = []
		# activate synapses
		for a in range(len(self.choose_syn_ampa)):
		    self.active.append([])
		    self.nc_ampa.append([])
		    self.nc_nmda.append([])
		    for b in range(len(self.choose_syn_ampa[a])):
		        if b in seg_idx[a]: # choose segments within each section
		            self.active[a].append(1)
	                self.nc_ampa[a].append( h.NetCon(self.stim[0],self.choose_syn_ampa[a][b],0,0,w_ampa)) 
            	    self.nc_nmda[a].append( h.NetCon(self.stim[0],self.choose_syn_nmda[a][b],0,0,w_nmda))
            	else:
		            self.active[a].append(0)
		            self.nc_ampa[a].append( h.NetCon(self.stim[0],self.choose_syn_nmda[a][b],0,0,0))
		            self.nc_nmda[a].append( h.NetCon(self.stim[0],self.choose_syn_nmda[a][b],0,0,0))

if __name__ == "__main__":
    pass
