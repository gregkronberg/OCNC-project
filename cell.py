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

class Cell_Migliore_2005():
	"""
	pyramidal neuron based on Migliore et al. 2005
	"""
	def __init__(self):
		#=============================================================================
		# Parameters from Migliore 2005 (signal propogation in oblique dendrites)
		#=============================================================================
		Vrest   = -65
		dt      = 0.1
		gna     =  .025
		ena = 55
		AXONM   = 5
		gkdr    = 0.01
		ek = -90
		h.celsius = 35.0  
		KMULT   =  0.03
		KMULTP  = 0.03
		ghd     = 0.00005
		ehd = -30
		vhalfl_prox  = -73
		vhalfl_dist = -81
		RaAll   = 150
		RaAx = 50
		RmAll   = 28000
		Cm =  1
		ka_grad = 1
		w_ampa = 4*.00018 # ampa weight (microsiemens or micro-omhos)
		w_nmda = 4*.00018 # nmda weight (microsiemens or micro-omhos)

		#==============================================================================
		# load morphology
		#==============================================================================
		h.load_file('geo5038804.hoc')   # load cell class geometry
		h.load_file('fixnseg.hoc')  # set segments based on dlambda rule

		self.soma = h.soma[0]
		self.axon = h.axon[0]
		self.dend_b = h.dendrite
		self.dend_a_trunk = h.user5
		self.dend_a_tuft = h.apical_dendrite

		# =============================================================================
		#%% axon
		# =============================================================================
		# passive biophysics
		self.axon.insert('pas')
		self.axon.g_pas = 1/RmAll		# passive conductance (S/cm^2)
		self.axon.e_pas = Vrest			# leak reversal potential (mV)
		self.axon.cm = Cm 				# specific capacitance (uf/cm^2)
		self.axon.Ra = RaAx 			# axial resistance (ohm cm) 

		# axon active biophysics (parmameters: Migliore 2005)
		self.axon.insert('nax')
		self.axon.gbar_nax = gna*AXONM
		self.axon.insert('kdr')
		self.axon.gkdrbar_kdr = gkdr
		self.axon.insert('kap')
		self.axon.gkabar_kap = KMULTP
		self.axon.ena = ena
		self.axon.ek = ek

		#==============================================================================
		#%% soma
		#==============================================================================
		# passive biophysics
		self.soma.insert('pas')
		self.soma.g_pas = 1/RmAll		# passive conductance (S/cm^2)
		self.soma.e_pas = Vrest			# leak reversal potential (mV) 
		self.soma.Ra = RaAll     
		self.soma.cm = Cm                           

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
		self.soma.ena = ena
		self.soma.ek = ek
		# self.soma.ehd_hd = ehd

		# set soma as origin for distance measurements
		h.distance(sec = self.soma)

		#==============================================================================
		#%% Basal Dendrites
		#==============================================================================    
		self.syn_b_ampa = []
		self.syn_b_nmda = []
		self.syn_b_clopath = []
		# loop over dendrite sections
		for a in range(0,len(self.dend_b)):      
		    # passive biophysics
		    self.dend_b[a].Ra = RaAll          # axial resistance (units?)
		    self.dend_b[a].insert('pas')       # insert passive mechanism
		    self.dend_b[a].g_pas = 1/RmAll   # passive conductance (S/cm^2)
		    self.dend_b[a].e_pas = Vrest        # leak reversal potential (mV)
		    self.dend_b[a].cm = Cm 				# specific capacitance (uf/cm^2)

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

		    self.dend_b[a].ena = ena
		    self.dend_b[a].ek = ek
		    # self.dend_b[a].ehd = ehd
		    
		    # channels that vary with distance from soma (parameters: Migliore 2005)
		    # loop over segments     
		    cnt=-1
		    self.syn_b_ampa.append([])
		    self.syn_b_nmda.append([])
		    self.syn_b_clopath.append([])
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
		            seg.gkabar_kap = KMULTP*(1+ka_grad*seg_dist/100)

			    # AMPA synapses (params: Kim 2015)
			    self.syn_b_ampa[a].append(h.Exp2Syn(self.dend_b[a](seg.x)))
			    self.syn_b_ampa[a][cnt].tau1 = 0.2
			    self.syn_b_ampa[a][cnt].tau2 = 2
			    self.syn_b_ampa[a][cnt].i = 0.18

			    # NMDA synapses (params: Kim 2015)
			    self.syn_b_nmda[a].append(h.Exp2SynNMDA(self.dend_b[a](seg.x)))
			    self.syn_b_nmda[a][cnt].tau1 = 1
			    self.syn_b_nmda[a][cnt].tau2 = 50

			     # Clopath learning rule
		        self.syn_b_clopath[a].append(h.STDPSynCCNon(self.dend_b[a](seg.x)))

		#==============================================================================
		#%% Apical Dendrites
		#==============================================================================  
		self.syn_a_trunk_ampa = []
		self.syn_a_trunk_nmda = []
		self.syn_a_trunk_clopath = []
		# loop over dendrite sections
		for a in range(0,len(self.dend_a_trunk)):      
		    # passive biophysics (params: Migliore 2005)
		    self.dend_a_trunk[a].Ra = RaAll          # axial resistance (units?)
		    self.dend_a_trunk[a].insert('pas')       # insert passive mechanism
		    self.dend_a_trunk[a].g_pas = 1/RmAll   # passive conductance (S/cm^2)
		    self.dend_a_trunk[a].e_pas = Vrest         # leak reversal potential (mV)
		    self.dend_a_trunk[a].cm = Cm 			# specific capacitance (uf/cm^2)

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
		    self.dend_a_trunk[a].gkabar_kad = 0    # set to zero, then change if the segment meets the "distal" criteria
		    
		    self.dend_a_trunk[a].ena = ena
		    self.dend_a_trunk[a].ek = ek
		    # self.dend_a_trunk[a].ehd = ehd
		    # channels that vary with distance from soma (params: Migliore 2005)
		    # loop over segments    
		    cnt=-1
		    self.syn_a_trunk_ampa.append([])
		    self.syn_a_trunk_nmda.append([])
		    self.syn_a_trunk_clopath.append([])
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
		            seg.gkabar_kap = KMULTP*(1+ka_grad*seg_dist/100)

				# AMPA synapses (params: Kim 2015)
		        self.syn_a_trunk_ampa[a].append(h.Exp2Syn(self.dend_a_trunk[a](seg.x)))
		        self.syn_a_trunk_ampa[a][cnt].tau1 = 0.2
		        self.syn_a_trunk_ampa[a][cnt].tau2 = 2
		        self.syn_a_trunk_ampa[a][cnt].i = 0.18
		    
		        # NMDA synapses (params: Kim 2015)
		        self.syn_a_trunk_nmda[a].append(h.Exp2SynNMDA(self.dend_a_trunk[a](seg.x)))
		        self.syn_a_trunk_nmda[a][cnt].tau1 = 1
		        self.syn_a_trunk_nmda[a][cnt].tau2 = 50

		        # Clopath learning rule
		        self.syn_a_trunk_clopath[a].append(h.STDPSynCCNon(self.dend_a_trunk[a](seg.x)))

		#==============================================================================
		#%% Apical tuft dendrites
		#==============================================================================  
		self.syn_a_tuft_ampa = []
		self.syn_a_tuft_nmda = []

		self.syn_a_tuft_clopath = []
		# loop over dendrite sections
		for a in range(0,len(self.dend_a_tuft)):      
		    # passive biophysics
		    self.dend_a_tuft[a].Ra = RaAll          # axial resistance (units?)
		    self.dend_a_tuft[a].insert('pas')       # insert passive mechanism
		    self.dend_a_tuft[a].g_pas = 1/RmAll   # passive conductance (S/cm^2)
		    self.dend_a_tuft[a].e_pas = Vrest         # leak reversal potential (mV)
		    self.dend_a_tuft[a].cm = Cm 			# specific capacitance (uf/cm^2)

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
		     
		    self.dend_a_tuft[a].ena = ena
		    self.dend_a_tuft[a].ek = ek
		    # self.dend_a_tuft[a].ehd = ehd   
		    
		    # loop over segments    
		    cnt=-1
		    self.syn_a_tuft_ampa.append([])
		    self.syn_a_tuft_nmda.append([])
		    self.syn_a_tuft_clopath.append([])
		    for seg in self.dend_a_tuft[a]:
		        cnt+=1
		        # channels that vary with distance from soma
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
		            seg.gkabar_kap = KMULTP*(1+ka_grad*seg_dist/100)

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

		        # Clopath learning rule
		        self.syn_a_tuft_clopath[a].append(h.STDPSynCCNon(self.dend_a_tuft[a](seg.x)))

# activate a specified subset of synapses with NetStim/NetCon
class Syn_act:
	"""
	Activate a specific subset of synpases with NetCon objects

	Arguments:

	netstim - NetStim object to be delivered to NetCon object

	subtree - list containing synapse point process objects of  given type on a given tree, e.g. ampa synapses on basal dendrites.  The list should be organized subtree_type[section][segment] 

	sec_idx - list of section indices to be activated within the specified subtree

	seg_idx  - list of segments to be activated within each section.  Should be organized as [sections][segments]

	w - weights for NetCon objects
	"""
	def __init__(self,netstim,subtree_ampa,subtree_nmda,subtree_clopath,sec_idx,seg_idx,w_ampa,w_nmda):
		
		# choose active segments
		# neuronal subtree_ampa (soma, apical, basal, etc.)
		self.subtree_ampa  = subtree_ampa 
		self.subtree_nmda = subtree_nmda
		self.subtree_clopath = subtree_clopath
		self.sec_idx = sec_idx	# list of section indeces with active synapses
		self.seg_idx = seg_idx	# nested list of active segment indices within each section, e.g [section][segment]
		self.w_ampa = w_ampa	# maximum ampa conductance (microsiemens or micro-omhos)
		self.w_nmda = w_nmda	# maximum ampa conductance (microsiemens or micro-omhos)
		# list sections with active synapses
		self.choose_syn_ampa = [self.subtree_ampa[a] for a in self.sec_idx]
		self.choose_syn_nmda = [self.subtree_nmda[a] for a in self.sec_idx]
		self.choose_syn_clopath = [self.subtree_clopath[a] for a in self.sec_idx]
		self.active = []

		# store NetCon objects
		self.nc_ampa = []
		self.nc_nmda = []
		self.nc_clopath = []
		# loop over active sections
		for a in range(len(self.choose_syn_ampa)):
			self.active.append([])
			self.nc_ampa.append([])
			self.nc_nmda.append([])
			self.nc_clopath.append([])
			# loop over segments
			for b in range(len(self.choose_syn_ampa[a])):
				# if segment is not in active segment list
				if b not in seg_idx[a]:
					# do not create NetCon object
					self.active[a].append(0)
					self.nc_ampa[a].append(0)
					self.nc_nmda[a].append(0)
					self.nc_clopath[a].append(0)

				# if segment is in active list
				if b in seg_idx[a]: 
					# create NetCon object
					self.active[a].append(1)
					self.nc_ampa[a].append( h.NetCon(netstim,self.choose_syn_ampa[a][b],0,0,w_ampa))
					self.nc_nmda[a].append( h.NetCon(netstim,self.choose_syn_nmda[a][b],0,0,w_nmda))
					self.nc_clopath[a].append( h.NetCon(netstim,self.choose_syn_clopath[a][b],0,0,w_ampa))

if __name__ == "__main__":
    pass
