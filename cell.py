"""
create cells and activate subsets of synapses
"""
# imports
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import stims

# create cell
class Cell_Migliore_2005:
	"""
	pyramidal neuron based on Migliore et al. 2005
	"""
	def __init__(self):
		self._init_parameters()
		self._init_geometry()
		self._init_axon()
		self._init_soma()
		self._init_basal()
		self._init_apical_trunk()
		self._init_apical_tuft()

	def _init_parameters(self):
		#=============================================================================
		# Parameters from Migliore 2005 (signal propogation in oblique dendrites)
		#=============================================================================
		# parameters are reported as (nS/um2) in paper, but need to be in (mho/cm2)
		# conversion 10*(nS/um2) = (mho/cm2)
		# units in paper are a typo, values are already reported in (mho/cm2)
		self.Vrest   = -65				# resting potential (mV)
		self.gna     =  .25				# peak sodium conductance (mho/cm2)
		self.ena = 55					# sodium reversal potential (mV)
		self.AXONM   = 5				# multiplicative factor for axonal conductance
		self.gkdr    = 0.1				# delayed rectifier potassium peak conductance (mho/cm2)
		self.ek = -90					# potassium reversal potential
		h.celsius = 35.0  				# temperature (degrees C)
		self.KMULT   =  0.3			# multiplicative factor for distal A-type potassium conductances
		self.KMULTP  = 0.3				# multiplicative factor for proximal A-type potassium conductances
		self.ghd     = 0.0005			# peak h-current conductance (mho/cm2)
		self.ehd = -30					# h-current reversal potential (mV)
		self.vhalfl_prox  = -73			# activation threshold for proximal a-type potassium (mV)
		self.vhalfl_dist = -81			# activation threshold for distal a-type potassium (mV)
		self.RaAll   = 150				# axial resistance, all compartments (ohm*cm)
		self.RaAx = 50					# axial resistance, axon (ohm*cm)					
		self.RmAll   = 28000			# specific membrane resistance (ohm/cm2)
		self.Cm =  1					# specific membrane capacitance (uf/cm2)
		self.ka_grad = 1				# slope of a-type potassium channel gradient with distance from soma 

		self.w_ampa = .00018 			# ampa peak conductance (microsiemens or micro-omhos)
		self.w_nmda = .00018 			# nmda peak conductance (microsiemens or micro-omhos)
		# self.syn_types = ['ampa','nmda','clopath']

	def _init_geometry(self):
		"""
		create cell geometry at hoc top level
		"""
		h.load_file('geo5038804.hoc')   # load cell class geometry from Migliore 2005
		h.load_file('fixnseg.hoc')  	# set discretization based on dlambda rule (set dlambda in hoc file)

		# self.geometry = {'soma':h.soma[0],
		# 'axon':h.axon[0],
		# 'dend_b':h.dendrite,
		# 'dend_a_trunk':h.user5,
		# 'dend_a_tuft':h.apical_dendrite}
		self.soma = h.soma[0]
		self.axon = h.axon[0]
		self.dend_b = h.dendrite
		self.dend_a_trunk = h.user5
		self.dend_a_tuft = h.apical_dendrite

	def _init_axon(self):
		"""
		insert axon biophysics
		"""
		# passive biophysics
		self.axon.insert('pas')
		self.axon.g_pas = 1/self.RmAll				# passive conductance (S/cm2)
		self.axon.e_pas = self.Vrest				# leak reversal potential (mV)
		self.axon.cm = self.Cm 						# specific capacitance (uf/cm2)
		self.axon.Ra = self.RaAx 					# axial resistance (ohm cm) 

		# active biophysics (parmameters: Migliore 2005)
		self.axon.insert('nax')						# voltage gated sodium
		self.axon.gbar_nax = self.gna*self.AXONM
		self.axon.insert('kdr')						# delayed rectifier potassium
		self.axon.gkdrbar_kdr = self.gkdr
		self.axon.insert('kap')						# a-type potassium
		self.axon.gkabar_kap = self.KMULTP
		self.axon.ena = self.ena					# sodium reversal potential (see _init_parameters)
		self.axon.ek = self.ek						# potassium reversal potential (see _init_parameters)

	def _init_soma(self):
		"""
		insert soma biophysics
		"""

		# passive biophysics
		self.soma.insert('pas')
		self.soma.g_pas = 1/self.RmAll				# passive conductance (S/cm^2)
		self.soma.e_pas = self.Vrest				# leak reversal potential (mV) 
		self.soma.Ra = self.RaAll					# axial resistance (ohm cm)     
		self.soma.cm = self.Cm                      # specific membrane capacitance (uf/cm2)    

		# soma biophysics (parmameters: Migliore 2005)
		self.soma.insert('na3')
		self.soma.gbar_na3 = self.gna				# voltage gated sodium
		self.soma.insert('hd')
		self.soma.ghdbar_hd = self.ghd				# h-current
		self.soma.vhalfl_hd = self.vhalfl_prox		# h-current activation threshold
		self.soma.insert('kdr')
		self.soma.gkdrbar_kdr = self.gkdr			# delayed rectifier potassium
		self.soma.insert('kap')
		self.soma.gkabar_kap = self.KMULTP			# scaling factor for a-type potassium current
		self.soma.ena = self.ena					# sodium reversal potential 
		self.soma.ek = self.ek						# potassium reversal potential 
		# self.soma.ehd_hd = ehd

		# set soma as origin for distance measurements
		h.distance(sec = self.soma)

	def _init_basal(self):
		"""
		insert basal dendrites biophysics and synapses
		"""  

		# preallocate synapse objects.  synapses will be indexed by [section number][segment number]
		self.syn_b_ampa = []						# ampa synapse object
		self.syn_b_nmda = []						# nmda_synapse object
		self.syn_b_clopath = []						# clopath synapse object (for implementing plasticity)
		# loop over basal dendrite sections
		for a in range(0,len(self.dend_b)):      
		    # passive biophysics
		    self.dend_b[a].Ra = self.RaAll          # axial resistance (units?)
		    self.dend_b[a].insert('pas')       		# insert passive mechanism
		    self.dend_b[a].g_pas = 1/self.RmAll   	# passive conductance (S/cm^2)
		    self.dend_b[a].e_pas = self.Vrest       # leak reversal potential (mV)
		    self.dend_b[a].cm = self.Cm 			# specific capacitance (uf/cm^2)

		    # active biophysics (parameters: Migliore 2005)
		    self.dend_b[a].insert('hd')
		    self.dend_b[a].ghdbar_hd = self.ghd		# h-current
		    self.dend_b[a].insert('na3')
		    self.dend_b[a].gbar_na3 = self.gna		# voltage gated sodium
		    self.dend_b[a].insert('kdr')
		    self.dend_b[a].gkdrbar_kdr = self.gkdr	# delayed rectifier potassium
		    self.dend_b[a].insert('kap')
		    self.dend_b[a].gkabar_kap = 0			# a-type potassium proximal
		    self.dend_b[a].insert('kad')
		    self.dend_b[a].gkabar_kad = 0			# a-type potassium distal
		    self.dend_b[a].ena = self.ena			# sodium reversal potential 
		    self.dend_b[a].ek = self.ek				# potassium reversal potential 
		    # self.dend_b[a].ehd = ehd
		    
		    # channels that vary with distance from soma (parameters: Migliore 2005)
		    # loop over segments     
		    cnt=-1									# counter for segment number
		    self.syn_b_ampa.append([])
		    self.syn_b_nmda.append([])
		    self.syn_b_clopath.append([])
		    for seg in self.dend_b[a]:
		        cnt+=1
		        # distance from soma
		        seg_dist = h.distance(seg.x,sec=self.dend_b[a])
		        
		        # h current
		        seg.ghdbar_hd = self.ghd*(1+3*seg_dist/100)
		        
		        # A-type potassium
		        if seg_dist > 100:	# distal
		            seg.vhalfl_hd = self.vhalfl_dist
		            seg.gkabar_kad = self.KMULT*(1+self.ka_grad*seg_dist/100)
		        else:	# proximal
		            seg.vhalfl_hd = self.vhalfl_prox
		            seg.gkabar_kap = self.KMULTP*(1+self.ka_grad*seg_dist/100)

			    # AMPA synapses (parameters: Kim 2015)
			    self.syn_b_ampa[a].append(h.Exp2Syn(self.dend_b[a](seg.x)))
			    self.syn_b_ampa[a][cnt].tau1 = 0.2
			    self.syn_b_ampa[a][cnt].tau2 = 2
			    self.syn_b_ampa[a][cnt].i = 0.18

			    # NMDA synapses (parameters: Kim 2015)
			    self.syn_b_nmda[a].append(h.Exp2SynNMDA(self.dend_b[a](seg.x)))
			    self.syn_b_nmda[a][cnt].tau1 = 1
			    self.syn_b_nmda[a][cnt].tau2 = 50

			     # Voltage based STDP learning rule (Clopath 2010)
		        self.syn_b_clopath[a].append(h.STDPSynCCNon(self.dend_b[a](seg.x)))

	def _init_apical_trunk(self):
		"""
		insert apical trunk dendrites biophysics and synapses
		""" 

		# preallocate synapse objects.  synapses will be indexed by [section number][segment number] 
		self.syn_a_trunk_ampa = []
		self.syn_a_trunk_nmda = []
		self.syn_a_trunk_clopath = []
		# loop over dendrite sections
		for a in range(0,len(self.dend_a_trunk)):      
		    # passive biophysics (params: Migliore 2005)
		    self.dend_a_trunk[a].Ra = self.RaAll          	# axial resistance (units?)
		    self.dend_a_trunk[a].insert('pas')       		# insert passive mechanism
		    self.dend_a_trunk[a].g_pas = 1/self.RmAll   	# passive conductance (S/cm^2)
		    self.dend_a_trunk[a].e_pas = self.Vrest         # leak reversal potential (mV)
		    self.dend_a_trunk[a].cm = self.Cm 				# specific capacitance (uf/cm^2)

		    # active biophysics (params Migliore 2005)
		    self.dend_a_trunk[a].insert('hd')
		    self.dend_a_trunk[a].ghdbar_hd = self.ghd		# h current
		    self.dend_a_trunk[a].insert('na3')
		    self.dend_a_trunk[a].gbar_na3 = self.gna		# voltage gated sodium current
		    self.dend_a_trunk[a].insert('kdr')
		    self.dend_a_trunk[a].gkdrbar_kdr = self.gkdr	# delayed rectifier potassium
		    # A-type potassium proximal
		    self.dend_a_trunk[a].insert('kap')
		    self.dend_a_trunk[a].gkabar_kap = 0    # set to zero, then change if the synapse meets the "proximal" criteria
		    # A-type potassium distal
		    self.dend_a_trunk[a].insert('kad')
		    self.dend_a_trunk[a].gkabar_kad = 0    # set to zero, then change if the segment meets the "distal" criteria
		    self.dend_a_trunk[a].ena = self.ena				# sodium reversal potential
		    self.dend_a_trunk[a].ek = self.ek				# sodium reversal potential 
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
		        seg.ghdbar_hd = self.ghd*(1+3*seg_dist/100)
		        
		        # A-type potassium 
		        if seg_dist > 100:	# distal
		            seg.vhalfl_hd = self.vhalfl_dist
		            seg.gkabar_kad = self.KMULT*(1+self.ka_grad*seg_dist/100)
		        else: # proximal
		            seg.vhalfl_hd = self.vhalfl_dist
		            seg.gkabar_kap = self.KMULTP*(1+self.ka_grad*seg_dist/100)

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

	def _init_apical_tuft(self):
		"""
		insert apical trunk dendrites biophysics and synapses
		""" 

		# preallocate synapse objects.  synapses will be indexed by [section number][segment number] 
		self.syn_a_tuft_ampa = []
		self.syn_a_tuft_nmda = []
		self.syn_a_tuft_clopath = []

		# loop over dendrite sections
		for a in range(0,len(self.dend_a_tuft)):  

		    # passive biophysics
		    self.dend_a_tuft[a].Ra = self.RaAll          # axial resistance (units?)
		    self.dend_a_tuft[a].insert('pas')       # insert passive mechanism
		    self.dend_a_tuft[a].g_pas = 1/self.RmAll   # passive conductance (S/cm^2)
		    self.dend_a_tuft[a].e_pas = self.Vrest         # leak reversal potential (mV)
		    self.dend_a_tuft[a].cm = self.Cm 			# specific capacitance (uf/cm^2)

		    # active biophysics (params: Migliore 2005)
		    # h current
		    self.dend_a_tuft[a].insert('hd')
		    self.dend_a_tuft[a].ghdbar_hd = self.ghd
		    #sodium current
		    self.dend_a_tuft[a].insert('na3')
		    self.dend_a_tuft[a].gbar_na3 = self.gna
		    # delayed rectifier potassium
		    self.dend_a_tuft[a].insert('kdr')
		    self.dend_a_tuft[a].gkdrbar_kdr = self.gkdr
		    # A-type potassium proximal
		    self.dend_a_tuft[a].insert('kap')
		    self.dend_a_tuft[a].gkabar_kap = 0    # set to zero, then change if the synapse meets the "proximal" criteria
		    # A-type potassium distal
		    self.dend_a_tuft[a].insert('kad')
		    self.dend_a_tuft[a].gkabar_kad = 0    # set to zero, then change if the synapse meets the "distal" criteria
		     
		    self.dend_a_tuft[a].ena = self.ena
		    self.dend_a_tuft[a].ek = self.ek
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
		        seg.ghdbar_hd = self.ghd*(1+3*seg_dist/100)
		        
		        # A-type potassium current
		        if seg_dist > 100:	# distal
		            seg.vhalfl_hd = self.vhalfl_dist
		            seg.gkabar_kad = self.KMULT*(1+self.ka_grad*seg_dist/100)
		        else:	# proximal
		            seg.vhalfl_hd = self.vhalfl_dist
		            seg.gkabar_kap = self.KMULTP*(1+self.ka_grad*seg_dist/100)

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

	subtree - list containing synapse point process objects of  given type on a given tree, e.g. ampa synapses on basal dendrites.  The list should be organized subtree_type[section][segment].  Create a new class instance for each subtree to be stimulated

	sec_idx - list of section indices to be activated within the specified subtree

	seg_idx  - list of segments to be activated within each section.  Should be organized as [sections][segments]

	w - weights for NetCon objects
	"""
	def __init__(self,netstim,subtree_ampa,subtree_nmda,subtree_clopath,sec_idx,seg_idx,w_ampa,w_nmda):
		
		# choose active segments
		# neuronal subtree_ampa (soma, apical, basal, etc.)
		self.subtree_ampa  = subtree_ampa 			# list all ampa synapses [section][segment]
		self.subtree_nmda = subtree_nmda			# list all nmda synapses [section][segment]		
		self.subtree_clopath = subtree_clopath		# list all clopath synapses [section][segment]
		self.sec_idx = sec_idx						# list sections with active synapses [section]
		self.seg_idx = seg_idx						# list active segments [section][segment]
		self.w_ampa = w_ampa						# maximum ampa conductance (uS or umho)
		self.w_nmda = w_nmda						# maximum nmda conductance (uS or umho)
		
		# list of only active synapse objects [section][segment]
		self.choose_syn_ampa = [self.subtree_ampa[a] for a in self.sec_idx]
		self.choose_syn_nmda = [self.subtree_nmda[a] for a in self.sec_idx]
		self.choose_syn_clopath = [self.subtree_clopath[a] for a in self.sec_idx]
		self.active = []	# keep track of which synapses are active [section][segment], 1 = active, 0 = inactive

		# store NetCon objects in lists [section][segment]. Inactive synapses will have a 0 as a placeholder 
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
	cell_1 = Cell_Migliore_2005()
