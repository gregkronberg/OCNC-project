"""
create cells and activate subsets of synapses
"""
# imports
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import stims
import param
import run_control

class PyramidalCell:
	""" 4 compartment pyramidal cell
	"""
	def __init__(self, p):
		self.geometry(p)
		self.mechanisms(p)

	def geometry(self, p):
		"""
		areas determined from cell geo5038804 from Migliore 2005
		basal: 19966.3598 um2
		apical: 23700.5664916 (tuft) + 11461.4440485 (trunk) = 35162.0105401 um2
		soma: 176.290723263 um2
		axon: 305.021056197 um2
		"""
		# list sections
		trees = ['basal', 'soma', 'apical_prox', 'apical_dist']
		areas = [19966.36, 176.29, 35162.01/2., 35162.01/2.]
		
		# store geometry
		self.geo = {}
		# store synapses
		self.syns = {}
		# create sections
		for tree_i, tree in enumerate(trees):
			self.geo[tree] = []
			self.syns[tree] = {
			'ampa' : [],
			'nmda' : [],
			'clopath' : []
			}
			for sec_i in range(p['nsec_'+tree]):
				self.geo[tree].append( h.Section( name=tree))

				# diameter basd on area of full morphology
				# diam = areas[tree_i]/(np.pi*p['L_'+tree])
				diam = p['diam_'+tree] 
				# p['diam_'+tree] = diam
				# create 3d specification, with cell arranged vertically
				h.pt3dadd(0, 0, 0, diam, sec=self.geo[tree][sec_i])
				if tree=='basal':
					h.pt3dadd(0, -p['L_'+tree], 0, diam, sec=self.geo[tree][sec_i])
				else:
					h.pt3dadd(0, p['L_'+tree], 0, diam, sec=self.geo[tree][sec_i])

				# add list to store synapses for each section
				for syn_key, syn in self.syns[tree].iteritems():
					syn.append([])

				self.geo[tree][sec_i].insert('pas')
				# passive conductance (S/cm2)
				self.geo[tree][sec_i].g_pas = 1/p['RmAll']			
				# leak reversal potential (mV)	
				self.geo[tree][sec_i].e_pas = p['Vrest']				
				# specific capacitance (uf/cm2)
				self.geo[tree][sec_i].cm = p['Cm'] 			
				# axial resistance (ohm cm) 		
				self.geo[tree][sec_i].Ra = p['RaAll'] 

				self.geo[tree][sec_i].L = p['L_'+tree]
				self.geo[tree][sec_i].diam = p['diam_'+tree]

		self.geo['basal'][0].connect(self.geo['soma'][0](0),0)
		self.geo['apical_prox'][0].connect(self.geo['soma'][0](1),0)
		self.geo['apical_dist'][0].connect(self.geo['apical_prox'][0](1),0)



		h.xopen('fixnseg.hoc')

		# set temperature in hoc
		h.celsius = p['celsius']
		# set soma as origin for distance measurements
		h.distance(sec=self.geo['soma'][0])
	
	# def measure_area(self, geo):
	# 	cell_temp = CellMigliore2005(p)
	# 	area={}
	# 	for tree_key, tree in geo.iteritems():
	# 		area[tree_key] = cell_temp.measure_area(tree)
			
	# 	return area

	def rotate(self, theta):
		"""Rotate the cell about the Z axis.
		"""
		for sec in h.allsec():
			for i in range(int(h.n3d(sec=sec))):
				x = h.x3d(i, sec=sec)
				y = h.y3d(i, sec=sec)
				c = np.cos(theta)
				s = np.sin(theta)
				xprime = x * c - y * s
				yprime = x * s + y * c
				h.pt3dchange(i, xprime, yprime, h.z3d(i, sec=sec), h.diam3d(i, sec=sec), sec=sec)

	def mechanisms(self, p):
		for tree_key, tree in self.geo.iteritems():
			for sec_i, sec in enumerate(tree): 

				if tree_key == 'soma':
					# voltage gated sodium
					sec.insert('na3')
					sec.gbar_na3 = p['gna']	
					# h-current			
					sec.insert('hd')
					sec.ghdbar_hd = p['ghd']				
					sec.vhalfl_hd = p['vhalfl_prox']
					# delayed rectifier potassium		
					sec.insert('kdr')
					sec.gkdrbar_kdr = p['gkdr']	
					# a-type potassium		
					sec.insert('kap')
					sec.gkabar_kap = p['KMULTP']
					# sodium reversal potential 
					sec.ena = p['ena']		
					# potassium reversal potential 
					sec.ek = p['ek']
					# dendrites active biophysics
				
				elif ((tree_key == 'basal') or 
				(tree_key == 'apical_prox') or 
				(tree_key == 'apical_dist')):

					# h-current
				    sec.insert('hd')
				    sec.ghdbar_hd = p['ghd']
				    
				    # voltage gated sodium		
				    sec.insert('na3')
				    sec.gbar_na3 = p['gna']

				    # delayed rectifier potassium	
				    sec.insert('kdr')
				    sec.gkdrbar_kdr = p['gkdr']	
				    # a-type potassium proximal
				    sec.insert('kap')
				    sec.gkabar_kap = 0	
				    # a-type potassium distal		
				    sec.insert('kad')
				    sec.gkabar_kad = 0	
				    # sodium reversal potential 
				    sec.ena = p['ena']
				    # potassium reversal potential
				    sec.ek = p['ek']		

				    # mechanisms that vary with distance from soma
				    # loop over segments
				    for seg_i,seg in enumerate(sec):
				    	# distance from soma
				    	seg_dist = h.distance(seg.x, sec=sec)
				    	
				    	# sodium
				    	seg.gbar_na3 = p['gna'] + p['dgna']*seg_dist
				    	# print seg_dist, seg.gbar_na3
				    	
				    	# h current
				    	seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*seg_dist/100)
				    	
				    	# A-type potassium
				        if seg_dist > 100:	# distal
				            seg.vhalfl_hd = p['vhalfl_dist']
				            seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*seg_dist/100)
				        else:	# proximal
				            seg.vhalfl_hd = p['vhalfl_prox']
				            seg.gkabar_kap = p['KMULTP']*(1+p['ka_grad']*seg_dist/100)

				        # loop over synapse types
				        for syn_key,syn in self.syns[tree_key].iteritems():
				        	if syn_key == 'ampa':
				        		syn[sec_i].append(h.Exp2Syn(sec(seg.x)))
				        		syn[sec_i][seg_i].tau1 = p['tau1_ampa']
				        		syn[sec_i][seg_i].tau2 = p['tau2_ampa']
				        		syn[sec_i][seg_i].i = p['i_ampa']
				        	elif syn_key == 'nmda':
				        		syn[sec_i].append(h.Exp2SynNMDA(sec(seg.x)))
				        		syn[sec_i][seg_i].tau1 = p['tau1_nmda']
				        		syn[sec_i][seg_i].tau2 = p['tau2_nmda']
				        	elif syn_key == 'clopath':
				        		syn[sec_i].append(h.STDPSynCCNon(sec(seg.x)))



# create cell
class CellMigliore2005:
	""" pyramidal neuron based on Migliore et al. 2005

	An instance of this object will creates a cell (hoc objects) at the top level of the hoc interpreter using the hoc files in _init_geometry.  The .geo attribute contains a python mapping to these hoc objects.  The geo object is organized as geo['section tree'][section](segment location)

	the syns attribute creates a container for synapse objects that are added to each segment in the hoc cell.  syns is organized as syns['section tree']['synapse type'][section][segment number]
	
	"""
	def __init__(self,p):
		# initialize geometry
		self.geometry(p)
		# insert membrane mechanisms
		self.mechanisms(p)

	def geometry(self,p):
		""" create cell geometry at hoc top level
		
		"""
		# load cell geometry into hoc interpreter
		h.load_file('geo5038804.hoc')  
		# h.load_file('geoc62564.hoc')
		# set discretization based on dlambda rule (set dlambda in hoc file) 
		h.load_file('fixnseg.hoc')  	
		# dictionary for storing geometry ['tree'][sec](seg location)
		self.geo = {}
		# dictionary for storing synapse objects ['tree']['type'][sec][seg]
		self.syns = {}
		# add section trees to geometry dictionary
		self.geo['soma'] = h.soma
		self.geo['axon'] =  h.axon
		self.geo['basal'] = h.dendrite
		self.geo['apical_trunk'] = h.user5
		self.geo['apical_tuft'] = h.apical_dendrite
		
		# set temperature in hoc
		h.celsius = p['celsius']
		# set soma as origin for distance measurements
		h.distance(sec = self.geo['soma'][0])

	def mechanisms(self,p):
		""" insert membrane mechanisms

		self.syns is updated to store an object for each synapse.  It is organized as ['tree']['synapse type'][section][segment].  Note that the last index will depend on how the cell is discretized as the number segments changes in each sections 

		the parameters for each membrane mechanism  are store in a dictionary called p.  See the param module for details.
		"""

		# loop over trees
		for tree_key,tree in self.geo.iteritems():
			
			# create sub-dictionary for different types of synapses
			self.syns[tree_key] = {
			'ampa' : [],
			'nmda' : [],
			'clopath' : []
			}

			# loop over sections
			for sec_i,sec in enumerate(tree):
				
				# add list to store synapses for each section
				for syn_key,syn in self.syns[tree_key].iteritems():
					syn.append([])

				# common passive biophysics for all sections
				sec.insert('pas')
				# passive conductance (S/cm2)
				sec.g_pas = 1/p['RmAll']			
				# leak reversal potential (mV)	
				sec.e_pas = p['Vrest']				
				# specific capacitance (uf/cm2)
				sec.cm = p['Cm'] 			
				# axial resistance (ohm cm) 		
				sec.Ra = p['RaAll'] 
										


				# axon active bipophysics
				if tree_key == 'axon':
					# voltage gated sodium
					sec.insert('nax')						
					sec.gbar_nax = p['gna']*p['AXONM']
					# delayed rectifier potassium
					sec.insert('kdr')						
					sec.gkdrbar_kdr = p['gkdr']
					# a-type potassium
					sec.insert('kap')						
					sec.gkabar_kap = p['KMULTP']
					# sodium reversal potential 
					sec.ena = p['ena']		
					# potassium reversal potential 
					sec.ek = p['ek']
					sec.Ra = p['RaAx']
					
				# soma active biophysics
				elif tree_key == 'soma':
					# voltage gated sodium
					sec.insert('na3')
					sec.gbar_na3 = p['gna']	
					# h-current			
					sec.insert('hd')
					sec.ghdbar_hd = p['ghd']				
					sec.vhalfl_hd = p['vhalfl_prox']
					# delayed rectifier potassium		
					sec.insert('kdr')
					sec.gkdrbar_kdr = p['gkdr']	
					# a-type potassium		
					sec.insert('kap')
					sec.gkabar_kap = p['KMULTP']
					# sodium reversal potential 
					sec.ena = p['ena']		
					# potassium reversal potential 
					sec.ek = p['ek']			

				# dendrites active biophysics
				elif ((tree_key == 'basal') or 
				(tree_key == 'apical_trunk') or 
				(tree_key == 'apical_tuft')):
					# h-current
				    sec.insert('hd')
				    sec.ghdbar_hd = p['ghd']
				    
				    # voltage gated sodium		
				    sec.insert('na3')
				    sec.gbar_na3 = p['gna']

				    # delayed rectifier potassium	
				    sec.insert('kdr')
				    sec.gkdrbar_kdr = p['gkdr']	
				    # a-type potassium proximal
				    sec.insert('kap')
				    sec.gkabar_kap = 0	
				    # a-type potassium distal		
				    sec.insert('kad')
				    sec.gkabar_kad = 0	
				    # sodium reversal potential 
				    sec.ena = p['ena']
				    # potassium reversal potential
				    sec.ek = p['ek']		

				    # mechanisms that vary with distance from soma
				    # loop over segments
				    for seg_i,seg in enumerate(sec):
				    	
				    	# distance from soma
				    	seg_dist = h.distance(seg.x,sec=sec)
				    	
				    	# sodium
				    	seg.gbar_na3 = p['gna'] + p['dgna']*seg_dist
				    	# print seg_dist, seg.gbar_na3
				    	
				    	# h current
				    	seg.ghdbar_hd = p['ghd']*(1+p['ghd_grad']*seg_dist/100)
				    	
				    	# A-type potassium
				        if seg_dist > 100:	# distal
				            seg.vhalfl_hd = p['vhalfl_dist']
				            seg.gkabar_kad = p['KMULT']*(1+p['ka_grad']*seg_dist/100)
				        else:	# proximal
				            seg.vhalfl_hd = p['vhalfl_prox']
				            seg.gkabar_kap = p['KMULTP']*(1+p['ka_grad']*seg_dist/100)

				        # loop over synapse types
				        for syn_key,syn in self.syns[tree_key].iteritems():
				        	if syn_key == 'ampa':
				        		syn[sec_i].append(h.Exp2Syn(sec(seg.x)))
				        		syn[sec_i][seg_i].tau1 = p['tau1_ampa']
				        		syn[sec_i][seg_i].tau2 = p['tau2_ampa']
				        		syn[sec_i][seg_i].i = p['i_ampa']
				        	elif syn_key == 'nmda':
				        		syn[sec_i].append(h.Exp2SynNMDA(sec(seg.x)))
				        		syn[sec_i][seg_i].tau1 = p['tau1_nmda']
				        		syn[sec_i][seg_i].tau2 = p['tau2_nmda']
				        	elif syn_key == 'clopath':
				        		syn[sec_i].append(h.STDPSynCCNon(sec(seg.x)))
	
	def measure_area(self, tree):
		"""
		given a tree measure the total area of all sections in the tree

		tree is a list of sections (hoc objects)
		"""
		area_all = []
		for sec_i, sec in enumerate(tree):
			# convert to um to cm (*.0001)
			L = .0001*sec.L
			a = .0001*sec.diam/2.
			rL = sec.Ra
			rm = 1/sec.g_pas
			area = 2*np.pi*a*L
			lam = np.sqrt(a*rm/(2*rL))
			area_all.append(area)

		return sum(area_all)



class CellKim2015:
	pass

class Syn_act:
	"""Activate a specific subset of synpases with NetCon objects
	
	arguments: 
	syns = nested dictionary containing synapse objects organized as ['section tree']['synapse type'][section number][segment number]

	p =  dictionary of parameters including sec_idx and seg_idx lists, including section and segment numbers to be activated.  p also contains 'w_list', which sets the weight for each of the activated synapses. w_list has the same dimensions as seg_idx

	stim = list of NetStim objects to be connected (via NetCon) to the synapses in syns that designated by sec_idx and seg_idx. The stim list will be iterated through and a distinct NetCon object will be created for each NetStim for each activated synapse.  The need for multiple NetStim objects arises from complicated stimulation patterns, like theta bursts, which are not easily programmed with a single NetStim

	The created NetCon objects are referenced by the nc object, which is organized the same way as syns, namely ['section tree']['synapse type'][section number][segment number][NetStim object]
	"""
	def __init__(self, syns, p, stim):
		# store netcon objects ['tree']['syn type'][section][segment][list of netstim objects]
		self.nc = {}
		
		# loop over synapse types
		for tree_key,tree in syns.iteritems():
			if tree_key in p['tree']:
				self.nc[tree_key] = {}

				for syntype_key,syn_type in syns[tree_key].iteritems():
					self.nc[tree_key][syntype_key] = []

					# loop over active sections
					for sec_i,sec in enumerate(p['sec_idx']):
						self.nc[tree_key][syntype_key].append([])
						
						# loop over active segments
						for seg_i,seg in enumerate(p['seg_idx'][sec_i]):
							self.nc[tree_key][syntype_key][sec_i].append([])

							# loop over stimulation bursts
							for syn_stim_i,syn_stim in enumerate(stim):
								self.nc[tree_key][syntype_key][sec_i][seg_i].append( h.NetCon(syn_stim, syns[tree_key][syntype_key][sec][seg], 0, 0, p['w_list'][sec_i][seg_i]))

class DendriteTransform:
	def __init__(self, p):
		cell1 = CellMigliore2005(p)
		apical_transform = self.dendrite_transform(geo=cell1.geo, python_tree=['apical_trunk','apical_tuft'], neuron_tree=['user5', 'apical_dendrite'])
		basal_transform = self.dendrite_transform(geo=cell1.geo, python_tree=['basal'], neuron_tree=['dendrite'])
		print 'apical:', apical_transform['a_cable'], apical_transform['L_cable']
		print 'basal:', basal_transform['a_cable'], basal_transform['L_cable']

	def measure_area(self, tree):
		"""
		given a tree measure the total area of all sections in the tree

		tree is a list of sections (hoc objects)
		"""
		area_all = []
		for sec_i, sec in enumerate(tree):
			# convert to um to cm (*.0001)
			L = .0001*sec.L
			a = .0001*sec.diam/2.
			rL = sec.Ra
			rm = 1/sec.g_pas
			area = 2*np.pi*a*L
			lam = np.sqrt(a*rm/(2*rL))
			area_all.append(area)

		return sum(area_all)

	def measure_length(self, geo):
		""" measure electrotonic length for each path along a cells dendritic tree
		"""
		# keep track of most recent section in each path [paths]
		secs = [geo['soma'][0]]

		# keep track of all sections in paths [paths][sections]
		# does not include soma
		paths=[[]]
		
		# iterate over paths (most recent section)
		for sec_i, sec in enumerate(secs):
			
			# current section
			current_sec_ref = h.SectionRef(sec=sec)
			# current children
			current_children = current_sec_ref.child
			
			# if there are children
			while len(current_children)>0:
				
				# iterate over current children 
				for child_i, child in enumerate(current_children):

					# add first child to current path
					if child_i==0:
						paths[sec_i].append(child)
						# update current section
						sec = child
					
					# if multiple children
					if child_i > 0:
						# create new path to copy previous path when tree splits
						new_path=[]
						for section in paths[sec_i]:
							new_path.append(section)

						# copy current path in list
						# if split occurs at soma, do not copy previous list
						if h.SectionRef(sec=child).parent.name()!='soma':
							paths.append(new_path)
						else:
							# create new list, not including soma
							paths.append([])
						
						# add corresponding child to the current section list
						secs.append(child)
						
						# add corresponding child to new path 
						paths[sec_i + child_i].append(child)
				
				# update current section and children		
				current_sec_ref = h.SectionRef(sec=sec)
				current_children = current_sec_ref.child


		# calculate electrotonic length of each path
		path_l = [] # [paths][electrotonic section length]
		sec_name = [] # [paths][section names]
		for path_i, path in enumerate(paths):
			path_l.append([])
			sec_name.append([])
			for sec_i, sec in enumerate(path):
				# convert all distances in cm
				# section length
				L = .0001*sec.L
				# section radius
				a = .0001*sec.diam/2
				# membrane resistivity
				rm = 1/sec.g_pas
				# axial resistivity
				rL = sec.Ra
				# space constant lambda
				lam = np.sqrt( (a*rm) / (2*rL) )
				# electrotonic length
				e_length = L/lam
				# electrotonic lengths for all paths and sections [paths][sections]
				path_l[path_i].append(e_length) 
				# keep track of section names [paths][sections]
				sec_name[path_i].append(sec.name())
		# print path_l[0]
		return {'path_l': path_l,
		'sec_name':sec_name}

	def dendrite_transform(self, geo, python_tree, neuron_tree):
		""" equivalent cable transform for dendritic tree
		"""
		# FIXME
		rL_cable = 150.
		rm_cable = 28000.
		# area
		A_full=0
		for tree_i, tree in enumerate(python_tree):
			A_full += self.measure_area(geo[tree])

		paths = self.measure_length(geo)
		e_lengths = []
		# iterate over all paths
		for path_i, path in enumerate(paths['path_l']):
			# only keep paths in the neuron_tree argument
			for tree in neuron_tree:
				for sec in paths['sec_name'][path_i]:
					if tree in sec:
						e_lengths.append(sum(path))
						break

						

		EL_full = np.mean(e_lengths) 

		# convert cm back to um (*10000)
		L_cable = (EL_full**(2./3)) * (A_full*rm_cable/(4.*np.pi*rL_cable))**(1./3)
		a_cable = A_full/(2.*np.pi*L_cable)

		return {'a_cable':10000*a_cable, 'L_cable':10000*L_cable}
# set procedure if called as a script
if __name__ == "__main__":
	args = run_control.Arguments('exp_1').kwargs
	p = param.Experiment(**args).p
	DendriteTransform(p)
