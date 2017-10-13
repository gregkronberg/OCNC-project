"""
stims_2
"""

import numpy as np
from neuron import h

class DCS:
	def __init__(self):
		self.dcs()

	def seg_location(self, sec):
		""" given a neuron section, output the 3d coordinates of segment in the section
		"""
		# number of 3d points in section
        n3d = int(h.n3d(sec=sec))
        seg_x = []
        seg_y = []
        seg_z = []
		for seg_i,seg in enumerate(sec):
                # relative position within section (0-1)
                seg_pos = seg.x            
                # preallocate 3d coordinates
                x = [None]*n3d
                y = [None]*n3d
                z = [None]*n3d
                position_3d =  [None]*n3d
                               
                # loop over 3d coordinates in each section
                for i in range(n3d):
                	# retrieve x,y,z
                    x[i] = h.x3d(i,sec=sec)
                    y[i] = h.y3d(i,sec=sec)
                    z[i] = h.z3d(i,sec=sec)

                    # calculate total distance of each 3d point from start of section
                    if i is 0:
                    	position_3d[i] = 0
                    else:
                    	position_3d[i] = position_3d[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + (z[i]-z[i-1])**2)
                
                # segment distance along section in 3D
                seg_dist = seg_pos*position_3d[-1]
                # find first 3D coordinate that contains the segment
                node_i = [dist_i for dist_i,dist in enumerate(position_3d) if dist >= position_seg]
                
                # if segement occurs exactly at a node set its location to the node location
                
                if position_3d[node_i[0]] == seg_dist:
                    seg_x,seg_y,seg_z = x[node_i[0]],z[node_i[0]],z[node_i[0]]
                # otherwise if segment falls between two coordinates, interpolate to get location
                else:
                    seg_x.append( np.mean( [x[ node_i[0]], x[ node_i[0]-1]]))
                    seg_y.append( np.mean( [y[ node_i[0]], y[ node_i[0]-1]]))
                    seg_z.append( np.mean( [z[ node_i[0]], z[ node_i[0]-1]]))
        return [seg_x,seg_y,seg_z]


	def dcs(self, cell=0):
        
        if cell == 0:
        	cell == h.allsec()

        # structure to store location and e_extracellular for each segment.  Organized as ['dimension'][section number][segment number]
        location = {'x':[], 'y':[],'z':[],'e':[]}
        # loop over sections
        for sec_i,sec in enumerate(cell):

        	# add list for each section to store data
        	for dim_key,dim in location.iteritems():
        		dim.append([])

            # insert extracellular mechanism
            sec.insert('extracellular')

            # number of 3d points in section
            n3d = int(h.n3d(sec=sec))
            
            # find 3d position of each segment
            for seg_i,seg in enumerate(sec):
                # relative position within section (0-1)
                seg_pos = seg.x            
                # preallocate 3d coordinates
                x = [None]*n3d
                y = [None]*n3d
                z = [None]*n3d
                position_3d =  [None]*n3d
                               
                # loop over 3d coordinates in each section
                for i in range(n3d):
                	# retrieve x,y,z
                    x[i] = h.x3d(i,sec=sec)
                    y[i] = h.y3d(i,sec=sec)
                    z[i] = h.z3d(i,sec=sec)

                    # calculate total distance of each 3d point from start of section
                    if i is 0:
                    	position_3d[i] = 0
                    else:
                    	position_3d[i] = position_3d[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + (z[i]-z[i-1])**2)
                
                # segment distance along section in 3D
                seg_dist = seg_pos*position_3d[-1]
                # find first 3D coordinate that contains the segment
                node_i = [dist_i for dist_i,dist in enumerate(position_3d) if dist >= position_seg]
                
                # if segement occurs exactly at a node set its location to the node location
                
                if position_3d[node_i[0]] == seg_dist:
                    seg_x,seg_y,seg_z = x[node_i[0]],z[node_i[0]],z[node_i[0]]
                # otherwise if segment falls between two coordinates, interpolate to get location
                else:
                    seg_x = np.mean([x[node_i[0]], x[node_i[0]-1]])
                    seg_y = np.mean([y[node_i[0]], y[node_i[0]-1]])
                    seg_z = np.mean([z[node_i[0]], z[node_i[0]-1]])
                 
                # angle of segment from somato-dendritic axis (neglect z axis)   
                if seg_y == 0:
                    angle = 0
                elif np.isnan(seg_x/float(seg_y)):
                    angle = 0
                else:
                    angle = np.arctan(seg_x/seg_y)

                # if y location is negative shift phase by pi
                if seg_y < 0:
                    angle = angle+np.pi
                
                # absolute distance of segment from (0,0) in um
                mag = np.sqrt(seg_x**2 + seg_y**2)
                
                # angle relative to electric field vector, zero angle means along somato-dendritic axis
                angle_field = angle + field_angle
                
                # convert um to mm
                conversion = .001 

                e = conversion*intensity*mag*np.cos(angle_field)
                # insert calculated extracellular potential in mV
                seg.e_extracellular = conversion*intensity*mag*np.cos(angle_field)

                # store in structure ['dimension'][section number][segment number]
                location['x'][sec_i].append(seg_x)
                location['y'][sec_i].append(seg_y)
                location['z'][sec_i].append(seg_z)
                location['e'][sec_i].append(e)

                return location