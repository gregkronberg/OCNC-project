"""
implement extracellular stimulation or presynaptic input patterns

Created on Wed Jun 28 03:06:39 2017

@author: Greg Kronberg
"""
import numpy as np
from neuron import h

# extracellular field
def dcs(field_angle,intensity,cell=0):
    """
    Apply DC extracellular field by inserting extracellular mechanism in each segment

    Arguments:

    field_angle - angle relative to the somato-dendritic axis

    intensity - stimulation intensity in V/m

    cell - which cell to stimulate. if 0, stimulate all segments in the top level of hoc
    """
    
    if cell == 0:
        n_sec = 0
        # loop over sections
        for section in h.allsec():
            
            section.insert('extracellular')
            n_sec += 1
            n3d = int(h.n3d(sec=section))
            
            # find 3d position of each segment
            for seg in section:
                position_seg = seg.x            # relative position within section (0-1)
                # preallocate 3d coordinates
                x = [None]*n3d
                y = [None]*n3d
                z = [None]*n3d
                position_3d =  [None]*n3d
                               
                # list all 3d coordinates in section
                for i in range(0,n3d):
                    x[i] = h.x3d(i,sec=section)
                    y[i] = h.y3d(i,sec=section)
                    z[i] = h.z3d(i,sec=section)
                    # find distance of each coordinate from first node (section range)
                    position_3d[i] =  np.sqrt((x[i]-x[0])**2 + (y[i]-y[0])**2 + (z[i]-z[0])**2)
                
                # segment location along section in 3D
                temp1 = position_seg*position_3d[-1]
                # find first 3D coordinate that contains the segment
                temp2 = [a for a in range(0,n3d) if  position_3d[a] >= temp1]
                # if segment falls between two coordinates, interpolate to get location
                if position_3d[temp2[0]] == temp1:
                    seg_position_x,seg_position_y,seg_position_z = x[temp2[0]],z[temp2[0]],z[temp2[0]]
                else:
                    seg_position_x = np.mean([x[temp2[0]],x[temp2[0]-1]])
                    seg_position_y = np.mean([y[temp2[0]],y[temp2[0]-1]])
                    seg_position_z = np.mean([z[temp2[0]],z[temp2[0]-1]])
                 
                # angle from somato-dendritic axis (neglect z axis)   
                if seg_position_y == 0:
                    angle = 0
                elif np.isnan(seg_position_x/float(seg_position_y)):
                    angle = 0
                else:
                    angle = np.arctan(seg_position_x/seg_position_y)
                    
                if seg_position_y < 0:
                    angle = angle+np.pi
                
                    # distance of segment from (0,0)
                mag = np.sqrt(seg_position_x**2 + seg_position_y**2)
                
                # angle relative to electric field vector
                angle_field = angle + field_angle
                
                # insert calculate extracellular potential
                seg.e_extracellular = .001*intensity*mag*np.cos(angle_field)
    else:        
        n_sec = 0
        #sec_list = [None]
        # loop over sections
        for section in cell.all:
            
            section.insert('extracellular')
            n_sec += 1
            n3d = int(h.n3d(sec=section))
            
            # find 3d position of each segment
            for seg in section:
                position_seg = seg.x            # relative position within section (0-1)
                # preallocate 3d coordinates
                x = [None]*n3d
                y = [None]*n3d
                z = [None]*n3d
                position_3d =  [None]*n3d
                               
                # list all 3d coordinates in section
                for i in range(0,n3d):
                    x[i] = h.x3d(i,sec=section)
                    y[i] = h.y3d(i,sec=section)
                    z[i] = h.z3d(i,sec=section)
                    # find distance of each coordinate from first node (section range)
                    position_3d[i] =  np.sqrt((x[i]-x[0])**2 + (y[i]-y[0])**2 + (z[i]-z[0])**2)
                
                # segment location along section in 3D
                temp1 = position_seg*position_3d[-1]
                # find first 3D coordinate that contains the segment
                temp2 = [a for a in range(0,n3d) if  position_3d[a] >= temp1]
                # if segment falls between two coordinates, interpolate to get location
                if position_3d[temp2[0]] == temp1:
                    seg_position_x,seg_position_y,seg_position_z = x[temp2[0]],z[temp2[0]],z[temp2[0]]
                else:
                    seg_position_x = np.mean([x[temp2[0]],x[temp2[0]-1]])
                    seg_position_y = np.mean([y[temp2[0]],y[temp2[0]-1]])
                    seg_position_z = np.mean([z[temp2[0]],z[temp2[0]-1]])
                 
                # angle from somato-dendritic axis (neglect z axis)   
                if seg_position_y == 0:
                    angle = 0
                elif np.isnan(seg_position_x/float(seg_position_y)):
                    angle = 0
                else:
                    angle = np.arctan(seg_position_x/seg_position_y)
                    
                if seg_position_y < 0:
                    angle = angle+np.pi
                
                    # distance of segment from (0,0)
                mag = np.sqrt(seg_position_x**2 + seg_position_y**2)
                
                # angle relative to electric field vector
                angle_field = angle + field_angle
                
                # insert calculate extracellular potential
                seg.e_extracellular = .001*intensity*mag*np.cos(angle_field)
            
class tbs:
    """
    creates NetStim object for delivering theta burst stimulation
    """
    def __init__(self,bursts=1,pulses=4,pulse_freq=100,burst_freq=5):
        self.warm_up = 60   # warm up time (ms)
        self.stim  = [] # list of stim objects
        for a in range(0,bursts): # create new object for each burst
            self.stim.append(h.NetStim())
            self.stim[a].start = self.warm_up + a*1000/burst_freq # start of burst
            self.stim[a].interval = 1000/pulse_freq
            self.stim[a].noise  = 0 
            self.stim[a].number = pulses