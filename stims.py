# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 03:06:39 2017

@author: Greg Kronberg
"""
import numpy as np
from neuron import h

def dcs(field_angle,intensity,cell=0):
    if cell == 0:
        n_sec = 0
    #sec_list = [None]
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
            
#%%
#def tbs(syns,bursts):
#    warm_up = 500   # warm up time (ms)
#    pulse_freq = 100
#    burst_freq = 5
##    bursts = 1
#    pulses = 4
#    stim = [None]*bursts
#    nc = [None]*bursts
#    for a in range(0,bursts):
#        stim[a] = h.NetStim()
#        stim[a].start = warm_up + a*1000/burst_freq
#        stim[a].interval = 1000/pulse_freq
#        stim[a].noise  = 0 
#        stim[a].number = pulses
#        w0 = 1
#        nc[a] = [None]*len(syns)
#        i = -1 
#        for act in subset_a:
#            i=i+1
#            nc[a][i] = h.NetCon(stim[a],syn_a[act],0,0,w0)