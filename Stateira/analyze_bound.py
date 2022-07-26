# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 14:33:51 2022

@author: Alex
"""

import numpy as np
from numpy.linalg import inv
import scipy 
from scipy import special
from pathlib import Path
import tables as tb
import os.path 
from .Stateira_config import Stateira_config as conf

"""
Maps [b1 < r < b2] to [-1, 1]
"""
def maprtox(b1, b2, r):
    x = (2*r - (b2 + b1)) / (b2 - b1)
    return x

"""
Maps [-1, 1] to [b1 < r < b2]
"""
def mapxtor(b1, b2, x):
    r = 0.5*(b2 - b1) * x + 0.5*(b2 + b1)
    return r


class bound_wf:
    def __init__(self, file_name, level_ind):
        f = tb.open_file(file_name, 'r')
        configArrayObject = f.get_node("/", "CONFIG")
        alphaInArrayObject = f.get_node("/WF_"+str(level_ind)+"/", "ALPHA_IN")
        alphaOutArrayObject = f.get_node("/WF_"+str(level_ind)+"/", "ALPHA_OUT")
        blArrayObject = f.get_node("/WF_"+str(level_ind)+"/", "BL")
        configWfArrayObject = f.get_node("/WF_"+str(level_ind)+"/", "CONFIG_WF")
        lvecArrayObject = f.get_node("/WF_"+str(level_ind)+"/", "L_VEC")
        rvecArrayObject = f.get_node("/WF_"+str(level_ind)+"/", "R_VEC")
        levelsArrayObject = f.get_node("/", "LEVELS")

        self.config_arr = configArrayObject.read()
        self.levels_arr = levelsArrayObject.read()
        self.alpha_in = alphaInArrayObject.read()
        self.alpha_out = alphaOutArrayObject.read()
        self.bl = blArrayObject.read()
        self.config_wf = configWfArrayObject.read()
        self.l_vec = lvecArrayObject.read()
        self.r_vec = rvecArrayObject.read()

        self.NC = int(self.config_arr[0])
        self.NP1 = int(self.config_arr[1])
        self.R_START = self.config_arr[2]
        self.R_STOP = self.config_arr[3]
        
        self.NUM_PART = int(self.config_wf[5])
        self.R_MID = int(self.config_wf[4])
        self.e_val = self.levels_arr[level_ind]
        self.norm = -1

        f.close()
    
    def get_psi_io(self, r, alpha_l):
        i = 0
        ret = np.full((self.NC,self.NC), 0.0)
        while i < self.NUM_PART:
            if self.bl[i+1] > r:
                for j in range(self.NP1):
                    ret = ret + alpha_l[i*self.NP1*self.NC + self.NC*j:i*self.NP1*self.NC + self.NC*j+self.NC, :] * special.eval_chebyt(j, maprtox(self.bl[i], self.bl[i+1], r))
                break
            i = i + 1
        return ret
    
    def get_psi(self, r):
        if r <= self.R_MID:
            psi = self.get_psi_io(r, self.alpha_out)
            psi = np.dot(psi, self.l_vec )
        else:
            psi = self.get_psi_io(r, self.alpha_in)
            psi = np.dot(psi, self.r_vec)
            
        return psi
    
    def get_psi_array(self, r_arr):
        psi_arr = np.full((len(r_arr), self.NC),0.0)
        for i in range(len(r_arr)):
            r = r_arr[i]
            psi_arr[i,:] = self.get_psi(r)
            
        return psi_arr
    
    def get_psi_norm(self):
        num_points = 600
        r_arr = np.linspace(self.R_START, self.R_STOP, num_points)
        psi_arr = self.get_psi_array(r_arr)
        delta_r = (self.R_STOP - self.R_START) / num_points
        s = 0
        for i in range(self.NC):
            s = s + np.sum(np.abs(psi_arr[:,i])**2)
        s = s * delta_r
        self.norm = s
        
        return s
    
    def get_psi_range(self):
        
        if self.norm < 0:
            self.get_psi_norm()
            
        num_points = 600
        r_arr = np.linspace(self.R_START, self.R_STOP, num_points)
        psi_arr = self.get_psi_array(r_arr)
        delta_r = (self.R_STOP - self.R_START) / num_points
        abs_sq = np.full((num_points),0.0)
        for i in range(self.NC):
            abs_sq = abs_sq + np.abs(psi_arr[:,i])**2
        abs_sq = abs_sq / self.norm * delta_r
        
        s0 = 0
        for i in range(num_points):
            s1 = s0 + abs_sq[i]
            
            if s1 > 0.001 and s0 <= 0.001:
                r_low = r_arr[i]
                
            if s1 > 0.999 and s0 <= 0.999:
                r_high = r_arr[i]
                
            s0 = s1
        
        return r_low, r_high
        
class bound:
    def __init__(self, file_name):
        self.file_name = file_name
        f = tb.open_file(file_name, 'r')
        
        configArrayObject = f.get_node("/", "CONFIG")
        self.config_arr = configArrayObject.read()
        
        levelArrayObject = f.get_node("/", "LEVELS")
        self.level_arr = levelArrayObject.read()
        self.num_level = len(self.level_arr)
        
        self.NC = int(self.config_arr[0])
        
        f.close()