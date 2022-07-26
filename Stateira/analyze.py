# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:02:31 2021

@author: alexa
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


class scat:
    def __init__(self, file_name):
        self.file_name = file_name
        f = tb.open_file(file_name, 'r')
        
        configArrayObject = f.get_node("/", "CONFIG")
        self.config_arr = configArrayObject.read()
        
        xmatArrayObject = f.get_node("/", "XMAT")
        kmatArrayObject = f.get_node("/", "KMAT")
        posindArrayObject = f.get_node("/", "POS_IND")
        
        kvecArrayObject = f.get_node("/", "KVEC")
                  
        if self.config_arr[7] == 1:
            alphaArrayObject = f.get_node("/", "ALPHA_ARR")
            self.alpha_arr = alphaArrayObject.read()
            blArrayObject = f.get_node("/", "B_L")
            self.bl_arr = blArrayObject.read()
            errorArrayObject = f.get_node("/", "ERROR")
            self.error_arr = errorArrayObject.read()
        
        
        self.kmat = kmatArrayObject.read()
        self.xmat = xmatArrayObject.read()
        self.pos_ind = posindArrayObject.read()
        
        self.kvec_arr = kvecArrayObject.read()
        
        self.NC = int(self.config_arr[0])
        self.NP1 = int(self.config_arr[1])
        self.R_START = self.config_arr[2]
        self.R_STOP = self.config_arr[3]
        self.E_IN = self.config_arr[4]
        self.NUM_PART = int(self.config_arr[5])
        self.NPOS = int(self.config_arr[6])
        
        f.close()
    
    def get_smat(self):
        tes = np.eye(self.NPOS) - 1.0j*self.kmat
        omega = np.dot(-2.0*self.kmat, inv(tes))
        smat = np.eye(self.NPOS) - 1.0j*omega
        return smat
    
    def get_scattering_length(self, in_chan, out_chan):
        k = self.kvec_arr[self.pos_ind[in_chan]]
        smat = self.get_smat()
        a = 1.0 /(1.0j * k) * (1-smat[in_chan, out_chan])/(1+smat[in_chan, out_chan])
        a = a * conf.r_scaling
        return a
    
    """
    Returns scattering cross section in cm^2
    """
    def get_cross_section_cm(self, in_chan, out_chan):
        
        smat = self.get_smat()
        s_ij = smat[in_chan, out_chan]
        
        if in_chan == out_chan:           
            t_ij = 1 - s_ij
        else:
            t_ij = - s_ij
            
        k = self.kvec_arr[self.pos_ind[in_chan]]
        sigma = np.pi/(k**2)*np.abs(t_ij)**2
        
        sigma = sigma*conf.r_scaling**2*10000
        
        return sigma
    
    def get_psi(self, r):
        i = 0
        ret = np.full((self.NC, self.NC), 0.0)
        while i < self.NUM_PART:
            if self.bl_arr[i+1] > r:
                for j in range(self.NP1):
                    ret = ret + self.alpha_arr[i*self.NP1*self.NC + self.NC*j:i*self.NP1*self.NC + self.NC*j+self.NC, :] * special.eval_chebyt(j, maprtox(self.bl_arr[i], self.bl_arr[i+1], r))
                break
            i = i + 1
        return ret
    
    def get_scat_sol(self, r, in_chan):
        c = 0
        k = self.kvec_arr[self.pos_ind[in_chan]]
        psimat = self.get_psi(r)
        psi_scat = np.full((self.NC), 0.0)
        for s in self.pos_ind:
            psi_scat = psi_scat + np.sqrt(k)*self.xmat[in_chan,c] * psimat[:,s]
            c = c + 1
            
        return psi_scat
    

    
    def get_scat_sol_arr(self, r_arr, in_chan):
        matrix_list = np.full((len(r_arr),self.NC),0.0)
        for i in range(len(r_arr)):
            matrix_list[i,:]= self.get_scat_sol(r_arr[i], in_chan)
        return matrix_list



    

