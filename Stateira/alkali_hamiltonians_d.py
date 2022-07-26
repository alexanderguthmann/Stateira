# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 15:38:04 2022
Generates the hamiltonian terms for distinguishable alkali atoms.

@author: alexa
"""

import numpy as np
from numpy.linalg import multi_dot
from numpy.linalg import inv
from numpy import linalg as lin
from sympy.physics.quantum.cg import CG
import tables as tb
from scipy import linalg as sclin
from .Stateira_config import Stateira_config as conf

g_e = 2.00231930436256
mu_b = 9.27400915E-24
h_planck = 6.62606896E-34



"""
Needed configuration data
Default Rb - Cs Values:
"""

s1 = 0.5        #Spin Atom 1
i1 = 1.5        #Nuclear Spin Atom 1
s2 = 0.5        #Spin Atom 2
i2 = 3.5        #Nuclear Spin Atom 2

Lmax = 0        #Maximum L included in Basis

g_nuc_1 = -0.000995141          #Nuclear g factor Atom 1
g_nuc_2 = -0.00039885395        #Nuclear g factor Atom 2

zeta1 = 3.417341305452145E09 * h_planck     #Hyperfine constant atom 1 [Joules]
zeta2 = 2.2981579425E09 * h_planck          #Hyperfine constant atom 2 [Joules]

Lparity = 'even'    #Parity of L. 'even': L = 0, 2, ... , Lmax ; 'odd': L = 1, 3, ... , Lmax 



def spher_int(L1, ML1, q, L3, ML3):
    res = (-1.0)**ML1 * (-1.0)**ML3 * np.sqrt(5*(2*L1+1) / (4*np.pi * (2*L3 + 1)))
    res = res * CG(L1, 0, 2, 0, L3, 0).doit().evalf() * CG(L1, -1.0*ML1, 2, q, L3, -1.0*ML3).doit().evalf()
    return res


def construct_basis(Mk):
    basis=list()
    #vec = np.full((6),0.0)
    
    if Lparity == 'even':
        Lrange = np.arange(0, Lmax+1, 2)
    else:
        Lrange = np.arange(1, Lmax+1, 2)
    for L in Lrange:
        for ML in np.arange(L, -L-1,-1):
            for ms1 in np.arange(s1, -s1-1,-1):
                for mi1 in np.arange(i1, -i1-1,-1):
                    for ms2 in np.arange(s2, -s2-1,-1):
                        for mi2 in np.arange(i2, -i2-1,-1):
                            
                            if ML + ms1 + mi1 + ms2 + mi2 == Mk:
                                vec = [L, ML, ms1, mi1, ms2, mi2]
                                basis.append(vec)
                                
    return  np.asarray(basis)

def construct_dipole_dipole_operator(basis):
    W = np.full((len(basis), len(basis)), 0.0)
    
    r=0
    c=0
    for bra in basis:
        L = bra[0]
        ML = bra[1]
        ms1 = bra[2]
        mi1 = bra[3]
        ms2 = bra[4]
        mi2 = bra[5]
        
        c=0
        for ket in basis:
            Lb = ket[0]
            MLb = ket[1]
            ms1b = ket[2]
            mi1b = ket[3]
            ms2b = ket[4]
            mi2b = ket[5]
            if mi1 == mi1b and mi2 == mi2b:
                if L == Lb or L == Lb + 2 or L == Lb - 2:
                    #T0:
                    if ML == MLb and ms1 == ms1b and ms2 == ms2b:
                         W[r, c] = W[r, c] + np.sqrt(16*np.pi/5.0) * ms1 * ms2 * spher_int(L, ML, 0, Lb, ML)
                         
                    #T0b:
    
                    if ML == MLb and ms1 == ms1b + 1 and ms2 == ms2b - 1:
                        a_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b + 1))
                        b_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b - 1))
                        W[r, c] = W[r, c] -0.25 * np.sqrt(16*np.pi/5.0) * a_s1 * b_s2 * spher_int(L, ML, 0, Lb, ML)
                    if ML == MLb and ms1 == ms1b - 1 and ms2 == ms2b + 1:
                        b_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b - 1))
                        a_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b + 1))
                        W[r, c] = W[r, c] -0.25 * np.sqrt(16*np.pi/5.0) * b_s1 * a_s2 * spher_int(L, ML, 0, Lb, ML)
                        
                    #Tp1:
                    if ML == MLb + 1 and ms1 == ms1b and ms2 == ms2b - 1:
                        b_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b - 1))
                        W[r, c] = W[r, c] - 1.5 * np.sqrt(8*np.pi/15.0) * ms1 * b_s2 * spher_int(L, ML, 1, Lb, MLb)
                    if ML == MLb + 1 and ms1 == ms1b - 1 and ms2 == ms2b:
                        b_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b - 1))
                        W[r, c] = W[r, c] - 1.5 * np.sqrt(8*np.pi/15.0) * ms2 * b_s1 * spher_int(L, ML, 1, Lb, MLb)
                        
                    #Tm1:
                    if ML == MLb - 1 and ms1 == ms1b and ms2 == ms2b + 1:
                        a_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b + 1))
                        W[r, c] = W[r, c] + 1.5 * np.sqrt(8*np.pi/15.0) * ms1 * a_s2 * spher_int(L, ML, -1, Lb, MLb)
                    if ML == MLb - 1 and ms1 == ms1b + 1 and ms2 == ms2b:
                        a_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b + 1))
                        W[r, c] = W[r, c] + 1.5 * np.sqrt(8*np.pi/15.0) * ms2 * a_s1 * spher_int(L, ML, -1, Lb, MLb)
                        
                    #Tp2:
                    if ML == MLb + 2 and ms1 == ms1b - 1 and ms2 == ms2b - 1:
                        b_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b - 1))
                        b_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b - 1))
                        W[r, c] = W[r, c] + 0.75 * np.sqrt(32*np.pi/15.0) * b_s1 * b_s2 * spher_int(L, ML, 2, Lb, MLb)
                        
                    
                    #Tm2:
                    if ML == MLb - 2 and ms1 == ms1b + 1 and ms2 == ms2b + 1:
                        a_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b + 1))
                        a_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b + 1))
                        W[r, c] = W[r, c] + 0.75 * np.sqrt(32*np.pi/15.0) * a_s1 * a_s2 * spher_int(L, ML, -2, Lb, MLb)
                    
            c = c + 1
        r = r + 1
        
    return W
def construct_L_square(basis):
    Lsq = np.full((len(basis), len(basis)), 0.0)

    r=0
    c=0
    for bra in basis:
        L = bra[0]
        ML = bra[1]
        ms1 = bra[2]
        mi1 = bra[3]
        ms2 = bra[4]
        mi2 = bra[5]
        Lsq[r,r] = L * (L + 1)
        r = r + 1
        
    return Lsq
def construct_projection_operators(basis):
    P_S0 = np.full((len(basis), len(basis)), 0.0)
    P_S1 = np.full((len(basis), len(basis)), 0.0)
    
    r=0
    c=0
    for bra in basis:
        L = bra[0]
        ML = bra[1]
        ms1 = bra[2]
        mi1 = bra[3]
        ms2 = bra[4]
        mi2 = bra[5]
        
        c=0
        for ket in basis:
            Lb = ket[0]
            MLb = ket[1]
            ms1b = ket[2]
            mi1b = ket[3]
            ms2b = ket[4]
            mi2b = ket[5]
            
            if L == Lb and ML == MLb and mi1 == mi1b and mi2 == mi2b:
                if ms1 + ms2 == 0 and ms1b + ms2b == 0:
                    P_S0[r,c] = CG(s1, ms1, s2, ms2, 0, ms1 + ms2).doit().evalf() * CG(s1, ms1b, s2, ms2b, 0, ms1b + ms2b).doit().evalf()
                if np.abs(ms1 + ms2) <= 1  and np.abs(ms1b + ms2b) <= 1 and ms1 + ms2 == ms1b + ms2b:
                    P_S1[r,c] = CG(s1, ms1, s2, ms2, 1, ms1 + ms2).doit().evalf() * CG(s1, ms1b, s2, ms2b, 1, ms1b + ms2b).doit().evalf()
                    
            c = c + 1
        r = r + 1

    return P_S0, P_S1


def construct_zeeman_operator(basis):
    H_zee = np.full((len(basis), len(basis)), 0.0)
    
    r=0
    c=0
    for bra in basis:
        L = bra[0]
        ML = bra[1]
        ms1 = bra[2]
        mi1 = bra[3]
        ms2 = bra[4]
        mi2 = bra[5]
        
        H_zee[r,r] = mu_b * (g_e * ms1 + g_nuc_1 * mi1 + g_e * ms2 + g_nuc_2 * mi2)
        
        r = r + 1
        
    return H_zee

def construct_hyperfine_operator(basis):
    H_hyp = np.full((len(basis), len(basis)), 0.0)
    
    r=0
    c=0
    for bra in basis:
        L = bra[0]
        ML = bra[1]
        ms1 = bra[2]
        mi1 = bra[3]
        ms2 = bra[4]
        mi2 = bra[5]
        
        c=0
        for ket in basis:
            Lb = ket[0]
            MLb = ket[1]
            ms1b = ket[2]
            mi1b = ket[3]
            ms2b = ket[4]
            mi2b = ket[5]
            
            if L == Lb and ML == MLb:
                if ms2 == ms2b and mi2 == mi2b:
                    if ms1 == ms1b + 1 and mi1 == mi1b - 1:
                        a_s1 = np.sqrt(s1*(s1+1) - ms1b*(ms1b + 1))
                        b_i1 = np.sqrt(i1*(i1+1) - mi1b*(mi1b - 1))
                        H_hyp[r, c] = H_hyp[r, c] + zeta1 * 0.5 * a_s1 * b_i1
                    if ms1 == ms1b - 1 and mi1 == mi1b + 1:
                        b_s1 = np.sqrt(s1*(s1+1) - ms1b*(ms1b - 1))
                        a_i1 = np.sqrt(i1*(i1+1) - mi1b*(mi1b + 1))
                        H_hyp[r, c] = H_hyp[r, c] + zeta1 * 0.5 * b_s1 * a_i1
                    if ms1 == ms1b and mi1 == mi1b:
                        H_hyp[r, c] = H_hyp[r, c] + zeta1 * ms1 * mi1
                        
                if ms1 == ms1b and mi1 == mi1b:
                    if ms2 == ms2b + 1 and mi2 == mi2b - 1:
                        a_s2 = np.sqrt(s2*(s2+1) - ms2b*(ms2b + 1))
                        b_i2 = np.sqrt(i2*(i2+1) - mi2b*(mi2b - 1))
                        H_hyp[r, c] = H_hyp[r, c] + zeta2 * 0.5 * a_s2 * b_i2
                    if ms2 == ms2b - 1 and mi2 == mi2b + 1:
                        b_s2 = np.sqrt(s2*(s2+1) - ms2b*(ms2b - 1))
                        a_i2 = np.sqrt(i2*(i2+1) - mi2b*(mi2b + 1))
                        H_hyp[r, c] = H_hyp[r, c] + zeta2 * 0.5 * b_s2 * a_i2
                    if ms2 == ms2b and mi2 == mi2b:
                        H_hyp[r, c] = H_hyp[r, c] + zeta2 * ms2 * mi2
                    
                    
            c = c + 1
        r = r + 1
        
    return H_hyp
        
    

"""
B_field: B-Field in Tesla
M_tot: Hamiltonian terms commute with total M = M_L + m_s1 + m_i1 + m_s2 + m_i2

Returns eigenvalues of limiting hamilitonian.
"""
def generate_ham(B_field, M_tot):
    
    basis = construct_basis(M_tot)
    basis_size = len(basis)
    
    conf.num_channels = basis_size
    
    P_S0, P_S1 = construct_projection_operators(basis)
    H_dip = construct_dipole_dipole_operator(basis)
    H_zee = construct_zeeman_operator(basis)
    H_hyp = construct_hyperfine_operator(basis)
    L_sq = construct_L_square(basis)
    
    H_zee = B_field * H_zee
    

    H_zee = 1.0/conf.E_scaling * H_zee
    H_hyp = 1.0/conf.E_scaling * H_hyp
    
    H_dip = 1.0/conf.E_scaling * H_dip
    P_S0 = 1.0/conf.E_scaling * P_S0
    P_S1 = 1.0/conf.E_scaling * P_S1
    
    H_infty = H_hyp + H_zee
    eig_val, Ut = lin.eigh(H_infty)
    
    H_infty_t = multi_dot([inv(Ut), H_infty, Ut])
    P_S0_t = multi_dot([inv(Ut), P_S0, Ut])
    P_S1_t = multi_dot([inv(Ut), P_S1, Ut])
    H_dip_t = multi_dot([inv(Ut), H_dip, Ut])
    H_zee_t = multi_dot([inv(Ut), H_zee, Ut])
    H_hyp_t = multi_dot([inv(Ut), H_hyp, Ut])
    L_sq_t = multi_dot([inv(Ut), L_sq, Ut])
    
    #print("Basis size: " + str(basis_size))
    
    #print("Lowest eigenvalue: " + str(np.amin(eig_val)) + " with index: " + str(np.argmin(eig_val)))
                 
    
    l_arr = list()
    for i in range(len(basis)):
        x = int(1.0*L_sq_t[i,i] + 0.1)
        ll = (-1 + np.sqrt(1+4*x))/2
        l_arr.append(int(ll))
    f = tb.open_file(conf.working_dir + '\\in\\l_list.h5', 'w')
    f.create_array('/', 'ARR', l_arr)
    f.close()
    
    f = tb.open_file(conf.working_dir + '\\in\\H_hyp.h5', 'w')
    f.create_array('/', 'ARR', H_hyp_t)
    f.close()
    
    f = tb.open_file(conf.working_dir + '\\in\\H_zee.h5', 'w')
    f.create_array('/', 'ARR', H_zee_t)
    f.close()
    
    f = tb.open_file(conf.working_dir + '\\in\\H_S0.h5', 'w')
    f.create_array('/', 'ARR', P_S0_t)
    f.close()
    
    f = tb.open_file(conf.working_dir + '\\in\\H_S1.h5', 'w')
    f.create_array('/', 'ARR', P_S1_t)
    f.close()
    
    f = tb.open_file(conf.working_dir + '\\in\\H_dip.h5', 'w')
    f.create_array('/', 'ARR', H_dip_t)
    f.close()

    conf.internal_states_energy = eig_val

    return eig_val


























