# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 21:38:21 2022
Generates the hamiltonian terms for non distinguishable alkali atoms.
@author: alexa
"""
import numpy as np
from sympy.physics.quantum.cg import CG
import tables as tb
from .Stateira_config import Stateira_config as conf

h_planck = 6.62606896E-34
g_e = 2.0023193043622
mu_b = 9.27400915E-24

"""
Needed configuration data
Default Li - Li Values:
"""

s1 = 0.5
i1 = 1.0
s2 = 0.5
i2 = 1.0

Lmax = 0            #Maximum L included in Basis
g_nuc_1 = -0.0004476540 #Nuclear g factor Atom 1
g_nuc_2 = -0.0004476540 #Nuclear g factor Atom 2

zeta1 = 152.1368407E06 * h_planck#* (i1 + 0.5) 
zeta2 = 152.1368407E06 * h_planck#* (i2 + 0.5)

Lparity = 'even'    #Parity of L. 'even': L = 0, 2, ... , Lmax ; 'odd': L = 1, 3, ... , Lmax 
boson_or_fermion = 'fermion'        #Particle type 'fermion' or 'boson'

def construct_symmetrized_basis(Mk):
    basic_basis=list()
    #vec = np.full((6),0.0)
    if Lparity == 'even':
        Lrange = np.arange(0, Lmax+1, 2)
    else:
        Lrange = np.arange(1, Lmax+1, 2)
        
    if boson_or_fermion == 'boson':
        epsilon = 1
    else:
        epsilon = -1
        
    for L in Lrange:
        for ML in np.arange(L, -L-1,-1):
            for ms1 in np.arange(s1, -s1-1,-1):
                for mi1 in np.arange(i1, -i1-1,-1):
                    for ms2 in np.arange(s2, -s2-1,-1):
                        for mi2 in np.arange(i2, -i2-1,-1):
                            
                            if ML + ms1 + mi1 + ms2 + mi2 == Mk:
                                vec = [L, ML, ms1, mi1, ms2, mi2]
                                basic_basis.append(vec)
                                
    sym_basis = list()
    
    for a_ket in basic_basis:
        L = a_ket[0]
        ML = a_ket[1]
        ms1 = a_ket[2]
        mi1 = a_ket[3]
        ms2 = a_ket[4]
        mi2 = a_ket[5]
        
        if ms1 == ms2 and mi1 == mi2:
            alpha = 0.5
        else:
            alpha = 1.0/np.sqrt(2)

        beta = epsilon * (-1.0)**L
        
        if ms1 == ms2 and mi1 == mi2:
            if beta == 1.0:
                vec = [alpha, beta, L, ML, ms1, mi1, ms2, mi2]
                vec2 = [alpha, beta, L, ML, ms2, mi2, ms1, mi1]
                if vec not in sym_basis and vec2 not in sym_basis:
                    sym_basis.append(vec)
        else:
            vec = [alpha, beta, L, ML, ms1, mi1, ms2, mi2]
            vec2 = [alpha, beta, L, ML, ms2, mi2, ms1, mi1]
            if vec not in sym_basis and vec2 not in sym_basis:
                sym_basis.append(vec)
              
    print(len(basic_basis))
    print(len(sym_basis))
    return  np.asarray(sym_basis)


def print_symmetrized_basis(sym_basis):
    for ele in sym_basis:
        alpha = ele[0]
        beta = ele[1]
        L = ele[2]
        ML = ele[3]
        ms1 = ele[4]
        mi1 = ele[5]
        ms2 = ele[6]
        mi2 = ele[7]
        
        st = "|" + str(ms1) + "," + str(mi1) + ";" + str(ms2) + "," + str(mi2) + ";" + str(L) + "," + str(ML) + ">"
        
        #st = st + " + (" + str(beta) + ")|" + str(ms2) + "," + str(mi2) + ";" + str(ms1) + "," + str(mi1) + ";" + str(L) + "," + str(ML) + ">"
        #st = str(2*ms1) + "  " + str(2*mi1) + "  " + str(2*ms2) + "  " + str(2*mi2)
        #st = "[" + str(2*ms1) + "," + str(2*mi1) + "," + str(2*ms2) + "," + str(2*mi2) + "," + str(L) + "],"
        
        print(st)
        
def spher_int(L1, ML1, q, L3, ML3):
    res = (-1.0)**ML1 * (-1.0)**ML3 * np.sqrt(5*(2*L1+1) / (4*np.pi * (2*L3 + 1)))
    res = res * CG(L1, 0, 2, 0, L3, 0).doit().evalf() * CG(L1, -1.0*ML1, 2, q, L3, -1.0*ML3).doit().evalf()
    return res
        
def construct_dipole_dipole_element(bra, ket):
    #[alpha, beta, L, ML, ms1, mi1, ms2, mi2]
    W = 0
    
    L = bra[0]
    ML = bra[1]
    ms1 = bra[2]
    mi1 = bra[3]
    ms2 = bra[4]
    mi2 = bra[5]
    
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
                W = W + np.sqrt(16*np.pi/5.0) * ms1 * ms2 * spher_int(L, ML, 0, Lb, ML)
                         
            #T0b:
    
            if ML == MLb and ms1 == ms1b + 1 and ms2 == ms2b - 1:
                a_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b + 1))
                b_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b - 1))
                W = W -0.25 * np.sqrt(16*np.pi/5.0) * a_s1 * b_s2 * spher_int(L, ML, 0, Lb, ML)
            if ML == MLb and ms1 == ms1b - 1 and ms2 == ms2b + 1:
                b_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b - 1))
                a_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b + 1))
                W = W -0.25 * np.sqrt(16*np.pi/5.0) * b_s1 * a_s2 * spher_int(L, ML, 0, Lb, ML)
                        
            #Tp1:
            if ML == MLb + 1 and ms1 == ms1b and ms2 == ms2b - 1:
                b_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b - 1))
                W = W - 1.5 * np.sqrt(8*np.pi/15.0) * ms1 * b_s2 * spher_int(L, ML, 1, Lb, MLb)
            if ML == MLb + 1 and ms1 == ms1b - 1 and ms2 == ms2b:
                b_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b - 1))
                W = W - 1.5 * np.sqrt(8*np.pi/15.0) * ms2 * b_s1 * spher_int(L, ML, 1, Lb, MLb)
                        
            #Tm1:
            if ML == MLb - 1 and ms1 == ms1b and ms2 == ms2b + 1:
                a_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b + 1))
                W = W + 1.5 * np.sqrt(8*np.pi/15.0) * ms1 * a_s2 * spher_int(L, ML, -1, Lb, MLb)
            if ML == MLb - 1 and ms1 == ms1b + 1 and ms2 == ms2b:
                a_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b + 1))
                W = W + 1.5 * np.sqrt(8*np.pi/15.0) * ms2 * a_s1 * spher_int(L, ML, -1, Lb, MLb)
                        
            #Tp2:
            if ML == MLb + 2 and ms1 == ms1b - 1 and ms2 == ms2b - 1:
                b_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b - 1))
                b_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b - 1))
                W = W + 0.75 * np.sqrt(32*np.pi/15.0) * b_s1 * b_s2 * spher_int(L, ML, 2, Lb, MLb)
                        
                    
            #Tm2:
            if ML == MLb - 2 and ms1 == ms1b + 1 and ms2 == ms2b + 1:
                a_s1 = np.sqrt(s1*(s1+1) - ms1b * (ms1b + 1))
                a_s2 = np.sqrt(s2*(s2+1) - ms2b * (ms2b + 1))
                W = W + 0.75 * np.sqrt(32*np.pi/15.0) * a_s1 * a_s2 * spher_int(L, ML, -2, Lb, MLb)
                    
        
    return W

def construct_P_S0_element(bra, ket):
    P_S0 = 0
    
    L = bra[0]
    ML = bra[1]
    ms1 = bra[2]
    mi1 = bra[3]
    ms2 = bra[4]
    mi2 = bra[5]
    
    Lb = ket[0]
    MLb = ket[1]
    ms1b = ket[2]
    mi1b = ket[3]
    ms2b = ket[4]
    mi2b = ket[5]
    

    if L == Lb and ML == MLb and mi1 == mi1b and mi2 == mi2b:
        if ms1 + ms2 == 0 and ms1b + ms2b == 0:
            P_S0 = CG(s1, ms1, s2, ms2, 0, ms1 + ms2).doit().evalf() * CG(s1, ms1b, s2, ms2b, 0, ms1b + ms2b).doit().evalf()
                    

    return P_S0

def construct_P_S1_element(bra, ket):
    P_S1 = 0
    
    L = bra[0]
    ML = bra[1]
    ms1 = bra[2]
    mi1 = bra[3]
    ms2 = bra[4]
    mi2 = bra[5]
    
    Lb = ket[0]
    MLb = ket[1]
    ms1b = ket[2]
    mi1b = ket[3]
    ms2b = ket[4]
    mi2b = ket[5]
    

    if L == Lb and ML == MLb and mi1 == mi1b and mi2 == mi2b:
        if np.abs(ms1 + ms2) <= 1  and np.abs(ms1b + ms2b) <= 1 and ms1 + ms2 == ms1b + ms2b:
            P_S1 = CG(s1, ms1, s2, ms2, 1, ms1 + ms2).doit().evalf() * CG(s1, ms1b, s2, ms2b, 1, ms1b + ms2b).doit().evalf()
                    

    return P_S1


def construct_zeeman_element(bra, ket):
    H_zee = 0
    
    L = bra[0]
    ML = bra[1]
    ms1 = bra[2]
    mi1 = bra[3]
    ms2 = bra[4]
    mi2 = bra[5]
    
    Lb = ket[0]
    MLb = ket[1]
    ms1b = ket[2]
    mi1b = ket[3]
    ms2b = ket[4]
    mi2b = ket[5]
        
    if bra == ket:
        H_zee = mu_b * (g_e * ms1 + g_nuc_1 * mi1 + g_e * ms2 + g_nuc_2 * mi2)
        #H_zee = mu_b * (g_e * ms1  + g_e * ms2 )
        
    return H_zee

def construct_L_square(bra, ket):
    Lsq = 0
    L = bra[0]
    ML = bra[1]
    ms1 = bra[2]
    mi1 = bra[3]
    ms2 = bra[4]
    mi2 = bra[5]
    
    Lb = ket[0]
    MLb = ket[1]
    ms1b = ket[2]
    mi1b = ket[3]
    ms2b = ket[4]
    mi2b = ket[5]
    
    if bra == ket:
        Lsq = L * (L + 1)
        
    return Lsq

def construct_hyperfine_element(bra, ket):
    H_hyp = 0
    
    L = bra[0]
    ML = bra[1]
    ms1 = bra[2]
    mi1 = bra[3]
    ms2 = bra[4]
    mi2 = bra[5]
    
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
                H_hyp = H_hyp + zeta1 * 0.5 * a_s1 * b_i1
            if ms1 == ms1b - 1 and mi1 == mi1b + 1:
                b_s1 = np.sqrt(s1*(s1+1) - ms1b*(ms1b - 1))
                a_i1 = np.sqrt(i1*(i1+1) - mi1b*(mi1b + 1))
                H_hyp = H_hyp + zeta1 * 0.5 * b_s1 * a_i1
            if ms1 == ms1b and mi1 == mi1b:
                H_hyp = H_hyp + zeta1 * ms1 * mi1
                        
        if ms1 == ms1b and mi1 == mi1b:
            if ms2 == ms2b + 1 and mi2 == mi2b - 1:
                a_s2 = np.sqrt(s2*(s2+1) - ms2b*(ms2b + 1))
                b_i2 = np.sqrt(i2*(i2+1) - mi2b*(mi2b - 1))
                H_hyp = H_hyp + zeta2 * 0.5 * a_s2 * b_i2
            if ms2 == ms2b - 1 and mi2 == mi2b + 1:
                b_s2 = np.sqrt(s2*(s2+1) - ms2b*(ms2b - 1))
                a_i2 = np.sqrt(i2*(i2+1) - mi2b*(mi2b + 1))
                H_hyp = H_hyp + zeta2 * 0.5 * b_s2 * a_i2
            if ms2 == ms2b and mi2 == mi2b:
                H_hyp = H_hyp + zeta2 * ms2 * mi2
                    
        
    return H_hyp

#[alpha, beta, L, ML, ms1, mi1, ms2, mi2]

def construct_operator(basis, constructor_function):
    O = np.full((len(basis),len(basis)), 0.0)
    
    for i in range(len(basis)):
        for j in range(len(basis)):
            alpha_i = basis[i,0]
            alpha_j = basis[j,0]
            beta_i = basis[i,1]
            beta_j = basis[j,1]
            
            bra = basis[i, 2:].copy()
            ket = basis[j, 2:].copy()
            L = bra[0]
            ML = bra[1]
            ms1 = bra[2]
            mi1 = bra[3]
            ms2 = bra[4]
            mi2 = bra[5]
            
            Lb = ket[0]
            MLb = ket[1]
            ms1b = ket[2]
            mi1b = ket[3]
            ms2b = ket[4]
            mi2b = ket[5]
            
            bra2 = [L, ML, ms1, mi1, ms2, mi2]
            ket2 = [Lb, MLb, ms1b, mi1b, ms2b, mi2b]
            
            O[i,j] = alpha_i * alpha_j * constructor_function(bra2, ket2)
            
            ket2 = [Lb, MLb, ms2b, mi2b, ms1b, mi1b]
            O[i,j] = O[i,j] + alpha_i * alpha_j * beta_j * constructor_function(bra2, ket2)
    
            bra2 = [L, ML, ms2, mi2, ms1, mi1]
            ket2 = [Lb, MLb, ms1b, mi1b, ms2b, mi2b]
            O[i,j] = O[i,j] + alpha_i * alpha_j * beta_i * constructor_function(bra2, ket2)
            
            ket2 = [Lb, MLb, ms2b, mi2b, ms1b, mi1b]
            O[i,j] = O[i,j] + alpha_i * alpha_j * beta_i * beta_j * constructor_function(bra2, ket2)
            
    return O
        
def generate_ham(B_field, M_tot):
    sym_bas = construct_symmetrized_basis(M_tot)
    basis_size = len(sym_bas)
    
    conf.num_channels = basis_size
    
    
    P_S0 = construct_operator(sym_bas, construct_P_S0_element)
    P_S1 = construct_operator(sym_bas, construct_P_S1_element)
    H_dip = construct_operator(sym_bas, construct_dipole_dipole_element)
    H_zee = construct_operator(sym_bas, construct_zeeman_element)
    H_hyp = construct_operator(sym_bas, construct_hyperfine_element)
    L_sq = construct_operator(sym_bas, construct_L_square)

    
    H_zee = B_field * H_zee
        
    H_zee = 1.0/conf.E_scaling * H_zee
    H_hyp = 1.0/conf.E_scaling * H_hyp
    
    H_dip = 1.0/conf.E_scaling * H_dip
    P_S0 = 1.0/conf.E_scaling * P_S0
    P_S1 = 1.0/conf.E_scaling * P_S1
    
    H_infty = H_hyp + H_zee
    eig_val, Ut = np.linalg.eigh(H_infty)
    

    L_sq_t = np.linalg.multi_dot([np.linalg.inv(Ut), L_sq, Ut])
    P_S0_t = np.linalg.multi_dot([np.linalg.inv(Ut), P_S0, Ut])
    P_S1_t = np.linalg.multi_dot([np.linalg.inv(Ut), P_S1, Ut])
    H_dip_t = np.linalg.multi_dot([np.linalg.inv(Ut), H_dip, Ut])
    H_zee_t = np.linalg.multi_dot([np.linalg.inv(Ut), H_zee, Ut])
    H_hyp_t = np.linalg.multi_dot([np.linalg.inv(Ut), H_hyp, Ut])
    

    l_arr = list()
    for i in range(len(sym_bas)):
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


#generate_ham(0.1)





















































        