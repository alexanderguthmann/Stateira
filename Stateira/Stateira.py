# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 14:12:01 2022

@author: Alex
"""
import scipy
import os
from .Stateira_config import Stateira_config as conf
#import Stateira_config as conf

"""
Run:
setx HDF5_DISABLE_VERSION_CHECK 2
"""

hbar = 1.0545718175E-34

"""
Sets and initializes the working directory
"""
def init_dir(directory):
    cmd_string = 'start /wait cmd /c "mkdir ' + str(directory) + '\\in"'
    os.system(cmd_string)
    conf.working_dir = directory
    
"""
Calculates and sets scaling energy
"""
def scaling_energy(m1, m2):
    mu = m1*m2 / (m1 + m2)
    conf.E_scaling = hbar**2 / (2*mu * conf.r_scaling**2)
    return conf.E_scaling

"""
Sets the energy for which the solution is calculated.
Ekin: Kinetic energy [Joules]
idx: Index of the eigenstate of the asymptotic hamiltonian. This sets the internal energy.
"""
def set_energy(Ekin, idx):
    conf.E_tot = Ekin/conf.E_scaling + conf.internal_states_energy[idx]
    
    

scattering_config = {
"num_threads":1,
"debug_mode":"false",
"num_chan":20,
"r_start":2.5,
"r_stop":400.0,
"stats_step":1,
"num_cheb":16,
"tolerance":1.0E-10,
"force_tolerance":"true",
"max_loops":10,
"e_in":0.0,
"output_file":"out/out.h5",
"pot_config_file":"lilipot_conf.config",
"v_infty_r":1.0E10,
"save_wavefunction":"false",
"max_step":10.0,
"boundary_conditions":0
}

bound_config = {
"type":"single",
"e_range_low":-1490,
"e_range_high":-1480,
"num_nodes":1,
"debug_mode":"false",
"num_chan":5,
"r_start":2.0,
"r_stop":3.6,
"r_mid":2.5,
"stats_step":1,
"num_cheb":16,
"tolerance":1.0E-7,
"energy_tolerance":1.0E-05,
"node_count_points":100,
"trace_points":1000,
"force_tolerance":"true",
"max_loops":30,
"e_in":-4961.150655306845,
"output_file":"examples\\LiLi_bound\\out\\out_s.h5",
"pot_config_file":"examples\\LiLi_bound\\lilipot_conf.config",
"v_infty_r":1.0E10,
"save_wavefunction":"true",
"max_step":10.0    
}

"""
Writes siem config file
"""
def write_conf():
    scattering_config["e_in"] = conf.E_tot
    scattering_config["num_chan"] = conf.num_channels
    with open(str(conf.working_dir) + '\\siem.config','w') as fout:
        for key in scattering_config:
            fout.write(key + "=" + str(scattering_config[key]) + "\n")
            
def run_siem(show_output):
    if show_output == True:
        os.system('start /wait cmd /c "set HDF5_DISABLE_VERSION_CHECK=2 && cd ' + conf.working_dir + ' && ' + conf.sp + 'SIEM.exe siem.config"')
    else:
        os.system('start /wait cmd /c "set HDF5_DISABLE_VERSION_CHECK=2 && cd ' + conf.working_dir + ' && ' + conf.sp + 'SIEM.exe siem.config"')

    
def write_conf_bound():
    bound_config["e_in"] = conf.E_tot
    bound_config["num_chan"] = conf.num_channels
    with open(str(conf.working_dir) + '\\siem_bound.config','w') as fout:
        for key in bound_config:
            fout.write(key + "=" + str(bound_config[key]) + "\n")
            
def run_siem_bound(show_output):
    if show_output == True:
        os.system('start /wait cmd /c "set HDF5_DISABLE_VERSION_CHECK=2 && cd ' + conf.working_dir + ' && ' + conf.sp + 'SIEM_BOUND.exe siem_bound.config"')
    else:
        os.system('start /wait cmd /c "set HDF5_DISABLE_VERSION_CHECK=2 && cd ' + conf.working_dir + ' && ' + conf.sp + 'SIEM_BOUND.exe siem_bound.config"')
