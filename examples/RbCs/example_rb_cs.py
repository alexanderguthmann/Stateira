# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 14:12:28 2022

@author: Alex
"""
import sys
sys.path.append('C:\\Users\\Alex\\Documents\\Stateira')
import Stateira as hm

amu = 1.66053906660E-27
kB = 1.380649E-23
mrb = 86.90918052 * amu
mcs = 132.905451931 * amu

hm.init_dir("C:\\Users\\Alex\\Documents\\Stateira\\examples\\RbCs")
hm.scaling_energy(mrb, mcs)
hm.alkali_hamiltonians_d.generate_ham(0.001, 2)

Tcoll = 100E-09
Ekin = Tcoll * kB

hm.set_energy(Ekin, 0)

hm.scattering_config["output_file"] = "out.h5"
hm.scattering_config["pot_config_file"] = "rbcspot_conf.config"
hm.scattering_config["save_wavefunction"] = "true"
hm.scattering_config["r_start"] = 2.5
hm.scattering_config["r_stop"] = 400
hm.write_conf()

hm.run_siem(True)
