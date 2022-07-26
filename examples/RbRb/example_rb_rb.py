# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:58:13 2022

@author: Alex
"""

import Stateira as hm

amu = 1.66053906660E-27
kB = 1.380649E-23
h_planck = 6.62606896E-34
mrb = 86.90918052 * amu

hm.init_dir("C:\\Users\\Alex\\Documents\\Stateira\\examples\\RbRb")
hm.scaling_energy(mrb, mrb)

hm.alkali_hamiltonians_nd.boson_or_fermion = 'boson'

hm.alkali_hamiltonians_nd.g_nuc_1 = -0.000995141
hm.alkali_hamiltonians_nd.g_nuc_2 = -0.000995141
hm.alkali_hamiltonians_nd.zeta1 = 3.417341305452145E09 * h_planck
hm.alkali_hamiltonians_nd.zeta2 = 3.417341305452145E09 * h_planck
hm.alkali_hamiltonians_nd.i1 = 1.5
hm.alkali_hamiltonians_nd.i2 = 1.5

hm.alkali_hamiltonians_nd.generate_ham(0.001, 0)

Tcoll = 100E-09
Ekin = Tcoll * kB

hm.set_energy(Ekin, 0)

hm.scattering_config["output_file"] = "out.h5"
hm.scattering_config["pot_config_file"] = "rb2pot.config"
hm.scattering_config["save_wavefunction"] = "true"
hm.scattering_config["r_start"] = 2.5
hm.scattering_config["r_stop"] = 400
hm.write_conf()

hm.run_siem(True)
