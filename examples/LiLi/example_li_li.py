# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:31:02 2022

@author: Alex
"""
import sys
sys.path.append('C:\\Users\\Alex\\Documents\\Stateira')
import Stateira as hm
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')

amu = 1.66053906660E-27
kB = 1.380649E-23
a0 = 5.29177210903E-11
mli = 6.01512279500000 * amu

hm.init_dir("C:\\Users\\Alex\\Documents\\Stateira\\examples\\LiLi")
hm.scaling_energy(mli, mli)
hm.alkali_hamiltonians_nd.generate_ham(0.08, 0)

Tcoll = 100E-09
Ekin = Tcoll * kB

hm.set_energy(Ekin, 0)

hm.scattering_config["output_file"] = "out.h5"
hm.scattering_config["pot_config_file"] = "lilipot_conf.config"
hm.scattering_config["save_wavefunction"] = "true"
hm.scattering_config["r_start"] = 1.4
hm.scattering_config["r_stop"] = 400
hm.write_conf()

hm.run_siem(True)


out = hm.analyze.scat("out.h5")
a = out.get_scattering_length(0, 0)
print("Scattering length: " + str(np.real(a/hm.conf.r_scaling)))

rarr = np.linspace(1.4,20,800)
vpsi = out.get_scat_sol_arr(rarr, 0)
plt.plot(rarr, vpsi)
plt.show()
