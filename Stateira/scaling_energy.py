# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 14:03:35 2022

@author: Alex
"""

import scipy

def scaling_energy(m1, m2, rs):
    mu = m1*m2 / (m1 + m2)
    Es = scipy.constants.atomic_mass**2 / (2*mu * rs**2)
    return Es

