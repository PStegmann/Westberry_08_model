# file: plot.py
# Description: plots the output of wb08.
#
# Copyright © 2018 Patrick Stegmann
#
# This file is part of Westberry_08_model.
#

import numpy as np
import matplotlib.pyplot as plt

phy = np.loadtxt("../build/output.txt", skiprows=1)

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 12,
        }

fig = plt.figure()
plt.semilogx(phy[:,1],-1.*phy[:,0])
plt.xlabel(r"Phytoplankton Carbon [mg $C \cdot m^{-3}$]")
plt.ylabel("Ocean Depth [m]")

plt.text(2, -45, "Mixed Layer", fontdict=font)
plt.text(4, -130, "Light-limited regime", fontdict=font)
plt.show()
