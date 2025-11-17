# -*- coding: utf-8 -*-

import uproot
import matplotlib.pyplot as plt
import numpy as np

# Ouvre les fichiers ROOT
Rdir = uproot.open("ggddir_onef.root")
Ronef = uproot.open("ggodir_onef.root")

fname = 'hp20;1'

# Recupere les histogrammes
hdir = Rdir[fname]
honef = Ronef[fname]

Dir = hdir.values()
Onef = honef.values()

# centre des bins en rapidite
edges = hdir.axes[0].edges()
centers = (edges[:-1] + edges[1:]) / 2.

# Calcul des fractions 

Tot = Dir+Onef
fdir  = Dir/Tot
fonef = Onef/Tot

# Prepare le plot
plt.plot(centers, fdir, label = r'$f_{dir}$')
plt.plot(centers, fonef, label = r'$f_{frag}$')

# Affiche tous les histogramme
plt.xlabel('y')
plt.ylabel(r"fractions , $p_\bot = 5 GeV$")
plt.legend(fontsize = 8)
# ~ plt.yscale('log')
plt.xlim(-6,6)
plt.ylim(0,1.1)
plt.grid(True)
plt.axhline(y=1, color='grey', linestyle='--')
plt.tight_layout()
plt.savefig('fractions.pdf')
plt.show()
