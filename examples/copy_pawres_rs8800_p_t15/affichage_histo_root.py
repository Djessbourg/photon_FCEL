# -*- coding: utf-8 -*-

import uproot
import matplotlib.pyplot as plt

# Ouvre le fichier ROOT
file = uproot.open("LO_frag_ratios.root")

# Recupere tous les histogrammes (par exemple, TH1F)
histos = {k: file[k] for k in file.keys() if file[k].classname.startswith("TH1")}

# Prepare le plot
plt.figure(figsize=(10, 6))

# Affiche tous les histogrammes
for name, hist in histos.items():
    values = hist.values()
    edges = hist.axes[0].edges()  # les bords des bins

    # Calcule les centres des bins pour le plot
    centers = (edges[:-1] + edges[1:]) / 2.

    plt.plot(centers, values, label=name)  # Ajoute chaque histo avec son nom

#plt.xlabel(r'$p_\bot$ (GeV)')
plt.xlabel(r'y')
plt.ylabel(r"Ratios, $p_\bot = 5 GeV$")
plt.legend(fontsize = 8)
# ~ plt.yscale('log')
#plt.xlim(2.5,15.5)
plt.grid(True)
plt.axhline(y=1, color='grey', linestyle='--')
plt.tight_layout()
plt.savefig('LO_frag_ratios.pdf')
plt.show()
