# -*- coding: utf-8 -*-

# =============================================================================
# File for computing R_pA of each process l
# with the R_pA of the entire direct photon processus.
# Last modified: 30/11/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
from src import sigma as sig
from src import Collision

### Initialisation of a sigma object ###
rs = 8800
# rs = 200
s = (rs)**2 # CM energy in Gev2

proton = ['NNPDF40_nlo_as_01180','MSHT20nlo_as118','CT18NLO']

Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

atom = 'Pb'
# atom = 'Au'
Atom = sig.Atom #dictionary where we stock important atomic aspect (now just Z and A) to plot and put into file names

Z = Atom[atom]['Z']
A = Atom[atom]['A']

PDF=[]
for prot in proton:
    PDF.append(sig.Sigma(prot,Pb,s,Z,A))

d = sig.Switch

#Kinematics
p_T = 5 #GeV
mu_f_2 = p_T**2											

# Plot characteristics
f_size = 17
alph = 0.5
X = np.logspace(-5,-1,100)

# Gluon pdf

Tot_G_p = []
for i,prot in enumerate(proton):
    sigma_p = PDF[i]
    p_set = sigma_p.p_set
    G_p=[]
    for num in range(p_set.size): 
        g_p = [sigma_p.Gluon_p(x,mu_f_2,num) for x in X]
        G_p.append(g_p)
    (G_p_min, G_p_max)=sig.min_max(G_p)
    Tot_G_p.append((G_p_min, G_p_max))

fig, ax = plt.subplots(constrained_layout=True)
for i,prot in enumerate(proton):
    plt.fill_between(X,Tot_G_p[i][0],Tot_G_p[i][1],alpha =alph, label= prot)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.title("gluon pdf set compared")
plt.text(0.1, 0.1, r'$Q^2 = 25$ GeV$^2$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.xlabel(r'x')
plt.ylabel(r'$F_2(x,Q^2)$')
plt.show()
plt.savefig(os.path.join(plots_dir, 'gluon_proton_pdf_set_comparison_rs'+str(rs)+'_P_T'+str(p_T)+'.pdf'))
plt.close()

# Gluon pdf

Tot_F_2 = []
for i,prot in enumerate(proton):
    sigma_p = PDF[i]
    p_set = sigma_p.p_set
    F_2=[]
    for num in range(p_set.size): 
        f_2 = [sigma_p.F2_p(x,mu_f_2,num) for x in X]
        F_2.append(f_2)
    (F_2_min, F_2_max)=sig.min_max(F_2)
    Tot_F_2.append((F_2_min, F_2_max))

fig, ax = plt.subplots(constrained_layout=True)
for i,prot in enumerate(proton):
    plt.fill_between(X,Tot_F_2[i][0],Tot_F_2[i][1],alpha =alph, label= prot)
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.title('F_2 pdf set compared')
plt.xlabel(r'x')
plt.ylabel(r'$F_2(x,Q^2)$')
plt.text(0.1, 0.1, r'$Q^2 = 25$ GeV$^2$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.show()
plt.savefig(os.path.join(plots_dir, 'F_2_proton_pdf_set_comparison_rs'+str(rs)+'_P_T'+str(p_T)+'.pdf'))
plt.close()