# -*- coding: utf-8 -*-

# =============================================================================
# Ploting R_pA without isosspin effect and comparing different RpA from
# different sets of PDF
# Last modified: 9/02/2026
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

proton = ['NNPDF40_nlo_as_01180'] #  ,'MSHT20nlo_as118','CT18NLO' Here put whatever PDF_set you want to compare 

Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

atom = 'Pb'
# atom = 'Au'
Atom = sig.Atom #dictionary where we stock important atomic aspect (now just Z and A) to plot and put into file names

Z = Atom[atom]['Z']
A = Atom[atom]['A']

#Kinematics
p_T = 5 #GeV
x_T = 2*p_T/rs
err=''

PDF=[]
RpA4 =[]
for i,prot in enumerate(proton):
    PDF.append(sig.Sigma(prot,Pb,s,Z,A))
    RpA4.append(PDF[i].Uncertainties_4R_iso_dy(x_T,var_err=err))

Y = sig.Y_list(x_T)

for R4 in RpA4:
    RpA, Rpp_iso, R_pA_iso, R_pA_wo_iso = R4
    plt.plot(Y,RpA[0],label='RpA clasic')
    plt.plot(Y,Rpp_iso[0], label='isospin effect')
    plt.plot(Y,R_pA_iso[0],label='iso+FCEL/g')
    plt.plot(Y,R_pA_wo_iso[0],label='FCEL iso corrected')
    plt.show()
d = sig.Switch

									

# Plot characteristics
# f_size = 17
# alph = 0.5
# X = np.logspace(-5,-1,100)
# def delta_ref(ref,data):
#     ref = np.array(ref)
#     data = np.array(data)
#     return (data)/ref

# colors = ['blue','orange','green']
