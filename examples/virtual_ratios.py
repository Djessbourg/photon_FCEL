# -*- coding: utf-8 -*-

# =============================================================================
# Study of the impact of adding a mass to the photon on the cross section
# Last modified: 12/12/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, ScalarFormatter
from src import sigma as sig
f_size = 15
a = 0
b= 0

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)
	
### Initialisation of a sigma object ###
rs = 8800
s = (rs)**2 # CM energy in Gev2
p = "nNNPDF30_nlo_as_0118_p"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"
Z = 82.
A = 208.
M = [5,10,15] # vitual photon mass in GeV
cen = 0

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

d = sig.Switch
convention = 'd2p_t'													# Define a cross section convention btw dp_t, dp_t2 ans d2p_t
col_type = 'pp'
Collision = sig.Collision

### y dependant ###
p_T = 10																# GeV
x_T = (2.*p_T)/rs 

for m in M:
    b = (m/p_T)**2
    Y = sig.Y_list_M(x_T,b)
    print((min(Y),max(Y)))
    R_qg, R_gq, R_qqbar, R_qbarq = [],[],[],[]
    for y in Y:
        R = pPb_cross_section.R_composition_dydpt_M(y,x_T,m,cen,is_pp=Collision(col_type),switch=convention)
        R_qg.append(R[0][0]);R_gq.append(R[1][0]);R_qqbar.append(R[2][0]);R_qbarq.append(R[3][0])
    plt.plot(Y,R_qg,label=r'qG/tot')
    plt.plot(Y,R_gq,label=r'Gq/tot')
    plt.plot(Y,R_qqbar,label=r'qqbar/tot')
    plt.plot(Y,R_qbarq,label=r'qbarq/tot')
    plt.axhline(y=1, color='grey', alpha=0.3)
    plot_usuals(n=2)
    plt.xlabel(r'$y$',fontsize=f_size)
    plt.ylabel(r'$R_i = \sigma_i/\sigma_{tot}$')
    plt.title(r'$\gamma^\star$ production for $M=$'+str(m)+r' GeV and $p_\perp = $'+str(p_T)+' GeV')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_component_virtual_ratios_'+str(rs)+'GeV_'+convention+'_y'+str(y)+'_M'+str(m)+'GeV'+'.pdf'))
    plt.show()
    plt.close()