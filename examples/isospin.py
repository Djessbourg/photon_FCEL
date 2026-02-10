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
RpA_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_dir')) # the direcory to save data from this file


import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
from src import sigma as sig
from src import Collision

### Initialisation of a sigma object ###
rs = 8800
# rs = 200
s = (rs)**2 # CM energy in Gev2

proton = ['NNPDF40_nlo_as_01180','MSHT20nlo_as118','CT18NLO' ] #  Here put whatever PDF_set you want to compare 

Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

atom = 'Pb'
# atom = 'Au'
Atom = sig.Atom #dictionary where we stock important atomic aspect (now just Z and A) to plot and put into file names

Z = Atom[atom]['Z']
A = Atom[atom]['A']

# Kinematics
p_T = 5 #GeV
x_T = 2*p_T/rs
err='' # works only for no errors 

# plot things

f_size = 15
a = 2
b= 3

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2) 

# Ploting 

PDF=[]
RpA4 =[]

for i,prot in enumerate(proton):
    f_name = prot+'_RpA4_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+str(p_T)+'GeV.txt'
    PDF.append(sig.Sigma(prot,Pb,s,Z,A))
    if os.path.exists(os.path.join(RpA_dir,f_name)):
        print(f" '{f_name}' already exists: loading file")
        RpA, Rpp_iso, R_pA_iso, R_pA_wo_iso = np.loadtxt(os.path.join(RpA_dir,f_name))
        RpA4.append([RpA, Rpp_iso, R_pA_iso, R_pA_wo_iso])
    else:
        print(f" '{f_name}' does not exists")
        R4 = PDF[i].Uncertainties_4R_iso_dy(x_T,var_err=err)
        r = [R4[0][0], R4[1][0], R4[2][0], R4[3][0]]
        RpA4.append(r)
        np.savetxt(os.path.join(RpA_dir,f_name),r)
        print(f" '{f_name}' has been created")

Y = sig.Y_list(x_T)
# Comparing for each pdf set the 4 ratios in the paper
# for i,prot in enumerate(proton):
#     fig, ax = plt.subplots(constrained_layout=True)
#     ax.xaxis.set_minor_locator(MultipleLocator(1))
#     ax.xaxis.set_major_locator(MultipleLocator(2))
#     RpA, Rpp_iso, R_pA_iso, R_pA_wo_iso = RpA4[i]
#     plt.plot(Y,RpA,label='RpA clasic')
#     plt.plot(Y,Rpp_iso, label='iso')
#     plt.plot(Y,R_pA_iso,label='iso+FCEL')
#     plt.plot(Y,R_pA_wo_iso,label='FCEL w/o iso')
#     plt.axhline(y=1, color='grey', alpha=0.3)
#     plot_usuals(n=2,loca='upper center',s1=f_size-a)
#     plt.xlabel(r'$y$',fontsize=f_size)
#     plt.ylim(0.8,1.1)
#     plt.tight_layout()
#     plt.text(0.6, 0.1, prot , horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-a)
#     plt.text(0.75, 0.2,r'$p_\bot =$ '+str(p_T)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-a) 
#     plt.text(0.75, 0.3,r'$\sqrt{s} =$ '+str(rs/1000)+' TeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-a)
#     ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
#     plt.savefig(os.path.join(plots_dir, prot+'_Riso_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_p_T'+str(p_T)+'GeV_M'+'.pdf'))
#     plt.show()
#     plt.close()

# Comparing the RpA without isospin effect for the different pdf set on the same figure	
# 						
fig, ax = plt.subplots(constrained_layout=True)
for i,prot in enumerate(proton):
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    RpA, Rpp_iso, R_pA_iso, R_pA_wo_iso = RpA4[i]
    plt.plot(Y,R_pA_wo_iso,label=prot)
plt.axhline(y=1, color='grey', alpha=0.3)
plot_usuals(n=1,loca='lower left',s1=f_size-a)
plt.xlabel(r'$y$',fontsize=f_size)
plt.ylabel(r'$R_{\text{pA}}$ without isospin effect',fontsize=f_size)
plt.ylim(0.8,1.1)
plt.text(0.75, 0.9,r'$p_\bot =$ '+str(p_T)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-a) 
plt.text(0.75, 0.8,r'$\sqrt{s} =$ '+str(rs/1000)+' TeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-a)
plt.tight_layout()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.savefig(os.path.join(plots_dir, 'pdf_comp_R_without_iso_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_p_T'+str(p_T)+'GeV_M'+'.pdf'))
plt.show()
plt.close()

fig, ax = plt.subplots(constrained_layout=True)
R_ref = RpA4[0][3]
for i,prot in enumerate(proton):
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(2))
    RpA, Rpp_iso, R_pA_iso, R_pA_wo_iso = RpA4[i]
    plt.plot(Y,R_pA_wo_iso/R_ref,label=prot)
plt.axhline(y=1, color='grey', alpha=0.3)
plot_usuals(n=1,loca='lower right',s1=f_size-a)
plt.xlabel(r'$y$',fontsize=f_size)
plt.ylabel(r'$R_{\text{pA}}/R_{[ref],\text{pA}}$ without isospin effect',fontsize=f_size)
plt.text(0.15, 0.9,r'$p_\bot =$ '+str(p_T)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-a) 
plt.text(0.15, 0.8,r'$\sqrt{s} =$ '+str(rs/1000)+' TeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-a)
plt.tight_layout()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.savefig(os.path.join(plots_dir, 'pdf_relative_comp_R_without_iso_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_p_T'+str(p_T)+'GeV_M'+'.pdf'))
plt.show()
plt.close()