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
proton = "NNPDF40_lo_as_01180"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

atom = 'Pb'
# atom = 'Au'

Atom = sig.Atom

Z = Atom[atom]['Z']
A = Atom[atom]['A']

pPb_cross_section = sig.Sigma(proton,Pb,s,Z,A)

d = sig.Switch

### first atempts to reproduce last results ###
p_T = 5 #GeV
# p_T = 2
x_T = (2.0*p_T)/rs
num = 0																			# The central set 
Y = sig.Y_list(x_T)														

# A plot for sigma tot and its components
f_size = 17
alph = 0.3
q0 =0.07
convention = 'dp_t2'
col_type = 'pp'
err ='q0,mu'

RpA_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_dir')) # the direcory to save data from this file
a=0		
b=2											

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)

# before RpA, cross section of each process and ratios of thoses
fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
(s1,er1),(s2,er2),(s3,er3),(s4,er4)=pPb_cross_section.dsimga_all_dy(x_T,num,is_pp=Collision(col_type),switch=convention,l=True)
stot = s1+s2+s3+s4

plt.plot(Y,s1,label=r'$\ell=1$')
plt.plot(Y,s2,label=r'$\ell=2$')
plt.plot(Y,s3,label=r'$\ell=3$')
plt.plot(Y,s4,label=r'$\ell=4$')

plt.ylabel(r'd$\sigma_{\ell}/$d$y$'+sig.Switch[convention][0]+' '+sig.Switch[convention][1],fontsize= f_size-a)
plt.yscale('log')
plt.ylim(0.01*min(s3),10*max(stot))
plt.xlabel('y',fontsize=f_size-a)
plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'best',n=2)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.text(0.5, 0.95, proton + " p-"+atom, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.text(0.5, 0.40, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.text(0.5, 0.50, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_all'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.show()
plt.close()

# fractions
fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
plt.plot(Y,s1/stot,label=r'$\ell=1$')
plt.plot(Y,s2/stot,label=r'$\ell=2$')
plt.plot(Y,s3/stot,label=r'$\ell=3$')
plt.plot(Y,s4/stot,label=r'$\ell=4$')
plt.axhline(y=1,color= 'grey',alpha=0.3)

plt.ylabel(r'd$f_{\ell}$',fontsize= f_size-a)
plt.ylim(0,1.1)
plt.xlabel('y',fontsize=f_size-a)
plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'best',n=2)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.text(0.5, 0.95, proton + " p-"+atom, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.text(0.5, 0.40, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.text(0.5, 0.50, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'f_all'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.show()
plt.close()


# # for the tot R_pA dir
# f_name = 'RpA_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.txt'
# if os.path.exists(os.path.join(RpA_dir,f_name)):
# 	print(f"The file '{f_name}' already exists. It is loaded.")
# 	Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
# else:
# 	print("The file does not exists")
# 	Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dy(x_T,switch = convention,var_err= err)
# 	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
# 	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
# 	r = [Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf]
# 	np.savetxt(os.path.join(RpA_dir,f_name), r)
# 	print(f" '{f_name}' has been created")

# # for the R_pA of each process
# f_name = 'R_all_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.txt'
# if os.path.exists(os.path.join(RpA_dir,f_name)):
# 	print(f"The file '{f_name}' already exists. It is loaded.")
# 	R_1,R_2,R_3,R_4= np.loadtxt(os.path.join(RpA_dir,f_name))
# else:
# 	print("The file does not exists")
# 	R_1,R_2,R_3,R_4 = pPb_cross_section.R_pp_All_dy(x_T,num,q0,switch = convention)
# 	r = [R_1,R_2,R_3,R_4]
# 	np.savetxt(os.path.join(RpA_dir,f_name), r)
# 	print(f" '{f_name}' has been created")


# fig, ax = plt.subplots(constrained_layout=True)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(2))
# plt.axhline(y=1, color='grey',alpha=alph)
# plt.plot(Y,R_1,label=r'$\ell=1$')
# plt.plot(Y,R_2,label=r'$\ell=2$')
# plt.plot(Y,R_3,label=r'$\ell=3$')
# plt.plot(Y,R_4,label=r'$\ell=4$')
# plt.plot(Y,Rpa,color='blue',linestyle= '--',label=r'$\ell=$ Sum')
# #plt.grid()
# # plt.ylim(bottom= 0.9,top = 1.05)
# plt.ylabel(r'$R_{\ell}=\sigma_{\ell}^\text{pA}/A\sigma_{\ell}^\text{pp}$',fontsize= f_size-a)
# plt.ylim(0.7,1.05)
# plt.xlabel('y',fontsize=f_size-a)
# plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'best',n=2)
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.text(0.15, 0.9, proton + " p-"+atom, horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.text(0.05, 0.40, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.text(0.05, 0.50, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'R_all'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
# plt.show()
# plt.close()

