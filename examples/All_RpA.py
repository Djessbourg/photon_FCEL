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
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
from src import sigma as sig
from src import Collision

### Initialisation of a sigma object ###
rs = 8800
# rs = 200
s = (rs)**2 # CM energy in Gev2
# proton = "nNNPDF30_nlo_as_0118_p"
# proton = "NNPDF40_lo_as_01180"
# proton = 'MSHT20lo_as130'
# proton = 'CT18LO'
proton = "NNPDF40_nlo_as_01180"
# proton = 'MSHT20nlo_as118'
# proton = 'CT18NLO'

Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

atom = 'Pb'
# atom = 'Au'

# y_abs = 6 # for Pb
# y_abs = 4.5 # for Au

Atom = sig.Atom #dictionary where we stock important atomic aspect (now just Z and A) to plot and put into file names

Z = Atom[atom]['Z']
A = Atom[atom]['A']

pPb_cross_section = sig.Sigma(proton,Pb,s,Z,A)

d = sig.Switch

### first atempts to reproduce last results ###
p_T = 5 #GeV
# p_T = 2
x_T = (2.0*p_T)/rs
num = 0																			# The central set 
Y = sig.Y_list(x_T,Z)														

# A plot for sigma tot and its components
f_size = 17
fig_size=(5,5)
alph = 0.3
q0 =0.07
convention = 'dp_t2'
col_type = 'pp'
err ='q0,mu'

RpA_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_dir')) # the direcory to save data from this file
a=0
b=0								

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)

# before RpA, cross section of each process and ratios of thoses
# fig, ax = plt.subplots(constrained_layout=True,figsize=fig_size)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(2))
# (s1,er1),(s2,er2),(s3,er3),(s4,er4)=pPb_cross_section.dsigma_all_dy(x_T,num,is_pp=Collision(col_type),switch=convention,l=False)
# stot = s1+s2+s3+s4

# plt.plot(Y,s1,label=r'$qg$')
# plt.plot(Y,s2,label=r'$gq$')
# plt.plot(Y,s3,label=r'$q\bar{q}$')
# plt.plot(Y,s4,label=r'$\bar{q}q$')

# plt.ylabel(r'd$\sigma_{\ell}/$d$y$'+sig.Switch[convention][0]+' '+sig.Switch[convention][1],fontsize= f_size-a)
# plt.yscale('log')
# # plt.ylim(5*10**1,7*10**3)
# plt.ylim(2*10**2,7*10**4)
# plt.xlabel(r'y',fontsize=f_size-a)
# # plt.xlim(-4,4)
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# # plt.text(0.5, 0.95, proton + " p-"+atom, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.text(0.5, 0.28, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# # plt.text(0.5, 0.38, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.text(0.5, 0.38, r'$\sqrt{s} =$'+str(rs) +r' GeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'lower center',n=2)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, proton+'sigma_all'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
# plt.show()
# plt.close()

# # fractions
# fig, ax = plt.subplots(constrained_layout=True)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(2))
# plt.plot(Y,s1/stot,label=r'$\ell=1$')
# plt.plot(Y,s2/stot,label=r'$\ell=2$')
# plt.plot(Y,s3/stot,label=r'$\ell=3$')
# plt.plot(Y,s4/stot,label=r'$\ell=4$')
# plt.axhline(y=1,color= 'grey',alpha=0.3)

# plt.ylabel(r'd$f_{\ell}$',fontsize= f_size-a)
# plt.ylim(0,1.1)
# plt.xlabel('y',fontsize=f_size-a)
# plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'best',n=2)
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.text(0.5, 0.95, proton + " p-"+atom, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.text(0.5, 0.40, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.text(0.5, 0.50, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'f_all'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
# plt.show()
# plt.close()


# for the tot R_pA dir (computed with NNPDF40_lo_as_01180)
f_name = proton+'_RpA_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_'+str(p_T)+'GeV.txt'
# f_name = proton+'_True_RpA_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_'+str(p_T)+'GeV.txt'
# f_name = proton+'_RpA_wo_iso_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_'+str(p_T)+'GeV.txt'
if os.path.exists(os.path.join(RpA_dir,f_name)):
	print(f"The file '{f_name}' already exists. It is loaded.")
	Y,Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
	# Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
	# Y = sig.Y_list(x_T,Z)
	# r = [Y,Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf]
	# # np.savetxt(os.path.join(RpA_dir,proton+'_'+f_name), r)
	# np.savetxt(os.path.join(RpA_dir,f_name), r)
else:
	print("The file does not exists")
	Y,Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dy(x_T,switch = convention,var_err= err)
	# Y,Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_True_RpA_dy(x_T,switch = convention,var_err= err)
	# Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_wo_iso_dy(x_T,switch=convention,var_err=err)
	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
	r = [Y,Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf]
	np.savetxt(os.path.join(RpA_dir,f_name), r)
	print(f" '{f_name}' has been created")

# for the R_pA of each process
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

fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
plt.axhline(y=1, color='grey',alpha=alph)
# plt.plot(Y,R_1,label=r'$\ell=1$')
# plt.plot(Y,R_2,label=r'$\ell=2$')
# plt.plot(Y,R_3,label=r'$\ell=3$')
# plt.plot(Y,R_4,label=r'$\ell=4$')
# plt.plot(Y,Rpa,color='blue',linestyle= '--',label=r'$\ell=$ Sum')

# plt.plot(Y,Rpa,color='blue',label=r'$R_{\text{pA}}$') 
# plt.fill_between(Y,Rpa-Rpa_minus,Rpa+Rpa_plus,color = 'blue',alpha=0.3)

xi_values = [0.05*i for i in range(1, 20)]
cmap = cm.plasma
norm = mcolors.Normalize(vmin=min(xi_values), vmax=max(xi_values))

R_pA_FCEL,R_pA_FCEG= pPb_cross_section.Rpp_FCELG_dy(x_T,num,q0,FCELG=True)
		 
for xi in xi_values:
	R_FCEL, R_FCEG, Rpa_taylor = pPb_cross_section.Rpp_Taylor_FCELG_dy(x_T,1-xi,num,q0,switch=convention,FCELG=True) # i put 1-xi to be coherent with the paper
	color = cmap(norm(xi))
	# plt.plot(Y, Rpa_taylor,alpha=0.6, color=color)
	plt.plot(Y, R_FCEL,alpha=0.6, color=color)
	plt.plot(Y, R_FCEG,alpha=0.6, color=color)

plt.plot(Y,R_pA_FCEL,color = 'blue',alpha=0.5,label='FCEL')			
plt.plot(Y,R_pA_FCEG,color = 'red',alpha= 0.5, label = 'FCEG')	

sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label(r'$\tilde{\xi}$', fontsize=f_size)

#plt.grid()
# plt.ylim(bottom= 0.9,top = 1.05)
# plt.ylabel(r'$R_{\ell}=\sigma_{\ell}^\text{pA}/A\sigma_{\ell}^\text{pp}$',fontsize= f_size-a)
plt.ylim(0.8,1.1)
plt.xlim(-6,6)
plt.xlabel('y',fontsize=f_size)
plot_usuals(s1=f_size,s2=f_size,loca = 'lower left',n=1)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.text(0.8, 0.9," p-"+atom, horizontalalignment='right', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.text(0.1, 0.9, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
# plt.text(0.1, 0.8, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
plt.text(0.1, 0.8, r'$\sqrt{s} =$'+str(rs) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'R_all'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.savefig(os.path.join(plots_dir, 'R_pA_vs_Taylor'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.savefig(os.path.join(plots_dir, 'R_pA_FCELG_vs_Taylor'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.show()
plt.close()

