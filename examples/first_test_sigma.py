# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity with the sigma class.
# Study on RpA and isospin effects.
# last modified: 13/02/2026
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
# rs = 200
rs = 8800
s = (rs)**2 # CM energy in Gev2
proton = "NNPDF40_nlo_as_01180"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

Atom = sig.Atom

atom = 'Pb'
# atom = 'Au'

Z = Atom[atom]['Z']
A = Atom[atom]['A']

pPb_cross_section = sig.Sigma(proton,Pb,s,Z,A)

d = sig.Switch

### first atempts to reproduce last results ###
# p_T = 2
p_T = 5
p_T2= 10
y = -4
y2= 0
y3= 4
x_T = (2.0*p_T)/rs
x_T2 = (2.0*p_T2)/rs
num = 0 																# The central set 

# A plot for sigma tot and its components
f_size = 17
alph = 0.3
convention = 'd2p_t'
col_type = 'pp'

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)


# Need to change here the mu into mu_factor
# ~ sigma_tot,err_tot= pPb_cross_section.dsigma_tot_dy(x_T,mu2,mu_f2,num,is_pp=Collision(col_type),switch = convention)
# ~ sigma_qG, err_qG  = pPb_cross_section.dsigma_qG_dy(x_T,mu2,mu_f2,num,is_pp=Collision(col_type),switch = convention)
# ~ sigma_Gq, err_Gq  = pPb_cross_section.dsigma_Gq_dy(x_T,mu2,mu_f2,num,is_pp=Collision(col_type),switch = convention)
# ~ sigma_qqbar, err_qqbar = pPb_cross_section.dsigma_qqbar_dy(x_T,mu2,mu_f2,num,is_pp=Collision(col_type),switch = convention)
# ~ sigma_qbarq, err_qbarq = pPb_cross_section.dsigma_qbarq_dy(x_T,mu2,mu_f2,num,is_pp=Collision(col_type),switch = convention)

# Ucen,(Uplus,Uminus), Umu , Updf = pPb_cross_section.Uncertaintites_dy(x_T,is_pp=Collision(col_type),switch = convention,var_err= err)
# Uplus_mu , Uminus_mu  = Umu[0],Umu[1]
# Uplus_pdf,Uminus_pdf = Updf[0],Updf[1]
# Y = sig.Y_list(x_T)

# ax = plt.subplot()
# plt.fill_between(Y,Ucen+Uplus_mu,Ucen-Uminus_mu,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# plt.fill_between(Y,Ucen+Uplus_pdf,Ucen-Uminus_pdf,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# plt.plot(Y,Ucen,color='blue',label=r'cental value')
# #plt.grid()
# plt.legend(loc='lower right',frameon =False,fontsize=10)
# plt.xlabel('y')
# plt.ylim(bottom=0)
# plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# plt.text(0.1, 0.1, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.text(0.1, 0.05, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# plt.show()
# plt.close()

# Ucen2,(Uplus2,Uminus2), Umu2 , Updf2 = pPb_cross_section.Uncertaintites_dy(x_T2,is_pp=Collision(col_type),switch = convention,var_err= err)
# Uplus_mu2 , Uminus_mu2  = Umu2[0],Umu2[1]
# Uplus_pdf2,Uminus_pdf2 = Updf2[0],Updf2[1]
# Y2 = sig.Y_list(x_T2)

# ax = plt.subplot()
# plt.fill_between(Y2,Ucen2+Uplus_mu2,Ucen2-Uminus_mu2,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# plt.fill_between(Y2,Ucen2+Uplus_pdf2,Ucen2-Uminus_pdf2,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# plt.plot(Y2,Ucen2,color='blue',label=r'cental value')
# #plt.grid()
# plt.legend(loc='lower right',frameon =False,fontsize=10)
# plt.xlabel('y')
# plt.ylim(bottom=0)
# plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# plt.text(0.1, 0.1, r'$p_\bot =$'+str(p_T2) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.text(0.1, 0.05, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T2)+'GeV.pdf'))
# plt.show()
# plt.close()

# ax = plt.subplot()
# plt.fill_between(Y,Ucen+Uplus_mu,Ucen-Uminus_mu,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# plt.fill_between(Y,Ucen+Uplus_pdf,Ucen-Uminus_pdf,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# plt.plot(Y,Ucen,color='blue',label=r'$p_\bot=$ '+str(p_T)+r' $GeV$')
# plt.fill_between(Y2,Ucen2+Uplus_mu2,Ucen2-Uminus_mu2,color = 'blue', alpha= 0.3 )
# plt.fill_between(Y2,Ucen2+Uplus_pdf2,Ucen2-Uminus_pdf2,color = 'green', alpha= 0.4)
# plt.plot(Y2,Ucen2,color='green',label=r'$p_\bot=$ '+str(p_T2)+r' $GeV$')
# #plt.grid()
# plt.legend(loc='lower right',frameon =False,fontsize=10)
# plt.xlabel('y')
# plt.ylim(bottom=1,top=1e4)
# plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# plt.yscale('log')
# plt.text(0.1, 0.05, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'_and_'+str(p_T2)+'GeV.pdf'))
# plt.show()
# plt.close()

# ~ plt.plot(Y,sigma_tot,label=r'$d\sigma_{tot}/dy$'+d[convention][0])
# ~ #plt.errorbar(Y,sigma_tot,yerr=err_tot,fmt='+')
# ~ plt.plot(Y,sigma_qG,label=r'$d\sigma_{qG}/dy$'+d[convention][0])
# ~ #plt.errorbar(Y,sigma_qG,yerr=err_qG,fmt='+')
# ~ plt.plot(Y,sigma_Gq,label=r'$d\sigma_{Gq}/dy$'+d[convention][0])
# ~ #plt.errorbar(Y,sigma_Gq,yerr=err_Gq,fmt='+')
# ~ plt.plot(Y,sigma_qqbar,label=r'$d\sigma_{qqbar}/dy$'+d[convention][0])
# ~ #plt.errorbar(Y,sigma_qqbar,yerr=err_qqbar,fmt='+')
# ~ plt.plot(Y,sigma_qbarq,label=r'$d\sigma_{qbarq}/dy$'+d[convention][0])
# ~ #plt.errorbar(Y,sigma_qbarq,yerr=err_qbarq,fmt='+')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,ncol=2)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'cross section $d\sigma/dy$'+d[convention][0]+ r' in '+d[convention][1])
# ~ #plt.yscale('log')
# ~ plt.title(r'$d\sigma/dy$'+d[convention][0]+ r'components for pp, with $p_\bot =$'+str(x_T*rs/2.)+r' $GeV$')
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_all_components_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))

#pt dependant
# Ucen, (Uplus,Uminus),Umu , Updf = pPb_cross_section.Uncertaintites_dpt(y,is_pp=Collision(col_type),switch = convention,var_err= err)
# Uplus_mu , Uminus_mu  = Umu[0],Umu[1]
# Uplus_pdf,Uminus_pdf = Updf[0],Updf[1]
# P_T = pPb_cross_section.P_T_list(y)

# ~ ax = plt.subplot()
# ~ plt.fill_between(P_T,Ucen+Uplus_mu,Ucen-Uminus_mu,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# ~ plt.fill_between(P_T,Ucen+Uplus_pdf,Ucen-Uminus_pdf,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# ~ plt.plot(P_T,Ucen,color='blue',label=r'cental value',linewidth=0.9)
# ~ #plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,fontsize=10)
# ~ plt.xlabel(r'$p_\bot$ ($GeV$)')
# ~ plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# ~ plt.yscale('log')
# ~ plt.text(0.1, 0.05, r'$y =$ '+str(y), horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.text(0.1, 0.1, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_y'+str(y)+'.pdf'))
# ~ plt.close()

# Ucen2,(Uplus2,Uminus2), Umu2 , Updf2 = pPb_cross_section.Uncertaintites_dpt(y2,is_pp=Collision(col_type),switch = convention,var_err= err)
# Uplus_mu2 , Uminus_mu2  = Umu2[0],Umu2[1]
# Uplus_pdf2,Uminus_pdf2 = Updf2[0],Updf2[1]
# P_T2 = pPb_cross_section.P_T_list(y2)

# ~ ax = plt.subplot()
# ~ plt.fill_between(P_T2,Ucen2+Uplus_mu2,Ucen2-Uminus_mu2,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# ~ plt.fill_between(P_T2,Ucen2+Uplus_pdf2,Ucen2-Uminus_pdf2,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# ~ plt.plot(P_T2,Ucen2,color='blue',label=r'cental value',linewidth=0.9)
# ~ #plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,fontsize=10)
# ~ plt.xlabel(r'$p_\bot$ ($GeV$)')
# ~ plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# ~ plt.yscale('log')
# ~ plt.text(0.1, 0.05, r'$y =$ '+str(y2), horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.text(0.1, 0.1, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_y'+str(y2)+'.pdf'))
# ~ plt.close()

# ~ Ucen2, Umu2 , Updf2 = pPb_cross_section.Uncertaintites_dpt(y3,is_pp=Collision(col_type),switch = convention,var_err= err)
# ~ Uplus_mu2 , Uminus_mu2  = Umu2[0],Umu2[1]
# ~ Uplus_pdf2,Uminus_pdf2 = Updf2[0],Updf2[1]
# ~ P_T2 = pPb_cross_section.P_T_list(y3)

# ~ ax = plt.subplot()
# ~ plt.fill_between(P_T2,Ucen2+Uplus_mu2,Ucen2-Uminus_mu2,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# ~ plt.fill_between(P_T2,Ucen2+Uplus_pdf2,Ucen2-Uminus_pdf2,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# ~ plt.plot(P_T2,Ucen2,color='blue',label=r'cental value',linewidth=0.9)
# ~ #plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,fontsize=10)
# ~ plt.xlabel(r'$p_\bot$ ($GeV$)')
# ~ plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# ~ plt.yscale('log')
# ~ plt.text(0.1, 0.05, r'$y =$ '+str(y3), horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.text(0.1, 0.1, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_y'+str(y3)+'.pdf'))
# ~ plt.close()

# p=0
# Delta = pow(10,p)																# To seperate the two different plots
# ax = plt.subplot()
# plt.fill_between(P_T,Ucen+Uplus_mu,Ucen-Uminus_mu,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# plt.fill_between(P_T,Ucen+Uplus_pdf,Ucen-Uminus_pdf,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# plt.plot(P_T,Ucen,color='blue',label=r'$y=$ '+str(y),linewidth=0.9)
# plt.fill_between(P_T2,(Ucen2+Uplus_mu2)*Delta,(Ucen2-Uminus_mu2)*Delta,color = 'blue', alpha= 0.3)
# plt.fill_between(P_T2,(Ucen2+Uplus_pdf2)*Delta,(Ucen2-Uminus_pdf2)*Delta,color = 'green', alpha= 0.4)
# plt.plot(P_T2,Ucen2*Delta,color='green',label=r'$y= $'+str(y2)+r', ($\times 10^{}$)'.format('{'+str(p)+'}'),linewidth=0.9)
# #plt.grid()
# plt.legend(loc='upper right',frameon =False,fontsize=10)
# plt.xlabel(r'$p_\bot$ ($GeV$)')
# plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# plt.yscale('log')
# plt.text(0.1, 0.05, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_y'+str(y)+'_and_'+str(y2)+'_Delta'+str(Delta)+'.pdf'))
# plt.show()
# plt.close()

# p = -1*2
# Delta = pow(10,p)																# To seperate the two different plots
# ax = plt.subplot()
# plt.fill_between(P_T,Ucen+Uplus_mu,Ucen-Uminus_mu,color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# plt.fill_between(P_T,Ucen+Uplus_pdf,Ucen-Uminus_pdf,color = 'green', alpha= 0.4 ,label=r'$\sigma_{pdf}$')
# plt.plot(P_T,Ucen,color='blue',label=r'$y=$ '+str(y),linewidth=0.9)
# plt.fill_between(P_T2,(Ucen2+Uplus_mu2)*Delta,(Ucen2-Uminus_mu2)*Delta,color = 'blue', alpha= 0.3)
# plt.fill_between(P_T2,(Ucen2+Uplus_pdf2)*Delta,(Ucen2-Uminus_pdf2)*Delta,color = 'green', alpha= 0.4)
# plt.plot(P_T2,Ucen2*Delta,color='green',label=r'$y= $'+str(y2)+r', ($\times 10^{}$)'.format('{'+str(p)+'}'),linewidth=0.9)
# #plt.grid()
# plt.legend(loc='upper right',frameon =False,fontsize=10)
# plt.xlabel(r'$p_\bot$ ($GeV$)')
# plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# plt.yscale('log')
# plt.text(0.1, 0.05, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_uncertainties_'+str(rs)+'GeV_y'+str(y)+'_and_'+str(y2)+'_Delta'+str(Delta)+'.pdf'))
# plt.show()
# plt.close()

### Ratios ###
# FCEL/G
# y dependant
# err ='q0,mu,pdf'
err = 'mu'
RpA_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_dir')) # the direcory to save data from this file
a=2

# f_name = 'RpA_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.txt'
f_name = proton+'_RpA_wo_iso_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+str(p_T)+'GeV.txt'
if os.path.exists(os.path.join(RpA_dir,f_name)):
	print(f"The file '{f_name}' already exists. It is loaded.")
	Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
else:
	print("The file does not exists")
	# Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dy(x_T,switch = convention,var_err= err)
	Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_wo_iso_dy(x_T,switch=convention,var_err=err)
	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
	r = [Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf]
	np.savetxt(os.path.join(RpA_dir,f_name), r)
	print(f" '{f_name}' has been created")

Y = sig.Y_list(x_T)

m_mu,M_mu = sig.min_max([Rpa+Rpa_plus_mu,Rpa-Rpa_minus_mu,Rpa])
m_q,M_q = sig.min_max([Rpa+Rpa_plus_q,Rpa-Rpa_minus_q,Rpa])
m_pdf, M_pdf = sig.min_max([Rpa+Rpa_plus_pdf,Rpa-Rpa_minus_pdf])

fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
plt.fill_between(Y,m_mu,M_mu,color = 'blue',alpha = alph,label=r'$\sigma_\mu $')
plt.fill_between(Y,m_q,M_q,color='red', alpha=alph,label=r'$\sigma_{\hat{q}_0}$')
plt.fill_between(Y,m_pdf,M_pdf,color='green', alpha=alph,label=r'$\sigma_{pdf}$')
plt.plot(Y,Rpa+Rpa_plus,color = 'blue', alpha = 0.5)
plt.plot(Y,Rpa-Rpa_minus,color = 'blue', alpha = 0.5)
plt.axhline(y=1, color='grey',alpha=alph)
plt.plot(Y,Rpa,color='blue')

plt.legend(loc='lower right',frameon =False,fontsize=f_size-a)
plt.xlabel('y',fontsize=f_size-a)
plt.ylabel(r'$R_{\text{pA}}^{\text{dir}}$',fontsize= f_size)
plt.ylim(0.8,1.1)
# plt.xlim(-4,4)
plot_usuals(s1=f_size,s2=f_size,loca = 'lower right')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.text(0.5, 0.9, proton + " p-"+atom, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
plt.text(0.2, 0.2, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
plt.text(0.2, 0.1, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'RpA_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.show()
plt.close()

# f_name = 'RpA_'+str(rs)+'GeV_'+convention+str(p_T2)+'GeV.txt'
# if os.path.exists(os.path.join(RpA_dir,f_name)):
# 	print(f"The file '{f_name}' already exists. It is loaded.")
# 	Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
# else:
# 	print("The file does not exists")
# 	Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dy(x_T2,switch = convention,var_err= err)
# 	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
# 	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
# 	r = [Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf]
# 	np.savetxt(os.path.join(RpA_dir,f_name), r)
# 	print(f" '{f_name}' has been created")

# Y = sig.Y_list(x_T2)

# m_mu,M_mu = sig.min_max([Rpa+Rpa_plus_mu,Rpa-Rpa_minus_mu,Rpa])
# m_q,M_q = sig.min_max([Rpa+Rpa_plus_q,Rpa-Rpa_minus_q,Rpa])
# fig, ax = plt.subplots(constrained_layout=True)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(2))
# plt.fill_between(Y,m_mu,M_mu,color = 'blue',alpha = alph,label=r'$\sigma_\mu $')
# plt.fill_between(Y,m_q,M_q,color='red', alpha=alph,label=r'$\sigma_{\hat{q}_0}$')
# plt.axhline(y=1, color='grey',alpha=alph)
# plt.plot(Y,Rpa,color='blue')
# #plt.grid()
# plt.legend(loc='lower right',frameon =False,fontsize=f_size-a)
# plt.xlabel('y')
# # plt.ylim(bottom= 0.9,top = 1.05)
# plt.ylabel(r'$R_{pA}^{dir}$',fontsize= f_size-a)
# plt.ylim(0.8,1.1)
# plt.xlabel('y',fontsize=f_size-a)
# plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'lower right')
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.text(0.1, 0.1, r'$p_\bot =$'+str(p_T2) +r' $GeV$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'RpA_'+str(rs)+'GeV_'+convention+str(p_T2)+'GeV.pdf'),bbox_inches="tight")
# plt.show()
# plt.close()


#pt dependant
# f_name = 'RpA_'+str(rs)+'GeV_'+convention+str(y)+'.txt'
# if os.path.exists(os.path.join(RpA_dir,f_name)):
# 	print(f"The file '{f_name}' already exists. It is loaded.")
# 	Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
# else:
# 	print("The file does not exists")
# 	Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dpt(y,switch = convention,var_err= err,)
# 	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
# 	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
# 	r = [Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf]
# 	np.savetxt(os.path.join(RpA_dir,f_name), r)
# 	print(f" '{f_name}' has been created")

# P_T = pPb_cross_section.P_T_list(y)

# m_mu,M_mu = sig.min_max([Rpa+Rpa_plus_mu,Rpa-Rpa_minus_mu,Rpa])
# m_q,M_q = sig.min_max([Rpa+Rpa_plus_q,Rpa-Rpa_minus_q,Rpa])
# fig, ax = plt.subplots(constrained_layout=True)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(5))
# plt.fill_between(P_T,m_mu,M_mu,color = 'blue',alpha = alph,label=r'$\sigma_\mu $')
# plt.fill_between(P_T,m_q,M_q,color='red', alpha=alph,label=r'$\sigma_{\hat{q}_0}$')
# plt.axhline(y=1, color='grey',alpha=alph)
# plt.plot(P_T,Rpa,color='blue')
# plt.legend(loc='lower right',frameon =False,fontsize=f_size-a)
# plt.ylabel(r'$R_{pA}^{dir}$',fontsize= f_size-a)
# plt.ylim(0.8,1.1)
# plt.xlabel(r'$p_\bot$ (GeV)',fontsize=f_size-a)
# plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'lower right')
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.text(0.1, 0.1, r'$y =$ '+str(y), horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'RpA_'+str(rs)+'GeV_'+'y'+str(y)+'.pdf'))
# plt.show()
# plt.close()


# f_name = 'RpA_'+str(rs)+'GeV_'+convention+str(y2)+'.txt'
# if os.path.exists(os.path.join(RpA_dir,f_name)):
# 	print(f"The file '{f_name}' already exists. It is loaded.")
# 	Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
# else:
# 	print("The file does not exists")
# 	Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dpt(y2,switch = convention,var_err= err,)
# 	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
# 	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
# 	r = [Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf]
# 	np.savetxt(os.path.join(RpA_dir,f_name), r)
# 	print(f" '{f_name}' has been created")
	
# P_T = pPb_cross_section.P_T_list(y2)

# m_mu,M_mu = sig.min_max([Rpa+Rpa_plus_mu,Rpa-Rpa_minus_mu,Rpa])
# m_q,M_q = sig.min_max([Rpa+Rpa_plus_q,Rpa-Rpa_minus_q,Rpa])
# fig, ax = plt.subplots(constrained_layout=True)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(5))
# plt.fill_between(P_T,m_mu,M_mu,color = 'blue',alpha = alph,label=r'$\sigma_\mu $')
# plt.fill_between(P_T,m_q,M_q,color='red', alpha=alph,label=r'$\sigma_{\hat{q}_0}$')
# plt.axhline(y=1, color='grey',alpha=alph)
# plt.plot(P_T,Rpa,color='blue')
# plt.legend(loc='lower right',frameon =False,fontsize=f_size-a)
# plt.ylabel(r'$R_{pA}^{dir}$',fontsize= f_size-a)
# plt.ylim(0.8,1.1)
# plt.xlabel(r'$p_\bot$ (GeV)',fontsize=f_size-a)
# plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'lower right')
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.text(0.1, 0.1, r'$y =$ '+str(y2), horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'RpA_'+str(rs)+'GeV_'+'y'+str(y2)+'.pdf'))
# plt.show()
# plt.close()

# f_name = 'RpA_'+str(rs)+'GeV_'+convention+str(y3)+'.txt'
# if os.path.exists(os.path.join(RpA_dir,f_name)):
# 	print(f"The file '{f_name}' already exists. It is loaded.")
# 	Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
# else:
# 	print("The file does not exists")
# 	Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dpt(y3,switch = convention,var_err= err,)
# 	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
# 	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
# 	r = [Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpas_minus_pdf]
# 	np.savetxt(os.path.join(RpA_dir,f_name), r)
# 	print(f" '{f_name}' has been created")
	
# P_T = pPb_cross_section.P_T_list(y3)

# m_mu,M_mu = sig.min_max([Rpa+Rpa_plus_mu,Rpa-Rpa_minus_mu,Rpa])
# m_q,M_q = sig.min_max([Rpa+Rpa_plus_q,Rpa-Rpa_minus_q,Rpa])
# fig, ax = plt.subplots(constrained_layout=True)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(5))
# plt.fill_between(P_T,m_mu,M_mu,color = 'blue',alpha = alph,label=r'$\sigma_\mu $')
# plt.fill_between(P_T,m_q,M_q,color='red', alpha=alph,label=r'$\sigma_{\hat{q}_0}$')
# plt.axhline(y=1, color='grey',alpha=alph)
# plt.plot(P_T,Rpa,color='blue')
# plt.legend(loc='lower right',frameon =False,fontsize=f_size-a)
# plt.ylabel(r'$R_{pA}^{dir}$',fontsize= f_size-a)
# plt.ylim(0.8,1.1)
# plt.xlabel(r'$p_\bot$ (GeV)',fontsize=f_size-a)
# plot_usuals(s1=f_size-a,s2=f_size-a,loca = 'lower right')
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.text(0.1, 0.1, r'$y =$ '+str(y3), horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'RpA_'+str(rs)+'GeV_'+'y'+str(y3)+'.pdf'))
# plt.show()
# plt.close()

# Delta Rpa q+/q-
# y

# ~ Rpa_plus,Rpa_minus = pPb_cross_section.Delta_RpA_plusminus_dy(x_T,switch = convention)
# ~ Y = sig.Y_list(x_T)

# ~ ax = plt.subplot()
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.plot(Y,Rpa_plus,color='red',label=r'$q^+$')
# ~ plt.plot(Y,Rpa_minus,color='blue',label=r'$q^-$')
# ~ #plt.grid()
# ~ plt.legend(loc='lower right',frameon =False,fontsize=8)
# ~ plt.xlabel(r'y')
# ~ plt.ylabel(r'$\Delta R_{pA}^{q^\pm}(y)$')
# ~ plt.text(0.90, 0.95, r'$p_\bot =$ '+str(p_T)+r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Delta_RpA_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()

# ~ Rpa_plus,Rpa_minus = pPb_cross_section.Delta_RpA_plusminus_dy(x_T2,switch = convention)
# ~ Y = sig.Y_list(x_T2)

# ~ ax = plt.subplot()
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.plot(Y,Rpa_plus,color='red',label=r'$q^+$')
# ~ plt.plot(Y,Rpa_minus,color='blue',label=r'$q^-$')
# ~ #plt.grid()
# ~ plt.legend(loc='lower right',frameon =False,fontsize=8)
# ~ plt.xlabel(r'y')
# ~ plt.ylabel(r'$\Delta R_{pA}^{q^\pm}(y)$')
# ~ plt.text(0.90, 0.95, r'$p_\bot =$ '+str(p_T2)+r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Delta_RpA_'+str(rs)+'GeV_'+convention+str(p_T2)+'GeV.pdf'))
# ~ plt.close()

# pt 

# ~ Rpa_plus,Rpa_minus = pPb_cross_section.Delta_RpA_plusminus_dpt(y,switch = convention)
# ~ P_T = pPb_cross_section.P_T_list(y)

# ~ ax = plt.subplot()
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.plot(P_T,Rpa_plus,color='red',label=r'$q^+$')
# ~ plt.plot(P_T,Rpa_minus,color='blue',label=r'$q^-$')
# ~ #plt.grid()
# ~ plt.legend(loc='lower right',frameon =False,fontsize=8)
# ~ plt.xlabel(r'$p_\bot$ ($GeV$)')
# ~ plt.ylabel(r'$\Delta R_{pA}^{q^\pm}(p_\bot)$')
# ~ plt.text(0.90, 0.95, r'$y =$ '+str(y), horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Delta_RpA_'+str(rs)+'GeV_'+'y'+str(y)+'.pdf'))
# ~ plt.close()

# ~ Rpa_plus,Rpa_minus = pPb_cross_section.Delta_RpA_plusminus_dpt(y2,switch = convention)
# ~ P_T = pPb_cross_section.P_T_list(y2)

# ~ ax = plt.subplot()
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.plot(P_T,Rpa_plus,color='red',label=r'$q^+$')
# ~ plt.plot(P_T,Rpa_minus,color='blue',label=r'$q^-$')
# ~ #plt.grid()
# ~ plt.legend(loc='lower right',frameon =False,fontsize=8)
# ~ plt.xlabel(r'$p_\bot$ ($GeV$)')
# ~ plt.ylabel(r'$\Delta R_{pA}^{q^\pm}(p_\bot)$')
# ~ plt.text(0.90, 0.95, r'$y =$ '+str(y2), horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Delta_RpA_'+str(rs)+'GeV_'+'y'+str(y2)+'.pdf'))
# ~ plt.close()

# ~ Rpa_plus,Rpa_minus = pPb_cross_section.Delta_RpA_plusminus_dpt(y3,switch = convention)
# ~ P_T = pPb_cross_section.P_T_list(y3)

# ~ ax = plt.subplot()
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.plot(P_T,Rpa_plus,color='red',label=r'$q^+$')
# ~ plt.plot(P_T,Rpa_minus,color='blue',label=r'$q^-$')
# ~ #plt.grid()
# ~ plt.legend(loc='lower right',frameon =False,fontsize=8)
# ~ plt.xlabel(r'$p_\bot$ ($GeV$)')
# ~ plt.ylabel(r'$\Delta R_{pA}^{q^\pm}(p_\bot)$')
# ~ plt.text(0.90, 0.95, r'$y =$ '+str(y3), horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Delta_RpA_'+str(rs)+'GeV_'+'y'+str(y3)+'.pdf'))
# ~ plt.close()

# Isospin

# U_Rpa, U_Riso, U_Riso_FCELG, U_Rwoiso = pPb_cross_section.Uncertainties_4R_iso_dy(x_T,switch = convention, var_err= 'q0,mu,pdf') 
# Rpa , Rpa_plus, Rpa_minus  = U_Rpa[0], U_Rpa[1], U_Rpa[2]
# Riso, Riso_plus, Riso_minus = U_Riso [0], U_Riso[1], U_Riso[2]
# Riso_FCELG, Riso_FCELG_plus, Riso_FCELG_minus = U_Riso_FCELG [0], U_Riso_FCELG[1], U_Riso_FCELG[2]
# Rwoiso, Rwoiso_plus, Rwoiso_minus = U_Rwoiso[0], U_Rwoiso[1], U_Rwoiso[2]

# Y = sig.Y_list(x_T)

# ax = plt.subplot()
# plt.fill_between(Y,Rpa+Rpa_plus,Rpa-Rpa_minus,color='blue', alpha = 0.2)
# plt.fill_between(Y,Riso+Riso_plus,Riso-Riso_minus,color='blue', alpha = 0.2)
# plt.fill_between(Y,Riso_FCELG+Riso_FCELG_plus,Riso_FCELG-Riso_FCELG_minus,color='blue', alpha = 0.2)
# # ~ plt.fill_between(Y,Rwoiso+Rwoiso_plus,Rwoiso-Rwoiso_minus,color='blue', alpha = 0.2)
# plt.plot(Y,Rpa,color = 'blue',label=r'$R_{pA}$')
# plt.plot(Y,Riso,color = 'orange',label=r'$R_{pA}^{iso}$')
# plt.plot(Y,Riso_FCELG,color = 'red',label=r'$R_{pA}^{FCELG+iso}$')
# plt.plot(Y,Rwoiso,color='green',label=r'$R_{pA}^{FCELG-iso}$')

# plt.legend(loc='lower right',frameon =False,fontsize=10)
# plt.axhline(y=1, color='grey', linestyle='--')
# plt.xlabel('y')
# plt.ylabel(r'$R$')
# plt.ylim(bottom=0.70,top=1.05)
# plt.text(0.90, 0.95, r'$p_\bot =$ '+str(p_T)+r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, '4R_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# plt.show()
# plt.close()

# ~ plt.plot(Y,R_all[0][0],label=r'qG/tot')
# ~ plt.plot(Y,R_all[1][0],label=r'Gq/tot')
# ~ plt.plot(Y,R_all[2][0],label=r'qqbar/tot')
# ~ plt.plot(Y,R_all[3][0],label=r'qbarq/tot')

# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,ncol=2)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'Ratios of different cross sections for pp')
# ~ #plt.yscale('log')
# ~ plt.title(r'Component ratios, with $p_\bot =$'+str(x_T*rs/2.)+r' $GeV$')
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_component_ratios_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()
