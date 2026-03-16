# -*- coding: utf-8 -*-

# =============================================================================
# Study of the impact of adding a mass to the photon on the cross section
# Last modified: 13/01/2026
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))
RpA_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_dir')) # the direcory to save data from this file

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, ScalarFormatter
from src import sigma as sig
f_size = 17
a = 0
b= 0

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)
	
### Initialisation of a sigma object ###
rs = 8800
s = (rs)**2 # CM energy in Gev2
proton = "NNPDF40_nlo_as_01180"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

M = [0,2,5,10,15] # vitual photon mass in GeV forgot the 2 GeV in the f_name of the True_RpA
cen = 0

Atom = sig.Atom

atom = 'Pb'
# atom = 'Au'

Z = Atom[atom]['Z']
A = Atom[atom]['A']

pPb_cross_section = sig.Sigma(proton,Pb,s,Z,A)

d = sig.Switch
convention = 'd2p_t'													# Define a cross section convention btw dp_t, dp_t2 ans d2p_t
col_type = 'pp'
Collision = sig.Collision

p_T = 5																# GeV
x_T = (2.*p_T)/rs 
### Xi function ###
Xi_dict = sig.Xi_dict
Xi_qg = Xi_dict["qg"]["M"]
Xi_gq = Xi_dict["gq"]["M"]
Xi_qqbar = Xi_dict["qqbar"]["M"]

y_space = [-6,6]
m = 0 # GeV
b = (m/p_T)**2

alph = 0.3
fig_size = (5,5)


#comparison btw two formulas (seems to be perfect at m = 0)
b = (m/p_T)**2
Y = sig.Y_list_M(x_T,b)

# qg_M, gq_M, qqbar_M, qbarq_M = [],[],[],[]
# qg, gq, qqbar, qbarq = [],[],[],[]
# R_M = pPb_cross_section.dsigma_all_dy_M(x_T,0,cen,iso='n',is_pp=Collision(col_type),switch=convention)
# qg_M=R_M[0][0];gq_M=R_M[1][0];qqbar_M=R_M[2][0];qbarq_M=R_M[3][0]
# R = pPb_cross_section.dsigma_all_dy(x_T,cen,iso='n',is_pp=Collision(col_type),switch=convention)
# qg=R[0][0];gq=R[1][0];qqbar=R[2][0];qbarq=R[3][0]
# plt.plot(Y,qg,label=r'qG')
# plt.plot(Y,gq,label=r'Gq')
# plt.plot(Y,qqbar,label=r'qqbar')
# plt.plot(Y,qbarq,label=r'qbarq')

# plt.plot(Y,qg_M,label=r'qG_M')
# plt.plot(Y,gq_M,label=r'Gq_M')
# plt.plot(Y,qqbar_M,label=r'qqbar_M')
# plt.plot(Y,qbarq_M,label=r'qbarq_M')

# C_np = pPb_cross_section.C_n_p_dy(x_T,cen,switch=convention)
# C_np_M = pPb_cross_section.C_n_p_dy_M(x_T,0,cen,switch=convention)

# plt.plot(Y,C_np,label = 'C_np')
# plt.plot(Y,C_np_M,label='C_np_M')
# plt.show()

# pp = pPb_cross_section.dsigma_tot_dy(x_T,cen,is_pp=Collision(col_type),switch=convention)[0]
# pn = pPb_cross_section.dsigma_tot_dy(x_T,cen,iso='n',is_pp=Collision(col_type),switch=convention)[0]

# pp_M = pPb_cross_section.dsigma_tot_dy_M(x_T,0,cen,is_pp=Collision(col_type),switch=convention)[0]
# pn_M = pPb_cross_section.dsigma_tot_dy_M(x_T,0,cen,iso='n',is_pp=Collision(col_type),switch=convention)[0]

# plt.plot(Y,pp,label='pp')
# plt.plot(Y,pp_M,label='pp_M')

# plt.plot(Y,pn,label='pn')
# plt.plot(Y,pn_M,label='pn_M')

# plt.legend()
# plt.show()

# plt.axhline(y=1, color='grey', alpha=0.3)
# plot_usuals(n=2)
# plt.xlabel(r'$y$',fontsize=f_size)
# plt.ylabel(r'$R_i = \sigma_i/\sigma_{tot}$')
# plt.title(r'$\gamma^\star$ production for $M=$'+str(m)+r' GeV and $p_\perp = $'+str(p_T)+' GeV')
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_component_virtual_ratios_'+str(rs)+'GeV_'+convention+'_p_T'+str(p_T)+'GeV_M'+str(m)+'GeV'+'.pdf'))
# plt.show()
# plt.close()


# for m in M:
#     b = (m/p_T)**2
#     epsilon = pow(10,-3)
#     Xi_min,Xi_max = sig.Xi_min_M(0,x_T,b)+epsilon,sig.Xi_max_M(0,x_T,b)-epsilon
#     Xi = np.linspace(Xi_min,Xi_max,100)
#     plt.plot(Xi,Xi_qg(Xi,b),label=r'$\Xi_{qg}$')
#     plt.plot(Xi,Xi_gq(Xi,b),label=r'$\Xi_{gq}$')
#     plt.plot(Xi,Xi_qqbar(Xi,b),label=r'$\Xi_{q\bar{q}}$')
#     plt.xlim(0,1)
#     plt.yscale('log')
#     plt.xlabel(r'$\xi$',fontsize=f_size)
#     plt.axvline(x=Xi_min, color='grey', alpha=0.3)
#     plt.axvline(x=Xi_max, color='grey', alpha=0.3)
#     plot_usuals(n=2)
#     # plt.ylabel(r'$d\sigma/d\xi$')
#     plt.title(r'$M=$'+str(m)+r' GeV, $p_\perp = $'+str(p_T)+r' GeV')
#     plt.tight_layout()
#     plt.savefig(os.path.join(plots_dir, 'Xi'+col_type+'_component_virtual_'+str(rs)+'GeV_'+convention+'_M'+str(m)+'GeV_pt'+str(p_T)+'GeV.pdf' ))
#     plt.show()
#     plt.close()


### y dependant ###

# fig, ax = plt.subplots(constrained_layout=True,figsize=fig_size)
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# # ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.xaxis.set_major_locator(MultipleLocator(2))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05))
# ax.yaxis.set_major_locator(MultipleLocator(0.1))
# plt.axhline(y=1, color='grey',alpha=alph)

for m in M:
    # b = (m/p_T)**2
    # Y = sig.Y_list_M(x_T,b)
    # # print((min(Y),max(Y)))
    # R_qg, R_gq, R_qqbar, R_qbarq = [],[],[],[]
    # for y in Y:
    #     R = pPb_cross_section.R_composition_dydpt_M(y,x_T,m,cen,is_pp=Collision(col_type),switch=convention)
    #     R_qg.append(R[0][0]);R_gq.append(R[1][0]);R_qqbar.append(R[2][0]);R_qbarq.append(R[3][0])
    # plt.plot(Y,R_qg,label=r'qG/tot')
    # plt.plot(Y,R_gq,label=r'Gq/tot')
    # plt.plot(Y,R_qqbar,label=r'qqbar/tot')
    # plt.plot(Y,R_qbarq,label=r'qbarq/tot')
    # plt.axhline(y=1, color='grey', alpha=0.3)
    # plot_usuals(n=2)
    # plt.xlabel(r'$y$',fontsize=f_size)
    # plt.ylabel(r'$R_i = \sigma_i/\sigma_{tot}$')
    # plt.title(r'$\gamma^\star$ production for $M=$'+str(m)+r' GeV and $p_\perp = $'+str(p_T)+' GeV')
    # plt.tight_layout()
    # plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_component_virtual_ratios_'+str(rs)+'GeV_'+convention+'_p_T'+str(p_T)+'GeV_M'+str(m)+'GeV'+'.pdf'))
    # plt.show()
    # plt.close()

    ### Ratios ###
    # FCEL/G
    # y dependant
    err ='q0,mu,pdf'
    # err = 'mu'
    # # err=''
    # a=0
    # if m == 0:
    #     f_name = proton+'_True_RpA_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_'+str(p_T)+'GeV.txt'
    # else:
    #     f_name = proton+'_True_RpA_M'+str(m)+'GeV_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_'+str(p_T)+'GeV.txt'
    # # f_name = proton+'_RpA_wo_iso_M_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_'+str(p_T)+'GeV.txt'
    # if os.path.exists(os.path.join(RpA_dir,f_name)):
    #     print(f"The file '{f_name}' already exists. It is loaded.")
    #     Y,Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
    # else:
    #     print("The file does not exists")
    #     Y,Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_True_RpA_dy_M(x_T,m,switch = convention,var_err= err)
    #     # Y, Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_wo_iso_dy_M(x_T,m,switch=convention,var_err=err)
    #     Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
    #     Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
    #     r = [Y,Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf]
    #     np.savetxt(os.path.join(RpA_dir,f_name), r)
    #     print(f" '{f_name}' has been created")

    # print(Rpa)
    # print(Rpa_minus)
    # print(Rpa_plus)
    # m_mu,M_mu = sig.min_max([Rpa+Rpa_plus_mu,Rpa-Rpa_minus_mu,Rpa])
    # m_q,M_q = sig.min_max([Rpa+Rpa_plus_q,Rpa-Rpa_minus_q,Rpa])
    # m_pdf, M_pdf = sig.min_max([Rpa+Rpa_plus_pdf,Rpa-Rpa_minus_pdf])

    # fig, ax = plt.subplots(constrained_layout=True,figsize=fig_size)
    # ax.xaxis.set_minor_locator(MultipleLocator(1))
    # # ax.xaxis.set_major_locator(MultipleLocator(1))
    # ax.xaxis.set_major_locator(MultipleLocator(2))
    # ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    # ax.yaxis.set_major_locator(MultipleLocator(0.1))
    # plt.fill_between(Y,m_mu,M_mu,color = 'blue',alpha = alph,label=r'$\sigma_\mu $')
    # plt.fill_between(Y,m_q,M_q,color='red', alpha=alph,label=r'$\sigma_{\hat{q}_0}$')
    # plt.fill_between(Y,m_pdf,M_pdf,color='green', alpha=alph,label=r'$\sigma_{pdf}$')
    # plt.plot(Y,Rpa+Rpa_plus,color = 'blue', alpha = 0.5)
    # plt.plot(Y,Rpa-Rpa_minus,color = 'blue', alpha = 0.5)
    # plt.axhline(y=1, color='grey',alpha=alph)
    # plt.plot(Y,Rpa,color='blue')

    # plt.fill_between(Y,Rpa-Rpa_minus,Rpa+Rpa_plus,alpha=alph,label= r'$M =$ '+str(m)+' GeV')

    # plt.legend(loc='lower right',frameon =False,fontsize=f_size-a)
    # plt.xlabel('y',fontsize=f_size-a)
    # plt.ylabel(r'$R_{\text{pA}}$',fontsize= f_size)
    # # plt.ylabel(r'$\tilde{R}_{\text{pA}}$',fontsize= f_size)
    # # plt.ylim(0.8,1.1)
    # plt.ylim(0.7,1.1)
    # plt.xlim(-4,4)
    # # plt.xlim(-3,3)
    # plot_usuals(s1=f_size,s2=f_size,loca = 'lower right')
    # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    # plt.text(0.8, 0.9, r'p'+atom, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
    # # plt.text(0.1, 0.9, r'$\sqrt{s} =$'+str(rs) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
    # plt.text(0.1, 0.9, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
    # plt.text(0.1, 0.8, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
    # plt.text(0.1, 0.7, r'$M =$'+str(m) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
    # plt.tight_layout()
    # plt.savefig(os.path.join(plots_dir, proton+'True_RpA_M_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
    # plt.savefig(os.path.join(plots_dir, proton+'RIpA_M_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
    # plt.show()

# plt.legend(loc='lower right',frameon =False,fontsize=f_size-a)
# plt.xlabel('y',fontsize=f_size-a)
# plt.ylabel(r'$R_{\text{pA}}$',fontsize= f_size)
# # plt.ylabel(r'$\tilde{R}_{\text{pA}}$',fontsize= f_size)
# # plt.ylim(0.8,1.1)
# plt.ylim(0.7,1.1)
# # plt.xlim(-4,4)
# plot_usuals(s1=f_size,s2=f_size,loca = 'lower right')
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.text(0.8, 0.9, r'p'+atom, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-a)
# # plt.text(0.1, 0.9, r'$\sqrt{s} =$'+str(rs) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
# plt.text(0.1, 0.9, r'$\sqrt{s} =$'+str(rs/1000) +r' TeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
# plt.text(0.1, 0.8, r'$p_\bot =$'+str(p_T) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
# # plt.text(0.1, 0.7, r'$M =$'+str(m) +r' GeV', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, proton+'True_RpA_M_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
# plt.show()
# plt.close()

#Mass depedancy test
y = -5
err = ''

f_name = 'Mass_test_repaired.txt'

if os.path.exists(os.path.join(RpA_dir,f_name)):
    print(f"The file '{f_name}' already exists. It is loaded.")
    Y,Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
else:
    print("The file does not exists")
    M_list,Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_True_RpA_dM(y,x_T,switch = convention,var_err= err)
    Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
    Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
    r = [M_list,Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf]
    np.savetxt(os.path.join(RpA_dir,f_name), r)
    print(f" '{f_name}' has been created")

plt.plot(M_list,Rpa)
plt.show()