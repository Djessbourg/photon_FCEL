# -*- coding: utf-8 -*-

# =============================================================================
# Study of the impact of adding a mass to the photon on the cross section
# Last modified: 13/01/2026
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
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
M = [0,2,5,10,50] # vitual photon mass in GeV
cen = 0

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

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

debug_mode = False

y_space = [-6,6]
m = 0 # GeV
b = (m/p_T)**2

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
### integrand test ###

# for y in y_space:
#     epsilon = pow(10,-3)
#     Xi_min,Xi_max = sig.Xi_min_M(y,x_T,b)+epsilon,sig.Xi_max_M(y,x_T,b)-epsilon
#     Xi = np.linspace(Xi_min,Xi_max,100)
#     qg,gq,qqbar,qbarq=[],[],[],[]
#     for xi in Xi:   
#         qg.append(pPb_cross_section.qG_M(y,x_T,xi,m,cen,is_pp=Collision(col_type),switch=convention,debug=debug_mode))
#         gq.append(pPb_cross_section.Gq_M(y,x_T,xi,m,cen,is_pp=Collision(col_type),switch=convention,debug=debug_mode))
#         qqbar.append(pPb_cross_section.qqbar_M(y,x_T,xi,m,cen,is_pp=Collision(col_type),switch=convention,debug=debug_mode))
#         qbarq.append(pPb_cross_section.qbarq_M(y,x_T,xi,m,cen,is_pp=Collision(col_type),switch=convention,debug=debug_mode))
#     plt.plot(Xi,qg,label='qg')
#     plt.plot(Xi,gq,label='gq')
#     plt.plot(Xi,qqbar,label='qqbar')
#     plt.plot(Xi,qbarq,label='qbarq')
#     plt.xlim(0,1)
#     plt.yscale('log')
#     plt.xlabel(r'$\xi$',fontsize=f_size)
#     plt.axvline(x=Xi_min, color='grey', alpha=0.3)
#     plt.axvline(x=Xi_max, color='grey', alpha=0.3)
#     plot_usuals(n=2)
#     plt.ylabel(r'$d\sigma/d\xi$')
#     # plt.title(r'$M=$'+str(m)+r' GeV, $y =$'+str(y)+r', $p_\perp = $'+str(p_T)+' GeV')
#     plt.title(r'DEBUG $M=$'+str(m)+r' GeV, $y =$'+str(y)+r', $p_\perp = $'+str(p_T)+' GeV')
#     plt.tight_layout()
#     # plt.savefig(os.path.join(plots_dir, 'sigma_integrand'+col_type+'_component_virtual_'+str(rs)+'GeV_'+convention+'_y'+str(y)+'_M'+str(m)+'GeV'+'.pdf'))
#     plt.savefig(os.path.join(plots_dir, 'DEBUG_sigma_integrand'+col_type+'_component_virtual_'+str(rs)+'GeV_'+convention+'_y'+str(y)+'_M'+str(m)+'GeV'+'.pdf'))
#     plt.show()
#     plt.close()
### y dependant ###


for m in M:
    b = (m/p_T)**2
    Y = sig.Y_list_M(x_T,b)
    # print((min(Y),max(Y)))
    R_qg, R_gq, R_qqbar, R_qbarq = [],[],[],[]
    for y in Y:
        R = pPb_cross_section.R_composition_dydpt_M(y,x_T,m,cen,is_pp=Collision(col_type),switch=convention,debug=debug_mode)
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
    # plt.title(r'DEBUG $\gamma^\star$ production for $M=$'+str(m)+r' GeV and $p_\perp = $'+str(p_T)+' GeV')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_component_virtual_ratios_'+str(rs)+'GeV_'+convention+'_p_T'+str(p_T)+'GeV_M'+str(m)+'GeV'+'.pdf'))
    # plt.savefig(os.path.join(plots_dir, 'DEBUG_sigma_'+col_type+'_component_virtual_ratios_'+str(rs)+'GeV_'+convention+'_p_T'+str(p_T)+'GeV_M'+str(m)+'GeV'+'.pdf'))
    plt.show()
    plt.close()

### plot pdf ###

# pdf_qg, pdf_gq, pdf_qqbar, pdf_qbarq = [],[],[],[]
# m = 0
# p_t = 5
# b = (m/p_t)**2
# mu_2 = p_t**2 
# x_t = 2*p_t/rs
# Y = sig.Y_list_M(x_t,b)
# f_corr = lambda xi: (b+1)/(b*(1-xi)+1)

# Gx1 = lambda xi,y: pPb_cross_section.Gluon_p(sig.x_1_M(y,x_t,xi,b),mu_2,0)
# F2x1 = lambda xi,y: pPb_cross_section.F2_p(sig.x_1_M(y,x_t,xi,b),mu_2,0)
# Gx2 = lambda xi,y: pPb_cross_section.Gluon_p(sig.x_2_M(y,x_t,xi,b),mu_2,0)
# F2x2 = lambda xi,y: pPb_cross_section.F2_p(sig.x_2_M(y,x_t,xi,b),mu_2,0)
# F_qqbar = lambda xi,y: pPb_cross_section.F_ij(sig.x_1_M(y,x_t,xi,b),sig.x_2_M(y,x_t,xi,b),mu_2,0,is_pp=Collision(col_type))
# F_qbarq = lambda xi,y: pPb_cross_section.F_ij(sig.x_1_M(y,x_t,xi,b),sig.x_2_M(y,x_t,xi,b),mu_2,0,direction='qbarq',is_pp=Collision(col_type))

# N_lim = 1000

# for y in Y:
#     Xi_min, Xi_max = sig.Xi_min_M(y,x_t,b),sig.Xi_max_M(y,x_t,b)
#     pdf_gq.append(quad(lambda xi:Gx1(xi,y)*F2x2(xi,y)*f_corr(xi),Xi_min, Xi_max,limit = N_lim)[0])
#     pdf_qg.append(quad(lambda xi:Gx2(xi,y)*F2x1(xi,y)*f_corr(xi),Xi_min, Xi_max,limit = N_lim)[0])
#     pdf_qqbar.append(quad(lambda xi:F_qqbar(xi,y)*f_corr(xi),Xi_min, Xi_max,limit = N_lim)[0])
#     pdf_qbarq.append(quad(lambda xi:F_qbarq(xi,y)*f_corr(xi),Xi_min, Xi_max,limit = N_lim)[0])

# plt.plot(Y,pdf_qg,label='qg')
# plt.plot(Y,pdf_gq,label='gq')
# plt.plot(Y,pdf_qqbar,label='qqbar')
# plt.plot(Y,pdf_qbarq,label='qbarq')
# plt.axvline(x = 0, linestyle = '--',color= 'grey', alpha = 0.5)
# plot_usuals(n=2)
# plt.xlabel(r'$y$',fontsize=f_size)
# plt.ylabel(r"$\int$ pdf's")
# # plt.title(r'$\gamma^\star$ production for $M=$'+str(m)+r' GeV and $p_\perp = $'+str(p_T)+' GeV')
# plt.title(r'DEBUG $M=$'+str(m)+r' GeV, $p_\perp = $'+str(p_t)+' GeV')
# plt.tight_layout()
# # plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_component_virtual_ratios_'+str(rs)+'GeV_'+convention+'_p_T'+str(p_T)+'GeV_M'+str(m)+'GeV'+'.pdf'))
# plt.savefig(os.path.join(plots_dir, 'DEBUG_pdf_'+col_type+'_component_virtual_ratios_'+str(rs)+'GeV_'+convention+'_p_T'+str(p_t)+'GeV_M'+str(m)+'GeV'+'.pdf'))
# plt.show() 
# plt.close()
