# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity with the sigma class.
# Plots of the different components of the cross section before FCEL/G effect.
# Last modified: 14/11/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from src import sigma as sig
from src import Collision

### Initialisation of a sigma object ###
rs = 8800
s = (rs)**2 # CM energy in Gev2
p = "nCTEQ15npFullNuc_1_1"
Pb = "nCTEQ15npFullNuc_208_82"
Z = 82.
A = 208.

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)


y = 4

mu2 = 1.																# you can put a random value when cteq = True, bc it's redefine for each p_T
mu_f2 = mu2																# The virtuality variables for alpha_s and the pdfs
num = 0 																# The central set 

d = sig.Switch
convention = 'd2p_t'													# Define a cross section convention btw dp_t, dp_t2 ans d2p_t
col_type = 'pPb'

### cross section pt_dependant with mu = mu_f = pt
sigma_tot,err_tot= pPb_cross_section.dsigma_tot_dpt(y,num,is_pp= Collision(col_type),switch = convention)
sigma_qG, err_qG  = pPb_cross_section.dsigma_qG_dpt(y,num,is_pp=Collision(col_type),switch = convention)
sigma_Gq, err_Gq  = pPb_cross_section.dsigma_Gq_dpt(y,num,is_pp=Collision(col_type),switch = convention)
sigma_qqbar, err_qqbar = pPb_cross_section.dsigma_qqbar_dpt(y,num,is_pp=Collision(col_type),switch = convention)
sigma_qbarq, err_qbarq = pPb_cross_section.dsigma_qbarq_dpt(y,num,is_pp=Collision(col_type),switch = convention)
P_T = pPb_cross_section.P_T_list(y)

plt.plot(P_T,sigma_tot,label=r'$d\sigma_{tot}/dy$'+d[convention][0])
#plt.errorbar(P_T,sigma_tot,yerr=err_tot,fmt='+')
plt.plot(P_T,sigma_qG,label=r'$d\sigma_{qG}/dy$'+d[convention][0])
#plt.errorbar(P_T,sigma_qG,yerr=err_qG,fmt='+')
plt.plot(P_T,sigma_Gq,label=r'$d\sigma_{Gq}/dy$'+d[convention][0])
#plt.errorbar(P_T,sigma_Gq,yerr=err_Gq,fmt='+')
plt.plot(P_T,sigma_qqbar,label=r'$d\sigma_{qqbar}/dy$'+d[convention][0])
#plt.errorbar(P_T,sigma_qqbar,yerr=err_qqbar,fmt='+')
plt.plot(P_T,sigma_qbarq,label=r'$d\sigma_{qbarq}/dy$'+d[convention][0])
#plt.errorbar(P_T,sigma_qbarq,yerr=err_qbarq,fmt='+')
plt.grid()
plt.legend(loc='upper right',frameon =False,ncol=2)
plt.xlabel(r'$p_\bot$ in GeV')
plt.ylabel(r'cross section $d\sigma/dy$'+d[convention][0]+' in '+d[convention][1])
plt.yscale('log')
plt.title(r'Cross sections for '+col_type+', with $y = $'+str(y)+r' at $\sqrt{s} = $ '+ str(rs)+r'$GeV$')
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_all_components_'+str(rs)+'GeV_'+convention+'_y'+str(y)+'.pdf'))
plt.show()

# Another one for ratios
R_all= pPb_cross_section.R_composition_dpt(y,num,is_pp = Collision(col_type),switch = convention,var='p')

plt.close()

plt.plot(P_T,R_all[0][0],label=r'qG/tot')
plt.plot(P_T,R_all[1][0],label=r'Gq/tot')
plt.plot(P_T,R_all[2][0],label=r'qqbar/tot')
plt.plot(P_T,R_all[3][0],label=r'qbarq/tot')
plt.axhline(y=1, color='grey', linestyle='--')

plt.grid()
plt.legend(loc='upper right',frameon =False,ncol =2)
plt.xlabel(r'$p_\bot$ in GeV')
plt.ylabel(r'Ratios of different cross sections for '+col_type)
#plt.yscale('log')
plt.title(r'Component ratios, with $y =$'+str(y))
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_component_ratios_'+str(rs)+'GeV_'+convention+'_y'+str(y)+'.pdf'))
plt.show()
plt.close()

