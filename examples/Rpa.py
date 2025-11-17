# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity with the sigma class.
# Here some plots comparing the cross sections from FCEL and FCEG, the shift
# effect and the ratios of it.
# Last modified: 13/11/2025
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
a = 2
b= 3

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

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

d = sig.Switch
convention = 'd2p_t'													# Define a cross section convention btw dp_t, dp_t2 ans d2p_t
col_type = 'pA'

### y dependant ###
p_T = 15																# GeV
x_T = (2.*p_T)/rs 
# ~ xi = 0.5
var_int = 'nu'

Y = sig.Y_list(x_T)
# ~ Y_xi = sig.Y_list_xi(x_T,xi)

mu2 = p_T**2
mu_f2= mu2
num = 0																	# The central set 

# FCEL, FCEG = pPb_cross_section.FCEL_G_integration_dy(x_T,num,switch = convention,var_int= var_int)

# sigma_FCEL , err_FCEL = FCEL[0], FCEL[1]
# sigma_FCEG, err_FCEG = FCEG[0],FCEG[1]

# sigma_tot_FCELG = np.add(sigma_FCEL,sigma_FCEG)
# err_tot_FCELG = np.sqrt(np.add(np.multiply(err_FCEL,err_FCEL),np.multiply(err_FCEG,err_FCEG)))

FCEL_tot, FCEG_tot = pPb_cross_section.dsigma_FCELG_dy(x_T,num,is_pp=True,switch = convention)
sigma_tot = pPb_cross_section.dsigma_tot_dy(x_T,num,is_pp=True,switch = convention)

# R_tot = np.divide(sigma_tot_FCELG,sigma_tot[0])
# R_FCEL = np.divide(sigma_FCEL,FCEL_tot[0])
# R_FCEG = np.divide(sigma_FCEG,FCEG_tot[0])
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 3))

fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_major_formatter(formatter)
# plt.plot(Y,sigma_FCEG,'r',label = r'$d\sigma$ FCEG')
# plt.plot(Y,sigma_FCEL,'b',label = r'$d\sigma$ FCEL')
plt.plot(Y,FCEG_tot[0],'r',label = r'$d\sigma_{tot}$ "FCEG"')
plt.plot(Y,FCEL_tot[0],'b',label = r'$d\sigma_{tot}$ "FCEL"')
# plt.plot(Y,sigma_tot_FCELG,'grey',label = r'$d\sigma$ tot')
plt.plot(Y,sigma_tot[0],'purple',label=r'$d\sigma_{pp}$')
plt.xlabel('y',fontsize= f_size-a)
plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')',fontsize= f_size-a)
plot_usuals(s1=f_size-b,s2=f_size-a,loca = 'lower center')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.tight_layout()
plt.text(0.05, 0.95, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-b)
plt.savefig(os.path.join(plots_dir, 'Sigma_FCEL_FCEG_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
plt.show()
plt.close()


# ax = plt.subplot()
# plt.plot(Y,R_FCEL,'b',label=r'$R^{FCEL}$')
# plt.plot(Y,R_FCEG,'r',label=r'$R^{FCEG}$')
# plt.plot(Y,R_tot,'purple',label = r'$R^{pA}$')
# plt.axhline(y=1, color='grey', linestyle='--')
# plt.grid(False)
# plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=10)
# plt.xlabel('rapidity y')
# plt.ylabel(r'ratios')
# #plt.yscale('log')
# plt.ylim(0.7,1.2)
# # ~ plt.title(r'FCEL and FCEG ratios at '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')
# plt.text(0.1, 0.1, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# plt.tight_layout()
# plt.savefig(os.path.join(plots_dir, 'Sigma_FCELG_ratios_'+col_type+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# plt.show()
# plt.close()

