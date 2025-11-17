# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity with the sigma class.
# Zoom on the integrand of FCEL and FCEG.
# Last modified: 30/10/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import matplotlib.pyplot as plt
from src import sigma as sig

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

p_T = 30 	# GeV
x_T = 2.*p_T/rs
xi =0.5
y = 0

mu2 = p_T**2
mu_f2= mu2
num = 0																	# The central set 

Nu = pPb_cross_section.Nu_list(y,xi,x_T,num,mu2)
alpha_s = pPb_cross_section.alpha_s_p(num,mu2)
chi = pPb_cross_section.FCEL_integrand(y,x_T,num,switch = convention)[3](xi)
sigma_hat = pPb_cross_section.FCEL_integrand(y,x_T,num,switch = convention)[4](xi)

print('alpha_s = '+str(alpha_s)+', chi = '+str(chi)+ ', sigma_hat = '+str(sigma_hat))

jacobian = [pPb_cross_section.FCEL_integrand(y,x_T,mu2,mu_f2,num,switch = convention)[0](nu,xi) for nu in Nu[0]]
P_FCEL = [pPb_cross_section.FCEL_integrand(y,x_T,mu2,mu_f2,num,switch = convention)[1](nu,xi) for nu in Nu[0]]
sigma_FCEL= [pPb_cross_section.FCEL_integrand(y,x_T,mu2,mu_f2,num,switch = convention)[2](nu,xi) for nu in Nu[0]]
P_FCEG = [pPb_cross_section.FCEG_integrand(y,x_T,mu2,mu_f2,num,switch = convention)[1](nu,xi) for nu in Nu[1]]
sigma_FCEG= [pPb_cross_section.FCEG_integrand(y,x_T,mu2,mu_f2,num,switch = convention)[2](nu,xi) for nu in Nu[1]]


# ~ plt.plot(Nu[0],FCEL,label= r'FCEL')
# ~ plt.plot(Nu[1],FCEG,label = r'FCEG')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
# ~ plt.xlabel(r'$\nu = \ln(u)$')
# ~ plt.ylabel(r'cross section integrand')
# ~ #plt.yscale('log')
# ~ plt.title(r'FCEL/G integrand for $y=$ '+str(y)+r', $p_\bot=$ '+str(p_T)+r' GeV, $\xi = $ '+str(xi) )
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'FCELG_integrand'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()

# to verify that the abscisse vector is in ascending order
# ~ plt.plot(Nu[0])
# ~ plt.grid()
# ~ plt.savefig(os.path.join(plots_dir, 'FCELG_integrand_abscisse_test.pdf'))
# ~ plt.close()

# Seperatly


plt.plot(Nu[0],P_FCEL,label= r'FCEL')
plt.plot(Nu[1],P_FCEG,label = r'FCEG')
plt.grid()
plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
plt.xlabel(r'$\nu = \ln(u)$')
plt.ylabel(r'Proba part')
#plt.yscale('log')
plt.title(r'$\hat{P}(e^\nu)$ for $y=$ '+str(y)+r', $p_\bot=$ '+str(p_T)+r' GeV, $\xi = $ '+str(xi) )
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'FCELG_P'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
plt.show()
plt.close()


plt.plot(Nu[0],sigma_FCEL,label= r'FCEL')
plt.plot(Nu[1],sigma_FCEG,label = r'FCEG')
plt.grid()
plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
plt.xlabel(r'$\nu = \ln(u)$')
plt.ylabel(r'$d\sigma_{pp}$ part')
#plt.yscale('log')
plt.title(r'$d\sigma_{pp}(y\pm \ln(1+\hat{\sigma}e^\nu)$ for $y=$ '+str(y)+r', $p_\bot=$ '+str(p_T)+r' GeV, $\xi = $ '+str(xi) )
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'FCELG_sigmapp'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
plt.show()
plt.close()


plt.plot(Nu[0],jacobian,label= r'FCEL')
plt.grid()
plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
plt.xlabel(r'$\nu = \ln(u)$')
plt.ylabel(r'jacobian part')
#plt.yscale('log')
plt.title(r'$\frac{1}{\hat{\sigma}+e^{-\nu}}$ for $y=$ '+str(y)+r', $p_\bot=$ '+str(p_T)+r' GeV, $\xi = $ '+str(xi) )
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'FCELG_jacobian'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
plt.show()
plt.close()
