# -*- coding: utf-8 -*-

# =============================================================================
# Finding if our results are
# compatible with the JetPHOX
# calculations for direct
# photons
# Last modified: 30/10/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))


import matplotlib.pyplot as plt
from src import sigma as sig
from src import Collision
from scipy.integrate import quad 

### Initialisation of a sigma object ###
rs = 8800
s = (rs)**2 # CM energy in Gev2
p = "nNNPDF30_nlo_as_0118_p"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"
Z = 82.
A = 208.

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

d = sig.Switch

convention = 'dp_t'
col_type = 'pp'
N_limit = sig.N_limit
num_cen = 0

def integrated_pt(pt_min,pt_max):
	Y = sig.Y_list(2*pt_min/rs)
	sigma = []
	for y in Y:
		print('y = '+str(y))
		sigma.append(quad(lambda x:pPb_cross_section.dsigma_tot_dydpt(y,2.*x/rs,num_cen,is_pp=Collision(col_type),switch = convention)[0],pt_min,pt_max,limit=N_limit)[0])
	return (Y,sigma)
	
pt_min = 3.
pt_max = 15.

pt = 100
xt = 2.*pt/rs
sigma, err_sigma = pPb_cross_section.dsigma_tot_dy(xt,num_cen,is_pp=Collision(col_type),switch = convention)
Y = sig.Y_list(xt)
print(Y)
print(err_sigma)
print(sigma)
ax = plt.subplot()
plt.plot(Y,sigma,color = 'blue' ,label=r'$\sigma_\mu$')
plt.grid()
# ~ plt.legend(loc='lower right',frameon =False,fontsize=10)
plt.xlabel('y')
# ~ plt.ylim(bottom=1,top=1e4)
plt.ylabel(r'$d\sigma/dy$'+ r'($pb$)')
# ~ plt.yscale('log')
# ~ plt.text(0.25, 0.05, r'collision: '+col_type+r' and $p_\bot \in \{3,15\} \ GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_integrated_pt_'+col_type+str(rs)+'GeV_'+convention+'GeV.pdf'))
plt.show()
plt.close()
