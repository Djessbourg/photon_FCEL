# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity with the sigma class
# Here is a special analysis for mu2 and mu_f2 variations on the cross section
# Last modified: 30/10/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from src import sigma as sig
from src import Collision
from scipy.integrate import quad

### Initialisation of a sigma object ###
rs = 8800.
s = (rs)**2 															# CM energy in Gev2
p = "nCTEQ15npFullNuc_1_1"
Pb = "nCTEQ15npFullNuc_208_82"
Z = 82.
A = 208.

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

p_T = 50. 																# GeV
x_T = 2.*p_T/rs
p_T2 = p_T**2
y = 0
num= 0 																# The central set 
N_limit = sig.N_limit

# A plot for sigma values depnding of mu2 and mu_f2
# ~ Mu2 = np.linspace(p_T2/4.,4*p_T2,15)
# ~ Mu_f2 = np.linspace(p_T2/4.,4*p_T2,15)
# ~ X, Y = np.meshgrid(Mu2, Mu_f2)
Mu_factor = np.linspace(1/4.,2.,50)
Mu_f_factor = np.linspace(1/4.,2.,50)
X, Y = np.meshgrid(Mu_factor, Mu_f_factor) 
Z = np.zeros(X.shape)													# for X and Y of same shapes 

convention = 'dp_t'
col_type = 'pp'
d = sig.Switch

def integrated_y(y_min,y_max,mu_factor=1.,mu_f_factor=1.):
	y_mean = (y_max+y_min)/2.
	P_T = pPb_cross_section.P_T_list(y_mean,p_t_min=7.5,p_t_max=200.)
	sigma = []
	for p_T in P_T:
		print('p_T = '+str(p_T))
		sigma.append(quad(lambda y:pPb_cross_section.dsigma_tot_dydpt(y,2.*p_T/rs,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,is_pp=Collision(col_type),switch = convention)[0],y_min,y_max,limit=N_limit)[0])
	return (P_T,sigma)

# ~ sig_tot = lambda mu2,mu_f2: pPb_cross_section.sigma_tot(y,x_T,mu2,mu_f2,num,is_pp=Collision(col_type),switch =convention)
sig_tot = lambda mu_factor,mu_f_factor: pPb_cross_section.dsigma_tot_dydpt(y,x_T,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,is_pp=Collision(col_type),switch =convention)

y_min, y_max = -0.67,0.67
# ~ sig_tot = lambda mu_factor,mu_f_factor: quad(lambda y:pPb_cross_section.dsigma_tot_dydpt(y,x_T,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,is_pp=Collision(col_type),switch = convention)[0],y_min,y_max,limit=N_limit)

for i in range(0,X[0].size -1): 										# loop on Y = mu_f2
	for j in range(0,X[0].size-1):										# loop on X = mu2
		Z[i][j] = sig_tot(X[i][j],Y[i][j])[0]							# only the value for[0]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Tracé de la nappe
ax.plot_surface(X, Y, Z, cmap='viridis')  # 'viridis' est un colormap

# Labels des axes
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3, 3))  									# Force la notation scientifique pour des valeurs petites
ax.zaxis.set_major_formatter(formatter)
ax.zaxis.offsetText.set_fontsize(10)
ax.set_xlabel(r'$a, \ \mu = a \times p_\bot$')
ax.set_ylabel(r'$b, \ \mu_f = b \times p_\bot$')
ax.set_zlabel(r'$d\sigma/dy$'+d[convention][0] +'(' +d[convention][1]+')')
plt.tight_layout()
# ~ plt.title(r'nappe for $p_\bot = $'+str(p_T)+ r' $GeV$ and $y = $'+str(y)+r' at $\sqrt{s} = $' +str(rs) +r'$GeV$')

# Affichage
plt.savefig(os.path.join(plots_dir, 'nappe_'+convention+str(p_T)+'GeV_y'+str(y)+'.pdf'))
plt.show()
# A second part with mu = mu_f, for y fixed, at pt 
