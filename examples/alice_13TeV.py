# -*- coding: utf-8 -*-

# =============================================================================
# Comparison of the cross
# section with Alice 13 TeV
# data:
# https://www.hepdata.net/record/ins2803487
# Last modified: 29/10/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import csv
import numpy as np
import matplotlib.pyplot as plt
from src import sigma as sig
from src import Collision
from scipy.integrate import quad

### Initialisation of a sigma object ###
rs = 13000
s = (rs)**2 # CM energy in Gev2
p = "nNNPDF30_nlo_as_0118_p"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"
Z = 82.
A = 208.

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

d = sig.Switch

y_range = [-0.67,0.67]
bin_width = y_range[1]-y_range[0]
num_cen = 0 																# The central set 

# A plot for sigma tot and its components
convention = 'dp_t'
col_type = 'pp'
N_limit = sig.N_limit
err = 'mu,pdf'															# variables to compute error from 

def integrated_yp_T(y_min,y_max,p_T_min,p_T_max,mu_factor=1.,mu_f_factor=1.):
	sigma=quad(lambda p_T:quad(lambda y:pPb_cross_section.dsigma_tot_dydpt(y,2.*p_T/rs,num_cen,mu_factor=mu_factor,mu_f_factor=mu_f_factor,is_pp=Collision(col_type),switch = convention)[0],y_min,y_max,limit=N_limit)[0],p_T_min,p_T_max,limit=N_limit)[0]
	return sigma
	
with open('HEPData-ins2803487-v1-table_1.csv', 'r',newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    rows = list(reader)

# find the line index 
for i, row in enumerate(rows):
    if '$p_' in row[0]:  # '$p_T'in line
        header_row_index = i
        break

# Read data from the net line
data_rows = rows[header_row_index + 1:]

pt_center =[]
pt_low = []
pt_high = []
cross_section = []
total_err = []
for row in data_rows:
    try:
        pt_center.append(float(row[0]))
        pt_low.append(float(row[1]))
        pt_high.append(float(row[2]))
        cross_section.append(float(row[3])*10**3) # It's in nanobarne per Gev !!! (that's why we *10^3) 
        stat_plus = float(row[4])
        stat_minus = float(row[5])
        sys_plus = float(row[6])
        sys_minus = float(row[7])

        # Mean of uncertainites up/down
        stat_err = (np.abs(stat_plus) + np.abs(stat_minus)) / 2.0
        sys_err = (np.abs(sys_plus) + np.abs(sys_minus)) / 2.0

        # Total quadratic uncertainty
        total_err.append(np.sqrt(stat_err**2 + sys_err**2))

    except:
        continue  # ingore incorect or empty lines

# Uncomment the following lines to create your own matrix to save (be aware not to transpose or nothing, unlike me)

# P_T = pPb_cross_section.P_T_list(0,p_t_min=7.5,p_t_max=15.)
# sigma = pPb_cross_section.dsigma_tot_dpt(0,num_cen,is_pp=Collision(col_type),switch = convention)[0]

def sigmai(mu_factor,mu_f_factor):
	sigma = []
	for i in range(len(pt_low)):
		sigma.append(integrated_yp_T(y_range[0],y_range[1],pt_low[i],pt_high[i],mu_factor=mu_factor,mu_f_factor=mu_f_factor))
	return sigma
	
# sigma = sigmai(1.,1.)
# sigma1 = sigmai(2.,2.)
# sigma2 = sigmai(0.5,2.)
# sigma3 = sigmai(2.,0.5)
# sigma4 = sigmai(0.5,0.5)
# (sigma_min,sigma_max) = sig.min_max([sigma1,sigma2,sigma3,sigma4])

# Sigma_matrix = np.array([pt_center,pt_low,pt_high,sigma,sigma_min,sigma_max])
# np.savetxt('Sigma_matrix_bin.txt',Sigma_matrix)

Mat = np.transpose(np.loadtxt('Sigma_matrix.txt'))						# I did a transposition for nothing
Mat2 = np.loadtxt('Sigma_matrix_bin.txt')

P_T = Mat[0]															# Integrated in y
sigma = Mat[1]
sigma_min = Mat[2]
sigma_max = Mat[3]

# ~ sigma_0 = [pPb_cross_section.dsigma_tot_dydpt(0,2.*p_T/rs,num_cen,mu_factor=1,mu_f_factor=1,is_pp=Collision(col_type),switch = convention)[0] for p_T in P_T] #y = 0 

P_T_Alice = Mat2[0]
P_T_min = Mat2[1]
P_T_max = Mat2[2]
sigma_bin = Mat2[3]
sigma_bin_min = Mat2[4]
sigma_bin_max = Mat2[5]

pt_bin_width= np.subtract(P_T_max,P_T_min)

sigma_bin = np.divide(sigma_bin,pt_bin_width)/bin_width
sigma_bin_min = np.divide(sigma_bin_min,pt_bin_width)/bin_width
sigma_bin_max = np.divide(sigma_bin_max,pt_bin_width)/bin_width

P_T_min = np.subtract(P_T_Alice,P_T_min)
P_T_max = np.subtract(P_T_max,P_T_Alice)

ax = plt.subplot()
plt.fill_between(P_T,np.divide(sigma_min,bin_width),np.divide(sigma_max,bin_width),color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
plt.plot(P_T,np.divide(sigma,bin_width),color='b',label=r'Python program (y integrated)')
plt.errorbar(P_T_Alice,sigma_bin,yerr=[sigma_bin-sigma_bin_min,sigma_bin_max-sigma_bin],xerr = [P_T_min,P_T_max],fmt='none',label = r'Python program (all integrated)')
plt.errorbar(pt_center,cross_section,yerr=total_err,xerr=[np.subtract(pt_center,pt_low),np.subtract(pt_high,pt_center)],fmt='.',label=r'Alice 13 TeV')
# ~ plt.plot(P_T,sigma_0,color = 'green',label =r'Python program for $y=0$')
#plt.grid()
plt.legend(frameon =False,fontsize=10)
plt.xlabel(r'$p_\bot$ ($GeV$)')
# ~ plt.ylim(bottom=0)
plt.xlim(6.,210.)
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
# ~ plt.text(0.1, 0.1, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.text(0.1, 0.05, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_Alice_'+str(rs)+'GeV.pdf'))
plt.show()
plt.close()

