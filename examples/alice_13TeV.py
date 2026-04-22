# -*- coding: utf-8 -*-

# =============================================================================
# Comparison of the cross
# section with Alice 13 TeV
# data:
# https://www.hepdata.net/record/ins2803487
# other section with Alice 8.16 TeV
# Last modified: 29/10/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from src import sigma as sig
from src import Collision
from scipy.integrate import quad

### Initialisation of a sigma object ###
# rs = 13000
# rs = 8800
rs = 8160
s = (rs)**2 # CM energy in Gev2
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

# y_range = [-0.67,0.67]
y_range = [-0.7,0.7]
bin_width = y_range[1]-y_range[0]
num_cen = 0 																# The central set 

# A plot for sigma tot and its components
convention = 'dp_t'
col_type = 'pp'
N_limit = sig.N_limit
err = 'mu,pdf'															# variables to compute error from 

fig_size= (5,5)
f_size= 17
a=0
b=0								

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)
      
def integrated_yp_T(y_min,y_max,p_T_min,p_T_max,mu_factor=1.,mu_f_factor=1.):
	sigma=quad(lambda p_T:quad(lambda y:pPb_cross_section.dsigma_tot_dydpt(y,2.*p_T/rs,num_cen,mu_factor=mu_factor,mu_f_factor=mu_f_factor,is_pp=Collision(col_type),switch = convention)[0],y_min,y_max,limit=N_limit)[0],p_T_min,p_T_max,limit=N_limit)[0]
	return sigma

# Do i need to do a integrated version of R_pA ?
n_13TeV = 'HEPData-ins2803487-v1-table_1.csv'
n_8_16TeV = 'HEPData-ins2895564-v1-RpA_8.16_TeV.csv' 
	
with open(os.path.join(data_dir,n_8_16TeV), 'r',newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    rows = list(reader)

# find the line index 
for i, row in enumerate(rows):
    if 'PT' in row[0]:  # '$p_T'in line
        header_row_index = i
        break

# Read data from the net line
data_rows = rows[header_row_index + 1:]

pt_center =[]
pt_low = []
pt_high = []
R_pA = []
total_err = []
for row in data_rows:
    try:
        pt_center.append(float(row[0]))
        pt_low.append(float(row[1]))
        pt_high.append(float(row[2]))
        R_pA.append(float(row[3])) # It's in nanobarne per Gev !!! (that's why we *10^3) 
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


# pt dependant R_pA
y=0
err= '' # there is barely no mu error in it 
RpA_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_dir'))
# f_name = proton+'_RpA_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_y'+str(y)+'.txt'
f_name = proton+'_True_RpA_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_y'+str(y)+'.txt'
# f_name = proton+'_RpA_wo_iso_'+str(rs)+'GeV_Z'+str(Z)+'_A'+str(A)+'_'+str(p_T)+'GeV.txt'
if os.path.exists(os.path.join(RpA_dir,f_name)):
	print(f"The file '{f_name}' already exists. It is loaded.")
	P_T,Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
	# Rpa ,Rpa_plus,Rpa_minus, Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf = np.loadtxt(os.path.join(RpA_dir,f_name))
	# Y = sig.Y_list(x_T,Z)
	# r = [Y,Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf]
	# # np.savetxt(os.path.join(RpA_dir,proton+'_'+f_name), r)
	# np.savetxt(os.path.join(RpA_dir,f_name), r)
else:
	print("The file does not exists")
	# P_T,Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_dpt(y,switch = convention,var_err= err)
	P_T,Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_True_RpA_dpt(y,switch = convention,var_err= err)
	# Rpa ,err_Rpa, err_var_Rpa = pPb_cross_section.Uncertainties_RpA_wo_iso_dy(x_T,switch=convention,var_err=err)
	Rpa_plus,Rpa_minus = err_Rpa[0],err_Rpa[1]
	Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf= err_var_Rpa[0],err_var_Rpa[1],err_var_Rpa[2],err_var_Rpa[3],err_var_Rpa[4],err_var_Rpa[5]
	r = [P_T,Rpa,Rpa_plus,Rpa_minus,Rpa_plus_q,Rpa_minus_q,Rpa_plus_mu,Rpa_minus_mu,Rpa_plus_pdf,Rpa_minus_pdf]
	np.savetxt(os.path.join(RpA_dir,f_name), r)
	print(f" '{f_name}' has been created")


#correction of integration with the bin size in the first approximation
# Rpa ,Rpa_plus,Rpa_minus = bin_width*Rpa,bin_width*Rpa_plus,bin_width*Rpa_minus

# Uncomment the following lines to create your own matrix to save (be aware not to transpose or nothing, unlike me)

# P_T = pPb_cross_section.P_T_list(0,p_t_min=7.5,p_t_max=15.)
# sigma = pPb_cross_section.dsigma_tot_dpt(0,num_cen,is_pp=Collision(col_type),switch = convention)[0]

# def sigmai(mu_factor,mu_f_factor):
# 	sigma = []
# 	for i in range(len(pt_low)):
# 		sigma.append(integrated_yp_T(y_range[0],y_range[1],pt_low[i],pt_high[i],mu_factor=mu_factor,mu_f_factor=mu_f_factor))
# 	return sigma
	
# sigma = sigmai(1.,1.)
# sigma1 = sigmai(2.,2.)
# sigma2 = sigmai(0.5,2.)
# sigma3 = sigmai(2.,0.5)
# sigma4 = sigmai(0.5,0.5)
# (sigma_min,sigma_max) = sig.min_max([sigma1,sigma2,sigma3,sigma4])

# Sigma_matrix = np.array([pt_center,pt_low,pt_high,sigma,sigma_min,sigma_max])
# np.savetxt('Sigma_matrix_bin.txt',Sigma_matrix)

# Mat = np.transpose(np.loadtxt('Sigma_matrix.txt'))						# I did a transposition for nothing
# Mat2 = np.loadtxt('Sigma_matrix_bin.txt')
# P_T = Mat[0]															# Integrated in y
# sigma = Mat[1]
# sigma_min = Mat[2]
# sigma_max = Mat[3]

# ~ sigma_0 = [pPb_cross_section.dsigma_tot_dydpt(0,2.*p_T/rs,num_cen,mu_factor=1,mu_f_factor=1,is_pp=Collision(col_type),switch = convention)[0] for p_T in P_T] #y = 0 

# P_T_Alice = Mat2[0]
# P_T_min = Mat2[1]
# P_T_max = Mat2[2]
# sigma_bin = Mat2[3]
# sigma_bin_min = Mat2[4]
# sigma_bin_max = Mat2[5]

# pt_bin_width= np.subtract(P_T_max,P_T_min)

# sigma_bin = np.divide(sigma_bin,pt_bin_width)/bin_width
# sigma_bin_min = np.divide(sigma_bin_min,pt_bin_width)/bin_width
# sigma_bin_max = np.divide(sigma_bin_max,pt_bin_width)/bin_width

# P_T_min = np.subtract(P_T_Alice,P_T_min)
# P_T_max = np.subtract(P_T_max,P_T_Alice)

fig, ax = plt.subplots(constrained_layout=True,figsize=fig_size)
# ax.xaxis.ticklabel_format(useMathText=True)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
# plt.fill_between(P_T,np.divide(sigma_min,bin_width),np.divide(sigma_max,bin_width),color = 'blue', alpha= 0.3 ,label=r'$\sigma_\mu$')
# plt.plot(P_T,np.divide(sigma,bin_width),color='b',label=r'Python program (y integrated)')
# plt.errorbar(P_T_Alice,sigma_bin,yerr=[sigma_bin-sigma_bin_min,sigma_bin_max-sigma_bin],xerr = [P_T_min,P_T_max],fmt='none',label = r'Python program (all integrated)')
plt.plot(P_T,Rpa,label = r'Our $R_{\text{pA}}$')
plt.fill_between(P_T,Rpa-Rpa_minus,Rpa+Rpa_plus,alpha=0.3,color='blue')
plt.errorbar(pt_center,R_pA,yerr=total_err,xerr=[np.subtract(pt_center,pt_low),np.subtract(pt_high,pt_center)],fmt='.',label=r'Alice pPb 8.16 TeV')
plt.axhline(y=1,color='grey',alpha=0.5)
# ~ plt.plot(P_T,sigma_0,color = 'green',label =r'Python program for $y=0$')
#plt.grid()
plot_usuals(loca='upper left')
# plt.legend(frameon =False,fontsize=10,loc='upper left')
plt.xlabel(r'$p_\bot$ GeV',fontsize=f_size)
# ~ plt.ylim(bottom=0)
# plt.xlim(5.,16.)
# plt.ylim(0.95,1.1)
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r'('+d[convention][1]+r')')
plt.ylabel(r'$R_{\text{pA}}$',fontsize=f_size)
plt.text(0.9, 0.1, r'$|y| \leq $'+ str(y_range[1]), horizontalalignment='right', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
plt.text(0.9, 0.2, r'$\sqrt{s} = $'+ str(rs/1000)+' TeV', horizontalalignment='right', verticalalignment='center',transform=ax.transAxes,fontsize=f_size)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# ~ plt.text(0.1, 0.05, r'collision: '+col_type, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'RpA_LHC_'+str(rs)+'GeV.pdf'))
plt.show()
plt.close()



