# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity with the sigma class.
# Plots of mu, pdf uncertainties.
# Last modified: 13/11/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import numpy as np
import matplotlib.pyplot as plt
from src import sigma as sig

def fill_between_n_curves(x, ys,label, color='gray', alpha=0.3, show_lines=True, labels=None):
	"""
	Fill the area between a number n of curves

	Parameters :
		x (array-like) : Abscisses vector.
		ys (list of array-like) : Curves list.
		label (str): Label of the filling.
		color (str) : Color of filling.
		alpha (float) : Transparency.
		show_lines (bool) : Show the n curves or not.
		labels (list of str) : Curves label (optional) .
	"""
	ys_array = np.vstack(ys)
	y_min = np.min(ys_array, axis=0)
	y_max = np.max(ys_array, axis=0)

	plt.figure(figsize=(10, 6))

	# Optionnal : draw lines
	if show_lines:
		for i, y in enumerate(ys):
			label_n = labels[i] if labels and i < len(labels) else 'courbe'.format(i+1)
			plt.plot(x, y, label=label_n)

	# filling between min and max
	plt.fill_between(x, y_min, y_max, color=color, alpha=alpha, label=label)


### Initialisation of a sigma object ###
rs = 8800
s = (rs)**2 # CM energy in Gev2
p = "nNNPDF30_nlo_as_0118_p"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"
Z = 82.
A = 208.

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

d = sig.Switch

### first atempts to reproduce last results ###
p_T = 5
x_T = (2.0*p_T)/rs														
mu2 = p_T**2
mu_f2 = mu2																# The virtuality variables for alpha_s and the pdfs

# A plot for sigma tot with its uncertainties
alph = 0.3
f_size = 15
convention = 'd2p_t'
col_type = 'pA'

Y = sig.Y_list(x_T)
num= 0
var_int ='nu'
# ~ Ucen, Uplus, Uminus, Usym = pPb_cross_section.Uncertainties_tot_dy(x_T,mu2,mu_f2,is_pp=Collision(col_type),switch = convention)
# ~ print(Ucen)
# ~ print (Uplus)
# ~ print(Uminus)

# ~ plt.plot(Y,Ucen,label = r'$d\sigma_{tot}/dy$'+d[convention][0])
# ~ plt.fill_between(Y,np.add(Ucen,Uplus),np.subtract(Ucen,Uminus),alpha= 0.5,label = r'CL at 1-$\sigma$')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'cross section $d\sigma/dy$'+d[convention][0]+ r' in '+d[convention][1])
# ~ #plt.yscale('log')
# ~ plt.title(r'$d\sigma/dy$'+d[convention][0]+ r'components for pp, with $p_\bot =$'+str(x_T*rs/2.)+r' $GeV$')
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_with_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()

# ~ etaplus = np.divide(Uplus,Ucen)
# ~ etaminus = np.divide(np.subtract(0,Uminus), Ucen)
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.fill_between(Y,etaplus,etaminus,alpha=0.5,label = r'relative CL at 1-$\sigma$')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'relative uncertainties $(Err-d\sigma)/d\sigma$')
# ~ #plt.yscale('log')
# ~ plt.title(r'$(Err-d\sigma)/d\sigma$components for '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_relative_'+col_type+'_with_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()


### mu uncertainties ###
# for y

# As it is coded that mu2 = (a*p_T)**2 ; mu_f2 = (b*p_T)**2, we have 5 possibilites to make: (1,1),(1/2,1/2),(1/2,2),(2,1/2),(2,2)
ab= [(1./2.,1./2.),(1./2.,2.),(2.,1./2.),(2.,2.)]
labels = [r'$\mu = $ '+str(couple[0])+r'$p_\bot \ \mu_f = $ '+str(couple[1])+r'$p_\bot$' for couple in ab]
# for the central value
FCEL_cen, FCEG_cen = pPb_cross_section.FCEL_G_integration_dy(x_T,num,switch = convention,var_int= var_int)

sigma_FCEL_cen , err_FCEL_cen = FCEL_cen[0], FCEL_cen[1]
sigma_FCEG_cen, err_FCEG_cen = FCEG_cen[0],FCEG_cen[1]

sigma_tot_FCELG_cen = np.add(sigma_FCEL_cen,sigma_FCEG_cen)
err_tot_FCELG_cen = np.sqrt(np.add(np.multiply(err_FCEL_cen,err_FCEL_cen),np.multiply(err_FCEG_cen,err_FCEG_cen)))

sigma_tot_cen = pPb_cross_section.dsigma_tot_dy(x_T,num,is_pp=True,switch = convention)

R_tot_cen = np.divide(sigma_tot_FCELG_cen,sigma_tot_cen[0])

sigma_FCELG = []
R = []
for couple in ab:
	mu_fact = couple[0]
	mu_f_fact =couple[1]
	FCEL, FCEG = pPb_cross_section.FCEL_G_integration_dy(x_T,num,mu_factor=mu_fact,mu_f_factor=mu_f_fact,switch = convention,var_int= var_int)
	sigma_FCEL , err_FCEL = FCEL[0], FCEL[1]
	sigma_FCEG, err_FCEG = FCEG[0],FCEG[1]
	sigma_tot_FCELG= np.add(sigma_FCEL,sigma_FCEG)
	sigma_FCELG.append(sigma_tot_FCELG)
	err_tot_FCELG = np.sqrt(np.add(np.multiply(err_FCEL,err_FCEL),np.multiply(err_FCEG,err_FCEG)))
	sigma_tot = pPb_cross_section.dsigma_tot_dy(x_T,num,mu_factor=mu_fact,mu_f_factor=mu_f_fact,is_pp=True,switch = convention)
	R.append(np.divide(sigma_tot_FCELG,sigma_tot[0]))

ax = plt.subplot()
fill_between_n_curves(Y,sigma_FCELG,r'$\sigma_\mu$',color='blue', alpha=0.4, show_lines= False,labels = labels)
plt.plot(Y,sigma_tot_FCELG_cen,label=r'$\mu = \mu_f = p_\bot$')
plt.grid()
plt.legend(loc='lower center',frameon =False)
plt.xlabel('y')
plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r' in '+d[convention][1])
#plt.yscale('log')
# ~ plt.title(r'$d\sigma/dy$'+d[convention][0]+ r' with uncertainties on $\mu, \mu_f$ at $p_\bot =$'+str(p_T)+r' $GeV$')
plt.text(0.5, 0.9, r'$p_\bot =$'+str(p_T) +r'$GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_with_mu_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
plt.show()
plt.close()

ax = plt.subplot()
fill_between_n_curves(Y,R,r'$\sigma_\mu$',color='blue', alpha=0.4, show_lines= True,labels = labels)
plt.plot(Y,R_tot_cen,label=r'$\mu = \mu_f = p_\bot$')
plt.axhline(y=1, color='grey', linestyle='--')
plt.grid()
plt.legend(loc='upper center',frameon =False)
plt.xlabel('y')
plt.ylabel(r'R ')
#plt.yscale('log')
# ~ plt.title(r'Ratios with $\mu, \mu_f$ uncertainties, with $p_\bot =$'+str(p_T)+r' $GeV$')
plt.text(0.5, 0.9, r'$p_\bot =$'+str(p_T) +r'$GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'ratios_'+col_type+'_with_mu_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
plt.show()
plt.close()

### All uncertainties ###

Ucen , Utot, Uall = pPb_cross_section.Uncertainties_RpA_dy(x_T,switch = convention,var_int= var_int,var_err='q0,mu')
Uplus, Uminus = Utot[0], Utot[1]
Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf=Uall[0],Uall[1],Uall[2],Uall[3],Uall[4],Uall[5]

fig, ax = plt.subplots(constrained_layout=True)
plt.axhline(y=1, color='grey', alpha=alph)

# plt.fill_between(Y,Ucen-Uminus,Ucen+Uplus,color='purple',alpha=0.4,label=r'$\sigma_{tot}$')
plt.fill_between(Y,Ucen-Uminus_q,Ucen+Uplus_q,color='green',alpha = 0.3,label = r'$\sigma_{\hat{q}_0}$')
plt.fill_between(Y,Ucen-Uminus_mu,Ucen+Uplus_mu,color = 'blue', alpha= 0.3,label = r'$\sigma_{\mu}$')
#plt.fill_between(Y,Ucen-Uminus_pdf,Ucen+Uplus_pdf,color = 'cyan',alpha=0.4,label=r' $\sigma_{pdf}$')
plt.plot(Y,Ucen)
plt.ylim(0.85,1.05)
plt.grid()
plt.legend(loc='upper right',fontsize=f_size-2, frameon = False,fontsize=8) 
plt.xlabel(r'y',fontsize=f_size)
#plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r' in '+d[convention][1])
#plt.title(r'pA cross section with uncertainties at $p_\bot =$'+str(p_T)+r' $GeV$')
plt.ylabel(r'$R_{pA}(y)$',fontsize=f_size)
 #plt.title(r'Cross section ratio with uncertainties at $p_\bot =$'+str(p_T)+r' $GeV$')
plt.text(0.5, 0.9, r'$p_\bot =$'+str(p_T) +r'$GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.tick_params(labelsize=f_size)
plt.savefig(os.path.join(plots_dir, 'RpA_'+col_type+'_with_all_uncertainties_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
plt.close()

