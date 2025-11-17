# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity with the sigma class
# Impact of mu on FCEL and FCEG cross section parts
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
p = "nNNPDF30_nlo_as_0118_p"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"
Z = 82.
A = 208.

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

d = sig.Switch
convention = 'd2p_t'													# Define a cross section convention btw dp_t, dp_t2 ans d2p_t
col_type = 'pp'

### y dependant ###
p_T = 15																# GeV
x_T = (2.*p_T)/rs 
# ~ xi = 0.5
var_int = 'nu'

# plot parameters
f_size = 15
a=1
b=2
alph = 0.3

Y = sig.Y_list(x_T)
# ~ Y_xi = sig.Y_list_xi(x_T,xi)

num = 0																	# The central set 

# ~ FCEL, FCEG = pPb_cross_section.FCEL_G_integration_dy(x_T,num,switch = convention,var_int= var_int)

# ~ sigma_FCEL , err_FCEL = FCEL[0], FCEL[1]
# ~ sigma_FCEG, err_FCEG = FCEG[0],FCEG[1]

# ~ sigma_tot_FCELG = np.add(sigma_FCEL,sigma_FCEG)
# ~ err_tot_FCELG = np.sqrt(np.add(np.multiply(err_FCEL,err_FCEL),np.multiply(err_FCEG,err_FCEG)))

# ~ FCEL_tot, FCEG_tot = pPb_cross_section.dsigma_FCELG_dy(x_T,num,is_pp=True,switch = convention)
# ~ sigma_tot = pPb_cross_section.dsigma_tot_dy(x_T,num,is_pp=True,switch = convention)

# ~ R_tot = np.divide(sigma_tot_FCELG,sigma_tot[0])
# ~ R_FCEL = np.divide(sigma_FCEL,FCEL_tot[0])
# ~ R_FCEG = np.divide(sigma_FCEG,FCEG_tot[0])

# ~ plt.plot(Y,sigma_FCEG,'r',label = r'$d\sigma$ FCEG')
# ~ plt.plot(Y,sigma_FCEL,'b',label = r'$d\sigma$ FCEL')
# ~ plt.plot(Y,FCEG_tot[0],'rx',label = r'$d\sigma_{tot}$ "FCEG"')
# ~ plt.plot(Y,FCEL_tot[0],'bx',label = r'$d\sigma_{tot}$ "FCEL"')
# ~ plt.plot(Y,sigma_tot_FCELG,'grey',label = r'$d\sigma$ tot')
# ~ plt.plot(Y,sigma_tot[0],'purple',label=r'$d\sigma_{pp}$')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'cross section $d\sigma/dy$'+d[convention][0]+ r' in '+d[convention][1])
# ~ #plt.yscale('log')
# ~ plt.title(r'$d\sigma/dy$'+d[convention][0]+ r'with FCEL/G at '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Test_sigma_FCEL_FCEG_'+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()

# ~ plt.plot(Y,R_FCEL,'b',label=r'$R^{FCEL}$')
# ~ plt.plot(Y,R_FCEG,'r',label=r'$R^{FCEG}$')
# ~ plt.plot(Y,R_tot,'purple',label = r'$R^{pA}$')
# ~ plt.axhline(y=1, color='grey', linestyle='--')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'ratios')
# ~ #plt.yscale('log')
# ~ plt.ylim(0.7,1.2)
# ~ plt.title(r'FCEL and FCEG ratios at '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Test_sigma_FCELG_ratios_'+col_type+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()

### Plot en p_t ###
# ~ y = 6
# ~ P_T = pPb_cross_section.P_T_list(y)


# ~ FCEL, FCEG = pPb_cross_section.FCEL_G_integration_dpt(y,num,switch = convention,var_int= var_int)

# ~ sigma_FCEL , err_FCEL = FCEL[0], FCEL[1]
# ~ sigma_FCEG, err_FCEG = FCEG[0],FCEG[1]

# ~ sigma_tot_FCELG = np.add(sigma_FCEL,sigma_FCEG)
# ~ err_tot_FCELG = np.sqrt(np.add(np.multiply(err_FCEL,err_FCEL),np.multiply(err_FCEG,err_FCEG)))

# ~ FCEL_tot, FCEG_tot = pPb_cross_section.dsigma_FCELG_dpt(y,num,is_pp=True,switch = convention)
# ~ sigma_tot = pPb_cross_section.dsigma_tot_dpt(y,num,is_pp=True,switch = convention)

# ~ R_tot = np.divide(sigma_tot_FCELG,sigma_tot[0])
# ~ R_FCEL = np.divide(sigma_FCEL,FCEL_tot[0])
# ~ R_FCEG = np.divide(sigma_FCEG,FCEG_tot[0])

# ~ plt.plot(P_T,sigma_FCEG,'r',label = r'$d\sigma$ FCEG')
# ~ plt.plot(P_T,sigma_FCEL,'b',label = r'$d\sigma$ FCEL')
# ~ plt.plot(P_T,FCEG_tot[0],'rx',label = r'$d\sigma_{tot}$ "FCEG"')
# ~ plt.plot(P_T,FCEL_tot[0],'bx',label = r'$d\sigma_{tot}$ "FCEL"')
# ~ plt.plot(P_T,sigma_tot_FCELG,'grey',label = r'$d\sigma$ tot')
# ~ plt.plot(P_T,sigma_tot[0],'purple',label=r'$d\sigma_{pp}$')
# ~ plt.grid()
# ~ plt.legend(frameon =False,ncol=2,fontsize=8)
# ~ plt.xlabel(r'transverse momentum $p_\bot$ in $GeV$')
# ~ plt.ylabel(r'$d\sigma/dy$'+d[convention][0]+ r' in '+d[convention][1])
# ~ #plt.yscale('log')
# ~ plt.title(r'$d\sigma/dy$'+d[convention][0]+ r'with FCEL/G at '+col_type+r', with $y =$'+str(y))
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_FCEL_FCEG_'+str(rs)+'GeV_'+convention+'_y'+str(y)+'.pdf'))
# ~ plt.close()

# ~ plt.plot(P_T,R_FCEL,'b',label=r'$R^{FCEL}$')
# ~ plt.plot(P_T,R_FCEG,'r',label=r'$R^{FCEG}$')
# ~ plt.plot(P_T,R_tot,'purple',label = r'$R^{pA}$')
# ~ plt.axhline(y=1, color='grey', linestyle='--')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
# ~ plt.xlabel(r'transverse momentum $p_\bot$ in $GeV$')
# ~ plt.ylabel(r'ratios')
# ~ #plt.yscale('log')
# ~ plt.ylim(0.7,1.2)
# ~ plt.title(r'FCEL and FCEG ratios at '+col_type+r', with $y =$'+str(y))
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_FCELG_ratios_'+col_type+str(rs)+'GeV_'+convention+'_y'+str(y)+'.pdf'))
# ~ plt.close()


# ~ etaplus = np.divide(np.subtract(sigma_tot[2],sigma_tot[1]),sigma_tot[1])
# ~ etaminus = np.divide(np.subtract(sigma_tot[0],sigma_tot[1]),sigma_tot[1])
# ~ Uetaplus = np.divide(Uplus,Ucen)
# ~ Uetaminus = np.divide(np.subtract(0,Uminus), Ucen)
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.fill_between(Y,etaplus,etaminus,color='blue',alpha=0.5,label = r'relative CL for $\mu \in [p_\bot/2,2p_\bot]$')
# ~ plt.fill_between(Y,Uetaplus,Uetaminus,color='pink',alpha=0.8,label = r'relavtive CL at 1-$\sigma$ (pdf)')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'relative uncertainties $(Err-d\sigma)/d\sigma$')
# ~ #plt.yscale('log')
# ~ plt.title(r'$(Err-d\sigma)/d\sigma$ for different $\mu = \mu_f$ and pdf members at '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_relative_mu_pdf_'+col_type+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()

### Ratios ###
R_FCEL_dy = []
R_FCEG_dy = []
mu = [0.5,2.]
for i in mu:
	for j in mu:
		print(str(i)+str(j))
		R_all = pPb_cross_section.R_composition_dy(x_T,num,mu_factor=i,mu_f_factor=j,is_pp=Collision(col_type),switch = convention)
		R_FCEG_dy.append(R_all[0][0])
		R_FCEL_dy.append(R_all[1][0]+R_all[2][0]+R_all[3][0])

R_FCEL_min, R_FCEL_max = sig.min_max(R_FCEL_dy)
R_FCEG_min, R_FCEG_max = sig.min_max(R_FCEG_dy)

R_all = pPb_cross_section.R_composition_dy(x_T,num,is_pp=Collision(col_type),switch = convention)
R_FCEG = R_all[0][0]
R_FCEL = R_all[1][0]+R_all[2][0]+R_all[3][0]

fig, ax = plt.subplots(constrained_layout=True)
ax.set_box_aspect(1)
plt.tick_params(labelsize=f_size-b)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
plt.plot(Y,R_FCEL,'b',label=r'$f_{FCEL}$')
plt.plot(Y,R_FCEG,'r',label=r'$f_{FCEG}$')
plt.axhline(y=1, color='grey',alpha=alph)
plt.fill_between(Y,R_FCEL_min,R_FCEL_max,color='b',alpha=0.5,edgecolor=None, linewidth=0)
plt.fill_between(Y,R_FCEG_min,R_FCEG_max,color='r',alpha=0.5,edgecolor=None, linewidth=0)
plt.grid(False)
plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=f_size-b)
plt.xlabel('y',fontsize=f_size-b)
plt.ylabel(r'$f_{FCEi}$',fontsize=f_size-b)
plt.ylim(0,1.1)
plt.text(0.20, 0.45, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-b)
#plt.yscale('log')
# ~ plt.title(r'FCEL and FCEG ratios for different $\mu = \mu_f$ at '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_mu_FCELG_ratios_'+col_type+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.show()
plt.close()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5),sharey = True)

# 	ax.set_anchor('W')
# 	ax.set_position(ax.get_position())

plt.subplots_adjust(wspace=0)

ax1.spines['right'].set_visible(False)
# ax2.spines['left'].set_visible(False)

plt.sca(ax1)
plt.tick_params(labelsize=f_size-b)
ax1.xaxis.set_minor_locator(MultipleLocator(1))
ax1.xaxis.set_major_locator(MultipleLocator(2))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.yaxis.set_major_locator(MultipleLocator(0.2))
plt.axhline(y=0, color='grey',alpha=alph)
plt.fill_between(Y,(R_FCEL_min-R_FCEL)/R_FCEL,(R_FCEL_max-R_FCEL)/R_FCEL,color='b',alpha=0.5,edgecolor=None, linewidth=0,label = r'$\mu$ CL FCEL')
plt.xlabel('y',fontsize=f_size-b)
plt.ylabel(r'rel $\Delta_\mu$ (%)',fontsize=f_size-b)
plt.legend(loc="upper left",frameon=False,fontsize=f_size-b)
plt.ylim(-0.5,0.5)
plt.text(0.20, 0.1, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax1.transAxes,fontsize= f_size-b)
#plt.yscale('log')
# ~ plt.title(r'Reltive CL for R_FCEL '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')

# Deuxième subplot
plt.sca(ax2)  # Définir ax2 comme subplot actif
plt.tick_params(labelsize=f_size-b)
ax2.xaxis.set_minor_locator(MultipleLocator(1))
ax2.xaxis.set_major_locator(MultipleLocator(2))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
ax2.yaxis.set_major_locator(MultipleLocator(0.2))
ax2.tick_params(right= True,labelright=False,which='both')
plt.axhline(y=0, color='grey',alpha = alph)
plt.fill_between(Y,(R_FCEG_min-R_FCEG)/R_FCEG,(R_FCEG_max-R_FCEG)/R_FCEG,color='r',alpha=0.5,edgecolor=None, linewidth=0,label = r'$\mu$ CL FCEG')
plt.grid(False)
plt.xlabel('y',fontsize=f_size-b)
plt.ylim(-0.5,0.5)
plt.legend(loc="upper right",frameon=False,fontsize=f_size-b)

#plt.yscale('log')
# ~ plt.title(r'Reltive CL for R_FCEG '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')

# Ajustement automatique de la mise en page

plt.savefig(os.path.join(plots_dir, 'sigma_mu_FCELG_relative_uncertainties_'+col_type+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'),bbox_inches="tight")
plt.show()
plt.close()

# ~ # Avec les pdf 

# ~ U_FCEL, U_FCEG= pPb_cross_section.Uncertainties_R_FCEl_FCEG_dy(x_T,is_pp=Collision(col_type),switch = convention)

# ~ Uetaplus_FCEL, Uetaminus_FCEL = np.divide(U_FCEL[1],U_FCEL[0]), np.divide(np.subtract(0,U_FCEL[2]), U_FCEL[0])
# ~ Uetaplus_FCEG, Uetaminus_FCEG = np.divide(U_FCEG[1],U_FCEG[0]), np.divide(np.subtract(0,U_FCEG[2]), U_FCEG[0])

# ~ ax = plt.subplot()
# ~ plt.plot(Y,U_FCEL[0],'b',label=r'$R_{FCEL}$')
# ~ plt.plot(Y,U_FCEG[0],'r',label=r'$R_{FCEG}$')
# ~ plt.fill_between(Y,np.add(U_FCEL[0],U_FCEL[1]),np.subtract(U_FCEL[0],U_FCEL[2]),color = 'blue',alpha=0.5,edgecolor=None, linewidth=0)
# ~ plt.fill_between(Y,np.add(U_FCEG[0],U_FCEG[1]),np.subtract(U_FCEG[0],U_FCEG[2]),color = 'red',alpha=0.5,edgecolor=None, linewidth=0)
# ~ plt.grid(False)
# ~ plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=10)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'Ratios')
# ~ #plt.yscale('log')
#plt.title(r'FCEL/G Ratios at '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')
# ~ plt.text(0.1, 0.15, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_pdf_ratios_FCELG_'+col_type+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))

# ~ plt.close()

# ~ fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# ~ plt.sca(ax1)
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.fill_between(Y,Uetaplus_FCEL,Uetaminus_FCEL,color='blue',alpha=0.5,edgecolor=None, linewidth=0,label=r'pdf CL (1$\sigma$) FCEL')
# ~ plt.grid(False)
# ~ plt.ylim(-0.1,0.1)
# ~ plt.legend(loc="upper left",frameon =False,fontsize=10)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylabel(r'relative uncertainties')
# ~ plt.text(0.15, 0.1, r'$p_\bot =$'+str(p_T) +r' $GeV$', horizontalalignment='center', verticalalignment='center',transform=ax1.transAxes)

# ~ #plt.yscale('log')
#plt.title(r'Reltive CL for FCEL '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')

# ~ plt.sca(ax2)  # Définir ax2 comme subplot actif
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.fill_between(Y,Uetaplus_FCEG,Uetaminus_FCEG,color='red',alpha=0.5,edgecolor=None, linewidth=0,label=r'pdf CL (1$\sigma$) FCEG')
# ~ plt.grid(False)
# ~ plt.legend(loc="upper right",frameon =False,fontsize=10)
# ~ plt.xlabel('rapidity y')
# ~ plt.ylim(-0.1,0.1)
# ~ #plt.yscale('log')
#plt.title(r'Reltive CL for FCEG '+col_type+r', with $p_\bot =$'+str(p_T)+r' $GeV$')

# ~ # Ajustement automatique de la mise en page
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'Relative_FCELG_'+col_type+str(rs)+'GeV_'+convention+str(p_T)+'GeV.pdf'))
# ~ plt.close()
# ~ ### p_t dependant ###
# ~ y = 4
# ~ mu_factor = [0.5,1.,2.]

# ~ sigma_FCEL, err_FCEL = [],[]
# ~ sigma_FCEG, err_FCEG = [],[]
# ~ sigma_tot, err_tot = [],[]

# ~ for i in mu_factor:															# move on the mu2 = mu_f2 line 
	# ~ sigma_qG, err_qG  = pPb_cross_section.dsigma_qG_dpt(y,num,is_pp=Collision(col_type),switch = convention,mu_factor=i,mu_f_factor=i)
	# ~ sigma_Gq, err_Gq  = pPb_cross_section.dsigma_Gq_dpt(x_T,num,is_pp=Collision(col_type),switch = convention,mu_factor=i,mu_f_factor=i)
	# ~ sigma_qqbar, err_qqbar = pPb_cross_section.dsigma_qqbar_dpt(x_T,num,is_pp=Collision(col_type), switch = convention,mu_factor=i,mu_f_factor=i)
	# ~ sigma_qbarq, err_qbarq = pPb_cross_section.dsigma_qbarq_dpt(x_T,num,is_pp=Collision(col_type), switch = convention,mu_factor=i,mu_f_factor=i)
	# ~ fcel = np.add(np.add(sigma_Gq,sigma_qqbar),sigma_qbarq)
	# ~ sigma_FCEL.append(fcel)
	# ~ err_fcel = np.add(np.add(err_Gq,err_qqbar),err_qbarq)
	# ~ err_FCEL.append(err_fcel)
	# ~ sigma_FCEG.append(sigma_qG)
	# ~ err_FCEG.append(err_qG)
	# ~ sigma_tot.append(np.add(fcel,sigma_qG))
	# ~ err_tot.append(np.add(err_fcel,err_qG))

# ~ P_T = pPb_cross_section.P_T_list(y)

# ~ plt.plot(P_T,sigma_FCEG[0],'r--',label = r'FCEG,$p_\bot/2$')
# ~ plt.plot(P_T,sigma_FCEG[1],'r',label = r'FCEG,$p_\bot$')
# ~ plt.plot(P_T,sigma_FCEG[2],'r-.',label = r'FCEG,$2p_\bot$')
# ~ plt.plot(P_T,sigma_FCEL[0],'b--',label = r'FCEL,$p_\bot/2$')
# ~ plt.plot(P_T,sigma_FCEL[1],'b',label = r'FCEL,$p_\bot$')
# ~ plt.plot(P_T,sigma_FCEL[2],'b-.',label = r'FCEL,$2p_\bot$')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False,ncol=2,fontsize=8)
# ~ plt.xlabel(r'transverse momentum $p_\bot$ in $GeV$')
# ~ plt.ylabel(r'cross section $d\sigma/dy$'+d[convention][0]+ r' in '+d[convention][1])
# ~ plt.yscale('log')
# ~ plt.title(r'$d\sigma/dy$'+d[convention][0]+ r'components for '+col_type+r', with $y =$ '+str(y))
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_'+col_type+'_FCEL_FCEG_'+str(rs)+'GeV_'+convention+'_y'+str(y)+'.pdf'))
# ~ plt.close()

# ~ etaplus = np.divide(np.subtract(sigma_tot[2],sigma_tot[1]),sigma_tot[1])
# ~ etaminus = np.divide(np.subtract(sigma_tot[0],sigma_tot[1]),sigma_tot[1])
# ~ plt.axhline(y=0, color='grey', linestyle='--')
# ~ plt.fill_between(Y,etaplus,etaminus,alpha=0.5,label =  r'relative CL for $\mu \in [ p_\bot/2,2p_\bot]$')
# ~ plt.grid()
# ~ plt.legend(loc='upper right',frameon =False)
# ~ plt.xlabel(r'transverse momentum $p_\bot$ in $GeV$')
# ~ plt.ylabel(r'relative uncertainties $(Err-d\sigma)/d\sigma$')
# ~ #plt.yscale('log')
# ~ plt.title(r'$(Err-d\sigma)/d\sigma$ for differents $\mu = \mu_f$  '+col_type+r', with $y =$ '+str(y))
# ~ plt.tight_layout()
# ~ plt.savefig(os.path.join(plots_dir, 'sigma_relative_mu'+col_type+str(rs)+'GeV_'+convention+'_y'+str(y)+'.pdf'))
# ~ plt.close()

