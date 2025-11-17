# -*- coding: utf-8 -*-

# =============================================================================
# Study the first impacts of compton QCD diagrams and annihilation 
# to parton luminosity.
# First attempt before sigma.py
# Last modified: 30/10/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import lhapdf as lha
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Ininitalisation of a patricle dictionary (i've not seen one in lhapdf documentation)
particle = {
    "d": {"id": 1, "charge": -1./3}, 
    "u": {"id": 2, "charge": 2./3}, 
    "s": {"id": 3, "charge": -1./3}, 
    "c": {"id": 4, "charge": 2./3}, 
    "b": {"id": 5, "charge": -1./3}, 
    "t": {"id": 6, "charge": 2./3},  # Quarks

    "dbar": {"id": -1, "charge": 1./3}, 
    "ubar": {"id": -2, "charge": -2./3}, 
    "sbar": {"id": -3, "charge": 1./3}, 
    "cbar": {"id": -4, "charge": -2./3}, 
    "bbar": {"id": -5, "charge": 1./3}, 
    "tbar": {"id": -6, "charge": -2./3},  # Antiquarks

    "g": {"id": 21, "charge": 0},  # Gluon
    "gamma": {"id": 22, "charge": 0},  # Photon
    "Z0": {"id": 23, "charge": 0},  # Boson Z
    "W+": {"id": 24, "charge": 1}, 
    "W-": {"id": -24, "charge": -1},  # Bosons W
    "h0": {"id": 25, "charge": 0},  # Boson de Higgs
}

id_to_charge = {v["id"]: v["charge"] for v in particle.values()}

def charge(particle_id):
	'''return the charge of a particle giving its id'''
	return(id_to_charge.get(particle_id,None))

# Initialisation of PDFs
p = lha.getPDFSet("nCTEQ15npFullNuc_1_1")
Pb = lha.getPDFSet("nCTEQ15npFullNuc_208_82")
p_pdf = p.mkPDF(0)
Pb_pdf = Pb.mkPDF(0)
N_c = 3					  # number of colors 
C_F = (N_c**2 -1)/(2*N_c) # color factor  

# Variables initiallisation
s = (8.8*pow(10,4))**2 #GeV^2, regarding the 2008 paper
alpha = 1./137
conv_fact = pow(197.,2)/100 # to convert Mev^-2 into barns

### Important functions ### 
def alpha_s(Q2,pdf):
	return pdf.alphasQ2(Q2)

def n_xfxQ2(pdf,particle_id,x,Q2):
	if np.abs(particle_id) >= 3:
		return pdf.xfxQ2(particle_id,x,Q2)
	else :
		return pdf.xfxQ2(2//particle_id,x,Q2)

def F2(pdf,x,Q2,iso = 'p',n_f=3,alpha_s_fun=True):
	'''gives back the F2 function for a x and a Q^2 given 
	(with an alpha_s object to extimates the number of possible flavours,
	 putting a random float would let n_f decides the number of flavours)'''
	F = 0.
	if alpha_s_fun:
		alphas = alpha_s(Q2,pdf)
	else:
		alphas = 0 # See for AlphaS Ojects
	if x == 0:
		return 0.
	else:
		if False:									# It doesnt interprete  if type(alphas) == lhapdf.Alpha_s
			n_f = alpha_s.numFlavoursQ2(Q2)
			if iso == 'p':
				for i in range(1,n_f+1):
					F += (charge(i)**2)*(pdf.xfxQ2(i,x,Q2)+pdf.xfxQ2(-i,x,Q2))
			elif iso == 'n':
				for i in range(1,n_f+1):
					F += (charge(i)**2)*(n_xfxQ2(pdf,i,x,Q2)+n_xfxQ2(pdf,-i,x,Q2))
		elif type(n_f) == int:
			if iso == 'p':
				for i in range(1,n_f+1):
					F += (charge(i)**2)*(pdf.xfxQ2(i,x,Q2)+pdf.xfxQ2(-i,x,Q2))
			elif iso == 'n':
				for i in range(1,n_f+1):
					F += (charge(i)**2)*(n_xfxQ2(pdf,i,x,Q2)+n_xfxQ2(pdf,-i,x,Q2))
	return F/x

def F_ij(xa,xb,Q2,pdf,npdf, direction = 'qqbar',n_f=3):
	'''return the F_ij function that appears in the annihilation process'''
	F = 0.
	if direction == 'qqbar':
		for i in range(1,n_f+1):
			F += (charge(i)**2)*(pdf.xfxQ2(i,xa,Q2))*(npdf.xfxQ2(-i,xb,Q2))
	elif direction == 'qbarq':
		for i in range(1,n_f+1):
			F += (charge(i)**2)*(pdf.xfxQ2(-i,xa,Q2))*(npdf.xfxQ2(i,xb,Q2))
	return F/(xa*xb)	

def Gluon(pdf,x,Q2):
	return pdf.xfxQ2(particle["g"]["id"],x,Q2)/x 

### Integrand functions and their plots ###
def qG(y,x_T,Xi,pdf,npdf,iso = 'p'):
	'''return the first part of compton process (replace iso p by n for a nX collision)'''
	Q2 = (x_T**2)*s
	x_proj = x_T*np.exp(y)/(2*Xi)
	x_targ = x_T*np.exp(-y)/(2*(1-Xi))
	F = F2(pdf,x_proj,Q2,iso)
	G = Gluon(npdf,x_targ,Q2)
	return[F*G*(1./3)*(1-Xi+(1./(1-Xi))),F*G,F,G]

def Gq(y,x_T,Xi,pdf,npdf,iso = 'p'):
	'''return the second part of compton process (replace iso p by n for a Xn collision)''' 
	Q2 = (x_T**2)*s
	x_proj = x_T*np.exp(y)/(2*Xi)
	x_targ = x_T*np.exp(-y)/(2*(1-Xi))
	F = F2(npdf,x_targ,Q2,iso)
	G = Gluon(pdf,x_proj,Q2)
	return[F*G*(1./3)*(Xi+1./Xi),F*G,F,G]
	
def qqbar(y,x_T,Xi,pdf,npdf):
	Q2 = (x_T**2)*s
	x_proj = x_T*np.exp(y)/(2*Xi)
	x_targ = x_T*np.exp(-y)/(2*(1-Xi))
	F_qqbar = F_ij(x_proj,x_targ,Q2,pdf,npdf)
	Xi_factor = Xi/(1-Xi) + (1-Xi)/Xi 
	return[F_qqbar*Xi_factor*(2*C_F/N_c),F_qqbar]
	
def qbarq(y,x_T,Xi,pdf,npdf):
	Q2 = (x_T**2)*s
	x_proj = x_T*np.exp(y)/(2*Xi)
	x_targ = x_T*np.exp(-y)/(2*(1-Xi))
	F_qbarq = F_ij(x_proj,x_targ,Q2,pdf,npdf,dirrection = 'qbarq')
	Xi_factor = Xi/(1-Xi) + (1-Xi)/Xi 
	return[F_qbarq*Xi_factor*(2*C_F/N_c),F_qbarq]

def Gq_qG_integrand(y,x_T,Xi,pdf,npdf,iso='p'):
	'''return the whole integrand part for the compton process'''
	Q2 = (x_T**2)*s
	if iso == 'n':
		return(((alpha*alpha_s(Q2,pdf))/pow(x_T*s/2,2))*(qG(y,x_T,Xi,pdf,npdf)[0] +Gq(y,x_T,Xi,pdf,npdf,iso)[0]))
	else:
		return(((alpha*alpha_s(Q2,pdf))/pow(x_T*s/2,2))*(qG(y,x_T,Xi,pdf,npdf)[0] +Gq(y,x_T,Xi,pdf,npdf)[0])) # verify alpha_s for pdf and npdf 

def Gq_integrand(y,x_T,Xi,pdf,npdf,iso='p'):
	'''return only the Gq canal of the compton process'''
	Q2 = (x_T**2)*s
	if iso == 'n':
		return(((alpha*alpha_s(Q2,pdf))/pow(x_T*s/2,2))*Gq(y,x_T,Xi,pdf,npdf,iso)[0])
	else:
		return(((alpha*alpha_s(Q2,pdf))/pow(x_T*s/2,2))*Gq(y,x_T,Xi,pdf,npdf)[0]) # verify alpha_s for pdf and npdf 

def qG_integrand(y,x_T,Xi,pdf,npdf,iso='p'):
	'''retrun only the qG canal of the compton process'''
	Q2 = (x_T**2)*s
	if iso == 'n':
		return(((alpha*alpha_s(Q2,pdf))/pow(x_T*s/2,2))*qG(y,x_T,Xi,pdf,npdf)[0])
	else:
		return(((alpha*alpha_s(Q2,pdf))/pow(x_T*s/2,2))*qG(y,x_T,Xi,pdf,npdf)[0]) # verify alpha_s for pdf and npdf 

def all_process_integrand(y,x_T,Xi,pdf,npdf,iso = 'p'):
	Q2 = (x_T**2)*s
	if iso == 'n':
		return(((alpha*alpha_s(Q2,pdf))/pow(s,2))*(2./pow(x_T,2) *(qG(y,x_T,Xi,pdf,npdf)[0]+Gq(y,x_T,Xi,pdf,npdf,iso)[0])+ qqbar(y,x_T,Xi,pdf,npdf)[0]+ qbarq(y,x_T,Xi,pdf,npdf)[0])) 
	else:
		return(((alpha*alpha_s(Q2,pdf))/pow(s,2))*(2./pow(x_T,2) *(qG(y,x_T,Xi,pdf,npdf)[0]+Gq(y,x_T,Xi,pdf,npdf)[0])+ qqbar(y,x_T,Xi,pdf,npdf)[0]+ qbarq(y,x_T,Xi,pdf,npdf)[0])) # verify alpha_s for pdf and npdf 
		

### Test part ###

# Verifying the neutron pdf 
X = np.logspace(-4,0,1000)*0.9999 # to stop before x = 1
Fn = [F2(p_pdf,x,Q2 = 100, iso ='n') for x in X]
Fp = [F2(p_pdf,x,Q2=100) for x in X]

plt.plot(X,Fn,label= r'Neutron')
plt.plot(X,Fp,label= r'Proton')
plt.legend()
plt.xlabel(r'parton\'s impulsion fraction x')
plt.xscale('log')
plt.ylabel(r'$F_2/x$')
plt.yscale('log')
plt.title(r'$F_2/x$ function with isospin parameter')
plt.savefig(os.path.join(plots_dir, 'isospin_comparison.pdf'))
plt.show()
plt.close()

### Integration part ###

# Not integration but some calculus
N_Xi = 3
Xi_test = [0.4,0.5,0.6]
x_t = 3.59*pow(10,-3) # for example (with Q2 = 100 GeV^2 and sqrt(s) = 8.8 TeV )
q2 = 100.
Y = []
# qG_list = []
# Gq_list = []
# Fp_list = []
# Fn_list=[]
# FPb_list= []
# Gp_list = []
# GPb_list= []
# for i in range(N_Xi):
	# #y_min = -np.log(2*(1-Xi_test[i])/x_t)
	# #y_max = np.log(2*Xi_test[i]/x_t)
	# qG_list.append([qG(y,x_t,Xi_test[i],p_pdf,Pb_pdf)[1] for y in Y[i]])
	# Gq_list.append([Gq(y,x_t,Xi_test[i],p_pdf,Pb_pdf)[1] for y in Y[i]])
	# Fp_list.append([qG(y,x_t,Xi_test[i],p_pdf,Pb_pdf)[2] for y in Y[i]])
	# Fn_list.append([Gq(y,x_t,Xi_test[i],p_pdf,p_pdf,'n')[2] for y in Y[i]])
	# FPb_list.append([Gq(y,x_t,Xi_test[i],p_pdf,Pb_pdf)[2] for y in Y[i]])
	# Gp_list.append([Gq(y,x_t,Xi_test[i],p_pdf,Pb_pdf)[3] for y in Y[i]])
	# GPb_list.append([qG(y,x_t,Xi_test[i],p_pdf,Pb_pdf)[3] for y in Y[i]])
y_min = -np.log(2./x_t)+0.1
y_max = np.log(2./x_t)-0.1
Y = np.linspace(y_min,y_max,1000)
	
# Now integration
i = 2 # Xi = 0.5 
epsilon = pow(10,-5) # to avoird unphysical values like x>1
Xi_min = lambda y,x_T: x_T*np.exp(y)/2 + epsilon
Xi_max = lambda y,x_T: 1-x_T*np.exp(-y)/2 - epsilon
Z = 82.
A = 208.
sigma_pn = []
sigma_pp = []
sigma_pp_qG=[]
sigma_pp_Gq=[]
sigma_pPb =[]
for y in Y: # We choose a well define y array
	Xi = np.linspace(Xi_min(y,x_t),Xi_max(y,x_t),100)
	Integrand_pn = [Gq_qG_integrand(y,x_t,xi,p_pdf,p_pdf,iso='n') for xi in Xi]
	Integrand_pp = [Gq_qG_integrand(y,x_t,xi,p_pdf,p_pdf) for xi in Xi]
	Integrand_pPb = [Gq_qG_integrand(y,x_t,xi,p_pdf,Pb_pdf) for xi in Xi]
	Integrand_pp_qG = [qG_integrand(y,x_t,xi,p_pdf,p_pdf) for xi in Xi ]
	Integrand_pp_Gq = [Gq_integrand(y,x_t,xi,p_pdf,p_pdf) for xi in Xi]
	sigma_pn.append(integrate.simpson(Integrand_pn,Xi))
	sigma_pp.append(integrate.simpson(Integrand_pp,Xi))
	sigma_pPb.append(integrate.simpson(Integrand_pPb,Xi))
	sigma_pp_qG.append(integrate.simpson(Integrand_pp_qG,Xi))
	sigma_pp_Gq.append(integrate.simpson(Integrand_pp_Gq,Xi))
	
# Test sur la forme de qG et Gq avec F2 et G

# plt.plot(Y[i],Fp_list[i],label = r'$F_p$ for $\xi =$'+str(Xi_test[i]))
# plt.plot(Y[i],Fn_list[i],label = r'$F_{n}$ for $\xi =$'+str(Xi_test[i]))
# #plt.plot(Y[i],Gp_list[i],label = r'$G_p$ for $\xi =$'+str(Xi_test[i]))
# #plt.plot(Y[i],GPb_list[i],label = r'$G_{Pb}$ for $\xi =$'+str(Xi_test[i]))
# plt.legend(loc='upper right',ncol=2,fontsize=5,frameon =False)
# plt.grid()
# plt.xlabel('rapidity y')
# plt.ylabel('$F_2/x$')
# plt.yscale('log')
# plt.title(r'$F_2/x$ and G for p and n at $\xi=$'+str(Xi_test[i]) + r' and $x_\bot=$'+str(x_t))
# plt.savefig(os.path.join(plots_dir, 'test_F2_p_and_n.pdf'))

# plt.close()	

# # changement de variable y -> x_1 (ou x_proj) et x_2 (ou x_targ)
# X_1 = np.exp(Y[i])*x_t/(2*Xi_test[i])
# X_2 = np.divide(x_t**2,X_1)

# plt.plot(X_1,Fp_list[i],label = r'$F_p$ for $\xi =$'+str(Xi_test[i]))
# plt.plot(X_2,FPb_list[i],label = r'$F_{Pb}$ for $\xi =$'+str(Xi_test[i]))
# plt.plot(X_1,Gp_list[i],label = r'$G_p$ for $\xi =$'+str(Xi_test[i]))
# plt.plot(X_2,GPb_list[i],label = r'$G_{Pb}$ for $\xi =$'+str(Xi_test[i]))
# plt.legend(loc='upper right',ncol=2,fontsize=5,frameon =False)
# plt.grid()
# plt.xlabel(r'parton impulsion $x_1(p)$ or $x_2(P_b)$')
# plt.ylabel('$F_2/x$ or G')
# plt.xscale('log')
# plt.title(r'$F_2/x$ and G for p and Pb at $\xi=$'+str(Xi_test[i]) + r' and $x_\bot=$'+str(x_t))
# plt.savefig(os.path.join(plots_dir, 'test2_G_F2.pdf'))

# plt.close()	

# ### Plot of different results ###
# # The first two for some values of Xi	
# for j in range(N_Xi):
	# plt.plot(Y[j],qG_list[j],label = r'qG for $\xi =$'+str(Xi_test[j]))

# for k in range(N_Xi):													# I do 2 loops for arranging the legend (it seems to not working for colours)
	# plt.plot(Y[k],Gq_list[k],label = r'Gq for $\xi =$'+str(Xi_test[k]))

# plt.legend(loc='upper right',ncol=2,fontsize=5,frameon =False)
# plt.grid()
# plt.xlabel('rapidity y')
# plt.ylabel('qG or gQ integrand')
# plt.title(r'qG or Gq for different $\xi$ values at $x_\bot=$'+str(x_t))
# plt.savefig(os.path.join(plots_dir, 'npdf_pdf_Qg_Gq.pdf'))


### A digression for ispospin comparisons ###
# the total plot
plt.close()
plt.plot(Y,np.multiply(sigma_pp,conv_fact),label=r'$d\sigma_{pp}/dy$')
#plt.plot(Y[i],np.multiply(sigma_pPb,conv_fact),label=r'$d\sigma_{pPb}/dy$')
#plt.plot(Y[i],np.multiply(sigma_pn,conv_fact), label= r'$d\sigma_{pn}/dy$')
plt.plot(Y,np.multiply(sigma_pp_Gq,conv_fact), label= r'$d\sigma_{pp,Gq}/dy$')
plt.plot(Y,np.multiply(sigma_pp_qG,conv_fact), label= r'$d\sigma_{pp,qG}/dy$')
plt.grid()
plt.legend(loc='upper right',frameon =False)
plt.xlabel('rapidity y')
#plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1e}'.format(x)))
plt.ylabel(r'cross section $d\sigma/dy$ in barn')
plt.title(r'$d\sigma/dy$ for pp and pPb, with $x_\bot =$'+str(x_t))
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'sigma_pp_compton_canal.pdf'))
plt.show()

R = np.divide(sigma_pPb,sigma_pp)
R_iso = (Z/A)+np.divide(sigma_pn,sigma_pp)*(1-Z/A)
R_no_iso = np.divide(R,R_iso)

# The ratio plots
# 1 : Cross section ratios
plt.close()
plt.plot(Y,R,label = r'$R_{pPb}$')
plt.plot(Y,R_iso, label= r'$R_{pn}$ (isospin effect)')
plt.plot(Y,R_no_iso,label = r'$R_{pPb}/R_{pn}$ (w/o isospsin effect)')
plt.grid()
plt.legend(loc='upper right',ncol=2,fontsize=5,frameon =False)
plt.xlabel('rapidity y')
plt.ylabel(r'Ratios')
#plt.yscale('log')
#plt.xlim(-6,6)
#plt.ylim(0,4)
plt.title(r'All Ratios with $x_\bot =$'+str(x_t))
plt.savefig(os.path.join(plots_dir, 'All_R.pdf'))
plt.show()
plt.close()

# 2 : Gluon npdf/pdf ratio
Gp = [Gluon(p_pdf,x,q2) for x in X]
GPb = [Gluon(Pb_pdf,x,q2) for x in X]
plt.plot(X,np.divide(GPb,Gp))
plt.grid()
plt.xlabel('parton impulsion ratio x')
plt.ylabel(r'$G_{Pb}/G_p$')
plt.xscale('log')
#plt.xlim(-6,6)
#plt.ylim(0,4)
plt.title(r'Gluon ratio for $x_\bot =$'+str(x_t))
plt.savefig(os.path.join(plots_dir, 'Gluon_ratio.pdf'))
plt.show()
plt.close()

### All of this is done on rapidity, but we have (to compare) to do it on x_T  (so at fixed rapidity but not p_T) ###
### A try for fixed rapidity, and a variable x_T (or p_T, same) ###
# plt.close()

# # with the same Xi_test, we have:
# y = 0 
# N_x_T = 1000
# x_T_min = 0 # maybe + epsilon, i feel the error comming
# x_T_max = [min(2*i*np.exp(-y),2*(1-i)*np.exp(y)) for i in Xi_test]
# x_T_list = [np.linspace(x_T_min,m,N_x_T) for m in x_T_max]
# qG_list_x_T = [[qG(y,x_T,Xi_test[i],p_pdf,Pb_pdf) for x_T in x_T_list[i]] for i in range(N_Xi)]
# Gq_list_x_T = [[Gq(y,x_T,Xi_test[i],p_pdf,Pb_pdf) for x_T in x_T_list[i]] for i in range(N_Xi)]

# sigma_pp = []
# sigma_pPb =[]
# for y in Y[0]: # We choose a well define y array
	# Xi = np.linspace(Xi_min(y,x_t),Xi_max(y,x_t),100)
	# Integrand_pp = [Gq_qG_integrand(y,x_t,xi,p_pdf,p_pdf) for xi in Xi]
	# Integrand_pPb = [Gq_qG_integrand(y,x_t,xi,p_pdf,Pb_pdf) for xi in Xi]
	# sigma_pp.append(integrate.simps(Integrand_pp,Xi))
	# sigma_pPb.append(integrate.simps(Integrand_pPb,Xi))
	
