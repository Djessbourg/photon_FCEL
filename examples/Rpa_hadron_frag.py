# -*- coding: utf-8 -*-

# =============================================================================
# File for computing R_pA tot 
# as a function of R_pA dir 
# and R_pA frag.
# Last modified: 12/11/2025
# =============================================================================
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))
examples_dir =  os.path.abspath(os.path.dirname(__file__))						# current directory
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_tot_dir_onef')) # the direcory to save data from this file

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from src import sigma as sig
from src import Probability as p
from scipy.integrate import quad

import uproot
from scipy.interpolate import CubicSpline
from scipy.interpolate import make_smoothing_spline

min_max = sig.min_max 															# to compute the min and max of a set of data


# Context data

rs = 8800
s = float((rs)**2) # CM energy in Gev^2
p_T_dispo = [3,5,10,15]
z = 0.6
Xi = 0.5

A = 208
Z = 82
B = 1
alpha_s = 0.5
m_photon = 0
m_pion = 0.134
m = m_photon

# plot variables
alph = 0.3
f_size = 17
tot_color, dir_color , frag_color = 'blue' , 'green' , 'orange'
ylims = {
	3:(0.2,1.10),
	5:(0.5,1.1),
	10:(0.4,1.3),
	15:(0.4,1.4)
}

# =============================================================================
# In any pt case if Ny == 41, Y_list is :
# >>> Y_list
# array([-6. , -5.7, -5.4, -5.1, -4.8, -4.5, -4.2, -3.9, -3.6, -3.3, -3. ,
#        -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3,  0. ,  0.3,
#         0.6,  0.9,  1.2,  1.5,  1.8,  2.1,  2.4,  2.7,  3. ,  3.3,  3.6,
#         3.9,  4.2,  4.5,  4.8,  5.1,  5.4,  5.7,  6. ])
# So to chose a fixed rapidity to plot R_pA pt dependant, pick it here.
# =============================================================================

fixed_y = [-3.9,0,3.9] #+/-3.9 would be rounded to +/-4

### Initialisation of a sigma object ###
proton = "nNNPDF30_nlo_as_0118_p"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

pPb_cross_section = sig.Sigma(proton,Pb,s,Z,A)

# Color informations on all processes

process ={
	"1":{"rho":{"3bar":2./3.,"6":1./3.},"Fc":{"3bar":4./3.,"6":2.}},
	"3":{"rho":{"1":8./9.,"8":1./9.},"Fc":{"1":0.,"8":3.}},
	"5":{"rho":{"3bar":4./22.,"6":19./22.},"Fc":{"3bar":4./3.,"6":2.}},
	"7":{"rho":{"1":0.,"8":1.},"Fc":{"1":0.,"8":3.}},
	"8":{"rho":{"1":80./105.,"8":25./105.},"Fc":{"1":0.,"8":3.}},
	"9":{"rho":{"1":0.,"8":1.},"Fc":{"1":0.,"8":3.}},
	"13a":{"rho":{"3":25./88.,"6bar":18./88.,"15":45./88.},"Fc":{"3":-1./3.,"6bar":5./3.,"15":11./3.}},
	"13b":{"rho":{"3":25./88.,"6bar":18./88.,"15":45./88.},"Fc":{"3":3.,"6bar":5.,"15":7.}},
	"14a":{"rho":{"3":25./88.,"6bar":18./88.,"15":45./88.},"Fc":{"3":-1./3.,"6bar":5./3.,"15":11./3.}},
	"14b":{"rho":{"3":25./88.,"6bar":18./88.,"15":45./88.},"Fc":{"3":3.,"6bar":5.,"15":7.}},
	"15":{"rho":{"1":2./7.,"8a":0.,"8s":5./7.},"Fc":{"1":0.,"8a":3.,"8s":3.}},
	"16":{"rho":{"1":1./6.,"8":1./3.,"27":1./2.},"Fc":{"1":0.,"8":3.,"27":8.}}
}

p_num_to_process ={
	"1":r"$qq'\rightarrow q(q'\rightarrow \gamma)$",
	"3":r"$q\bar{q'} \rightarrow q(\bar{q'}\rightarrow \gamma)$",
	"5":r"$qq \rightarrow q(q\rightarrow \gamma)$",
	"7":r"$q\bar{q} \rightarrow q'(\bar{q'}\rightarrow \gamma)$",
	"8":r"$q\bar{q} \rightarrow q(\bar{q}\rightarrow \gamma)$",
	"9":r"$q\bar{q} \rightarrow g(g\rightarrow \gamma)$",
	"13a":r"$q_1g_2\rightarrow q(g\rightarrow \gamma)$",
	"13b":r"$g_1q_2\rightarrow q(g\rightarrow \gamma)$",
	"13":r"$qg\rightarrow q(g\rightarrow \gamma)$",
	"14a":r"$q_1g_2\rightarrow g(q\rightarrow \gamma)$",
	"14b":r"$g_1q_2 \rightarrow g(q\rightarrow \gamma)$",
	"14":r"$qg\rightarrow g(q\rightarrow \gamma)$",
	"15":r"$gg\rightarrow \bar{q}(q\rightarrow \gamma)$",
	"16":r"$gg\rightarrow g(g\rightarrow \gamma)$"
}

# Extraction of f_alpha as functions of y for p_T = 5 GeV

# =============================================================================
# Precisions: to compute 14a/14b/13a/13b you have to modify the 
# htermtwo_jet.f file in the directory /jetphox_1.3.1_4/src/onef/hcoeff .
# At the very end of the file you will have an object called SF. To compute 14a
# just comment with 'C's the lines for SF(14+J0MAX) and write next line :
# SF(14+J0MAX)=0. You do the inverse for 14b and the process is the same for 
# channels 13a and 13b. And in the param.indat file just do the frag loop over 
# 14,14 or 13,13.
# =============================================================================

# open ROOT file 

def prefix(string):																# useful for normalizing root histograms
	'''return the prefix of a string with a ; separator'''
	s = ''
	for i in string:
		if i ==';':
			break
		s = s + i
	return s

f_alphas = {}
f_dirs={}
Processes = {}

for pt in p_T_dispo:
	pawres_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'copy_pawres_'+str(pt))) # put here your jetphox data
	file = uproot.open(os.path.join(pawres_dir,"LO_frag_ratios.root"))
	# retrieve all fragmentaion histograms
	histo = {k: file[k] for k in file.keys() if file[k].classname.startswith("TH1")} #labels are the same as the process dict but with a ';1' in addition
	f_alpha_cs = {}
	f_alpha_smooth = {}
	# Create a dictionary of cubic sline fits of the f_alpha
	for k in histo.keys():
		values = histo[k].values()
		edges = histo[k].axes[0].edges()
		centers = (edges[:-1] + edges[1:]) / 2.
		cs = CubicSpline(centers, values)
		smooth = make_smoothing_spline(centers, values)
		lil_Y = np.linspace(-6, 6,31)
		smooth2 = make_smoothing_spline(lil_Y, smooth(lil_Y))
		f_alpha_cs[prefix(k)]=cs 													# cs gives an array for a float or an array of float given
		f_alpha_smooth[prefix(k)]=smooth2 											# smooth also, but smoother
	f_alphas[pt]=f_alpha_smooth
	# Same steps for f_dir
	Rdir = uproot.open(os.path.join(pawres_dir,"ggddir_onef.root"))
	Ronef = uproot.open(os.path.join(pawres_dir,"ggodir_onef.root"))
	fname = 'hp20;1'
	hdir = Rdir[fname]
	honef = Ronef[fname]
	Dir = hdir.values()
	Onef = honef.values()
	edges = hdir.axes[0].edges()
	centers = (edges[:-1] + edges[1:]) / 2.
	Tot = Dir+Onef
	fdir  = Dir/Tot
	f_dir = CubicSpline(centers,fdir)
	f_dirs[pt]=f_dir
	Processus = {} 																	# dict for all cross sections interpolated from jetphox data
	for p_num in list(process.keys()):
		fn = 'ggoLOdiag_'+p_num+'.root'
		file = uproot.open(os.path.join(pawres_dir,fn))
		histo = file[fname]
		Values = histo.values()
		edges = histo.axes[0].edges()
		rapidity = (edges[:-1] + edges[1:]) / 2.
		# 	cs = CubicSpline(rapidity,Values)
		smooth = make_smoothing_spline(rapidity,Values)
		lil_Y = np.linspace(-6, 6,31)
		smooth2 = make_smoothing_spline(lil_Y, smooth(lil_Y))
		Processus[p_num] = smooth2
	Processes[pt]=Processus


def f_alpha_pt(pt):
	'''Retrun the dict of f_aphas for a given pt'''
	return f_alphas[pt]

def f_dir_pt(pt):
	return f_dir[pt]

def Processus_pt(pt):
	return Processes[pt]

def cross_section(y,pt,p_num):
	return Processus_pt(pt)[p_num](y)

def sign(a):
	'''Gives back the sign of a real number a'''
	if a > 0:
		return 1
	if a < 0:
		return -1
	if a == 0:
		raise ValueError("a = 0")

def f_xi(Fc,xi):
	'''The f_xi function to differentiate FCEL to FCEG in the integration bounds'''
	if Fc > 0:
		return xi
	if Fc < 0:
		return 1-xi
	elif Fc == 0:
		raise ValueError("Fc = 0")

def Integrand(y,pt,p_num,Fc,z,q):
	'''Gives back the integrand for a process p_num with one of its color
	factor Fc'''
	prob = p.proba(A, B, rs, pt, y, alpha_s, Fc, m,q0=q,z=z)
	sigma_hat = lambda xi : prob.sigma_hat(xi)
	chi = lambda xi:prob.Chi(xi)
	g_FCEL = lambda u,xi: p.g_u(u,chi(xi),Fc,alpha_s)
	P = lambda nu,xi,nu_min: p.p_tilde_u(np.exp(nu),chi(xi),Fc,alpha=alpha_s)/(1-np.exp(-1*g_FCEL(np.exp(nu_min),xi)))
	jacobian = lambda nu,xi: np.exp(nu)*pow(sigma_hat(xi)*np.exp(nu)+1,-1*sign(Fc))
	delta = lambda nu,xi: np.log(1+sigma_hat(xi)*np.exp(nu))
	dsigma = lambda nu,xi: cross_section(y+sign(Fc)*delta(nu,xi),pt,p_num)		# if FCEL: y+delta ; if FCEG: y-delta
	return lambda nu,xi,nu_min:jacobian(nu,xi)*P(nu,xi,nu_min)*dsigma(nu,xi)

def Integration_xi_fixed(y,pt,p_num,Fc,xi,z,q,eps=10**(-10)):
	'''Integration over nu at fixed xi (bc the cross section is not xi dependant)
	for a process p_num with one of its color factor Fc'''
	prob = p.proba(A,B,rs,pt,y,alpha_s,Fc,m,q0=q,z=z)
	sigma_hat = prob.sigma_hat(xi)
	k = rs/pt
	delta_y_max = min(np.log(k*f_xi(Fc, xi))-sign(Fc)*y,np.log(2))
	delta_y_min = eps
	nu_max = np.log(np.exp(delta_y_max-1.)/sigma_hat)
	if delta_y_min > 1e-15:
		nu_min =  np.log((np.exp(delta_y_min)-1)/sigma_hat)
	elif delta_y_min > 0:
		nu_min = np.log(delta_y_min/sigma_hat)
	I = Integrand(y,pt,p_num,Fc,z,q)
	return quad(lambda nu:I(nu,xi,nu_min), nu_min, nu_max)

def RpA_xi(y,pt,p_num,Fc,xi,z,q):
	'''Ratio of sigma pA over sigma pp for a process p_num with one of its
	color factor Fc'''
	I = Integration_xi_fixed(y,pt,p_num,Fc, xi,z,q)
	sigma_pA = I[0]
	sigma_pp = cross_section(y,pt,p_num)
	return sigma_pA/sigma_pp

def all_alpha_colors_xi(Y,pt,p_num,xi,z,q):
	'''Ratios of sigma pA over sigma pp for all color representation of the
	process p_num at a fixed xi'''
	alpha = process[p_num]
	Rho = alpha['rho']
	Fc = alpha['Fc']
	colors = list(Rho.keys())
	R_tot = 0
	Result = {}
	for c in colors:
		R =[]
		rho , fc = Rho[c], Fc[c]
		if rho == 0:
			continue
		if fc == 0:
			for y in Y:
				R.append(1.)
		else:
			for y in Y:
				R.append(RpA_xi(y,pt,p_num,fc,xi,z,q))
		R = np.array(R)															# To do basic array operations
		Result[c]=R
		R_tot = R_tot+rho*R
	Result['tot'] = R_tot
	return Result

p_num_plot_list = ['14a','14b','1','3','15']
p_num_color_list = ['blue','orange','green','red','purple']
z = [0.6,0.4,0.8]																# with the norm u = [u_0,u_-,u_+]
q = [0.07,0.05,0.09]
xi = [0.5,0.25,0.75]

def R_frag_alpha(Y,pt,p_num,z,q,xi):
	'''Return the RpA for a given process p_num'''
	R = all_alpha_colors_xi(Y,pt,p_num,xi, z,q)
	R_tot = R['tot']
	return R_tot

def R_frag(Y,pt,z,q,xi,plot = False):
	'''Return the RpA of all fragmentation processes'''
	keys = list(process.keys())
	R = np.zeros_like(Y)
	for p_num in keys:
		f = f_alpha_pt(pt)[p_num]
		R_alpha = R_frag_alpha(Y,pt, p_num, z, q, xi,)
		if plot and (p_num in p_num_plot_list):
			plt.plot(Y,R_alpha,label=p_num_to_process[p_num])
		R = R+f(Y)*R_alpha
	return R

# here fo the uncertainty functions

def R_frag_uncertainties(Y,pt,z,q,xi):
	'''Return the RpA of all fragmentation processes (i=0) 
	with z,xi uncertainties (i=1,2) 
	and then z,q,xi uncertainties (i=3,4)
	with z, q, xi  should be lists of [v0,v-,v+]'''
	R_cen = R_frag(Y,pt,z[0],q[0],xi[0])
	R_Z = [R_cen,R_frag(Y,pt,z[1],q[0],xi[0]),R_frag(Y,pt,z[2],q[0],xi[0])]
	R_xi = [R_cen,R_frag(Y,pt,z[0],q[0],xi[1]),R_frag(Y,pt,z[0],q[0],xi[2])]
	R_q = [R_cen,R_frag(Y,pt,z[0],q[1],xi[0]),R_frag(Y,pt,z[0],q[2],xi[0])]
	U_Z = min_max(R_Z)
	U_xi = min_max(R_xi)
	U_q = min_max(R_q)
	U_min = np.sqrt((R_cen-U_Z[0])**2+(R_cen-U_xi[0])**2)
	U_max = np.sqrt((R_cen-U_Z[1])**2+(R_cen-U_xi[1])**2)
	U_tot_min = np.sqrt((R_cen-U_Z[0])**2+(R_cen-U_xi[0])**2+(R_cen-U_q[0])**2)
	U_tot_max = np.sqrt((R_cen-U_Z[1])**2+(R_cen-U_xi[1])**2+(R_cen-U_q[1])**2)
	return [R_cen,U_min,U_max,U_tot_min,U_tot_max]

# here do the tot RpA with dir photons

def RpA_dir_uncertainties(pt,q):
	'''Return the R_pA for direct photons from the sigma.py file, with a variable q'''
	x_T = 2*pt/rs
	R_cen = pPb_cross_section.Rpp_FCELG_dy(x_T, 0, q[0])
	R_q = [R_cen, pPb_cross_section.Rpp_FCELG_dy(x_T, 0, q[1]),pPb_cross_section.Rpp_FCELG_dy(x_T, 0, q[2])]
	U_q = min_max(R_q)
	return [R_cen, abs(R_cen-U_q[0]),abs(U_q[1]-R_cen)]

def RpA_tot(Y,pt,z,q,xi):
	'''Return the total R_pA = f_dir*R_pA_dir + (1-f_dir)*R_pA_frag with given:
	- Y a list of rapidities (corresponding to the range of frag fits)
	- pt the trnsverse momentum (available in p_T_dispo)
	- z the fragmentation momentum ratio
	- q the transport coefficient
	- xi the parton fraction p+_3/p+_1 = -\hat{u}/\hat{s}'''
	x_T = 2*pt/rs
	Rpa_frag = R_frag(Y,pt, z, q, xi)
	Rpa_dir = pPb_cross_section.Rpp_FCELG_dy(x_T, 0, q)
	f = f_dir_pt(pt)(Y)
	return f*Rpa_dir+(1-f)*Rpa_frag

def RpA_tot_uncertainties(Y,pt,z,q,xi):
	'''Return the total R_pA = f_dir*R_pA_dir + (1-f_dir)*R_pA_frag and its uncertainties with given:
	- pt the trnsverse momentum (available in p_T_dispo)
	- z the fragmentation momentum ratio
	- q the transport coefficient
	- xi the parton fraction p+_3/p+_1 = -\hat{u}/\hat{s}
	In the result list you will find:
	- i = 0 : Y (because the file is saved so if you want to use it elsewhere, you can)
	- i = 1,2,3: R_tot, R_tot_minus, R_tot_plus
	- i = 4,5,6: R_dir, R_dir_minus, R_dir_plus
	- i = 7,8,9: R_frag, R_frag_min, R_frag_max'''
	Ny = len(Y)
	file_name = 'RpA_tot_dir_onef_Ny'+str(Ny)+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.txt'
	if os.path.exists(os.path.join(data_dir,file_name)):
		print(f"The file '{file_name}' already exists. It is loaded.")
		r = np.loadtxt(os.path.join(data_dir,file_name))
	else:
		print("The file does not exists")
		R_frag_unc = R_frag_uncertainties(Y,pt, z, q, xi)
		R_frag = R_frag_unc[0]
		D_all_min, D_all_max = R_frag_unc[1],R_frag_unc[2]
		R_dir_unc = RpA_dir_uncertainties(pt,q)
		R_dir = R_dir_unc[0]
		f = f_dir_pt(pt)(Y)
		R_cen = f*R_dir + (1-f)*R_frag
		R_q = [R_cen, RpA_tot(Y,pt, z[0], q[1], xi[0]), RpA_tot(Y,pt, z[0], q[2], xi[0])]
		U_q = min_max(R_q)
		U_q_min , U_q_max = R_cen - U_q[0], U_q[1]-R_cen
		U_min = np.sqrt(U_q_min**2+((1-f)**2)*D_all_min**2)
		U_max = np.sqrt(U_q_max**2+((1-f)**2)*D_all_max**2)
		r = np.array([Y,R_cen, R_cen - U_min ,R_cen + U_max,R_dir,R_dir-R_dir_unc[1],R_dir+R_dir_unc[2],R_frag,R_frag-R_frag_unc[3],R_frag+R_frag_unc[4]])
		np.savetxt(os.path.join(data_dir,file_name), r)
		print(f" '{file_name}' has been created")
	return r

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)
	
def plot_sigma_spline(Y,pt):
	keys = list(Processus.keys())
	for p_num in keys:
		plt.plot(Y,Processus_pt(pt)[p_num](Y),label = p_num)
		plt.title(p_num + ' = '+ p_num_to_process[p_num]+' at pt ='+str(pt))
		plt.show()

def plot_splines(Y,pt,L,C,n=2,All = True):
	'''Return the plots of cross section splines in a Y array giving pt 
	and L the list of the p_num wanted. If all == True, then all the cross sections 
	are ploted on the same figure with n colluns for legend'''
	a=0
	b=0
	if All:
		fig, ax = plt.subplots(constrained_layout=True)
		ax.xaxis.set_minor_locator(MultipleLocator(1))
		for i,p_num in enumerate(L):
			cross = np.array([cross_section(y, pt, p_num) for y in Y])
			if (p_num =='14a') | (p_num =='14b') | (p_num =='13a')| (p_num =='13b'):
				style = '--'
			else:
				style = '-'
			plt.plot(Y,cross,label= p_num_to_process[p_num],linestyle=style,color = C[i])
		if ('14a' in L) and ('14b' in L):
			cross1 = np.array([cross_section(y, pt, '14a') for y in Y])
			cross2 = np.array([cross_section(y, pt, '14b') for y in Y])
			plt.plot(Y,cross1+cross2,label=p_num_to_process['14'] ,color='brown')
		if ('13a' in L) and ('13b' in L):
			cross1 = np.array([cross_section(y, pt, '13a') for y in Y])
			cross2 = np.array([cross_section(y, pt, '13b') for y in Y])
			plt.plot(Y,cross1+cross2,label=p_num_to_process['13'] ,color='brown') # change the color here if you want to plot both 13 and 14 process
		plt.xlabel(r'$y$',fontsize=f_size-a)
		plt.ylabel(r'$d\sigma/dy$ (pb)',fontsize=f_size-a)
		plt.yscale('log')
		plt.text(0.5, 0.45,r'$p_\bot =$ '+str(pt)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b) #bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=alph)
		plot_usuals(n=n,s1=f_size-b,s2=f_size-a)
		n_fig = 'all_splines_Ny'+str(len(Y))+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
		ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
		plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
		plt.show()
	elif not All: 																# do the same here
		for p_num in L:
			cross = np.array([cross_section(y, pt, p_num) for y in Y])
			fig, ax = plt.subplots(constrained_layout=True)
			plt.plot(Y,cross,label= p_num_to_process[p_num])
			plot_usuals(s=f_size)
			plt.xlabel(r'$y$',fontsize=f_size)
			plt.text(0.25, 0.1,r'$p_\bot =$ '+str(pt),bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=alph),horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
			n_fig = p_num+'_spline_Ny'+str(len(Y))+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
			ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
			plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
			plt.show()

def plot_f_alpha(Y,pt,L,C,n=2):
	a=0
	b=0
	fig, ax = plt.subplots(constrained_layout=True)
	ax.xaxis.set_minor_locator(MultipleLocator(1))
	ax.yaxis.set_minor_locator(MultipleLocator(0.1))
	for i,p_num in enumerate(L):
		f_alpha = f_alpha_pt(pt)[p_num](Y)
		if (p_num =='14a') | (p_num =='14b') | (p_num =='13a')| (p_num =='13b'):
			style = '--'
		else:
			style = '-'
		plt.plot(Y,f_alpha,label= p_num_to_process[p_num],linestyle=style,color = C[i])
	if ('14a' in L) and ('14b' in L):
		f1 = f_alpha_pt(pt)['14a'](Y)
		f2 = f_alpha_pt(pt)['14b'](Y)
		plt.plot(Y,f1+f2,label=p_num_to_process['14'] ,color='brown')
	if ('13a' in L) and ('13b' in L):
		f1 = f_alpha_pt(pt)['13a'](Y)
		f2 = f_alpha_pt(pt)['13b'](Y)
		plt.plot(Y,f1+f2,label=p_num_to_process['13'] ,color='brown') # change the color here if you want to plot both 13 and 14 process
	plt.axhline(y=1,color='grey',alpha=alph)
	plt.xlabel(r'$y$',fontsize=f_size-a)
	plt.ylabel(r'$f_\alpha$',fontsize=f_size-a)
	plt.ylim(0,1.1)
	plt.text(0.8, 0.35,r'$p_\bot =$ '+str(pt)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b) #bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=alph)
	plot_usuals(n=n,s1=f_size-b,s2=f_size-a,loca='upper center')
	n_fig = 'all_f_alpha_Ny'+str(len(Y))+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
	ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
	plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
	plt.show()

def pt_plots(Y,y_list,z,q,xi):
	'''Y a vector of rapidities and y_list some rapidities inside it'''
	a=0
	b=0
	Yr = np.round(Y,1)															# adapt the round parameters folowing your Y vector and your y_list
	for y in y_list:
		if not(y in Yr):
			raise ValueError(str(y)+' is not in the Y argument you gave')
		else:
			i = np.where(Yr == y)
			R_tot , R_tot_minus, R_tot_plus = [],[],[]
			R_dir , R_dir_minus, R_dir_plus = [],[],[]
			R_frag, R_frag_minus, R_frag_plus = [],[],[]
			for pt in p_T_dispo:													# Supposing that the coputation is already done, otherwise there will be a lot of calculation
				r = RpA_tot_uncertainties(Y,pt ,z, q, xi)
				r_tot , r_tot_minus, r_tot_plus = r[1:4]
				r_dir , r_dir_minus, r_dir_plus = r[4:7]
				r_frag, r_frag_minus, r_frag_plus = r[7:]
				R_tot.append(r_tot[i][0]);R_tot_minus.append(r_tot_minus[i][0]);R_tot_plus.append(r_tot_plus[i][0])
				R_dir.append(r_dir[i][0]);R_dir_minus.append(r_dir_minus[i][0]);R_dir_plus.append(r_dir_plus[i][0])
				R_frag.append(r_frag[i][0]);R_frag_minus.append(r_frag_minus[i][0]);R_frag_plus.append(r_frag_plus[i][0])
			fig, ax = plt.subplots(constrained_layout=True)
			ax.xaxis.set_minor_locator(MultipleLocator(1))
			ax.xaxis.set_major_locator(MultipleLocator(5))
			plt.axhline(y=1, color='grey', alpha=alph)
			plt.plot(p_T_dispo,R_dir,color = dir_color, label = r'$R_{pA}^{dir}$')
			plt.plot(p_T_dispo,R_frag,color = frag_color,label = r'$R_{pA}^{frag}$')
			plt.plot(p_T_dispo,R_tot,color = tot_color,label = r'$R_{pA}^{tot}$')
			plt.fill_between(p_T_dispo, R_tot_minus, R_tot_plus, color=tot_color, alpha=alph)
			plt.fill_between(p_T_dispo, R_dir_minus, R_dir_plus, color=dir_color, alpha=alph)
			plt.fill_between(p_T_dispo, R_frag_minus, R_frag_plus, color=frag_color, alpha=alph)
			plt.ylim(0.8,1.1)
			plot_usuals(s1=f_size-b,s2=f_size-a,n=2,loca='lower left')
			plt.xlabel(r'$p_\bot$ (GeV)',fontsize=f_size-a)
			plt.ylabel(r'$R_{pA}$',fontsize=f_size-a)
			plt.text(0.75, 0.9,r'$y =$ '+str(int(np.round(y))),horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b)
			plt.text(0.75, 0.8,r'$\sqrt{s} =$ '+str(rs/1000)+' TeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b)
			n_fig = 'tot_dir_frag_Npt'+str(len(p_T_dispo))+'_y'+str(y)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
			ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
			plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
			plt.show()
			fig, ax = plt.subplots(constrained_layout=True)
			ax.xaxis.set_minor_locator(MultipleLocator(1))
			ax.xaxis.set_major_locator(MultipleLocator(5))
			plt.axhline(y=1, color='grey', alpha=alph)
			plt.plot(p_T_dispo,R_dir,color = dir_color, label = r'$R_{pA}^{dir}$')
			plt.plot(p_T_dispo,R_tot,color = tot_color,label = r'$R_{pA}^{tot}$')
			plt.fill_between(p_T_dispo, R_tot_minus, R_tot_plus, color=tot_color, alpha=alph)
			plt.ylim(0.8,1.1)
			plot_usuals(s1=f_size-b,s2=f_size-a,loca='lower left')
			plt.xlabel(r'$p_\bot$ (GeV)',fontsize=f_size-a)
			plt.ylabel(r'$R_{pA}$',fontsize=f_size-a)
			plt.text(0.8, 0.85,r'$y =$'+str(int(np.round(y))),horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b)
			n_fig = 'tot_dir_Npt'+str(len(p_T_dispo))+'_y'+str(y)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
			ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
			plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
			plt.show()

def plot_mains(Y,pt,z,q,xi):
	Ny = len(Y)
	fig, ax = plt.subplots(constrained_layout=True)
	ax.xaxis.set_minor_locator(MultipleLocator(1))
	plt.axhline(y=1, color='grey', linestyle='--')
	NY2 = 100
	Y2 = np.linspace(-6,6, NY2)
	R = R_frag(Y2,pt, z[0], q[0], xi[0],plot=True)
	plt.plot(Y2,R,linestyle='--',label = r'$R_{pA}^{frag}$')
	plot_usuals(n=2)
	plt.xlabel('y',fontsize=f_size)
	ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
	n_fig = 'frag_and_mains_Ny'+str(Ny)+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
	plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")
	plt.show()
	# all the main_components with their color representations
	for p_num in p_num_plot_list:
		R_alpha = all_alpha_colors_xi(Y2, pt,p_num, xi[0], z[0],q[0])
		keys = list(R_alpha.keys())
		fig, ax = plt.subplots(constrained_layout=True)
		ax.xaxis.set_minor_locator(MultipleLocator(1))
		for color in keys:
			if color == 'tot':
				plt.plot(Y2,R_alpha[color],linestyle='--',color = 'red',label = color)
			else:
				plt.plot(Y2,R_alpha[color],label = color)
		plt.xlabel('y',fontsize=f_size)
		plt.ylabel(r'$R_{pA}^{frag,R_\alpha}$ for $\alpha=$ '+p_num_to_process[p_num],fontsize=f_size)
		plt.axhline(y=1, color='grey', linestyle='--')
		plot_usuals()
		ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
		n_fig = 'frag_and_colors_alpha'+str(p_num)+'_Ny'+str(Ny)+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
		plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")
		plt.show()
 
def all_plots(Y,pt,z,q,xi):
	a=0
	b=0
	Ny = len(Y)
	r = RpA_tot_uncertainties(Y,pt, z, q, xi)
	r_tot , r_tot_minus, r_tot_plus = r[1:4]
	r_dir , r_dir_minus, r_dir_plus = r[4:7]
	r_frag, r_frag_minus, r_frag_plus = r[7:]
	# first plot comparing R_frag, R_dir and R_tot, with all uncertainties
	fig, ax = plt.subplots(constrained_layout=True)
	ax.xaxis.set_minor_locator(MultipleLocator(1))
	plt.axhline(y=1, color='grey', alpha=alph)
	plt.plot(Y,r_dir,color = dir_color, label = r'$R_{pA}^{dir}$')
	plt.plot(Y,r_frag,color = frag_color,label = r'$R_{pA}^{frag}$')
	plt.plot(Y,r_tot,color = tot_color,label = r'$R_{pA}^{tot}$')
	plt.fill_between(Y, r_tot_minus, r_tot_plus, color=tot_color, alpha=alph)
	plt.fill_between(Y, r_dir_minus, r_dir_plus, color=dir_color, alpha=alph)
	plt.fill_between(Y, r_frag_minus, r_frag_plus, color=frag_color, alpha=alph)
	plot_usuals(n=2,s1=f_size-b,s2=f_size-a,loca = 'lower left')
	plt.ylim(0.8,1.1)
	plt.xlabel('y',fontsize=f_size-a)
	plt.ylabel(r'$R_{pA}$', fontsize=f_size-a)
	plt.text(0.75, 0.9,r'$p_\bot =$ '+str(pt)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b) #bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=alph)
	plt.text(0.75, 0.8,r'$\sqrt{s} =$ '+str(rs/1000)+' TeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b)
	ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
	n_fig = 'tot_dir_frag_Ny'+str(Ny)+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
	plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")
	plt.show()
	# RpA_tot unc only + dir in dash
	fig, ax = plt.subplots(constrained_layout=True)
	ax.xaxis.set_minor_locator(MultipleLocator(1))
	plt.axhline(y=1, color='grey', alpha=alph)
	plt.plot(Y,r_dir,color= dir_color,linestyle='--',label= r'$R_{pA}^{dir}$')
	plt.plot(Y,r_tot,color = tot_color,label = r'$R_{pA}^{tot}$')
	plt.fill_between(Y, r_tot_minus, r_tot_plus,color=tot_color,alpha=alph)
	plot_usuals(s1=f_size-b,s2=f_size-a,loca = 'lower left')
	plt.ylim(0.8,1.1)
	plt.xlabel('y',fontsize=f_size-a)
	plt.ylabel(r'$R_{pA}$',fontsize=f_size-a)
	plt.text(0.75, 0.9,r'$p_\bot =$ '+str(pt)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b) #bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=alph)
	plt.text(0.75, 0.8,r'$\sqrt{s} =$ '+str(rs/1000)+' TeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b)
	ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
	n_fig = 'tot_vs_dir_Ny'+str(Ny)+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
	plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")
	plt.show()
# 	# the frag Rpa with the main components
# 	fig, ax = plt.subplots(constrained_layout=True)
# 	ax.xaxis.set_minor_locator(MultipleLocator(1))
# 	plt.axhline(y=1, color='grey',alpha=alph)
# 	NY2 = 100
# 	Y2 = np.linspace(min(Y),max(Y), NY2)
# 	R = R_frag(Y2,pt, z[0], q[0], xi[0],plot=True)
# 	plt.plot(Y2,R,linestyle='--',label = r'$R_{pA}^{frag}$')
# 	plot_usuals(n=1,s1=f_size-b,s2=f_size-a,loca='lower left')
# 	plt.ylim(ylims[pt])
# 	plt.xlabel('y',fontsize=f_size-a)
# 	plt.text(0.85, 0.95,r'$p_\bot =$ '+str(pt)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b)
# 	ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# 	n_fig = 'frag_and_mains_Ny'+str(Ny)+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
# 	plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")
# 	plt.show()
# 	# all the main_components with their color representations
# 	for p_num in p_num_plot_list:
# 		R_alpha = all_alpha_colors_xi(Y2,pt, p_num, xi[0], z[0],q[0])
# 		keys = list(R_alpha.keys())
# 		fig, ax = plt.subplots(constrained_layout=True)
# 		ax.xaxis.set_minor_locator(MultipleLocator(1))
# 		for color in keys:
# 			if color == 'tot':
# 				plt.plot(Y2,R_alpha[color],linestyle='--',color = 'red',label = 'Sum')
# 			else:
# 				plt.plot(Y2,R_alpha[color],label = color)
# 		plt.xlabel('y',fontsize=f_size-a)
# 		plt.ylim(ylims[pt])
# 		plt.ylabel(r'$R_{pA}^{frag,R}$',fontsize=f_size-a)
# 		plt.text(0.1, 0.5, p_num_to_process[p_num],horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize=f_size-b )
# 		plt.text(0.85, 0.95,r'$p_\bot =$ '+str(pt)+' GeV',horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size-b)
# 		plt.axhline(y=1, color='grey',alpha=alph)
# 		plot_usuals(s1 = f_size-b,s2=f_size-a)
# 		ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# 		n_fig = 'frag_and_colors_alpha'+str(p_num)+'_Ny'+str(Ny)+'_pt'+str(pt)+'_rs'+str(rs)+'_A'+str(A)+'_Z'+str(Z)+'.pdf'
# 		plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")
# 		plt.show()

pt = p_T_dispo[1]
x_T = 2*pt/rs
Y_list = sig.Y_list(x_T)
pt_plots(Y_list, fixed_y, z, q, xi)
all_plots(Y_list,pt, z, q, xi)
for p_t in p_T_dispo:
	x_T = 2*p_t/rs
	Y_list = sig.Y_list(x_T)
# 	plot_splines(np.linspace(-6,6,100), p_t, p_num_plot_list,p_num_color_list,n=1)
# 	plot_f_alpha(np.linspace(-6,6,100), p_t, p_num_plot_list,p_num_color_list,n=2)
# 	all_plots(Y_list,p_t, z, q, xi)
# plot_mains(np.linspace(-6, 6,100),z, q, xi, n)
# plot_sigma_spline(np.linspace(-6, 6,100))