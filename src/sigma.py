# -*- coding: utf-8 -*-

# =============================================================================
# This code is made to compute direct photon production in pp collision with 
# the implication of FCEL/G effect to compute the pA photon production.
# It shows for different arguments the nuclear ratio R_pA, with the possibility 
# to add uncertainties on \mu, \mu_f, \hat{q} and pdf uncertainites.
# Credits : BOURGEAIS Djessy; ARLEO Francois.
# Subatech (IN2P3,CNRS)
# Nantes, France
# Last modified: 12/12/2025
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import lhapdf as lha
import numpy as np 
from scipy import integrate
try :
	from . import Probability as p
except ImportError:
	print('Caution! : running from src directory, using local import')
	import Probability as p 

# Ininitalisation of a patricle dictionary
particle = {
	"d": {"id": 1, "charge": -1./3}, 
	"u": {"id": 2, "charge": 2./3}, 
	"s": {"id": 3, "charge": -1./3}, 
	"c": {"id": 4, "charge": 2./3}, 
	"b": {"id": 5, "charge": -1./3}, 
	"t": {"id": 6, "charge": 2./3},												# Quarks

	"dbar": {"id": -1, "charge": 1./3}, 
	"ubar": {"id": -2, "charge": -2./3}, 
	"sbar": {"id": -3, "charge": 1./3}, 
	"cbar": {"id": -4, "charge": -2./3}, 
	"bbar": {"id": -5, "charge": 1./3}, 
	"tbar": {"id": -6, "charge": -2./3},										# Antiquarks

	"g": {"id": 21, "charge": 0},												# Gluon
	"gamma": {"id": 22, "charge": 0},											# Photon
	"Z0": {"id": 23, "charge": 0},												# Z boson
	"W+": {"id": 24, "charge": 1},
	"W-": {"id": -24, "charge": -1},											# W bosons 
	"h0": {"id": 25, "charge": 0},												# Higgs boson
}

id_to_charge = {v["id"]: v["charge"] for v in particle.values()}

def charge(particle_id):
	'''return the charge of a particle giving its id'''
	return(id_to_charge.get(particle_id,None))

# Variables initiallisation
N_c = 3. 							  											# number of colors 
C_F = (N_c**2 -1)/(2*N_c) 														# color factor  
alpha = 1./137.																	# alpha EM
conv_fact = pow(0.197,2)*pow(10,10)												# to convert Gev^-2 into barns
pi = np.pi

# Integration parameters
epsilon = pow(10,-1)															# to avoid phase space borders 
N_pt = 41
N_Xi = 1000
N_y= 41
N_limit = 1000																	# number of subdivisions for quad integration

# Extern functions
def x_1_tilde(y,x_T):
	return x_T*np.exp(y)/2

def x_2_tilde(y,x_T):
	return x_T*np.exp(-1*y)/2

def x_1(y,x_T,Xi):
	return x_1_tilde(y,x_T)/Xi

def x_2(y,x_T,Xi):
	return x_2_tilde(y,x_T)/(1-Xi)

Xi_min = lambda y,x_T: x_1_tilde(y,x_T)
Xi_max = lambda y,x_T: 1-x_2_tilde(y,x_T)

# for virtual photons, p_T -> M_T, so M_T/sqrt(s) = (x_T/2)*sqrt(b+1), b = (M/p_T)**2 and x_T = 2*p_T/rs
def tau(Xi,b):
	return b*Xi*(1-Xi)/(b*(1-Xi)+1)

def x_1_tilde_M(y,x_T,b):
	return x_T*np.sqrt(b+1)*np.exp(y)/2

def x_2_tilde_M(y,x_T,b):
	return x_T*np.sqrt(b+1)*np.exp(-1*y)/2

def x_1_M(y,x_T,Xi,b):
	return x_1_tilde_M(y,x_T,b)/Xi

def x_2_M(y,x_T,Xi,b):
	return x_2_tilde_M(y,x_T,b)/(1-Xi+tau(Xi,b))

Xi_min_M = lambda y,x_T,b: x_1_tilde_M(y,x_T,b)
Xi_max_M = lambda y,x_T,b: (1-x_2_tilde_M(y,x_T,b))/(1-(b/(b+1))*x_2_tilde_M(y,x_T,b))

def Y_list(x_T):
	'''Return a numpy linespace array of N_y rapidities given x_T'''
	y_min = max(-np.log(2./x_T),-6)
	y_max = min(np.log(2./x_T),6)
	Y = np.linspace(y_min,y_max,N_y)
	return Y

def Y_list_M(x_T,b):
	'''Return a numpy linespace array of N_y rapidities given x_T'''
	eps_M = np.log(b+1)/2.
	y_min = max(-np.log(2./x_T)+eps_M + epsilon,-6)
	y_max = min(np.log(2./x_T)-eps_M - epsilon,6)
	Y = np.linspace(y_min,y_max,N_y)
	return Y

def Y_list_xi(x_T,xi):
	'''Debuging version of Y_list when the integration over xi didnt work well '''
	y_min = -np.log((1-xi)/x_T)
	y_max = np.log(xi/x_T)
	return np.linspace(y_min,y_max,N_y)

def min_max(ys):
	'''return the min list and max list of a list of lists'''
	ys_array = np.vstack(ys)
	y_min = np.min(ys_array, axis=0)
	y_max = np.max(ys_array, axis=0)
	return (y_min,y_max)

# Conveniance objects
Switch = {																		# To convert in LaTeX for plots
	'dp_t' : (r'$dp_\bot$',r'$pb/GeV$'),
	'd2p_t' : (r'$d^2p_\bot$',r'$pb/GeV^2$'),
	'dp_t2' : (r'$dp_\bot^2$',r'$pb/GeV^2$'),
	}

def Collision(string):
	'''Return True(bool) if the argument is "pp", False if it is "pA"'''
	if string == 'pp':
		return True
	elif (string == 'pA') or (string == 'Ap') :
		return False
	else :
		raise ValueError("The argument must be 'pp' or 'pA'")

class Sigma:
	'''To initialize the sigma object :
		- p_set (string) is the name of the proton pdf in the lhapdf list
		- A_set (string) is the name of the nuclei pdf in the lhapdf list (put a random one for pp collisions)
		- s (float) is the Mandelstam variable for CM collisions  (in GeV^2)
		- Z (integer) is the Atomic number of the chosen nuclei
		- A (integer) is the total nucleon number of the chosen nuclei'''
	def __init__(self, p_set, A_set, s, Z, A):
		self.p_set = lha.getPDFSet(p_set)										# set of pdf for the projectile (with the pdf name as a string) 
		self.A_set = lha.getPDFSet(A_set)				 						# same for target
		self.s = float(s) 														# center of mass energy
		self.rs =np.sqrt(self.s)
		self.Z = float(Z) 														# the number of proton of the target element
		self.A = float(A)														# the number of nucleus of the target element
		self.p_list = self.p_set.mkPDFs()										# create the list of all members for the pdf set
		self.A_list = self.A_set.mkPDFs()
		self.current_p = self.p_list[0]											# Initialise the pdf member to 0
		self.current_A = self.A_list[0] 
		self.current_p_num = 0
		self.current_A_num = 0
		self.tot_dy_parameters = []												# Those next ones are juste for saving time in computation,
		self.Gq_dy_parameters = []												# for not redoing the same calculations 2 times in a row
		self.qG_dy_parameters = []
		self.qqbar_dy_parameters = []
		self.qbarq_dy_parameters = []
		self.sigma_tot_dy = ([],[])
		self.sigma_Gq_dy = ([],[])
		self.sigma_qG_dy = ([],[])
		self.sigma_qqbar_dy = ([],[])
		self.sigma_qbarq_dy = ([],[])
		self.tot_dpt_parameters = []
		self.Gq_dpt_parameters = []
		self.qG_dpt_parameters = []
		self.qqbar_dpt_parameters = []
		self.qbarq_dpt_parameters = []
		self.sigma_tot_dpt = ([],[])
		self.sigma_Gq_dpt = ([],[])
		self.sigma_qG_dpt = ([],[])
		self.sigma_qqbar_dpt = ([],[])
		self.sigma_qbarq_dpt = ([],[])
		self.tot_dy_xi_parameters = []
		self.Gq_dy_xi_parameters = []
		self.qG_dy_xi_parameters = []
		self.qqbar_dy_xi_parameters = []
		self.qbarq_dy_xi_parameters = []
		self.sigma_tot_dy_xi = []
		self.sigma_Gq_dy_xi = []
		self.sigma_qG_dy_xi = []
		self.sigma_qqbar_dy_xi = []
		self.sigma_qbarq_dy_xi = []
		
	### init functions ###
	
	def pdf_p(self,num):
		'''Redefine the current pdf of the projectile (proton) for a certain 
		member number of the set (if it is different from the current one) 
		and return the new current proton pdf set'''
		Num = self.current_p_num
		size = self.p_set.size
		if num < 0 or num > size:
			raise ValueError("num should be between 0 and "+str(size))
		if num != Num:
			p = self.p_list														# verify that the member belongs to the list
			self.current_p_num = num
			self.current_p = p[num]
			return self.current_p
		else:
			return self.current_p
	
	def pdf_A(self,num):
		'''Redefine the current pdf of the target (nuclei) for a certain member
		number of the set (if it is different from the current one) 
		and return the new current nuclei pdf set'''
		Num = self.current_A_num
		size = self.p_set.size
		if num < 0 or num > size:
			raise ValueError("num should be between 0 and "+str(size))
		if num != Num:
			A = self.A_list 													# verify that the member belongs to the list 
			self.current_A_num = num
			self.current_A = A[num]
			return self.current_A
		else:	
			return self.current_A
	
	def alpha_s_A(self,num,mu2):
		'''return alpha_s for the nuclei at mu2'''
		A_num = self.pdf_A(num)
		return A_num.alphasQ2(mu2)
	
	def alpha_s_p(self,num,mu2):
		'''return alpha_s for the proton at mu2'''
		p_num = self.pdf_p(num)
		return p_num.alphasQ2(mu2)
		
	def n_xfxQ2(self,particle_id,x,mu_f2,num):
		'''return the pdf of the neutron with the mu_f2 parameter'''
		p = self.pdf_p(num)
		if np.abs(particle_id) >= 3:
			return p.xfxQ2(particle_id,x,mu_f2)
		else :
			return p.xfxQ2(2//particle_id,x,mu_f2)								# Math trick using euclide division to make 1->2 and 2->1
	
	def F2_p(self,x,mu_f2,num,iso = 'p',n_f=3):
		'''gives back the F2 function of the proton for an x and a mu_f2 given.
		- the iso argument stands for isospin, i.e chose if it uses proton
		or neutron pdf 
		- the n_f argument stands for the number of flavours taken in count
		(by default, it is 3: u,d,s)'''
		F = 0.
		if iso == 'p':
			p = self.pdf_p(num)
			for i in range(1,n_f+1):
				F += (charge(i)**2)*(p.xfxQ2(i,x,mu_f2)+p.xfxQ2(-i,x,mu_f2))
		elif iso == 'n':
			for i in range(1,n_f+1):
				F += (charge(i)**2)*(self.n_xfxQ2(i,x,mu_f2,num)+self.n_xfxQ2(-i,x,mu_f2,num))
		return F/x
	
	def F2_A(self,x,mu_f2,num,n_f=3):
		'''gives back the F2 function of the target nuclide for an x and a mu_f2 given
		- the n_f argument stands for the number of flavours taken in count
		(by default, it is 3: u,d,s)'''
		F = 0.
		A = self.pdf_A(num)
		for i in range(1,n_f+1):
			F += (charge(i)**2)*(A.xfxQ2(i,x,mu_f2)+A.xfxQ2(-i,x,mu_f2))
		return F/x
	
	def F_ij(self,xa,xb,mu_f2,num, direction = 'qqbar',iso='p',n_f=3,is_pp = False):
		'''return the F_ij function that appears in the annihilation process'''
		F = 0.
		p = self.pdf_p(num)
		A = self.pdf_A(num)
		if iso == 'n':
			if direction == 'qqbar':
				for i in range(1,n_f+1):
					val1 = p.xfxQ2(i, xa, mu_f2)
					val2 = self.n_xfxQ2(-i, xb, mu_f2,num)
					if val1 is None or val2 is None:
						raise ValueError("xfxQ2 returned None for i={}, xa={}, xb={}, mu_f2={}".format(i, xa, xb, mu_f2))
					F += (charge(i)**2) * val1 * val2
				if xa == 0 or xb == 0:
					raise ValueError("xa or xb is zero, leading to division by zero.")
				return F / (xa * xb)
			elif direction == 'qbarq':
				for i in range(1,n_f+1):
					val1 = p.xfxQ2(-i, xa, mu_f2)
					val2 = self.n_xfxQ2(i, xb, mu_f2,num)
					if val1 is None or val2 is None:
						raise ValueError("xfxQ2 returned None for i={}, xa={}, xb={}, mu_f2={}".format(i, xa, xb, mu_f2))
					F += (charge(i)**2) * val1 * val2
				if xa == 0 or xb == 0:
					raise ValueError("xa or xb is zero, leading to division by zero.")
				return F / (xa * xb)
		elif iso == 'p':
			if not is_pp:
				if direction == 'qqbar':
					for i in range(1,n_f+1):
						val1 = p.xfxQ2(i, xa, mu_f2)
						val2 = A.xfxQ2(-i, xb, mu_f2)
						if val1 is None or val2 is None:
							raise ValueError("xfxQ2 returned None for i={}, xa={}, xb={}, mu_f2={}".format(i, xa, xb, mu_f2))
						F += (charge(i)**2) * val1 * val2
				elif direction == 'qbarq':
					for i in range(1,n_f+1):
						val1 = p.xfxQ2(-i, xa, mu_f2)
						val2 = A.xfxQ2(i, xb, mu_f2)
						if val1 is None or val2 is None:
							raise ValueError("xfxQ2 returned None for i={}, xa={}, xb={}, mu_f2={}".format(i, xa, xb, mu_f2))
						F += (charge(i)**2) * val1 * val2
				if xa == 0 or xb == 0:
					raise ValueError("xa or xb is zero, leading to division by zero.")
				return F / (xa * xb) 
			elif is_pp:
				A = self.pdf_p(num)
				if direction == 'qqbar':
					for i in range(1,n_f+1):
						val1 = p.xfxQ2(i, xa, mu_f2)
						val2 = A.xfxQ2(-i, xb, mu_f2)
						if val1 is None or val2 is None:
							raise ValueError("xfxQ2 returned None for i={}, xa={}, xb={}, mu_f2={}".format(i, xa, xb, mu_f2))
						F += (charge(i)**2) * val1 * val2
				elif direction == 'qbarq':
					for i in range(1,n_f+1):
						val1 = p.xfxQ2(-i, xa, mu_f2)
						val2 = A.xfxQ2(i, xb, mu_f2)
						if val1 is None or val2 is None:
							raise ValueError("xfxQ2 returned None for i={}, xa={}, xb={}, mu_f2={}".format(i, xa, xb, mu_f2))
						F += (charge(i)**2) * val1 * val2
				if xa == 0 or xb == 0:
					raise ValueError("xa or xb is zero, leading to division by zero.")
				return F / (xa * xb)
		else:
			raise ValueError("Iso type '{}' is not recognized. Expected 'p' or 'n'.".format(iso))
	
	def Gluon_p(self,x,mu_f2,num):
		'''Return the gluon pdf of the proton'''
		p = self.pdf_p(num)
		return p.xfxQ2(particle["g"]["id"],x,mu_f2)/x
	
	def Gluon_A(self,x,mu_f2,num):
		'''Return the gluon pdf of the nuclei'''
		A = self.pdf_A(num)
		return A.xfxQ2(particle["g"]["id"],x,mu_f2)/x
	
	def jac_xi(self,Xi,b):
		return 1/(self.rs*Xi*(1-Xi+tau(Xi,b)))
	
	### integrand functions ###
	
	def qG(self,y,x_T,Xi,num,mu_factor=1,mu_f_factor=1,n_f=3,is_pp = False,switch = 'dp_t'): 				# not affected by isospin because G_p  = G_n 
		'''Return the q(p)G(A)-> gamma q integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		x_proj = x_1(y,x_T,Xi)
		x_targ = x_2(y,x_T,Xi)
		Xi_factor = 1-Xi+(1./(1-Xi))
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		mu2 = (p_t*mu_factor)**2
		mu_f2 = (p_t*mu_f_factor)**2
		alpha_s = self.alpha_s_p(num,mu2)
		if switch == 'dp_t':
			prefactor = 4*pi*alpha*alpha_s/(pow(s,1.5)*x_T)
		elif switch == 'd2p_t':
			prefactor = 4*alpha*alpha_s/((s*x_T)**2)
		elif switch == 'dp_t2':
			prefactor = 4*pi*alpha*alpha_s/((s*x_T)**2)
		if not is_pp:
			F = self.F2_p(x_proj,mu_f2,num,iso ='p',n_f=n_f)
			G = self.Gluon_A(x_targ,mu_f2,num)
		elif is_pp:
			F = self.F2_p(x_proj,mu_f2,num,iso='p',n_f=n_f)
			G = self.Gluon_p(x_targ,mu_f2,num)
		return F*G*Xi_factor*prefactor/N_c
	
	def qG_M(self,y,x_T,Xi,M,num,mu_factor=1,mu_f_factor=1,n_f=3,is_pp = False,switch = 'dp_t'): 				# not affected by isospin because G_p  = G_n 
		'''Return the q(p)G(A)-> gamma^\star q integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		M_t = np.sqrt(M**2+p_t**2)
		b = (M/p_t)**2 #if M = 0, b = 0 and then you get the same formula as for real photons
		x_proj = x_1_M(y,x_T,Xi,b)
		x_targ = x_2_M(y,x_T,Xi,b)
		hat_s = s*x_proj*x_targ
		Xi_factor = 1-Xi+(1./(1-Xi))-2*tau(Xi,b)*(Xi-tau(Xi,b))/(1-Xi)
		mu2 = (M_t*mu_factor)**2
		mu_f2 = (M_t*mu_f_factor)**2
		dsigma_dxi = Xi_factor/hat_s
		alpha_s = self.alpha_s_p(num,mu2)
		if switch == 'dp_t':
			prefactor = pi*alpha*alpha_s*x_T*rs
		elif switch == 'd2p_t':
			prefactor = alpha*alpha_s
		elif switch == 'dp_t2':
			prefactor = pi*alpha*alpha_s
		if not is_pp:
			F = self.F2_p(x_proj,mu_f2,num,iso ='p',n_f=n_f)
			G = self.Gluon_A(x_targ,mu_f2,num)
		elif is_pp:
			F = self.F2_p(x_proj,mu_f2,num,iso='p',n_f=n_f)
			G = self.Gluon_p(x_targ,mu_f2,num)
		return F*G*dsigma_dxi*prefactor*self.jac_xi(Xi,b)/N_c
		
	def Gq(self,y,x_T,Xi,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch ='dp_t'):  
		'''Return the G(p)q(A)-> gamma q integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		x_proj = x_1(y,x_T,Xi)
		x_targ = x_2(y,x_T,Xi)
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		mu2 = (p_t*mu_factor)**2
		mu_f2 = (p_t*mu_f_factor)**2
		alpha_s = self.alpha_s_p(num,mu2)
		Xi_factor = Xi+(1./Xi)
		if switch == 'dp_t':
			prefactor = 4*pi*alpha*alpha_s/(pow(s,1.5)*x_T)
		elif switch == 'd2p_t':
			prefactor = 4*alpha*alpha_s/((s*x_T)**2)
		elif switch == 'dp_t2':
			prefactor = 4*pi*alpha*alpha_s/((s*x_T)**2)
		if not is_pp:
			F = self.F2_A(x_targ,mu_f2,num,n_f)
			G = self.Gluon_p(x_proj,mu_f2,num)
		elif is_pp:
			F = self.F2_p(x_targ,mu_f2,num,iso,n_f) 					# here, the target could be a neutron, so we put an iso parameter
			G = self.Gluon_p(x_proj,mu_f2,num)
		return F*G*Xi_factor*prefactor/N_c
	
	def Gq_M(self,y,x_T,Xi,M,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch ='dp_t'):  
		'''Return the G(p)q(A)-> gamma^\star q integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		M_t = np.sqrt(M**2+p_t**2)
		b = (M/p_t)**2 #if M = 0, b = 0 and then you get the same formula as for real photons
		x_proj = x_1_M(y,x_T,Xi,b)
		x_targ = x_2_M(y,x_T,Xi,b)
		hat_s = s*x_proj*x_targ
		mu2 = (M_t*mu_factor)**2
		mu_f2 = (M_t*mu_f_factor)**2
		Xi_factor = Xi - tau(Xi,b) + 1/(Xi-tau(Xi,b)) - 2*tau(Xi,b)*(1-Xi)/(Xi-tau(Xi,b))
		dsigma_dxi = Xi_factor/hat_s
		alpha_s = self.alpha_s_p(num,mu2)
		if switch == 'dp_t':
			prefactor = pi*alpha*alpha_s*x_T*rs
		elif switch == 'd2p_t':
			prefactor = alpha*alpha_s
		elif switch == 'dp_t2':
			prefactor = pi*alpha*alpha_s
		if not is_pp:
			F = self.F2_p(x_proj,mu_f2,num,iso ='p',n_f=n_f)
			G = self.Gluon_A(x_targ,mu_f2,num)
		elif is_pp:
			F = self.F2_p(x_proj,mu_f2,num,iso='p',n_f=n_f)
			G = self.Gluon_p(x_targ,mu_f2,num)
		return F*G*dsigma_dxi*self.jac_xi(Xi,b)*prefactor*N_c
		
	def qqbar(self,y,x_T,Xi,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp = False,switch = 'dp_t'):
		'''Return the q(p)qbar(A)-> gamma G integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		x_proj = x_1(y,x_T,Xi)
		x_targ = x_2(y,x_T,Xi)
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		mu2 = (p_t*mu_factor)**2
		mu_f2 = (p_t*mu_f_factor)**2
		alpha_s = self.alpha_s_p(num,mu2)
		if switch == 'dp_t':
			prefactor = 4*pi*alpha*alpha_s/(pow(s,1.5)*x_T)
		elif switch == 'd2p_t':
			prefactor = 4*alpha*alpha_s/((s*x_T)**2)
		elif switch == 'dp_t2':
			prefactor = 4*pi*alpha*alpha_s/((s*x_T)**2)
		Xi_factor = (Xi/(1-Xi)) + ((1-Xi)/Xi)
		F_qqbar = self.F_ij(x_proj,x_targ,mu_f2,num,direction='qqbar',iso=iso,n_f=n_f,is_pp=is_pp)
		return F_qqbar*Xi_factor*(2*C_F/N_c)*prefactor
	
	def qqbar_M(self,y,x_T,Xi,M,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp = False,switch = 'dp_t'):
		'''Return the q(p)qbar(A)-> gamma G integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		M_t = np.sqrt(M**2+p_t**2)
		b = (M/p_t)**2 #if M = 0, b = 0 and then you get the same formula as for real photons
		tau = b*Xi*(1-Xi)/(b*(1-Xi)+1)
		x_proj = (x_T*np.sqrt(b+1)*np.exp(y)/2)/Xi
		x_targ = (x_T*np.sqrt(b+1)*np.exp(-1*y)/2)/(1-Xi+tau)
		hat_s = s*x_proj*x_targ
		mu2 = (M_t*mu_factor)**2
		mu_f2 = (M_t*mu_f_factor)**2
		Xi_factor = (1-Xi)/(Xi-tau) +(Xi-tau)/(1-Xi) + 2*tau/((Xi-tau)*(1-Xi))
		dsigma_dxi = Xi_factor/hat_s
		alpha_s = self.alpha_s_p(num,mu2)
		if switch == 'dp_t':
			prefactor = pi*alpha*alpha_s*x_T*rs
		elif switch == 'd2p_t':
			prefactor = alpha*alpha_s
		elif switch == 'dp_t2':
			prefactor = pi*alpha*alpha_s
		F_qqbar = self.F_ij(x_proj,x_targ,mu_f2,num,direction='qqbar',iso=iso,n_f=n_f,is_pp=is_pp)
		return F_qqbar*dsigma_dxi*self.jac_xi(Xi,b)*(2*C_F/N_c)*prefactor
	
	def qbarq(self,y,x_T,Xi,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp = False, switch = 'dp_t'): # y,x_T,xi,mu_f2,num,n_f, is_pp
		'''Return the qbar(p)q(A)-> gamma G integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		x_proj = x_1(y,x_T,Xi)
		x_targ = x_2(y,x_T,Xi)
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		mu2 = (p_t*mu_factor)**2
		mu_f2 = (p_t*mu_f_factor)**2
		alpha_s = self.alpha_s_p(num,mu2)
		if switch == 'dp_t':
			prefactor = 4*pi*alpha*alpha_s/(pow(s,1.5)*x_T)
		elif switch == 'd2p_t':
			prefactor = 4*alpha*alpha_s/((s*x_T)**2)
		elif switch == 'dp_t2':
			prefactor = 4*pi*alpha*alpha_s/((s*x_T)**2)
		Xi_factor = (Xi/(1-Xi)) + ((1-Xi)/Xi)
		F_qbarq = self.F_ij(x_proj,x_targ,mu_f2,num,direction='qbarq',iso=iso,n_f=n_f,is_pp=is_pp)
		return F_qbarq*Xi_factor*(2*C_F/N_c)*prefactor
	
	def qbarq_M(self,y,x_T,Xi,M,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp = False, switch = 'dp_t'): # y,x_T,xi,mu_f2,num,n_f, is_pp
		'''Return the q(p)qbar(A)-> gamma G integrand with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		s = self.s
		rs = self.rs 
		p_t = x_T*rs/2.
		M_t = np.sqrt(M**2+p_t**2)
		b = (M/p_t)**2 #if M = 0, b = 0 and then you get the same formula as for real photons
		tau = b*Xi*(1-Xi)/(b*(1-Xi)+1)
		x_proj = (x_T*np.sqrt(b+1)*np.exp(y)/2)/Xi
		x_targ = (x_T*np.sqrt(b+1)*np.exp(-1*y)/2)/(1-Xi+tau)
		hat_s = s*x_proj*x_targ
		mu2 = (M_t*mu_factor)**2
		mu_f2 = (M_t*mu_f_factor)**2
		Xi_factor = (1-Xi)/(Xi-tau) +(Xi-tau)/(1-Xi) + 2*tau/((Xi-tau)*(1-Xi))
		dsigma_dxi = Xi_factor/hat_s
		alpha_s = self.alpha_s_p(num,mu2)
		if switch == 'dp_t':
			prefactor = pi*alpha*alpha_s*x_T*rs
		elif switch == 'd2p_t':
			prefactor = alpha*alpha_s
		elif switch == 'dp_t2':
			prefactor = pi*alpha*alpha_s
		F_qqbar = self.F_ij(x_proj,x_targ,mu_f2,num,direction='qbarq',iso=iso,n_f=n_f,is_pp=is_pp)
		return F_qqbar*dsigma_dxi*self.jac_xi(Xi,b)*(2*C_F/N_c)*prefactor
	
	def all_process_integrand(self,y,x_T,Xi,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the total pA (or pn, or pp if is_pp = True) collision 
		integrand for single gamma production. But be carefull for the pn 
		collisions, (curently,) you have to init your sigma with 2 proton pdf :
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		I_Gq = self.Gq(y,x_T,Xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		I_qG = self.qG(y,x_T,Xi,num,mu_factor,mu_f_factor,n_f,is_pp,switch) 					# not concerned for isospin effects 
		I_qqbar = self.qqbar(y,x_T,Xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		I_qbarq = self.qbarq(y,x_T,Xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		return (I_Gq + I_qG + I_qqbar + I_qbarq)
	
	def all_process_integrand_M(self,y,x_T,Xi,M,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the total pA (or pn, or pp if is_pp = True) collision 
		integrand for single gamma^\star production. But be carefull for the pn 
		collisions, (curently,) you have to init your sigma with 2 proton pdf :
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- Xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		I_Gq = self.Gq_M(y,x_T,Xi,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		I_qG = self.qG_M(y,x_T,Xi,M,num,mu_factor,mu_f_factor,n_f,is_pp,switch) 					# not concerned for isospin effects 
		I_qqbar = self.qqbar_M(y,x_T,Xi,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		I_qbarq = self.qbarq_M(y,x_T,Xi,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		return (I_Gq + I_qG + I_qqbar + I_qbarq)
		
	def dsigma_tot_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the total cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		Integrand = lambda xi: self.all_process_integrand(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min(y,x_T),Xi_max(y,x_T), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_tot_dydpt_M(self,y,x_T,M,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the total cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		rs = self.rs
		p_t = rs*x_T/2
		b = (M/p_t)**2
		Integrand = lambda xi: self.all_process_integrand_M(y,x_T,xi,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min_M(y,x_T,b),Xi_max_M(y,x_T,b), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_qG_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor=1,n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the qG cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		Integrand = lambda xi: self.qG(y,x_T,xi,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min(y,x_T),Xi_max(y,x_T), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_qG_dydpt_M(self,y,x_T,M,num,mu_factor=1,mu_f_factor=1,n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the qG cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		rs = self.rs
		p_t = rs*x_T/2
		b = (M/p_t)**2
		Integrand = lambda xi: self.qG_M(y,x_T,xi,M,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		# print('(xi_min,xi_max) = '+ str((Xi_min_M(y,x_T,b),Xi_max_M(y,x_T,b))))
		sigma, err = integrate.quad(Integrand,Xi_min_M(y,x_T,b),Xi_max_M(y,x_T,b), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_Gq_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the Gq cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		Integrand = lambda xi: self.Gq(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min(y,x_T),Xi_max(y,x_T), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_Gq_dydpt_M(self,y,x_T,M,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the Gq cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		rs = self.rs
		p_t = rs*x_T/2
		b = (M/p_t)**2
		Integrand = lambda xi: self.Gq_M(y,x_T,xi,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min_M(y,x_T,b),Xi_max_M(y,x_T,b), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_qqbar_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the qqbar cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		Integrand = lambda xi: self.qqbar(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min(y,x_T),Xi_max(y,x_T), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_qqbar_dydpt_M(self,y,x_T,M,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the qqbar cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		rs = self.rs
		p_t = rs*x_T/2
		b = (M/p_t)**2
		Integrand = lambda xi: self.qqbar_M(y,x_T,xi,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min_M(y,x_T,b),Xi_max_M(y,x_T,b), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_qbarq_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the qbarq cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		Integrand = lambda xi: self.qbarq(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min(y,x_T),Xi_max(y,x_T), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	def dsigma_qbarq_dydpt_M(self,y,x_T,M,num,mu_factor=1,mu_f_factor=1,iso = 'p',n_f = 3, is_pp = False,switch = 'dp_t'):
		'''Return the qbarq cross section for a point of the phase space (y,x_T)
		integrated over xi with its integration uncertainites, with:
		- y, the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		rs = self.rs
		p_t = rs*x_T/2
		b = (M/p_t)**2
		Integrand = lambda xi: self.qbarq_M(y,x_T,xi,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		sigma, err = integrate.quad(Integrand,Xi_min_M(y,x_T,b),Xi_max_M(y,x_T,b), limit = N_limit)
		return (conv_fact*sigma, conv_fact*err)
	
	### integration part ###
	# Rapidity integration
	# no urgent need here to do the virtual photon 
	
	def dsigma_tot_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch = 'dp_t'):
		'''Return the list of the total cross section rapidity dependent and its
		integration uncertainties
		(an amelioration would be to add those on pdf) with :
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.tot_dy_parameters
		if (sP == []) or ( P != sP):											# Re-compute only if it's the first time, or if parameters are differents from the last ones
			self.tot_dy_parameters = P
			Y = Y_list(x_T)
			sigma = []
			err_sigma = []
			for y in Y:
				Int = self.dsigma_tot_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_tot_dy = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_tot_dy
	
	def dsigma_qG_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,n_f=3,is_pp = False,switch = 'dp_t'):
		'''Return the list of the qG cross section rapidity dependent
		and its integration uncertainties, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,num,mu_factor,mu_f_factor,n_f,is_pp,switch]
		sP = self.qG_dy_parameters
		if (sP == []) or (P != sP):
			self.qG_dy_parameters = P
			Y = Y_list(x_T)
			sigma = []
			err_sigma = []
			for y in Y:
				Int = self.dsigma_qG_dydpt(y,x_T,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_qG_dy = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_qG_dy
	
	def dsigma_Gq_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch= 'dp_t'):
		'''Return the list of the Gq cross section rapidity dependent
		and its integration uncertainties, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.Gq_dy_parameters
		if (sP == []) or (P != sP):
			self.Gq_dy_parameters = P
			Y = Y_list(x_T)
			sigma = []
			err_sigma = []
			for y in Y:
				Int = self.dsigma_Gq_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_Gq_dy = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_Gq_dy
	
	def dsigma_qqbar_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp = False,switch='dp_t'):
		'''Return the list of the qqbar cross section rapidity dependent
		and its integration uncertainties, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.qqbar_dy_parameters
		if (sP == []) or (P != sP):
			self.qqbar_dy_parameters = P
			Y = Y_list(x_T)
			sigma = []
			err_sigma = []
			for y in Y:
				Int = self.dsigma_qqbar_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_qqbar_dy = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_qqbar_dy
	
	def dsigma_qbarq_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch='dp_t'):
		'''Return the list of the qbarq cross section rapidity dependent
		and its integration uncertainties, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.qbarq_dy_parameters
		if (sP == []) or (P != sP):
			self.qbarq_dy_parameters = P
			Y = Y_list(x_T)
			sigma = []
			err_sigma = []
			for y in Y:
				Int = self.dsigma_qbarq_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_qbarq_dy = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_qbarq_dy
	
	def dsimga_all_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch= 'dp_t'):
		'''Return all cross section rapidity dependant components seperately,
		and their integration uncertainties with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as :
		[qG,Gq,qqbar,qbarq] with e.g. qG = (sigma_qG_dy, I_uncertainties_qG_dy)'''
		qG = self.dsigma_qG_dy(x_T,num,mu_factor,mu_f_factor,n_f,is_pp, switch)
		Gq = self.dsigma_Gq_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qqbar = self.dsigma_qqbar_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qbarq = self.dsigma_qbarq_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		return[qG,Gq,qqbar,qbarq]
		
	def dsigma_FCELG_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch= 'dp_t'):
		'''Return the components of the total cross section that contribute to
		FCEL and FCEG effects and their integration uncertainites, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as :
		[FCEL,FCEG] with FCEL = Gq+ qqbar+qbarq and FCEG = qG'''
		FCEG = self.dsigma_qG_dy(x_T,num,mu_factor,mu_f_factor,n_f,is_pp, switch)
		Gq = self.dsigma_Gq_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qqbar = self.dsigma_qqbar_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qbarq = self.dsigma_qbarq_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		FCEL = (Gq[0]+qqbar[0]+qbarq[0],Gq[1]+qqbar[1]+qbarq[1])
		return [FCEL,FCEG]
		
	# debug par at fixed xi
	
	def dsigma_tot_dy_xi(self,x_T,xi,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch = 'dp_t'):
		'''Return the list of the total cross section rapidity dependent at a 
		fixed xi, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.tot_dy_xi_parameters
		if (sP == []) or ( P != sP):											# Re-compute only if it's the first time, or if parameters are differents from the last ones
			self.tot_dy_xi_parameters = P
			Y = Y_list_xi(x_T,xi)
			sigma = []
			for y in Y:
				Integrand = self.all_process_integrand(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
				sigma.append(np.multiply(conv_fact,Integrand))
			self.sigma_tot_dy_xi = np.array(sigma)
			return np.array(sigma)
		else:
			return self.sigma_tot_dy
	
	def dsigma_qG_dy_xi(self,x_T,xi,num,mu_factor=1,mu_f_factor=1,n_f=3,is_pp = False,switch = 'dp_t'):
		'''Return the list of the qG cross section rapidity dependent at a
		fixed xi, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,xi,num,mu_factor,mu_f_factor,n_f,is_pp,switch]
		sP = self.qG_dy_xi_parameters
		if (sP == []) or (P != sP):
			self.qG_dy_xi_parameters = P
			Y = Y_list_xi(x_T,xi)
			sigma = []
			for y in Y:
				Integrand = self.qG(y,x_T,xi,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
				sigma.append(conv_fact*Integrand)
			self.sigma_qG_dy_xi = np.array(sigma)
			return np.array(sigma)
		else:
			return self.sigma_qG_dy_xi
	
	def dsigma_Gq_dy_xi(self,x_T,xi,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch= 'dp_t'):
		'''Return the list of the Gq cross section rapidity dependent at a
		fixed xi, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.Gq_dy_xi_parameters
		if (sP == []) or (P != sP):
			self.Gq_dy_xi_parameters = P
			Y = Y_list_xi(x_T,xi)
			sigma = []
			for y in Y:
				Integrand = self.Gq(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
				sigma.append(conv_fact*Integrand)
			self.sigma_Gq_dy_xi = np.array(sigma)
			return np.array(sigma)
		else:
			return self.sigma_Gq_dy_xi
	
	def dsigma_qqbar_dy_xi(self,x_T,xi,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp = False,switch='dp_t'):
		'''return the list of the qqbar cross section rapidity dependent at a
		fixed xi, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.qqbar_dy_xi_parameters
		if (sP == []) or (P != sP):
			self.qqbar_dy_xi_parameters = P
			Y = Y_list_xi(x_T,xi)
			sigma = []
			for y in Y:
				Integrand = self.qqbar(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
				sigma.append(conv_fact*Integrand)
			self.sigma_qqbar_dy_xi = np.array(sigma)
			return np.array(sigma)
		else:
			return self.sigma_qqbar_dy_xi
	
	def dsigma_qbarq_dy_xi(self,x_T,xi,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch='dp_t'):
		'''return the list of the qbarq cross section rapidity dependent at a 
		fixed xi, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.qbarq_dy_xi_parameters
		if (sP == []) or (P != sP):
			self.qbarq_dy_xi_parameters = P
			Y = Y_list_xi(x_T,xi)
			sigma = []
			for y in Y:
				Integrand = self.qbarq(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f, is_pp,switch)
				sigma.append(conv_fact*Integrand)
			self.sigma_qbarq_dy_xi = np.array(sigma)
			return np.array(sigma)
		else:
			return self.sigma_qbarq_dy_xi
	
	def dsimga_all_dy_xi(self,x_T,xi,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch= 'dp_t'):
		'''Return all cross section rapidity dependant components seperately,
		with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as :
		[qG,Gq,qqbar,qbarq]'''
		qG = self.dsigma_qG_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,n_f,is_pp, switch)
		Gq = self.dsigma_Gq_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qqbar = self.dsigma_qqbar_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qbarq = self.dsigma_qbarq_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		return[qG,Gq,qqbar,qbarq]
		
	def dsigma_FCELG_dy_xi(self,x_T,xi,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch= 'dp_t'):
		''''Return the components of the total cross section that contribute to
		FCEL and FCEG effects, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as :
		[FCEL,FCEG] with FCEL = Gq+ qqbar+qbarq and FCEG = qG'''
		FCEG = self.dsigma_qG_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,n_f,is_pp, switch)
		Gq = self.dsigma_Gq_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qqbar = self.dsigma_qqbar_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		qbarq = self.dsigma_qbarq_dy_xi(x_T,xi,num,mu_factor,mu_f_factor,iso,n_f,is_pp ,switch)
		FCEL = Gq+qqbar+qbarq
		return[FCEL,FCEG]
	
	# p_T integration
	
	def dsigma_tot_dpt(self,y,num,mu_factor = 1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch ='dp_t'):
		'''Return the list of the total cross section transverse momentum dependent
		and its integration uncertainties, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.tot_dpt_parameters 
		if (sP == []) or ( P != sP):											# Re-compute only if it's the first time, or if parameters are differents from the last ones
			self.tot_pt_parameters = P
			rs = self.rs
			P_T = self.P_T_list(y)
			sigma = []
			err_sigma = []
			for p_t in P_T:
				x_T = 2.*p_t/rs
				Int = self.dsigma_tot_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_tot_dpt = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_tot_dpt
	
	def dsigma_qG_dpt(self,y,num,mu_factor = 1,mu_f_factor=1,n_f=3,is_pp = False,switch='dp_t'):
		'''Return the list of the qG cross section transverse momentum dependent
		and its integration uncertainties, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [y,num,mu_factor,mu_f_factor,n_f,is_pp,switch]
		sP = self.qG_dpt_parameters
		if (sP == []) or (P != sP):
			self.qG_dpt_parameters = P
			s = self.s
			P_T = self.P_T_list(y)
			sigma = []
			err_sigma = []
			for p_t in P_T:
				x_T = 2.*p_t/pow(s,0.5)
				Int = self.dsigma_qG_dydpt(y,x_T,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_qG_dpt = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_qG_dpt
	
	def dsigma_Gq_dpt(self,y,num,mu_factor = 1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch='dp_t'):
		'''Return the list of the Gq cross section transverse momentum dependent
		and its integration uncertainties, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.Gq_dpt_parameters
		if (sP == []) or (P != sP):
			self.Gq_dpt_parameters = P
			s = self.s
			P_T = self.P_T_list(y)
			sigma = []
			err_sigma = []
			for p_t in P_T:
				x_T = 2.*p_t/pow(s,0.5)
				Int = self.dsigma_Gq_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_Gq_dpt = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_Gq_dpt
	
	def dsigma_qqbar_dpt(self,y,num,mu_factor = 1,mu_f_factor=1,iso='p',n_f=3,is_pp = False ,switch='dp_t'):
		'''Return the list of the qqbar cross section transverse momentum dependent
		and its integration uncertainties, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.qqbar_dpt_parameters
		if (sP == []) or (P != sP):
			self.qqbar_dpt_parameters = P
			s = self.s
			P_T = self.P_T_list(y)
			sigma = []
			err_sigma = []
			for p_t in P_T:
				x_T = 2.*p_t/pow(s,0.5)
				Int = self.dsigma_qqbar_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_qqbar_dpt = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_qqbar_dpt
	
	def dsigma_qbarq_dpt(self,y,num,mu_factor = 1,mu_f_factor=1,iso ='p',n_f=3,is_pp = False,switch ='dp_t'):
		'''Return the list of the qbarq cross section transverse momentum dependent
		and its integration uncertainties, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		P = [y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch]
		sP = self.qbarq_dpt_parameters
		if (sP == []) or (P != sP):
			self.qbarq_dpt_parameters = P
			s = self.s
			P_T = self.P_T_list(y)
			sigma = []
			err_sigma = []
			for p_t in P_T:
				x_T = 2.*p_t/pow(s,0.5)
				Int = self.dsigma_qbarq_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso ,n_f,is_pp,switch)
				sigma.append(Int[0])
				err_sigma.append(Int[1])
			self.sigma_qbarq_dpt = (np.array(sigma),np.array(err_sigma))
			return (np.array(sigma),np.array(err_sigma))
		else:
			return self.sigma_qbarq_dpt
	
	def sigmal_all_dpt(self,y,num,iso ='p',n_f=3,is_pp = False,switch ='dp_t',mu_factor = 1,mu_f_factor=1):
		'''Return all cross section transverse momentum dependant components seperately,
		and their integration uncertainties with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as :
		[qG,Gq,qqbar,qbarq] with e.g. qG = (sigma_qG_dpt, I_uncertainties_qG_dpt)'''
		qG = self.dsigma_qG_dpt(y,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		Gq = self.dsigma_Gq_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qqbar = self.dsigma_qqbar_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qbarq = self.dsigma_qbarq_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		return[qG,Gq,qqbar,qbarq]
		
	def dsigma_FCELG_dpt(self,y,num,iso ='p',n_f=3,is_pp = False,switch ='dp_t',mu_factor = 1,mu_f_factor=1):
		'''Return the components of the total cross section transeverse momentum 
		dependant that contribute to FCEL and FCEG effects and their 
		integration uncertainites, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as :
		[FCEL,FCEG] with FCEL = Gq+ qqbar+qbarq and FCEG = qG'''
		FCEG = self.dsigma_qG_dpt(y,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		Gq = self.dsigma_Gq_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qqbar = self.dsigma_qqbar_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qbarq = self.dsigma_qbarq_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		FCEL = (Gq[0]+qqbar[0]+qbarq[0],Gq[1]+qqbar[1]+qbarq[1])
		return[FCEL,FCEG]
		
	# Integration in one point (y,p_T) of the phase space 
	
	def sigma_tot(self,y,x_T,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp= False,switch='dp_t'):
		'''Return the cross section for a single point of the phase space (y,p_T)
		and its integration uncertainties, with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		I_tot = lambda xi: self.all_process_integrand(y,x_T,xi,num,mu_factor,mu_f_factor,iso,n_f, is_pp ,switch)
		return integrate.quad(I_tot,Xi_min(y,x_T),Xi_max(y,x_T))
	
	# Uncertainties part (for Bayesian sets)
	
	def Uncertaintites_dy(self,x_T,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp=False,switch='dp_t',var_err= 'mu,pdf'):
		'''Return the total cross section rapidity dependant with various uncertainties.
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		- varr_err a string that contains "mu","pdf" or both to compute uncertainties
		on mu and mu_f (as 1/2 and 2 times mu_factors) and/or the pdf set
		Returned as:
		[Ucen,(Uplus,Uminus),(Uplus_mu,Uminus_mu),(Uplus_pdf,Uminus_pdf)]
		Uplus and Uminus are the quadratic sum of mu and pdf uncertainties'''
		p_set = self.p_set
		Ucen = []
		Uplus_pdf = []
		Uminus_pdf = []
		Uplus_mu = []
		Uminus_mu= []
		num_cen= 0
		Y = Y_list(x_T)
		factor = [1./2.,2]
		mu_factor_list = [(i,j) for i in factor for j in factor]
		mu_val =[]
		print('central value')
		sigma_cen =  self.dsigma_tot_dy(x_T,num_cen,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp=is_pp ,switch=switch)[0]
		if 'mu' in var_err:
			for (a,b) in mu_factor_list:
				print('(mu_factor,mu_f_factor) = '+str((a,b)))
				sigma_mu , err_mu = self.dsigma_tot_dy(x_T,num_cen,mu_factor=a*mu_factor,mu_f_factor=b*mu_f_factor,iso=iso,n_f=n_f,is_pp=is_pp ,switch=switch)
				mu_val.append(sigma_mu)
			Uminus_mu, Uplus_mu = min_max(mu_val)
			Uminus_mu , Uplus_mu = sigma_cen-Uminus_mu, Uplus_mu-sigma_cen
			Ucen = sigma_cen
		if 'pdf' in var_err:
			Ucen_pdf = []
			for j,y in enumerate(Y):
				print('y = '+str(y) +', '+str(j+1)+'/'+str(N_y))
				dsigma = []
				for i in range(p_set.size):
					dsigma.append(self.dsigma_tot_dydpt(y,x_T,i,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp=is_pp ,switch=switch)[0])
				unc = p_set.uncertainty(dsigma)
				Ucen_pdf.append(unc.central)
				Uplus_pdf.append(unc.errplus)
				Uminus_pdf.append(unc.errminus)
			Ucen_pdf, Uplus_pdf, Uminus_pdf = np.array(Ucen_pdf),np.array(Uplus_pdf),np.array(Uminus_pdf)
		if not('mu' in var_err) and ('pdf'in var_err):
			Ucen = Ucen_pdf
		elif not('mu' in var_err) and not('pdf'in var_err):
			raise ValueError("no error variable adequate to this function")
		Uplus ,Uminus= (Uplus_pdf**2+Uplus_mu**2)**0.5,(Uminus_pdf**2+Uminus_mu**2)**0.5
		return[Ucen,(Uplus,Uminus),(Uplus_mu,Uminus_mu),(Uplus_pdf,Uminus_pdf)]
		
	def Uncertaintites_dpt(self,y,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp=False,switch='dp_t',var_err= 'mu,pdf'):
		'''Return the total cross section transverse momentum dependant with various uncertainties.
		- y the rapidity
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		- varr_err a string that contains "mu","pdf" or both to compute uncertainties
		on mu and mu_f (as 1/2 and 2 times mu_factors) and/or the pdf set
		Returned as:
		[Ucen,(Uplus,Uminus),(Uplus_mu,Uminus_mu),(Uplus_pdf,Uminus_pdf)]
		Uplus and Uminus are the quadratic sum of mu and pdf uncertainties'''
		p_set = self.p_set
		Ucen = []
		Uplus_pdf = []
		Uminus_pdf = []
		Uplus_mu = []
		Uminus_mu= []
		num_cen= 0
		P_T = self.P_T_list(y)
		rs = self.rs
		factor = [1./2.,2]
		mu_factor_list = [(i,j) for i in factor for j in factor]
		mu_val =[]
		print('central value')
		sigma_cen = self.dsigma_tot_dpt(y,num_cen,iso=iso,n_f=n_f,is_pp=is_pp ,switch=switch,mu_factor=mu_factor,mu_f_factor=mu_f_factor)[0]
		if 'mu' in var_err:
			for (a,b) in mu_factor_list:
				print('(mu_factor,mu_f_factor) = '+str((a,b)))
				sigma_mu , err_mu = self.dsigma_tot_dpt(y,num_cen,iso=iso,n_f=n_f,is_pp=is_pp ,switch=switch,mu_factor=a*mu_factor,mu_f_factor=b*mu_f_factor)
				mu_val.append(sigma_mu)
			Uminus_mu, Uplus_mu = min_max(mu_val)
			Uminus_mu , Uplus_mu = sigma_cen-Uminus_mu, Uplus_mu-sigma_cen
			Ucen = sigma_cen
		if 'pdf' in var_err:
			Ucen_pdf = []
			for j,p_t in enumerate(P_T):
				print('p_T = '+str(p_t) +', '+str(j+1)+'/'+str(N_pt))
				dsigma = []
				x_T = 2.*p_t/rs
				for i in range(p_set.size):
					dsigma.append(self.dsigma_tot_dydpt(y,x_T,i,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp=is_pp ,switch=switch)[0])
				unc = p_set.uncertainty(dsigma)
				Ucen_pdf.append(unc.central)
				Uplus_pdf.append(unc.errplus)
				Uminus_pdf.append(unc.errminus)
			Ucen_pdf, Uplus_pdf, Uminus_pdf = np.array(Ucen_pdf),np.array(Uplus_pdf),np.array(Uminus_pdf)
		if not('mu' in var_err) and ('pdf'in var_err):
			Ucen = Ucen_pdf
		elif not('mu' in var_err) and not('pdf'in var_err):
			raise ValueError("no error variable adequate to this function")
		Uplus ,Uminus= (Uplus_pdf**2+Uplus_mu**2)**0.5,(Uminus_pdf**2+Uminus_mu**2)**0.5
		return[Ucen,(Uplus,Uminus),(Uplus_mu,Uminus_mu),(Uplus_pdf,Uminus_pdf)]
	
	### Ratios part ###
	
	def R_pA_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,n_f=3,switch='dp_t'):
		'''Return the ratio of the pA (no FCEL/G) cross section over the pp one
		and its integration uncertainties, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		sigma_pA, err_pA = self.dsigma_tot_dy(x_T,num,mu_factor,mu_f_factor,n_f=n_f,switch=switch )
		sigma_pp, err_pp = self.dsigma_tot_dy(x_T,num,mu_factor,mu_f_factor,n_f=n_f,is_pp = True,switch=switch)
		R_pA =  sigma_pA/sigma_pp
		err_R_pA = R_pA*((err_pA/sigma_pA)**2+(err_pp/sigma_pp)**2)**0.5
		return (R_pA,err_R_pA)
		
	def R_pA_dpt(self,y,num,mu_factor = 1,mu_f_factor = 1,n_f=3,switch = 'dp_t'):
		'''Return the ratio of the pA (no FCEL/G) cross section over the pp one
		and its integration uncertainties, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		sigma_pA, err_pA = self.dsigma_tot_dpt(y,num,n_f=n_f ,switch=switch,mu_factor=mu_factor,mu_f_factor = mu_f_factor)
		sigma_pp, err_pp = self.dsigma_tot_dpt(y,num,n_f=n_f ,is_pp = True,switch =switch,mu_factor=mu_factor,mu_f_factor = mu_f_factor)
		R_pA = sigma_pA/sigma_pp
		err_R_pA = R_pA*((err_pA/sigma_pA)**2+(err_pp/sigma_pp)**2)**0.5
		return (R_pA,err_R_pA)
	
	def R_pA_dydpt(self,y,x_T,num,mu_factor = 1,mu_f_factor = 1,n_f=3,switch = 'dp_t'):
		'''Return the ratio of the pA (no FCEL/G) cross section over the pp one
		and its integration uncertainties, with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		sigma_pA, err_pA = self.dsigma_tot_dydpt(y,x_T,num,n_f=n_f ,switch=switch,mu_factor=mu_factor,mu_f_factor = mu_f_factor)
		sigma_pp, err_pp = self.dsigma_tot_dydpt(y,x_T,num,n_f=n_f ,is_pp = True,switch =switch,mu_factor=mu_factor,mu_f_factor = mu_f_factor)
		R_pA = sigma_pA/sigma_pp
		err_R_pA = R_pA*((err_pA/sigma_pA)**2+(err_pp/sigma_pp)**2)**0.5
		return (R_pA,err_R_pA)
		
	def R_iso_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor=1,n_f=3,switch='dp_t'):
		'''Return the ratio of the reconstruction of the nucleide target as 
		Z*sigma_pp + (A-Z)*sigma_pn over A*sigma_pp and its integration
		uncertainties, with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		sigma_pp , err_pp = self.dsigma_tot_dydpt(y,x_T,num,mu_factor,mu_f_factor,'p',n_f=n_f,switch=switch ,is_pp = True)
		sigma_pn , err_pn= self.dsigma_tot_dydpt(y,x_T,num,mu_factor,mu_f_factor,'n',n_f=n_f, switch=switch  ,is_pp = True)
		Z = float(self.Z)
		A = float(self.A)
		R_iso = (sigma_pn/sigma_pp)*(1-Z/A)
		err_R_iso = R_iso*((err_pp/sigma_pp)**2+(err_pn/sigma_pn)**2)**0.5
		return(R_iso+(Z/A),err_R_iso)
	
	def R_iso_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,n_f=3,switch='dp_t'):
		'''Return the ratio of the reconstruction of the nucleide target as 
		Z*sigma_pp + (A-Z)*sigma_pn over A*sigma_pp and its integration
		uncertainties, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)''' 
		sigma_pp, err_pp = self.dsigma_tot_dy(x_T,num,mu_factor,mu_f_factor,'p',n_f=n_f,switch=switch ,is_pp = True)
		sigma_pn, err_pn = self.dsigma_tot_dy(x_T,num,mu_factor,mu_f_factor,'n',n_f=n_f, switch=switch  ,is_pp = True)
		Z = float(self.Z)
		A = float(self.A)
		R_iso = (sigma_pn/sigma_pp)*(1-Z/A)
		err_R_iso = R_iso*((err_pp/sigma_pp)**2+(err_pn/sigma_pn)**2)**0.5
		return(R_iso+(Z/A),err_R_iso)
		
	def R_iso_dpt(self,y,num,mu_factor = 1,mu_f_factor = 1,n_f=3,switch='dp_t'):
		'''Return the ratio of the reconstruction of the nucleide target as 
		Z*sigma_pp + (A-Z)*sigma_pn over A*sigma_pp and its integration
		uncertainties, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		sigma_pp, err_pp = self.dsigma_tot_dpt(y,num,iso='p',n_f=n_f,switch = switch ,is_pp = True,mu_factor=mu_factor,mu_f_factor= mu_f_factor)
		sigma_pn, err_pn = self.dsigma_tot_dpt(y,num,iso='n',n_f=n_f,switch=switch ,is_pp = True,mu_factor=mu_factor,mu_f_factor= mu_f_factor)
		Z = float(self.Z)
		A = float(self.A)
		R_iso = sigma_pn/sigma_pp*(1-Z/A)
		err_R_iso = R_iso*((err_pp/sigma_pp)**2+(err_pn/sigma_pn)**2)**0.5
		return(R_iso+(Z/A),err_R_iso)
	
	def R_wo_iso_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,n_f = 3,switch = 'dp_t'):
		'''Return the ratio of R_pA over R_iso, rapidity dependent, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		R_pA, err_R_pA = self.R_pA_dy(x_T,num,mu_factor,mu_f_factor,n_f,switch)
		R_iso, err_R_iso = self.R_iso_dy(x_T,num,mu_factor,mu_f_factor,n_f,switch)
		R_wo_iso = R_pA/R_iso
		err_R_wo_iso = R_wo_iso*((err_R_pA/R_pA)**2+(err_R_iso/R_iso)**2)**0.5
		return (R_wo_iso, err_R_wo_iso)
		
	def R_wo_iso_dpt(self,y,num,mu_factor = 1,mu_f_factor = 1,n_f=3,switch='dp_t'):
		'''Return the ratio of R_pA over R_iso, transverse momentum dependent,
		with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		R_pA, err_R_pA = self.R_pA_dpt(y,num,mu_factor,mu_f_factor,n_f,switch)
		R_iso, err_R_iso = self.R_iso_dpt(y,num,mu_factor,mu_f_factor,n_f,switch)
		R_wo_iso = R_pA/R_iso
		err_R_wo_iso = R_wo_iso*((err_R_pA/R_pA)**2+(err_R_iso/R_iso)**2)**0.5
		return (R_wo_iso, err_R_wo_iso)
	
	def R_wo_iso_dydpt(self,y,x_T,num,mu_factor = 1,mu_f_factor = 1,n_f=3,switch='dp_t'):
		'''Return the ratio of R_pA over R_iso, transverse momentum dependent,
		with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)'''
		R_pA, err_R_pA = self.R_pA_dydpt(y,x_T,num,mu_factor,mu_f_factor,n_f,switch)
		R_iso, err_R_iso = self.R_iso_dydpt(y,x_T,num,mu_factor,mu_f_factor,n_f,switch)
		R_wo_iso = R_pA/R_iso
		err_R_wo_iso = R_wo_iso*((err_R_pA/R_pA)**2+(err_R_iso/R_iso)**2)**0.5
		return (R_wo_iso, err_R_wo_iso)
		
	def R3_dy(self,x_T,num,mu_factor = 1,mu_f_factor = 1,n_f = 3,switch = 'dp_t'):
		'''Return simultaneously R_pA_dy, R_iso_dy and R_wo_iso_dy. Gains in 
		computational time, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : [(R_pA,err_R_pA),(R_iso,err_R_iso),(R_wo_iso,err_R_wo_iso)]'''
		sigma_pA, err_pA = self.dsigma_tot_dy(x_T,num,mu_factor,mu_f_factor,n_f=n_f,switch=switch )
		sigma_pp, err_pp = self.dsigma_tot_dy(x_T,num,mu_factor,mu_f_factor,n_f=n_f,is_pp = True,switch=switch)
		sigma_pn, err_pn = self.dsigma_tot_dy(x_T,num,mu_factor,mu_f_factor,iso='n',n_f=n_f, switch=switch  ,is_pp = True)
		Z = float(self.Z)
		A = float(self.A)
		R_pA = sigma_pA/sigma_pp
		err_R_pA = R_pA*((err_pA/sigma_pA)**2+(err_pp/sigma_pp)**2)**0.5
		R_iso = sigma_pn/sigma_pp*(1-Z/A)
		err_R_iso = R_iso*((err_pp/sigma_pp)**2+(err_pn/sigma_pn)**2)**0.5
		R_wo_iso = R_pA/R_iso
		err_R_wo_iso = R_wo_iso*((err_R_pA/R_pA)**2+(err_R_iso/R_iso)**2)**0.5
		return[(R_pA,err_R_pA),(R_iso,err_R_iso),(R_wo_iso,err_R_wo_iso)]
		
	def R3_dpt(self,y,num,mu_factor=1,mu_f_factor=1,n_f=3,switch='dp_t'):
		'''Return simultaneously R_pA_dpt, R_iso_dpt and R_wo_iso_dpt. Gains in
		computational time, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : [(R_pA,err_R_pA),(R_iso,err_R_iso),(R_wo_iso,err_R_wo_iso)]'''
		sigma_pA, err_pA = self.dsigma_tot_dpt(y,num,n_f=n_f,switch=switch,mu_factor=mu_factor,mu_f_factor=mu_f_factor )
		sigma_pp, err_pp = self.dsigma_tot_dpt(y,num,n_f=n_f,is_pp = True,switch=switch,mu_factor=mu_factor,mu_f_factor=mu_f_factor)
		sigma_pn, err_pn = self.dsigma_tot_dpt(y,num,iso='n',n_f=n_f, switch=switch,is_pp = True,mu_factor=mu_factor,mu_f_factor=mu_f_factor)
		Z = float(self.Z)
		A = float(self.A)
		R_pA = sigma_pA/sigma_pp
		err_R_pA = R_pA*((err_pA/sigma_pA)**2+(err_pp/sigma_pp)**2)**0.5
		R_iso = sigma_pn/sigma_pp*(1-Z/A)
		err_R_iso = R_iso*((err_pp/sigma_pp)**2+(err_pn/sigma_pn)**2)**0.5
		R_wo_iso = R_pA/R_iso
		err_R_wo_iso = R_wo_iso*((err_R_pA/R_pA)**2+(err_R_iso/R_iso)**2)**0.5
		return[(R_pA,err_R_pA),(R_iso,err_R_iso),(R_wo_iso,err_R_wo_iso)]
		
		
	def R_composition_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp=False,switch='dp_t'):
		'''Return the list of all 4 ratios for the 4 composants of the cross 
		section over the total cross section and their integration uncertainties
		with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : [(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]'''
		qG = self.dsigma_qG_dy(x_T,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		Gq = self.dsigma_Gq_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qqbar = self.dsigma_qqbar_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qbarq = self.dsigma_qbarq_dy(x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		tot = (qG[0]+Gq[0]+qqbar[0]+qbarq[0],qG[1]+Gq[1]+qqbar[1]+qbarq[1])
		R_qG , R_Gq , R_qqbar, R_qbarq = qG[0]/tot[0],Gq[0]/tot[0],qqbar[0]/tot[0],qbarq[0]/tot[0]
		err_R_qG = R_qG*((qG[1]/qG[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_Gq = R_Gq*((Gq[1]/Gq[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qqbar = R_qqbar*((qqbar[1]/qqbar[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qbarq = R_qbarq*((qbarq[1]/qbarq[0])**2+(tot[1]/tot[0])**2)**0.5
		return[(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]
		
	def R_composition_dpt(self,y,num,mu_factor=1,mu_f_factor =1,iso = 'p',n_f=3,is_pp=False, switch='dp_t'):
		'''Return the list of all 4 ratios for the 4 composants of the cross 
		section over the total cross section and their integration uncertainties
		with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : [(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]'''
		qG = self.dsigma_qG_dpt(y,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		Gq = self.dsigma_Gq_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qqbar = self.dsigma_qqbar_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qbarq = self.dsigma_qbarq_dpt(y,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		tot = (qG[0]+Gq[0]+qqbar[0]+qbarq[0],qG[1]+Gq[1]+qqbar[1]+qbarq[1])
		R_qG , R_Gq , R_qqbar, R_qbarq = qG[0]/tot[0],Gq[0]/tot[0],qqbar[0]/tot[0],qbarq[0]/tot[0]
		err_R_qG = R_qG*((qG[1]/qG[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_Gq = R_Gq*((Gq[1]/Gq[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qqbar = R_qqbar*((qqbar[1]/qqbar[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qbarq = R_qbarq*((qbarq[1]/qbarq[0])**2+(tot[1]/tot[0])**2)**0.5
		return[(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]
	
	def R_composition_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor =1,iso = 'p',n_f=3,is_pp=False, switch='dp_t'):
		'''Return the list of all 4 ratios for the 4 composants of the cross 
		section over the total cross section and their integration uncertainties
		with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : [(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]'''
		qG = self.dsigma_qG_dydpt(y,x_T,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		Gq = self.dsigma_Gq_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qqbar = self.dsigma_qqbar_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qbarq = self.dsigma_qbarq_dydpt(y,x_T,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		tot = (qG[0]+Gq[0]+qqbar[0]+qbarq[0],qG[1]+Gq[1]+qqbar[1]+qbarq[1])
		R_qG , R_Gq , R_qqbar, R_qbarq = qG[0]/tot[0],Gq[0]/tot[0],qqbar[0]/tot[0],qbarq[0]/tot[0]
		err_R_qG = R_qG*((qG[1]/qG[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_Gq = R_Gq*((Gq[1]/Gq[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qqbar = R_qqbar*((qqbar[1]/qqbar[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qbarq = R_qbarq*((qbarq[1]/qbarq[0])**2+(tot[1]/tot[0])**2)**0.5
		return[(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]
	
	def R_composition_dydpt_M(self,y,x_T,M,num,mu_factor=1,mu_f_factor =1,iso = 'p',n_f=3,is_pp=False, switch='dp_t'):
		'''Return the list of all 4 ratios for the 4 composants of the cross 
		section over the total cross section and their integration uncertainties
		with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- M the mass of the virtual photon
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : [(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]'''
		qG = self.dsigma_qG_dydpt_M(y,x_T,M,num,mu_factor,mu_f_factor,n_f,is_pp,switch)
		Gq = self.dsigma_Gq_dydpt_M(y,x_T,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qqbar = self.dsigma_qqbar_dydpt_M(y,x_T,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		qbarq = self.dsigma_qbarq_dydpt_M(y,x_T,M,num,mu_factor,mu_f_factor,iso,n_f,is_pp,switch)
		tot = (qG[0]+Gq[0]+qqbar[0]+qbarq[0],qG[1]+Gq[1]+qqbar[1]+qbarq[1])
		R_qG , R_Gq , R_qqbar, R_qbarq = qG[0]/tot[0],Gq[0]/tot[0],qqbar[0]/tot[0],qbarq[0]/tot[0]
		err_R_qG = R_qG*((qG[1]/qG[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_Gq = R_Gq*((Gq[1]/Gq[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qqbar = R_qqbar*((qqbar[1]/qqbar[0])**2+(tot[1]/tot[0])**2)**0.5
		err_R_qbarq = R_qbarq*((qbarq[1]/qbarq[0])**2+(tot[1]/tot[0])**2)**0.5
		return[(R_qG,err_R_qG),(R_Gq,err_R_Gq),(R_qqbar,err_R_qqbar),(R_qbarq,err_R_qbarq)]
		
	# With Uncertainties
	
	def Uncertainties_R_composition_dy(self,x_T,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp=False,switch='dp_t'):
		'''Return the pdf set uncertainites on the ratio of componentns over 
		the total, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : 
			[(Ucen_qG, Uplus_qG, Uminus_qG, Usym_qG),
			(Ucen_Gq, Uplus_Gq, Uminus_Gq, Usym_Gq),
			(Ucen_qqbar, Uplus_qqbar, Uminus_qqbar, Usym_qqbar),
			(Ucen_qbarq, Uplus_qbarq, Uminus_qbarq, Usym_qbarq)]'''
		p = self.p_set
		Ucen_qG, Uplus_qG, Uminus_qG, Usym_qG = [],[],[],[]
		Ucen_Gq, Uplus_Gq, Uminus_Gq, Usym_Gq = [],[],[],[]
		Ucen_qqbar, Uplus_qqbar, Uminus_qqbar, Usym_qqbar = [],[],[],[]
		Ucen_qbarq, Uplus_qbarq, Uminus_qbarq, Usym_qbarq = [],[],[],[]
		Y = Y_list(x_T)
		for j in range(len(Y)):
			R_qG, R_Gq, R_qqbar, R_qbarq= [],[],[],[]
			for i in range(p.size):
				print('y_idx = '+str(j+1)+' over ' + str(N_y)+ ' and mem_idx = '+ str(i+1) +' over '+ str(p.size)) 
				R = self.R_composition_dydpt(Y[j],x_T,i,mu_factor,mu_f_factor,iso='p',n_f=3,is_pp=False,switch='dp_t')
				R_qG.append(R[0][0])
				R_Gq.append(R[1][0])
				R_qqbar.append(R[2][0])
				R_qbarq.append(R[3][0])
			unc_qG, unc_Gq, unc_qqbar, unc_qbarq = p.uncertainty(R_qG), p.uncertainty(R_Gq), p.uncertainty(R_qqbar), p.uncertainty(R_qbarq)
			Ucen_qG.append(unc_qG.central); Uplus_qG.append(unc_qG.errplus); Uminus_qG.append(unc_qG.errminus); Usym_qG.append(unc_qG.errsymm)
			Ucen_Gq.append(unc_Gq.central); Uplus_Gq.append(unc_Gq.errplus); Uminus_Gq.append(unc_Gq.errminus); Usym_Gq.append(unc_Gq.errsymm)
			Ucen_qqbar.append(unc_qqbar.central); Uplus_qqbar.append(unc_qqbar.errplus); Uminus_qqbar.append(unc_qqbar.errminus); Usym_qqbar.append(unc_qqbar.errsymm)
			Ucen_qbarq.append(unc_qbarq.central); Uplus_qbarq.append(unc_qbarq.errplus); Uminus_qbarq.append(unc_qbarq.errminus); Usym_qbarq.append(unc_qbarq.errsymm)
		return[(Ucen_qG, Uplus_qG, Uminus_qG, Usym_qG),(Ucen_Gq, Uplus_Gq, Uminus_Gq, Usym_Gq),(Ucen_qqbar, Uplus_qqbar, Uminus_qqbar, Usym_qqbar),(Ucen_qbarq, Uplus_qbarq, Uminus_qbarq, Usym_qbarq)]
	
	def Uncertainties_R_composition_dpt(self,y,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp=False,switch='dp_t'):
		'''Return the pdf set uncertainites on the ratio of componentns over 
		the total, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : 
			[(Ucen_qG, Uplus_qG, Uminus_qG, Usym_qG),
			(Ucen_Gq, Uplus_Gq, Uminus_Gq, Usym_Gq),
			(Ucen_qqbar, Uplus_qqbar, Uminus_qqbar, Usym_qqbar),
			(Ucen_qbarq, Uplus_qbarq, Uminus_qbarq, Usym_qbarq)]'''
		p = self.p_set
		Ucen_qG, Uplus_qG, Uminus_qG, Usym_qG = [],[],[],[]
		Ucen_Gq, Uplus_Gq, Uminus_Gq, Usym_Gq = [],[],[],[]
		Ucen_qqbar, Uplus_qqbar, Uminus_qqbar, Usym_qqbar = [],[],[],[]
		Ucen_qbarq, Uplus_qbarq, Uminus_qbarq, Usym_qbarq = [],[],[],[]
		P_T = self.P_T_list(y)
		rs = self.rs
		for j in range(len(P_T)):
			R_qG, R_Gq, R_qqbar, R_qbarq= [],[],[],[]
			for i in range(p.size):
				x_T = 2*P_T[j]/rs
				print('pt_idx = '+str(j+1)+' over ' + str(N_pt)+ ' and mem_idx = '+ str(i+1) +' over '+ str(p.size)) 
				R = self.R_composition_dydpt(y,x_T,i,mu_factor,mu_f_factor,iso='p',n_f=3,is_pp=False,switch='dp_t')
				R_qG.append(R[0][0])
				R_Gq.append(R[1][0])
				R_qqbar.append(R[2][0])
				R_qbarq.append(R[3][0])
			unc_qG, unc_Gq, unc_qqbar, unc_qbarq = p.uncertainty(R_qG), p.uncertainty(R_Gq), p.uncertainty(R_qqbar), p.uncertainty(R_qbarq)
			Ucen_qG.append(unc_qG.central); Uplus_qG.append(unc_qG.errplus); Uminus_qG.append(unc_qG.errminus); Usym_qG.append(unc_qG.errsymm)
			Ucen_Gq.append(unc_Gq.central); Uplus_Gq.append(unc_Gq.errplus); Uminus_Gq.append(unc_Gq.errminus); Usym_Gq.append(unc_Gq.errsymm)
			Ucen_qqbar.append(unc_qqbar.central); Uplus_qqbar.append(unc_qqbar.errplus); Uminus_qqbar.append(unc_qqbar.errminus); Usym_qqbar.append(unc_qqbar.errsymm)
			Ucen_qbarq.append(unc_qbarq.central); Uplus_qbarq.append(unc_qbarq.errplus); Uminus_qbarq.append(unc_qbarq.errminus); Usym_qbarq.append(unc_qbarq.errsymm)
		return[(Ucen_qG, Uplus_qG, Uminus_qG, Usym_qG),(Ucen_Gq, Uplus_Gq, Uminus_Gq, Usym_Gq),(Ucen_qqbar, Uplus_qqbar, Uminus_qqbar, Usym_qqbar),(Ucen_qbarq, Uplus_qbarq, Uminus_qbarq, Usym_qbarq)]
	
	def Uncertainties_R_FCEl_FCEG_dy(self,x_T,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp=False,switch ='dp_t'):
		'''Return the pdf set uncertainites on the ratio of FCEL/G components
		over the total, with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : 
			[(Ucen_FCEL, Uplus_FCEL, Uminus_FCEL, Usym_FCEL),
			(Ucen_FCEG, Uplus_FCEG, Uminus_FCEG, Usym_FCEG)]'''
		p = self.p_set
		Ucen_FCEL, Uplus_FCEL, Uminus_FCEL, Usym_FCEL = [],[],[],[]
		Ucen_FCEG, Uplus_FCEG, Uminus_FCEG, Usym_FCEG = [],[],[],[]
		Y = Y_list(x_T)
		for j in range(len(Y)):
			R_qG, R_Gq, R_qqbar, R_qbarq= [],[],[],[]
			for i in range(p.size):
				print('y_idx = '+str(j+1)+' over ' + str(N_y)+ ' and mem_idx = '+ str(i+1) +' over '+ str(p.size)) 
				R = self.R_composition_dydpt(Y[j],x_T,i,mu_factor,mu_f_factor,iso='p',n_f=3,is_pp=False,switch='dp_t')
				R_qG.append(R[0][0])
				R_Gq.append(R[1][0])
				R_qqbar.append(R[2][0])
				R_qbarq.append(R[3][0])
			FCEL, FCEG = np.add(np.add(R_Gq,R_qqbar),R_qbarq), np.array(R_qG)
			unc_FCEL, unc_FCEG = p.uncertainty(FCEL), p.uncertainty(FCEG)
			Ucen_FCEL.append(unc_FCEL.central); Uplus_FCEL.append(unc_FCEL.errplus); Uminus_FCEL.append(unc_FCEL.errminus); Usym_FCEL.append(unc_FCEL.errsymm)
			Ucen_FCEG.append(unc_FCEG.central); Uplus_FCEG.append(unc_FCEG.errplus); Uminus_FCEG.append(unc_FCEG.errminus); Usym_FCEG.append(unc_FCEG.errsymm)
		return [(Ucen_FCEL, Uplus_FCEL, Uminus_FCEL, Usym_FCEL),(Ucen_FCEG, Uplus_FCEG, Uminus_FCEG, Usym_FCEG)]
	
	def Uncertainties_R_FCEl_FCEG_dpt(self,y,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,is_pp=False,switch ='dp_t'):
		'''Return the pdf set uncertainites on the ratio of FCEL/G components
		over the total, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		as : 
			[(Ucen_FCEL, Uplus_FCEL, Uminus_FCEL, Usym_FCEL),
			(Ucen_FCEG, Uplus_FCEG, Uminus_FCEG, Usym_FCEG)]'''
		p = self.p_set
		Ucen_FCEL, Uplus_FCEL, Uminus_FCEL, Usym_FCEL = [],[],[],[]
		Ucen_FCEG, Uplus_FCEG, Uminus_FCEG, Usym_FCEG = [],[],[],[]
		P_T = self.P_T_list(y)
		rs = self.rs
		for j in range(len(P_T)):
			R_qG, R_Gq, R_qqbar, R_qbarq= [],[],[],[]
			for i in range(p.size):
				x_T = 2*P_T[j]/rs
				print('pt_idx = '+str(j+1)+' over ' + str(N_pt)+ ' and mem_idx = '+ str(i+1) +' over '+ str(p.size)) 
				R = self.R_composition_dydpt(y,x_T,i,mu_factor,mu_f_factor,iso='p',n_f=3,is_pp=False,switch='dp_t')
				R_qG.append(R[0][0])
				R_Gq.append(R[1][0])
				R_qqbar.append(R[2][0])
				R_qbarq.append(R[3][0])
			FCEL, FCEG = np.add(np.add(R_Gq,R_qqbar),R_qbarq), np.array(R_qG)
			unc_FCEL, unc_FCEG = p.uncertainty(FCEL), p.uncertainty(FCEG)
			Ucen_FCEL.append(unc_FCEL.central); Uplus_FCEL.append(unc_FCEL.errplus); Uminus_FCEL.append(unc_FCEL.errminus); Usym_FCEL.append(unc_FCEL.errsymm)
			Ucen_FCEG.append(unc_FCEG.central); Uplus_FCEG.append(unc_FCEG.errplus); Uminus_FCEG.append(unc_FCEG.errminus); Usym_FCEG.append(unc_FCEG.errsymm)
		return [(Ucen_FCEL, Uplus_FCEL, Uminus_FCEL, Usym_FCEL),(Ucen_FCEG, Uplus_FCEG, Uminus_FCEG, Usym_FCEG)]
		
	### FCEL/FCEG implementation ###
	
	def FCEL_integrand(self,y,x_T,num,mu_factor=1,mu_f_factor=1,iso ='p',n_f=3,switch ='dp_t',var_int='nu',q0=0.07):
		'''return the FCEL integrand as a function (lambda class) of nu and xi,
		with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		- q0 the transport coefficient
		as :
			lambda nu,xi,nu_min:jacobian(nu,xi)*P(nu,xi,nu_min)*dsigma(nu,xi)
			- nu the default integration variable, which can be set in var_int:
				* "nu" = ln(u)
				* "delta" = ln(1+x)
				* "nu2" = -ln(u)
				(for more info see the reedme)
			- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
			- nu_min refers to the minimum of the integration bounds for a good
			normalisation of P (same, see reedme for more info)'''
		alpha_s = 0.5
		rs = self.rs
		p_t = rs*x_T/2.
		Fc = N_c
		A = self.A
		m = 0 																	# masse of the photon
		proba_FCEL = p.proba(A,1.,rs,p_t,y,alpha_s,Fc,m,q0)
		sigma_hat = lambda xi: proba_FCEL.sigma_hat(xi)
		chi = lambda xi:proba_FCEL.Chi(xi)
		g_FCEL = lambda u,xi: p.g_u(u,chi(xi),Fc,alpha_s)						#functions for a better normalisation of p_tilde_u
		if var_int == 'nu':														# u = exp(nu)
			P = lambda nu,xi,nu_min: p.p_tilde_u(np.exp(nu),chi(xi),Fc,alpha=alpha_s)/(1-np.exp(-1*g_FCEL(np.exp(nu_min),xi)))
			jacobian = lambda nu,xi: np.exp(nu)/(sigma_hat(xi)*np.exp(nu)+1)
			delta = lambda nu,xi: np.log(1+sigma_hat(xi)*np.exp(nu))
			dsigma =lambda nu,xi: self.Gq(y+delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso =iso,n_f=n_f,is_pp = True,switch =switch) + self.qqbar(y+delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp = True,switch = switch)+self.qbarq(y+delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp = True,switch = switch)
			return lambda nu,xi,nu_min:jacobian(nu,xi)*P(nu,xi,nu_min)*dsigma(nu,xi)
		elif var_int == 'delta':												# delta = ln(1+x)
			P = lambda delta,xi,delta_min: p.p_tilde_u((np.exp(delta)-1)/sigma_hat(xi),chi(xi),Fc,alpha=alpha_s)/(1-np.exp(-1*g_FCEL((np.exp(delta_min)-1)/sigma_hat(xi),xi)))
			jacobian = lambda delta,xi: 1.
			dsigma =lambda delta,xi: self.Gq(y+delta,x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso =iso,n_f=n_f,is_pp = True,switch =switch) + self.qqbar(y+delta,x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp = True,switch = switch)+self.qbarq(y+delta,x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp = True,switch = switch)
			return lambda delta,xi,delta_min:jacobian(delta,xi)*P(delta,xi,delta_min)*dsigma(delta,xi)
		elif var_int == 'nu2':													# u = exp(-nu)
			P = lambda nu,xi,nu_min: p.p_tilde_u(np.exp(-1*nu),chi(xi),Fc,alpha=alpha_s)/(1-np.exp(-1*g_FCEL(np.exp(-1*nu_min),xi)))
			jacobian = lambda nu,xi: np.exp(-nu)/(sigma_hat(xi)*np.exp(-1*nu))
			delta = lambda nu,xi: np.log(1+sigma_hat(xi)*np.exp(-1*nu))
			dsigma =lambda nu,xi: self.Gq(y+delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso =iso,n_f=n_f,is_pp = True,switch =switch) + self.qqbar(y+delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp = True,switch = switch)+self.qbarq(y+delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp = True,switch = switch)
			return lambda nu,xi,nu_min:jacobian(nu,xi)*P(nu,xi,nu_min)*dsigma(nu,xi)
		
	def FCEG_integrand(self,y,x_T,num,mu_factor=1,mu_f_factor=1,n_f=3,switch ='dp_t',var_int = 'nu',q0=0.07):
		'''return the FCEG integrand as a function (lambda class) of nu and xi,
		with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		- q0 the transport coefficient
		as :
			lambda nu,xi,nu_min:jacobian(nu,xi)*P(nu,xi,nu_min)*dsigma(nu,xi)
			- nu the default integration variable, which can be set in var_int:
				* "nu" = ln(u)
				* "delta" = ln(1+x)
				* "nu2" = -ln(u)
				(for more info see the reedme)
			- xi the partonic fraction p+_3/p+_1 = -\hat{u}/\hat{s}
			- nu_min refers to the minimum of the integration bounds for a good
			normalisation of P (same, see reedme for more info)'''
		alpha_s = 0.5
		rs = self.rs
		p_t = rs*x_T/2.
		Fc = -1./N_c
		A = self.A
		m=0
		proba_FCEG = p.proba(A,1.,rs,p_t,y,alpha_s,Fc,m,q0)
		sigma_hat = lambda xi: proba_FCEG.sigma_hat(xi)
		chi = lambda xi:proba_FCEG.Chi(xi)
		g_FCEG = lambda u,xi: p.g_u(u,chi(xi),Fc,alpha_s)						# functions for a better normalisation of p_tilde_u
		if var_int == 'nu':
			P = lambda nu,xi,nu_min: p.p_tilde_u(np.exp(nu),chi(xi),Fc,alpha=alpha_s)/(1-np.exp(-1*g_FCEG(np.exp(nu_min),xi)))
			jacobian = lambda nu,xi: np.exp(nu)*(1+np.exp(nu)*sigma_hat(xi))
			delta = lambda nu,xi: np.log(1+sigma_hat(xi)*np.exp(nu))
			dsigma = lambda nu,xi: self.qG(y-delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,is_pp = True,switch = switch)
			return lambda nu,xi,nu_min:jacobian(nu,xi)*P(nu,xi,nu_min)*dsigma(nu,xi)
		elif var_int == 'delta':
			P = lambda delta,xi,delta_min: p.p_tilde_u((np.exp(delta)-1)/sigma_hat(xi),chi(xi),Fc,alpha=alpha_s)/(1-np.exp(-1*g_FCEG((np.exp(delta_min)-1)/sigma_hat(xi),xi)))
			jacobian = lambda delta: np.exp(2*delta)
			dsigma = lambda delta,xi: self.qG(y-delta,x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,is_pp = True,switch = switch)
			return lambda delta,xi,delta_min:jacobian(delta)*P(delta,xi,delta_min)*dsigma(delta,xi)
		elif var_int == 'nu2':
			P = lambda nu,xi,nu_min: p.p_tilde_u(np.exp(-1*nu),chi(xi),Fc,alpha=alpha_s)/(1-np.exp(-1*g_FCEG(np.exp(-1*nu_min),xi)))
			jacobian = lambda nu,xi: np.exp(-1*nu)*(1+np.exp(-1*nu)*sigma_hat(xi))
			delta = lambda nu,xi: np.log(1+sigma_hat(xi)*np.exp(-1*nu))
			dsigma = lambda nu,xi: self.qG(y-delta(nu,xi),x_T,xi,num,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,is_pp = True,switch = switch)
			return lambda nu,xi,nu_min:jacobian(nu,xi)*P(nu,xi,nu_min)*dsigma(nu,xi)
		# replace by [lambda nu,xi:jacobian(nu,xi), lambda nu,xi:P(nu,xi,nu_min),lambda nu,xi:dsigma(nu,xi),chi,sigma_hat] for the file integrand_FCELG.py
		
	def FCEL_G_integration_dy(self,x_T,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0 = 0.07):
		'''Return the integrated FCEL and FCEG componnents rapidity dependent,
		with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0 the transport coefficient
		as:
			[(sigma_FCEL,err_sigma_FCEL),
			(sigma_FCEG,err_sigma_FCEG)]'''
		if eps<=0:
			raise ValueError("eps <= 0, leading to integration error.")
		else:
			Y = Y_list(x_T)
			A = self.A
			rs = self.rs
			p_T = rs*x_T/2.
			alpha_s = 0.5														# don't use the pdf one, bc the energy scale here is sqrt(q_hat*L)
			Fc_FCEG = -1./N_c
			sigma_FCEL = []
			err_sigma_FCEL = []
			sigma_FCEG = []
			err_sigma_FCEG = []
			for i,y in enumerate(Y):
				print('y = '+str(y)+'| '+str(i+1)+'/'+str(N_y))
				prob_FCEG = p.proba(A,1.,rs,p_T,y,alpha_s,Fc_FCEG,0,q0)
				sigma_hat = lambda xi: prob_FCEG.sigma_hat(xi)
				k = 2./x_T
				delta_y_FCEL_max = lambda xi: min(np.log(k*xi)-y,np.log(2))
				delta_y_FCEG_max = lambda xi: min(y+np.log((1-xi)*k),np.log(2))
				delta_y_min = eps
				if var_int == 'nu':
					max_FCEL = lambda xi: np.log((np.exp(delta_y_FCEL_max(xi))-1.)/sigma_hat(xi))
					max_FCEG = lambda xi: np.log((np.exp(delta_y_FCEG_max(xi))-1.)/sigma_hat(xi))
					if delta_y_min > 1e-15:
						min_FCEL = lambda xi: np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
						min_FCEG = lambda xi: np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
					elif delta_y_min > 0:
						min_FCEL = lambda xi: np.log(delta_y_min/sigma_hat(xi))
						min_FCEG = lambda xi: np.log(delta_y_min/sigma_hat(xi))
				elif var_int == 'nu2':
					min_FCEL = lambda xi: -1*np.log((np.exp(delta_y_FCEL_max(xi))-1.)/sigma_hat(xi))
					min_FCEG = lambda xi: -1*np.log((np.exp(delta_y_FCEG_max(xi))-1.)/sigma_hat(xi))
					if delta_y_min > 1e-15:
						max_FCEL = lambda xi: -1*np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
						max_FCEG = lambda xi: -1*np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
					elif delta_y_min > 0:
						max_FCEL = lambda xi: -1*np.log(delta_y_min/sigma_hat(xi))
						max_FCEG = lambda xi: -1*np.log(delta_y_min/sigma_hat(xi))
				elif var_int == 'delta':
					max_FCEL = delta_y_FCEL_max
					max_FCEG = delta_y_FCEG_max
					min_FCEL = lambda xi:delta_y_min
					min_FCEG = lambda xi:delta_y_min
				Integrand_FCEL = self.FCEL_integrand(y,x_T,num,mu_factor,mu_f_factor,iso,n_f,switch,var_int=var_int,q0=q0)
				Int_FCEL = integrate.quad(lambda xi: integrate.quad(lambda nu: Integrand_FCEL(nu,xi,min_FCEL(xi)),min_FCEL(xi),max_FCEL(xi))[0],Xi_min(y,x_T),Xi_max(y,x_T))
				Integrand_FCEG = self.FCEG_integrand(y,x_T,num,mu_factor,mu_f_factor,n_f,switch,var_int=var_int,q0=q0)
				Int_FCEG = integrate.quad(lambda xi: integrate.quad(lambda nu: Integrand_FCEG(nu,xi,min_FCEG(xi)),min_FCEG(xi),max_FCEG(xi))[0],Xi_min(y,x_T),Xi_max(y,x_T))
				sigma_FCEL.append(conv_fact*Int_FCEL[0])
				err_sigma_FCEL.append(conv_fact*Int_FCEL[1])
				sigma_FCEG.append(conv_fact*Int_FCEG[0])
				err_sigma_FCEG.append(conv_fact*Int_FCEG[1])
			return [(np.array(sigma_FCEL),np.array(err_sigma_FCEL)),(np.array(sigma_FCEG),np.array(err_sigma_FCEG))]
	
	def FCEL_G_integration_dpt(self,y,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0 = 0.07):
		'''Return the integrated FCEL and FCEG componnents transverse momentum
		dependent, with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0 the transport coefficient
		as:
			[(sigma_FCEL,err_sigma_FCEL),
			(sigma_FCEG,err_sigma_FCEG)]'''
		if eps<=0:
			raise ValueError("eps <= 0, leading to integration error.")
		else:
			P_T = self.P_T_list(y)
			A = self.A
			rs = self.rs
			alpha_s = 0.5														# don't use the pdf one, bc the energy scale here is sqrt(q_hat*L)
			Fc_FCEG = -1./N_c
			sigma_FCEL = []
			err_sigma_FCEL = []
			sigma_FCEG = []
			err_sigma_FCEG = []
			for i,pt in enumerate(P_T):
				x_T = 2*pt/rs
				print('pt = '+str(pt)+'| '+str(i+1)+'/'+str(N_pt))
				prob_FCEG = p.proba(A,1.,rs,pt,y,alpha_s,Fc_FCEG,0,q0)
				sigma_hat = lambda xi: prob_FCEG.sigma_hat(xi)
				k = 2./x_T
				delta_y_FCEL_max = lambda xi: min(np.log(k*xi)-y,np.log(2))
				delta_y_FCEG_max = lambda xi: min(y+np.log((1-xi)*k),np.log(2))
				delta_y_min = eps
				if var_int == 'nu':
					max_FCEL = lambda xi: np.log((np.exp(delta_y_FCEL_max(xi))-1.)/sigma_hat(xi))
					max_FCEG = lambda xi: np.log((np.exp(delta_y_FCEG_max(xi))-1.)/sigma_hat(xi))
					if delta_y_min > 1e-15:
						min_FCEL = lambda xi: np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
						min_FCEG = lambda xi: np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
					elif delta_y_min > 0:
						min_FCEL = lambda xi: np.log(delta_y_min/sigma_hat(xi))
						min_FCEG = lambda xi: np.log(delta_y_min/sigma_hat(xi))
				elif var_int == 'nu2':
					min_FCEL = lambda xi: -1*np.log((np.exp(delta_y_FCEL_max(xi))-1.)/sigma_hat(xi))
					min_FCEG = lambda xi: -1*np.log((np.exp(delta_y_FCEG_max(xi))-1.)/sigma_hat(xi))
					if delta_y_min > 1e-15:
						max_FCEL = lambda xi: -1*np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
						max_FCEG = lambda xi: -1*np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
					elif delta_y_min > 0:
						max_FCEL = lambda xi: -1*np.log(delta_y_min/sigma_hat(xi))
						max_FCEG = lambda xi: -1*np.log(delta_y_min/sigma_hat(xi))
				elif var_int == 'delta':
					max_FCEL = delta_y_FCEL_max
					max_FCEG = delta_y_FCEG_max
					min_FCEL = lambda xi:delta_y_min
					min_FCEG = lambda xi:delta_y_min
				Integrand_FCEL = self.FCEL_integrand(y,x_T,num,mu_factor,mu_f_factor,iso,n_f,switch,var_int=var_int,q0=q0)
				Int_FCEL = integrate.quad(lambda xi: integrate.quad(lambda nu: Integrand_FCEL(nu,xi,min_FCEL(xi)),min_FCEL(xi),max_FCEL(xi))[0],Xi_min(y,x_T),Xi_max(y,x_T))
				Integrand_FCEG = self.FCEG_integrand(y,x_T,num,mu_factor,mu_f_factor,n_f,switch,var_int=var_int,q0=q0)
				Int_FCEG = integrate.quad(lambda xi: integrate.quad(lambda nu: Integrand_FCEG(nu,xi,min_FCEG(xi)),min_FCEG(xi),max_FCEG(xi))[0],Xi_min(y,x_T),Xi_max(y,x_T))
				sigma_FCEL.append(conv_fact*Int_FCEL[0])
				err_sigma_FCEL.append(conv_fact*Int_FCEL[1])
				sigma_FCEG.append(conv_fact*Int_FCEG[0])
				err_sigma_FCEG.append(conv_fact*Int_FCEG[1])
			return [(np.array(sigma_FCEL),np.array(err_sigma_FCEL)),(np.array(sigma_FCEG),np.array(err_sigma_FCEG))]
	
	def FCEL_G_integration_dydpt(self,y,x_T,num,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0 = 0.07):
		'''Return the integrated FCEL and FCEG componnents for (y,p_T) of the phase space,
		with:
		- y the rapidity
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- is_pp a booleen (=False by default) to tell the collision type
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0 the transport coefficient
		as:
			[(sigma_FCEL,err_sigma_FCEL),
			(sigma_FCEG,err_sigma_FCEG)]'''
		if eps<=0:
			raise ValueError("eps <= 0, leading to integration error.")
		else:
			A = self.A
			rs = self.rs
			p_T = rs*x_T/2.
			alpha_s = 0.5														# don't use the pdf one, bc the energy scale here is sqrt(q_hat*L)
			Fc_FCEG = -1./N_c
			prob_FCEG = p.proba(A,1.,rs,p_T,y,alpha_s,Fc_FCEG,0,q0)
			sigma_hat = lambda xi: prob_FCEG.sigma_hat(xi)
			k = 2./x_T
			delta_y_FCEL_max = lambda xi: min(np.log(k*xi)-y,np.log(2))
			delta_y_FCEG_max = lambda xi: min(y+np.log((1-xi)*k),np.log(2))
			delta_y_min = eps
			if var_int == 'nu':
				max_FCEL = lambda xi: np.log((np.exp(delta_y_FCEL_max(xi))-1.)/sigma_hat(xi))
				max_FCEG = lambda xi: np.log((np.exp(delta_y_FCEG_max(xi))-1.)/sigma_hat(xi))
				if delta_y_min > 1e-15:
					min_FCEL = lambda xi: np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
					min_FCEG = lambda xi: np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
				elif delta_y_min > 0:
					min_FCEL = lambda xi: np.log(delta_y_min/sigma_hat(xi))
					min_FCEG = lambda xi: np.log(delta_y_min/sigma_hat(xi))
			elif var_int == 'nu2':
				min_FCEL = lambda xi: -1*np.log((np.exp(delta_y_FCEL_max(xi))-1.)/sigma_hat(xi))
				min_FCEG = lambda xi: -1*np.log((np.exp(delta_y_FCEG_max(xi))-1.)/sigma_hat(xi))
				if delta_y_min > 1e-15:
					max_FCEL = lambda xi: -1*np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
					max_FCEG = lambda xi: -1*np.log((np.exp(delta_y_min)-1)/sigma_hat(xi))
				elif delta_y_min > 0:
					max_FCEL = lambda xi: -1*np.log(delta_y_min/sigma_hat(xi))
					max_FCEG = lambda xi: -1*np.log(delta_y_min/sigma_hat(xi))
			elif var_int == 'delta':
				max_FCEL = delta_y_FCEL_max
				max_FCEG = delta_y_FCEG_max
				min_FCEL = lambda xi:delta_y_min
				min_FCEG = lambda xi:delta_y_min
			Integrand_FCEL = self.FCEL_integrand(y,x_T,num,mu_factor,mu_f_factor,iso,n_f,switch,var_int=var_int,q0=q0)
			Int_FCEL = integrate.quad(lambda xi: integrate.quad(lambda nu: Integrand_FCEL(nu,xi,min_FCEL(xi)),min_FCEL(xi),max_FCEL(xi))[0],Xi_min(y,x_T),Xi_max(y,x_T))
			Integrand_FCEG = self.FCEG_integrand(y,x_T,num,mu_factor,mu_f_factor,n_f,switch,var_int=var_int,q0=q0)
			Int_FCEG = integrate.quad(lambda xi: integrate.quad(lambda nu: Integrand_FCEG(nu,xi,min_FCEG(xi)),min_FCEG(xi),max_FCEG(xi))[0],Xi_min(y,x_T),Xi_max(y,x_T))
			sigma_FCEL = conv_fact*Int_FCEL[0]
			err_sigma_FCEL = conv_fact*Int_FCEL[1]
			sigma_FCEG = conv_fact*Int_FCEG[0]
			err_sigma_FCEG = conv_fact*Int_FCEG[1]
		return [(sigma_FCEL,err_sigma_FCEL),(sigma_FCEG,err_sigma_FCEG)]
	
	def Uncertainties_pA_dy(self,x_T,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0_list=[0.07,0.05,0.09],var_err='q0,mu,pdf'):
		'''Return the uncertinites on the final sigma_pA (i.e taking in count FCEL
		and FCEG effects), with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0_list the transport coefficient with its uncertaintites as [q0,q0-,q0+]
		- var_err (str) containing the error on variable wanted like: "q0","mu","pdf" (can contains multiple ones like:"q0,pdf")
		as:
			[Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
		so :
		- for the 0 index the central value
		- for the 1 index the + and - total uncertainties
		- for the 2 index the + and - uncertainties for (in order) q0, mu, and pdfs'''
		p_set = self.p_set
		Uplus_pdf, Uminus_pdf = [],[]
		Y = Y_list(x_T)
		q_list = [q0_list[1],q0_list[2]]
		q_val= []
		factor = [1./2.,2]
		mu_factor_list = [(i,j) for i in factor for j in factor]
		mu_val =[]
		q0_cen =q0_list[0]
		num_cen = 0
		pdf_val = []
		print('Central value')
		FCEL_cen, FCEG_cen = self.FCEL_G_integration_dy(x_T,num_cen,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
		Ucen = np.array(FCEL_cen[0]+FCEG_cen[0])
		Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf = np.zeros_like(Y),np.zeros_like(Y),np.zeros_like(Y),np.zeros_like(Y),np.zeros_like(Y),np.zeros_like(Y)
		pdf_val.append(FCEL_cen[0]+FCEG_cen[0])
		if 'q0' in var_err:
			for q in q_list:
				print('q0 = '+str(q))
				FCEL_q , FCEG_q = self.FCEL_G_integration_dy(x_T,num_cen,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q)
				q_val.append(FCEL_q[0]+FCEG_q[0])
			Uminus_q , Uplus_q = min_max(q_val)
			Uminus_q, Uplus_q = Ucen - Uminus_q, Uplus_q - Ucen
		if 'mu' in var_err:
			for (a,b) in mu_factor_list:
				print('(mu_factor,mu_f_factor) = '+str((a,b)))
				FCEL_mu , FCEG_mu = self.FCEL_G_integration_dy(x_T,num_cen,mu_factor=a*mu_factor,mu_f_factor=b*mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
				mu_val.append(FCEL_mu[0]+FCEG_mu[0])
			Uminus_mu, Uplus_mu = min_max(mu_val)
			Uminus_mu , Uplus_mu = Ucen-Uminus_mu, Uplus_mu-Ucen 
		if 'pdf' in var_err:
			print('Computation on pdf set')
			for j,y in enumerate(Y):
				print('y = '+str(y)+'| '+str(j+1)+'/'+str(N_y))
				y_pdf_val =[]
				for i in range(1,p_set.size):
					FCEL_pdf, FCEG_pdf = self.FCEL_G_integration_dydpt(y,x_T,i,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
					y_pdf_val.append(FCEL_pdf[0]+FCEG_pdf[0])
				unc_pdf = p_set.uncertainty(y_pdf_val)
				Uplus_pdf.append(unc_pdf.errplus); Uminus_pdf.append(unc_pdf.errminus)
			Uplus_pdf , Uminus_pdf = np.array(Uplus_pdf), np.array(Uminus_pdf)
		if 'pdf' not in var_err:
			Uplus_pdf, Uminus_pdf = np.zeros_like(Y),np.zeros_like(Y)
		if 'mu' not in var_err:
			Uplus_mu, Uminus_mu = np.zeros_like(Y),np.zeros_like(Y)
		if 'q0' not in var_err:
			Uplus_q, Uminus_q = np.zeros_like(Y),np.zeros_like(Y)
		Uplus = np.sqrt(Uplus_q**2+Uplus_mu**2+Uplus_pdf**2)
		Uminus = np.sqrt(Uminus_q**2+Uminus_mu**2+Uminus_pdf**2)
		return [Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
	
	def Uncertainties_pA_dpt(self,y,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0_list=[0.07,0.05,0.09],var_err='q0,mu,pdf'):
		'''Return the uncertinites on the final sigma_pA (i.e taking in count FCEL
		and FCEG effects), with:
		- y the rapidity
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0_list the transport coefficient with its uncertaintites as [q0,q0-,q0+]
		- var_err (str) containing the error on variable wanted like: "q0","mu","pdf" (can contains multiple ones like:"q0,pdf")
		as:
			[Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
		so :
		- for the 0 index the central value
		- for the 1 index the + and - total uncertainties
		- for the 2 index the + and - uncertainties for (in order) q0, mu, and pdfs'''
		p_set = self.p_set
		Uplus_pdf, Uminus_pdf = [],[]
		P_T = self.P_T_list(y)
		rs = self.rs
		q_list = [q0_list[1],q0_list[2]]
		q_val= []
		factor = [1./2.,2]
		mu_factor_list = [(i,j) for i in factor for j in factor]
		mu_val =[]
		q0_cen =q0_list[0]
		num_cen = 0
		pdf_val = []
		print('Central value')
		FCEL_cen, FCEG_cen = self.FCEL_G_integration_dpt(y,num_cen,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
		Ucen = np.array(FCEL_cen[0]+FCEG_cen[0])
		Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf = np.zeros_like(P_T),np.zeros_like(P_T),np.zeros_like(P_T),np.zeros_like(P_T),np.zeros_like(P_T),np.zeros_like(P_T)
		pdf_val.append(FCEL_cen[0]+FCEG_cen[0])
		if 'q0' in var_err:
			for q in q_list:
				print('q0 = '+str(q))
				FCEL_q , FCEG_q = self.FCEL_G_integration_dpt(y,num_cen,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q)
				q_val.append(FCEL_q[0]+FCEG_q[0])
			Uminus_q , Uplus_q = min_max(q_val)
			Uminus_q, Uplus_q = Ucen - Uminus_q, Uplus_q - Ucen
		if 'mu' in var_err:
			for (a,b) in mu_factor_list:
				print('(mu_factor,mu_f_factor) = '+str((a,b)))
				FCEL_mu , FCEG_mu = self.FCEL_G_integration_dpt(y,num_cen,mu_factor=a*mu_factor,mu_f_factor=b*mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
				mu_val.append(FCEL_mu[0]+FCEG_mu[0])
			Uminus_mu, Uplus_mu = min_max(mu_val)
			Uminus_mu , Uplus_mu = Ucen-Uminus_mu, Uplus_mu-Ucen 
		if 'pdf' in var_err:
			print('Computation on pdf set')
			for j,pt in enumerate(P_T):
				print('pt = '+str(pt)+'| '+str(j+1)+'/'+str(N_pt))
				pt_pdf_val =[]
				x_T = 2*pt/rs
				for i in range(1,p_set.size):
					FCEL_pdf, FCEG_pdf = self.FCEL_G_integration_dydpt(y,x_T,i,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
					pt_pdf_val.append(FCEL_pdf[0]+FCEG_pdf[0])
				unc_pdf = p_set.uncertainty(pt_pdf_val)
				Uplus_pdf.append(unc_pdf.errplus); Uminus_pdf.append(unc_pdf.errminus)
			Uplus_pdf , Uminus_pdf = np.array(Uplus_pdf), np.array(Uminus_pdf)
		if 'pdf' not in var_err:
			Uplus_pdf, Uminus_pdf = np.zeros_like(P_T),np.zeros_like(P_T)
		if 'mu' not in var_err:
			Uplus_mu, Uminus_mu = np.zeros_like(P_T),np.zeros_like(P_T)
		if 'q0' not in var_err:
			Uplus_q, Uminus_q = np.zeros_like(P_T),np.zeros_like(P_T)
		Uplus = np.sqrt(Uplus_q**2+Uplus_mu**2+Uplus_pdf**2)
		Uminus = np.sqrt(Uminus_q**2+Uminus_mu**2+Uminus_pdf**2)
		return [Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
	
	def Uncertainties_RpA_dy(self,x_T,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0_list=[0.07,0.05,0.09],var_err='q0,mu,pdf'):
		'''Return the uncertinites on the final R_pA (i.e taking in count FCEL
		and FCEG effects), with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0_list the transport coefficient with its uncertaintites as [q0,q0-,q0+]
		- var_err (str) containing the error on variable wanted like: "q0","mu","pdf" (can contains multiple ones like:"q0,pdf")
		as:
			[Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
		so :
		- for the 0 index the central value
		- for the 1 index the + and - total uncertainties
		- for the 2 index the + and - uncertainties for (in order) q0, mu, and pdfs'''
		p_set = self.p_set
		Uplus_pdf, Uminus_pdf = [],[]
		Y = Y_list(x_T)
		q_list = [q0_list[1],q0_list[2]]
		q_val= []
		factor = [1./2.,2]
		mu_factor_list = [(i,j) for i in factor for j in factor]
		mu_val =[]
		q0_cen =q0_list[0]
		num_cen = 0
		pdf_val = []
		print('Central value')
		FCEL_cen, FCEG_cen = self.FCEL_G_integration_dy(x_T,num_cen,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
		sigma_pp_cen,err_pp_cen = self.dsigma_tot_dy(x_T,num_cen, mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso =iso,n_f=n_f,is_pp=True,switch =switch)
		Ucen = np.array((FCEL_cen[0]+FCEG_cen[0])/sigma_pp_cen)
		Uplus_q,Uminus_q,Uplus_mu,Uminus_mu= np.zeros_like(Y),np.zeros_like(Y), np.zeros_like(Y),np.zeros_like(Y)
		pdf_val.append((FCEL_cen[0]+FCEG_cen[0])/sigma_pp_cen)
		if 'q0' in var_err:
			for q in q_list:
				print('q0 = '+str(q))
				FCEL_q , FCEG_q = self.FCEL_G_integration_dy(x_T,num_cen,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q)
				q_val.append((FCEL_q[0]+FCEG_q[0])/sigma_pp_cen)
			Uminus_q , Uplus_q = min_max(q_val)
			Uminus_q, Uplus_q = Ucen - Uminus_q, Uplus_q - Ucen
		if 'mu' in var_err:
			for (a,b) in mu_factor_list:
				print('(mu_factor,mu_f_factor) = '+str((a,b)))
				FCEL_mu , FCEG_mu = self.FCEL_G_integration_dy(x_T,num_cen,mu_factor= a*mu_factor,mu_f_factor=b*mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
				sigma_pp_mu , err_pp_mu = self.dsigma_tot_dy(x_T,num_cen,mu_factor= a*mu_factor,mu_f_factor=b*mu_f_factor,iso = iso,n_f=n_f,is_pp=True,switch =switch)
				mu_val.append((FCEL_mu[0]+FCEG_mu[0])/sigma_pp_mu)
			Uminus_mu, Uplus_mu = min_max(mu_val)
			Uminus_mu , Uplus_mu = Ucen-Uminus_mu, Uplus_mu-Ucen
		if 'pdf' in var_err:
			print('Computation on pdf set')
			for j,y in enumerate(Y):
				print('y = '+str(y)+'| '+str(j+1)+'/'+str(N_y))
				y_pdf_val =[]
				for i in range(1,p_set.size):
					FCEL_pdf, FCEG_pdf = self.FCEL_G_integration_dydpt(y,x_T,i,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
					sigma_pp_pdf , err_pp_pdf = self.dsigma_tot_dydpt(y,x_T,i,mu_factor,mu_f_factor,iso = iso,n_f=n_f,is_pp=True,switch =switch)
					y_pdf_val.append((FCEL_pdf[0]+FCEG_pdf[0])/sigma_pp_pdf)
				unc_pdf = p_set.uncertainty(y_pdf_val)
				Uplus_pdf.append(unc_pdf.errplus); Uminus_pdf.append(unc_pdf.errminus)
			Uplus_pdf , Uminus_pdf = np.array(Uplus_pdf), np.array(Uminus_pdf)
		if 'pdf' not in var_err:
			Uplus_pdf, Uminus_pdf = np.zeros_like(Y),np.zeros_like(Y)
		if 'mu' not in var_err:
			Uplus_mu, Uminus_mu = np.zeros_like(Y),np.zeros_like(Y)
		if 'q0' not in var_err:
			Uplus_q, Uminus_q = np.zeros_like(Y),np.zeros_like(Y)
		Uplus = np.sqrt(Uplus_q**2+Uplus_mu**2+Uplus_pdf**2)
		Uminus = np.sqrt(Uminus_q**2+Uminus_mu**2+Uminus_pdf**2)
		return [Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
	
	def Uncertainties_RpA_dpt(self,y,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0_list=[0.07,0.05,0.09],var_err='q0,mu,pdf'):
		'''Return the uncertinites on the final R_pA (i.e taking in count FCEL
		and FCEG effects), with:
		- y the rapidity
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0_list the transport coefficient with its uncertaintites as [q0,q0-,q0+]
		- var_err (str) containing the error on variable wanted like: "q0","mu","pdf" (can contains multiple ones like:"q0,pdf")
		as:
			[Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
		so :
		- for the 0 index the central value
		- for the 1 index the + and - total uncertainties
		- for the 2 index the + and - uncertainties for (in order) q0, mu, and pdfs'''
		p_set = self.p_set
		Uplus_pdf, Uminus_pdf = [],[]
		P_T = self.P_T_list(y)
		rs = self.rs
		q_list = [q0_list[1],q0_list[2]]
		q_val= []
		factor = [1./2.,2]
		mu_factor_list = [(i,j) for i in factor for j in factor]
		mu_val =[]
		q0_cen =q0_list[0]
		num_cen = 0
		pdf_val = []
		print('Central value')
		FCEL_cen, FCEG_cen = self.FCEL_G_integration_dpt(y,num_cen,iso=iso,n_f=n_f,switch =switch,eps = eps,mu_factor=mu_factor,mu_f_factor=mu_f_factor,var_int=var_int,q0=q0_cen)
		sigma_pp_cen,err_pp_cen = self.dsigma_tot_dpt(y,num_cen,iso=iso,n_f=n_f,is_pp=True,switch=switch,mu_factor=mu_factor,mu_f_factor=mu_f_factor)
		Ucen = np.array((FCEL_cen[0]+FCEG_cen[0])/sigma_pp_cen)
		pdf_val.append((FCEL_cen[0]+FCEG_cen[0])/sigma_pp_cen)
		if 'q0' in var_err:
			for q in q_list:
				print('q0 = '+str(q))
				FCEL_q , FCEG_q = self.FCEL_G_integration_dpt(y,num_cen,iso=iso,n_f=n_f,switch =switch,eps = eps,mu_factor=mu_factor,mu_f_factor=mu_f_factor,var_int=var_int,q0=q)
				q_val.append((FCEL_q[0]+FCEG_q[0])/sigma_pp_cen)
			Uminus_q , Uplus_q = min_max(q_val)
			Uminus_q, Uplus_q = Ucen - Uminus_q, Uplus_q - Ucen
		if 'mu' in var_err:
			for (a,b) in mu_factor_list:
				print('(mu_factor,mu_f_factor) = '+str((a,b)))
				FCEL_mu , FCEG_mu = self.FCEL_G_integration_dpt(y,num_cen,iso=iso,n_f=n_f,switch =switch,eps = eps,mu_factor=a*mu_factor,mu_f_factor=b*mu_f_factor,var_int=var_int,q0=q0_cen)
				sigma_pp_mu , err_pp_mu = self.dsigma_tot_dpt(y,num_cen,iso = iso,n_f=n_f,is_pp=True,switch =switch,mu_factor=a*mu_factor,mu_f_factor=b*mu_f_factor)
				mu_val.append((FCEL_mu[0]+FCEG_mu[0])/sigma_pp_mu)
			Uminus_mu, Uplus_mu = min_max(mu_val)
			Uminus_mu , Uplus_mu = Ucen-Uminus_mu, Uplus_mu-Ucen
		if 'pdf' in var_err:
			print('Computation on pdf set')
			for j,pt in enumerate(P_T):
				print('pt = '+str(pt)+'| '+str(j+1)+'/'+str(N_pt))
				pt_pdf_val =[]
				x_T = 2*pt/rs
				for i in range(1,p_set.size):
					FCEL_pdf, FCEG_pdf = self.FCEL_G_integration_dydpt(y,x_T,i,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0_cen)
					sigma_pp_pdf , err_pp_pdf = self.dsigma_tot_dydpt(y,x_T,i,mu_factor,mu_f_factor,iso = iso,n_f=n_f,is_pp=True,switch =switch)
					pt_pdf_val.append((FCEL_pdf[0]+FCEG_pdf[0])/sigma_pp_pdf)
				unc_pdf = p_set.uncertainty(pt_pdf_val)
				Uplus_pdf.append(unc_pdf.errplus); Uminus_pdf.append(unc_pdf.errminus)
			Uplus_pdf , Uminus_pdf = np.array(Uplus_pdf), np.array(Uminus_pdf)
		if 'pdf' not in var_err:
			Uplus_pdf, Uminus_pdf = np.zeros_like(P_T),np.zeros_like(P_T)
		if 'mu' not in var_err:
			Uplus_mu, Uminus_mu = np.zeros_like(P_T),np.zeros_like(P_T)
		if 'q0' not in var_err:
			Uplus_q, Uminus_q = np.zeros_like(P_T),np.zeros_like(P_T)
		Uplus = np.sqrt(Uplus_q**2+Uplus_mu**2+Uplus_pdf**2)
		Uminus = np.sqrt(Uminus_q**2+Uminus_mu**2+Uminus_pdf**2)
		return [Ucen,[Uplus,Uminus],[Uplus_q,Uminus_q,Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
	
	def Rpp_FCELG_dy(self,x_T,num,q0,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu'):
		'''Return the ratio of the shifted sigma_pp (or pn if iso = 'n') 
		over sigma_pp (resp pn), with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- q0 the transport coefficient
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)'''
		FCEL_cen, FCEG_cen = self.FCEL_G_integration_dy(x_T,num,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0)
		sigma_pp_cen,err_pp_cen = self.dsigma_tot_dy(x_T,num, mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,is_pp=True,switch =switch)
		return (FCEL_cen[0]+FCEG_cen[0])/sigma_pp_cen
	
	def Delta_RpA_plusminus_dy(self,x_T,num=0,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',q0_list = [0.07,0.05,0.09]):
		'''Return the relative RpA variation on q0 by icresing or reducing it,
		with:
		- x_T = 2*p_T/√s ,p_T the transverse momentum
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0_list the transport coefficient with its uncertaintites as [q0,q0-,q0+]
		as:
		- ((Rpa_plus-Rpa)/Rpa,(Rpa_minus-Rpa)/Rpa)'''
		q0_cen , q0_minus, q0_plus = q0_list
		Rpa = self.Rpp_FCELG_dy(x_T,num,q0_cen,mu_factor,mu_f_factor,iso,n_f,switch ,eps ,var_int)
		Rpa_plus = self.Rpp_FCELG_dy(x_T,num,q0_plus,mu_factor,mu_f_factor,iso,n_f,switch ,eps ,var_int)
		Rpa_minus = self.Rpp_FCELG_dy(x_T,num,q0_minus,mu_factor,mu_f_factor,iso,n_f,switch ,eps ,var_int)
		return ((Rpa_plus-Rpa)/Rpa,(Rpa_minus-Rpa)/Rpa)
	
	def Rpp_FCELG_dpt(self,y,num,q0,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu'):
		'''Return the ratio of the shifted sigma_pp (or pn if iso = 'n') 
		over sigma_pp (resp pn), with:
		- y the rapidity
		- num the member of the pdf set
		- q0 the transport coefficient 
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)'''
		FCEL_cen, FCEG_cen = self.FCEL_G_integration_dpt(y,num,mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso=iso,n_f=n_f,switch =switch,eps = eps,var_int=var_int,q0=q0)
		sigma_pp_cen,err_pp_cen = self.dsigma_tot_dpt(y,num, mu_factor= mu_factor,mu_f_factor=mu_f_factor,iso =iso,n_f=n_f,is_pp=True,switch =switch)
		return (FCEL_cen[0]+FCEG_cen[0])/sigma_pp_cen
	
	def Delta_RpA_plusminus_dpt(self,y,num=0,mu_factor=1,mu_f_factor=1,iso='p',n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu'):
		'''Return the relative RpA variation on q0 by icresing or reducing it,
		with:
		- y the rapidity
		- num the member of the pdf set
		- mu_factor that describes mu = p_T*mu_factor (same for mu_f_factor)
		- iso the isospin variable (='p' by default)
		- n_f the number of flavours (=3 by default)
		- switch the convention of p_T integration (= 'dp_t' by default)
		- eps the minimun delta_y can value
		- var_int (str) the integration variable (="nu" by default)
		- q0_list the transport coefficient with its uncertaintites as [q0,q0-,q0+]
		as:
		- ((Rpa_plus-Rpa)/Rpa,(Rpa_minus-Rpa)/Rpa)'''
		q0_plus = 0.09
		q0_minus = 0.05
		q0_cen = 0.07
		Rpa = self.Rpp_FCELG_dpt(y,num,q0_cen,mu_factor,mu_f_factor,iso,n_f,switch ,eps ,var_int)
		Rpa_plus = self.Rpp_FCELG_dpt(y,num,q0_plus,mu_factor,mu_f_factor,iso,n_f,switch ,eps ,var_int)
		Rpa_minus = self.Rpp_FCELG_dpt(y,num,q0_minus,mu_factor,mu_f_factor,iso,n_f,switch ,eps ,var_int)
		return ((Rpa_plus-Rpa)/Rpa,(Rpa_minus-Rpa)/Rpa)
		
	# Isospin effect with FCEL/G
	
	def Uncertainties_R_iso_dy(self,x_T,mu_factor=1,mu_f_factor=1,n_f=3,switch ='dp_t',var_err='mu,pdf'):
		''''''
		p_set = self.p_set
		Ucen = []
		Uplus_pdf = []
		Uminus_pdf = []
		Uplus_mu = []
		Uminus_mu= []
		num_cen= 0
		Y = Y_list(x_T)
		factor = [1./2.,2]
		mu_factor_list = [(i,j) for i in factor for j in factor]
		mu_val =[]
		Uplus_mu,Uminus_mu= np.zeros_like(Y),np.zeros_like(Y)
		print('central value')
		R_cen =  self.R_iso_dy(x_T,num_cen,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,switch=switch)[0]
		if 'mu' in var_err:
			for (a,b) in mu_factor_list:
				print('(mu_factor,mu_f_factor) = '+str((a,b)))
				R_mu , err_mu = self.R_iso_dy(x_T,num_cen,mu_factor=a*mu_factor,mu_f_factor=b*mu_f_factor,n_f=n_f,switch=switch)
				mu_val.append(R_mu)
			Uminus_mu, Uplus_mu = min_max(mu_val)
			Uminus_mu , Uplus_mu = R_cen-Uminus_mu, Uplus_mu-R_cen
			Uplus_mu,Uminus_mu = np.array(Uplus_mu), np.array(Uminus_mu)
			Ucen = R_cen
		if 'pdf' in var_err:
			Ucen_pdf = []
			for j,y in enumerate(Y):
				print('y = '+str(y) +', '+str(j)+'/'+str(N_y-1))
				dsigma = []
				for i in range(p_set.size):
					dsigma.append(self.R_iso_dydpt(y,x_T,i,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,switch=switch)[0])
				unc = p_set.uncertainty(dsigma)
				Ucen_pdf.append(unc.central)
				Uplus_pdf.append(unc.errplus)
				Uminus_pdf.append(unc.errminus)
			Ucen_pdf, Uplus_pdf, Uminus_pdf = np.array(Ucen_pdf),np.array(Uplus_pdf),np.array(Uminus_pdf)
		if not('pdf' in var_err):
			Uplus_pdf, Uminus_pdf = np.zeros_like(Y),np.zeros_like(Y)
		if not('mu' in var_err) and ('pdf'in var_err):
			Ucen = Ucen_pdf
		elif not('mu' in var_err) and not('pdf'in var_err):
			raise ValueError("no error variable adequate to this function")
		Uplus = np.sqrt(Uplus_mu**2+Uplus_pdf**2)
		Uminus = np.sqrt(Uminus_mu**2+Uminus_pdf**2)
		return[Ucen,(Uplus,Uminus),[Uplus_mu,Uminus_mu,Uplus_pdf,Uminus_pdf]]
	
	def Uncertainties_4R_iso_dy(self,x_T,mu_factor=1,mu_f_factor=1,n_f=3,switch ='dp_t',eps = 1e-15,var_int='nu',var_err='q0,mu,pdf'):
		'''Return the 4 ratios for the isospin study:
		- RpA (or called Rpp^FCEL/G) = sigma_pp^FCEL/sigma_pp
		- RpA^Iso = sigma_pA^Iso/sigma_pp ( with sigma_pA^Iso = Z*sigma_pp+(A-Z)*sigma_pn)
		- RpA^(FCEL/G+Iso) = sigma_pA^(Iso + FCEL/G)/sigma_pp
		- RpA^(FCEL/G - Iso) = RpA^(FCEL/G+Iso)/RpA^Iso
		, implementing the FCEL/G effect, with their uncertainties as :
		 [(Rpp,Rpp_plus,Rpp_minus),(Riso,Riso_plus,Riso_minus),(Riso_FCELG,Riso_FCELG_plus,Riso_FCELG_minus),(Rwoiso,Rwoiso_plus,Rwoiso_minus)]'''
		Z = self.Z
		A = self.A
		Riso ,err_Riso, err_var_Riso = self.Uncertainties_R_iso_dy(x_T,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,switch =switch,var_err=var_err)
		Riso_plus,Riso_minus = err_Riso[0],err_Riso[1]
		Fnp =  (Riso-Z/A)*A/(A-Z)
		Fnp_plus, Fnp_minus = A/(A-Z)*Riso_plus, A/(A-Z)*Riso_minus 
		Rpp ,err_Rpp, err_var_Rpp = self.Uncertainties_RpA_dy(x_T,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,switch =switch,eps = eps,var_int=var_int,var_err=var_err)
		Rpp_plus,Rpp_minus = err_Rpp[0],err_Rpp[1]
		Rpn ,err_Rpn, err_var_Rpn = self.Uncertainties_RpA_dy(x_T,mu_factor=mu_factor,mu_f_factor=mu_f_factor,n_f=n_f,switch =switch,eps = eps,var_int=var_int,var_err=var_err)
		Rpn_plus,Rpn_minus = err_Rpn[0],err_Rpn[1]
		Riso_FCELG = (Z/A)*(Rpp-Fnp*Rpn)+Fnp*Rpn
		Riso_FCELG_plus,Riso_FCELG_minus=  ((Rpp_plus*Z/A)**2+(Fnp_plus**2+Rpn_plus**2)*(1-Z/A)**2)**0.5, ((Rpp_minus*Z/A)**2+(Fnp_minus**2+Rpn_minus**2)*(1-Z/A)**2)**0.5
		Rwoiso = Riso_FCELG/Riso
		Rwoiso_plus,Rwoiso_minus = Rwoiso*((Riso_FCELG_plus/Riso_FCELG)**2+(Riso_plus/Riso)**2)**0.5,Rwoiso*((Riso_minus/Riso)**2+(Riso_FCELG_minus/Riso_FCELG)**2)**0.5
		return [(Rpp,Rpp_plus,Rpp_minus),(Riso,Riso_plus,Riso_minus),(Riso_FCELG,Riso_FCELG_plus,Riso_FCELG_minus),(Rwoiso,Rwoiso_plus,Rwoiso_minus)]
		
	# debug for xi fixed
	
	def FCEL_G_integration_dy_dxi(self,x_T,xi,mu2,mu_f2,num,n_f=3,switch ='dp_t'):
		''''''
		Y = Y_list_xi(x_T,xi)
		sigma_FCEL = []
		err_sigma_FCEL = []
		sigma_FCEG = []
		err_sigma_FCEG = []
		iso ='p'
		for y in Y:
			print('y = '+str(y))
			delta_y_FCEL_max = min(max(Y)-y,np.log(2))
			delta_y_FCEG_max = min(y-min(Y),np.log(2))
			A = (delta_y_FCEL_max == 0.)
			B = (delta_y_FCEG_max == 0.)
			if A and B :
				sigma_FCEL.append(np.nan)
				err_sigma_FCEL.append(np.nan)
				sigma_FCEG.append(np.nan)
				err_sigma_FCEG.append(np.nan)
			elif A and not(B):
				sigma_FCEL.append(np.nan)
				err_sigma_FCEL.append(np.nan)
				Integrand_FCEG = lambda delta_y: self.FCEG_integrand(y,x_T,mu2,mu_f2,num,n_f,switch)(delta_y,xi)
				Int_FCEG = integrate.quad(Integrand_FCEG,0,delta_y_FCEG_max)
				sigma_FCEG.append(conv_fact*Int_FCEG[0])
				err_sigma_FCEG.append(conv_fact*Int_FCEG[1])
			elif not(A) and B:
				Integrand_FCEL = lambda delta_y : self.FCEL_integrand(y,x_T,mu2,mu_f2,num,iso,n_f,switch)(delta_y,xi)
				Int_FCEL = integrate.quad(Integrand_FCEL,0,delta_y_FCEL_max)
				sigma_FCEL.append(conv_fact*Int_FCEL[0])
				err_sigma_FCEL.append(conv_fact*Int_FCEL[1])
				sigma_FCEG.append(np.nan)
				err_sigma_FCEG.append(np.nan)
			elif not(A) and not(B):
				Integrand_FCEL = lambda delta_y : self.FCEL_integrand(y,x_T,mu2,mu_f2,num,iso,n_f,switch)(delta_y,xi)
				Int_FCEL = integrate.quad(Integrand_FCEL,0,delta_y_FCEL_max)
				Integrand_FCEG = lambda delta_y: self.FCEG_integrand(y,x_T,mu2,mu_f2,num,n_f,switch)(delta_y,xi)
				Int_FCEG = integrate.quad(Integrand_FCEG,0,delta_y_FCEG_max)
				sigma_FCEL.append(conv_fact*Int_FCEL[0])
				err_sigma_FCEL.append(conv_fact*Int_FCEL[1])
				sigma_FCEG.append(conv_fact*Int_FCEG[0])
				err_sigma_FCEG.append(conv_fact*Int_FCEG[1])
		return [(sigma_FCEL,err_sigma_FCEL),(sigma_FCEG,err_sigma_FCEG)]
	
	### Other functions ###
	
	def P_T_list(self,y,p_t_min = 3,p_t_max= 15):
		'''The numpy linspace of transverse momentum in GeV, with:
		- y the rapidity
		- p_t_min the minimum wanted if the phase space is too large
		- p_t_max same for the maximum wanted'''
		rs = self.rs
		p_max = min(p_t_max,min(rs*np.exp(y),rs*np.exp(-y)))
		P_T = np.linspace(p_t_min,p_max,N_pt)
		return P_T
	
	def Nu_list(self,y,xi,x_T,num,mu2,N_nu = 100,eps = 1e-10):
		''''''
		A = self.A
		rs = self.rs
		p_T = rs*x_T/2.
		alpha_s = self.alpha_s_p(num,mu2)
		Fc_FCEG = -1./N_c
		Fc_FCEL = N_c
		prob_FCEG = p.proba(A,1,rs,p_T,y,alpha_s,Fc_FCEG,0)
		prob_FCEL = p.proba(A,1,rs,p_T,y,alpha_s,Fc_FCEL,0)
		sigma_hat_FCEG = prob_FCEG.sigma_hat(xi)
		sigma_hat_FCEL = prob_FCEL.sigma_hat(xi)
		k = 2./x_T
		delta_y_FCEL_max = min(np.log(k*xi)-y,np.log(2))
		delta_y_FCEG_max = min(y+np.log((1-xi)*k),np.log(2))
		delta_y_min = eps
		nu_max_FCEL = np.log((np.exp(delta_y_FCEL_max)-1.)/sigma_hat_FCEL)
		nu_max_FCEG = np.log((np.exp(delta_y_FCEG_max)-1.)/sigma_hat_FCEG)
		if delta_y_min > 0:
			nu_min_FCEL = np.log((np.exp(delta_y_min)-1)/sigma_hat_FCEL)
			nu_min_FCEG = np.log((np.exp(delta_y_min)-1)/sigma_hat_FCEG)
		else:
			nu_min_FCEL = -1e-60
			nu_min_FCEG = -1e-60
		FCEL = np.linspace(nu_min_FCEL,nu_max_FCEL,N_nu)
		FCEG = np.linspace(nu_min_FCEG,nu_max_FCEG,N_nu)
		return (FCEL,FCEG)
