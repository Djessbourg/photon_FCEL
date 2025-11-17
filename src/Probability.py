import numpy as np
from scipy.special import spence

# Based on the articles of Francois Arleo about FCEL effect in QCD diagrams

# attention, la fonction spence(z) ici est defini comme l'integrale de 1 a z de ln(t)/(1-t)dt, alors que le dilogarithme est defini comme etant l'integrale de 1 a 1-z du meme integrand
# on a donc Li_2(z) = spence(1-z)
def dilog(z):
    '''fonction dilogarythme'''
    return spence(1-z)

###On code les fonctions generales ici###
def heav(x):
    """ retourne la fonction de heaviside pour x la fraction des elements a dissocier"""
    if x > 1:
        return 1
    else:
        return 0

def dI_u(u,chi,Fc,alpha =0.5):
    '''retourne dI/du avec u = omega/omega_hat'''
    return (np.abs(Fc)*alpha)/(u*np.pi)*(np.log(1+pow(1./u,2))-np.log(1+pow(1./u,2)*chi))

def u_dI_u(u,chi,Fc,alpha=0.5):
    '''retourne udI/du avec u = omega/omega_hat'''
    return (np.abs(Fc)*alpha)/(np.pi)*(np.log(1+pow(1./u,2))-np.log(1+pow(1./u,2)*chi))

def first_dilog_u(u,chi):
    """renvoie la premiere evaluation du dilogarithme dans l'expression de P tilde"""
    return dilog(-pow(1./u,2)*chi)
    
def second_dilog_u(u):
    """renvoie la seconde evaluation du dilogarithme dans l'expression de P tilde"""
    return dilog(-pow(1./u,2))
    
def g_u(u,chi,Fc,alpha=0.5):
    """revoie la version approchee"""
    return alpha*np.abs(Fc)*np.log(chi)*(2*np.log(u)-0.5*np.log(chi))/(2*np.pi)

def p_tilde_u(u,chi,Fc,alpha = 0.5):
    """renvoie la densitee de probabilite en fonction des variables sans dimensions u et chi """
    return dI_u(u,chi,Fc,alpha)*np.exp(-(np.abs(Fc)*alpha)/(2*np.pi)*(first_dilog_u(u,chi)-second_dilog_u(u)))
        
def p_tilde_u_approx(u,chi,Fc,alpha = 0.5):
    """renvoie la densitee de probabilite en fonction des variables sans dimensions u et chi """
    return dI_u(u,chi,Fc,alpha)*np.exp(-1*g_u(u,chi,Fc,alpha))

def p_tilde_u_advanced(u,chi,Fc,umin,alpha=0.5):
    if u > umin:
        return p_tilde_u(u,chi,Fc,alpha)
    elif u <= umin:
        return p_tilde_u_approx(u,chi,Fc,alpha)

def Exp_u(u,chi,Fc,alpha=0.5):
    return np.exp(-(np.abs(Fc)*alpha)/(2*np.pi)*(first_dilog_u(u,chi)-second_dilog_u(u)))
    
def dI_nu(nu,chi,alpha):
    '''retourne dI/du avec u = omega/omega_hat'''
    return (Nc*alpha)/(nu*np.pi)*(np.log(1+pow(nu,2))-np.log(1+pow(nu,2)*chi))*heav(1-chi)

def nu_dI_nu(nu,chi,alpha):
    '''retourne udI/du avec u = omega/omega_hat'''
    return (Nc*alpha)/(np.pi)*(np.log(1+pow(nu,2))-np.log(1+pow(nu,2)*chi))*heav(1-chi)

def first_dilog_nu(nu,chi):
    """renvoie la premiere evaluation du dilogarithme dans l'expression de P tilde"""
    return dilog(-pow(nu,2)*chi)
    
def second_dilog_nu(nu):
    """renvoie la seconde evaluation du dilogarithme dans l'expression de P tilde"""
    return dilog(-pow(nu,2))

def p_tilde_nu(nu,chi,alpha):
    """renvoie la densitee de probabilite en fonction des variables sans dimensions u et chi """
    return dI_nu(nu,chi)*np.exp(-(Nc*alpha)/(2*np.pi)*(first_dilog_nu(nu,chi)-second_dilog_nu(nu)))

rn = 1.12 #fm
hbar_c = 0.197 #GeV.fm
mp = 0.938 #GeV
Lambda_QCD = 0.25 #GeV 
Lp = 1.5 #fm
Nc = 3 #Nombre de couleurs
#M(Mbot) = 3GeV pour cc_bar = 9 GeV pour bb_bar

def L(A):
    if A >1:
        return (3./2.)*rn*pow(A,1./3.)
    else :
        return Lp

class proba():
    def __init__(self,A,B,rs,p_t,y,alpha_s,Fc,m,q0=0.07,z=0):
        self.A = A          											#Nombre de masses du noyau cible 
        self.B = B          											#Nombre de masse de la particule incidente
        self.p_t = p_t          										#transverse momentum
        self.rs = rs													#enegie du centre de masse
        self.Ep = self.rs**2/(2*mp)										#approximative proton energy
        self.y = y         											 	#rapidite
        self.alpha_s = alpha_s											#Running alpha_s (depend de l'echelle)
        self.Fc = Fc													#Le facteur de couleur intervenant sous la forme C1+Cr-C2
        self.m = m 														#particle mass
        self.q0 = q0													#transport coefficient constant in GeV^2/fm
        self.m_t = np.sqrt(m**2+p_t**2)
        self.z = z
        self.x2 = lambda xi: np.exp(-1*self.y)*self.M_xi(xi)/self.rs
        self.x0A = hbar_c/(2*mp*L(self.A))
        self.x0B = hbar_c/(2*mp*L(self.B))
        self.Lambda_B = lambda xi: max(self.Lambda(self.B,xi),Lambda_QCD)
        self.Lambda_A = lambda xi: self.Lambda(self.A,xi)

    def xA(self,xi):
        return min(self.x0A,self.x2(xi))
        
    def xB(self,xi):
        return min(self.x0B,self.x2(xi))

    def q(self,A,xi):
        q0 = self.q0
        if A == self.A:
            return q0*pow(0.01/self.xA(xi),0.3)
        elif A == self.B:
            return q0*pow(0.01/self.xB(xi),0.3)

    def Lambda(self,A,xi):
        return np.sqrt(self.q(A,xi)*L(A))

    def M_xi(self,xi):
        '''retrourne la masse invariante en fonction de chi, (voir pour le cas particulier du photon)'''
        m_t = float(self.m_t)
        m = self.m
        z = self.z
        if z != 0:
            m_t = m_t/z
        if m == 0:
            return(m_t/(1-xi))
        elif m >0:
            return(m_t/np.sqrt(xi*(1-xi)))

    def sigma_hat(self,xi):
        return self.Lambda(self.A,xi)/self.M_xi(xi)
        
    def Chi(self,xi):
        return pow(self.Lambda(self.B,xi)/self.Lambda(self.A,xi),2)

    def first_dilog(self,x,xi):
        """renvoie la premiere evaluation du dilogarithme dans l'expression de P"""
        return dilog(-1.*pow((self.Lambda_B(xi))/(self.M_xi(xi)*x),2))
    
    def second_dilog(self,x,xi):
        """renvoie la seconde evaluation du dilogarithme dans l'expression de P"""
        return dilog(-1.*pow((self.Lambda_A(xi))/(self.M_xi(xi)*x),2))
    
    def heaviside(self,xi):
        """fonction 'step' ou marche de Heaviside pour l'operation lbotA^2-lambdaB^2"""
        if ((pow(self.Lambda_A(xi),2)-pow(self.Lambda_B(xi),2)) > 0):
            return 1
        else:
            return 0

    def dI(self,x,xi):
        """retourne dI/dx"""
        return (np.abs(self.Fc)*self.alpha_s)/(x*np.pi)*(np.log(1+pow((self.Lambda_A(xi))/(self.M_xi(xi)*x),2))-np.log(1+pow((self.Lambda_B(xi))/(self.M_xi(xi)*x),2)))*self.heaviside(xi)

    def x_dI(self,x,xi):
        """retourne x * dI/dx"""
        return (np.abs(self.Fc)*self.alpha_s)/(np.pi)*(np.log(1+pow((self.Lambda_A(xi))/(self.M_xi(xi)*x),2))-np.log(1+pow((self.Lambda_B(xi))/(self.M_xi(xi)*x),2)))*self.heaviside(xi)
    
    def Pbar_x(self,x,xi):
        """retourne la fonction densite de probabilite utile pour l'integrale du FCEL"""
        return self.dI(x,xi)*np.exp(-1*np.abs(self.Fc)*self.alpha_s/(2*np.pi)*(self.first_dilog(x,xi)-self.second_dilog(x,xi)))
        
    def Exp_x(self,x,xi):
        return np.exp(-1.*np.abs(self.Fc)*self.alpha_s/(2.*np.pi)*(self.first_dilog(x,xi)-self.second_dilog(x,xi)))
        
    def Arg_x(self,x,xi):
        return (-1.*np.abs(self.Fc)*self.alpha_s)/(2.*np.pi)*(self.first_dilog(x,xi)-self.second_dilog(x,xi))

    def delta_y_max(self):
        return min(np.log(2),np.log(self.rs/self.p_t)-self.y)

    # E = M^2/(2*m_p*x2)
