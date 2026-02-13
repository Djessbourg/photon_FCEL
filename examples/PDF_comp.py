# -*- coding: utf-8 -*-

# =============================================================================
# Annex file to understand the choce of the PDF set 
# last modified : 13/02/2026
# =============================================================================

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
from src import sigma as sig
from src import Collision

### Initialisation of a sigma object ###
rs = 8800
# rs = 200
s = (rs)**2 # CM energy in Gev2

proton = ['NNPDF40_nlo_as_01180'] #,'MSHT20nlo_as118','CT18NLO' Here put whatever PDF_set you want to compare 

Pb = "nNNPDF30_nlo_as_0118_A208_Z82"

atom = 'Pb'
# atom = 'Au'
Atom = sig.Atom #dictionary where we stock important atomic aspect (now just Z and A) to plot and put into file names

Z = Atom[atom]['Z']
A = Atom[atom]['A']

PDF=[]
for prot in proton:
    PDF.append(sig.Sigma(prot,Pb,s,Z,A))

d = sig.Switch

#Kinematics
p_T = 5 #GeV
# p_T = 2 
mu_f_2 = p_T**2											

# Plot characteristics
f_size = 17
alph = 0.5
X = np.logspace(-5,-1,100)
def delta_ref(ref,data):
    ref = np.array(ref)
    data = np.array(data)
    return (data)/ref

colors = ['blue','orange','green']
# Gluon pdf

Tot_G_p = []
cen_G_p = []
for i,prot in enumerate(proton):
    sigma_p = PDF[i]
    p_set = sigma_p.p_set
    G_p=[]
    for num in range(p_set.size):
        if num == 0:
            cen_G_p.append([sigma_p.Gluon_p(x,mu_f_2,num) for x in X])
            continue
        g_p = [sigma_p.Gluon_p(x,mu_f_2,num) for x in X]
        G_p.append(g_p)
    (G_p_min, G_p_max)=sig.min_max(G_p)
    Tot_G_p.append((G_p_min, G_p_max))



fig, ax = plt.subplots()
ref  = cen_G_p[0]
for i,prot in enumerate(proton):
    cen =  cen_G_p[i]
    (p_min,p_max) = Tot_G_p[i]
    plt.fill_between(X,delta_ref(ref,p_min),delta_ref(ref,p_max),alpha =alph, color=colors[i], label= prot)
    # plt.plot(X,delta_ref(ref,cen))
plt.legend()
plt.xscale('log')
# plt.yscale('log')
# plt.title("gluon pdf set compared")
plt.text(0.8, 0.1, r'$Q^2 =$'+str(mu_f_2)+r' GeV$^2$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
plt.xlabel(r'x')
plt.ylabel(r'$G(x,Q^2)/G_{[ref]}(x,Q^2)$')
plt.savefig(os.path.join(plots_dir, 'gluon_proton_set_comparison_rs'+str(rs)+'_P_T'+str(p_T)+'.pdf'))
plt.show()
plt.close()

# F_2 pdf

Tot_F_2 = []
cen_F_2 = []
for i,prot in enumerate(proton):
    sigma_p = PDF[i]
    p_set = sigma_p.p_set
    F_2=[]
    for num in range(p_set.size):
        if num == 0:
            cen_F_2.append([sigma_p.F2_p(x,mu_f_2,num) for x in X])
            continue
        f_2 = [sigma_p.F2_p(x,mu_f_2,num) for x in X]
        F_2.append(f_2)
    (F_2_min, F_2_max)=sig.min_max(F_2)
    Tot_F_2.append((F_2_min, F_2_max))

fig, ax = plt.subplots()
ref  = cen_F_2[0]
for i,prot in enumerate(proton):
    cen =  cen_F_2[i]
    (p_min,p_max) = Tot_F_2[i]
    plt.fill_between(X,delta_ref(ref,p_min),delta_ref(ref,p_max),alpha =alph, color=colors[i], label= prot)
plt.legend()
plt.xscale('log') 
# plt.yscale('log')
# plt.title('F_2 pdf set compared')
plt.xlabel(r'x')
plt.ylabel(r'$F_2(x,Q^2)/F_{2,[ref]}(x,Q^2)$')
plt.text(0.8, 0.1, r'$Q^2 =$ '+str(mu_f_2)+' GeV$^2$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
plt.savefig(os.path.join(plots_dir, 'F_2_proton_set_comparison_rs'+str(rs)+'_P_T'+str(p_T)+'.pdf'))
plt.show()
plt.close()

# u+ubar pdf

Tot_U = []
cen_U = []
for i,prot in enumerate(proton):
    sigma_p = PDF[i]
    p_set = sigma_p.p_set
    U=[]
    for num in range(p_set.size):
        p_pdf = sigma_p.pdf_p(num)
        if num == 0:
            cen_U.append([p_pdf.xfxQ2(-sig.particle['u']['id'],x,mu_f_2)+p_pdf.xfxQ2(sig.particle['u']['id'],x,mu_f_2) for x in X])
            continue
        u = [p_pdf.xfxQ2(-sig.particle['u']['id'],x,mu_f_2)+p_pdf.xfxQ2(sig.particle['u']['id'],x,mu_f_2)for x in X]
        U.append(u)
    (U_min, U_max)=sig.min_max(U)
    Tot_U.append((U_min, U_max))

fig, ax = plt.subplots()
ref  = cen_U[0]
for i,prot in enumerate(proton):
    cen =  cen_U[i]
    (p_min,p_max) = Tot_U[i]
    plt.fill_between(X,delta_ref(ref,p_min),delta_ref(ref,p_max),alpha =alph, color=colors[i], label= prot)
plt.legend()
plt.xscale('log') 
# plt.yscale('log')
# plt.title('F_2 pdf set compared')
plt.xlabel(r'x')
plt.ylabel(r'$(u+\bar{u})(x,Q^2)/(u+\bar{u})_{[ref]}(x,Q^2)$')
plt.text(0.8, 0.1, r'$Q^2 =$ '+str(mu_f_2)+' GeV$^2$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
plt.savefig(os.path.join(plots_dir, 'u_u_bar_proton_set_comparison_rs'+str(rs)+'_P_T'+str(p_T)+'.pdf'))
plt.show()
plt.close()

# d+d bar pdf

Tot_D = []
cen_D = []
for i,prot in enumerate(proton):
    sigma_p = PDF[i]
    p_set = sigma_p.p_set
    D=[]
    for num in range(p_set.size):
        p_pdf = sigma_p.pdf_p(num)
        if num == 0:
            cen_D.append([p_pdf.xfxQ2(-sig.particle['d']['id'],x,mu_f_2)+p_pdf.xfxQ2(sig.particle['d']['id'],x,mu_f_2) for x in X])
            continue
        d = [p_pdf.xfxQ2(-sig.particle['d']['id'],x,mu_f_2)+p_pdf.xfxQ2(sig.particle['d']['id'],x,mu_f_2)for x in X]
        D.append(d)
    (D_min, D_max)=sig.min_max(D)
    Tot_D.append((D_min, D_max))

fig, ax = plt.subplots()
ref  = cen_D[0]
for i,prot in enumerate(proton):
    cen =  cen_D[i]
    (p_min,p_max) = Tot_D[i]
    plt.fill_between(X,delta_ref(ref,p_min),delta_ref(ref,p_max),alpha =alph, color=colors[i], label= prot)
plt.legend()
plt.xscale('log') 
# plt.yscale('log')
# plt.title('F_2 pdf set compared')
plt.xlabel(r'x')
plt.ylabel(r'$(d+\bar{d})(x,Q^2)/(d+ \bar{d})_{[ref]}(x,Q^2)$')
plt.text(0.8, 0.1, r'$Q^2 =$ '+str(mu_f_2)+' GeV$^2$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
plt.savefig(os.path.join(plots_dir, 'd_d_bar_proton_set_comparison_rs'+str(rs)+'_P_T'+str(p_T)+'.pdf'))
plt.show()
plt.close()

# s pdf

Tot_S = []
cen_S = []
for i,prot in enumerate(proton):
    sigma_p = PDF[i]
    p_set = sigma_p.p_set
    S = []
    for num in range(p_set.size):
        p_pdf = sigma_p.pdf_p(num)
        if num == 0:
            cen_S.append([p_pdf.xfxQ2(sig.particle['s']['id'],x,mu_f_2) for x in X])
            continue
        s = [p_pdf.xfxQ2(sig.particle['s']['id'],x,mu_f_2)for x in X]
        S.append(s)
    (S_min, S_max)=sig.min_max(S)
    Tot_S.append((S_min, S_max))

fig, ax = plt.subplots()
ref  = cen_S[0]
for i,prot in enumerate(proton):
    cen =  cen_S[i]
    (p_min,p_max) = Tot_S[i]
    plt.fill_between(X,delta_ref(ref,p_min),delta_ref(ref,p_max),alpha =alph, color=colors[i], label= prot)
plt.legend()
plt.xscale('log') 
# plt.yscale('log')
# plt.title('F_2 pdf set compared')
plt.xlabel(r'x')
plt.ylabel(r'$s(x,Q^2)/s_{[ref]}(x,Q^2)$')
plt.text(0.8, 0.1, r'$Q^2 =$ '+str(mu_f_2)+' GeV$^2$', horizontalalignment='center', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
plt.savefig(os.path.join(plots_dir, 's_proton_set_comparison_rs'+str(rs)+'_P_T'+str(p_T)+'.pdf'))
plt.show()
plt.close()