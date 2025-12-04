# -*- coding: utf-8 -*-


###--------------------------###
# Do a plot of proton pdf 
###--------------------------###
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
#demander quand meme qu'on installe Tkinter pour l'interface
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator
import lhapdf as lha

### Initialisation ###
p = "nNNPDF30_nlo_as_0118_p"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"
Z = 82.
A = 208.

p_pdf = lha.mkPDF(p,0)
A_pdf = lha.mkPDF(Pb,0)

Q = 91 # scale in GeV , 91 is Z mass

#proton

d = lambda x: p_pdf.xfxQ(1,x,Q)
db = lambda x: p_pdf.xfxQ(-1,x,Q)
u = lambda x: p_pdf.xfxQ(2,x,Q)
ub = lambda x: p_pdf.xfxQ(-2,x,Q)
s = lambda x: p_pdf.xfxQ(3,x,Q)
sb = lambda x: p_pdf.xfxQ(-3,x,Q)
g = lambda x: p_pdf.xfxQ(21,x,Q)/4
u_val = lambda x: u(x)-ub(x)
F_2p = lambda x: ((2./3.)**2)*(u(x)+ub(x))+((1./3.)**2)*(d(x)+db(x)+s(x)+sb(x))

#nuclei

dA = lambda x: A_pdf.xfxQ(1,x,Q)
dbA = lambda x: A_pdf.xfxQ(-1,x,Q)
uA = lambda x: A_pdf.xfxQ(2,x,Q)
ubA = lambda x: A_pdf.xfxQ(-2,x,Q)
sA = lambda x: A_pdf.xfxQ(3,x,Q)
sbA = lambda x: A_pdf.xfxQ(-3,x,Q)
gA = lambda x: A_pdf.xfxQ(21,x,Q)/4
u_Aval = lambda x: uA(x)-ubA(x)
F_2A = lambda x: ((2./3.)**2)*(uA(x)+ubA(x))+((1./3.)**2)*(dA(x)+dbA(x)+s(x)+sbA(x))

F_2pA = lambda x: F_2A(x)/(F_2p(x))
X = np.logspace(-5,0,100)
g_pA = lambda x: gA(x)/g(x)
u_pA = lambda x: (u_Aval(x))/(u_val(x))
d_pA = lambda x: (dA(x)+dbA(x))/(d(x)+db(x))
# plot
f_size = 17

def plot_usuals(n=1,s1=f_size,s2=f_size,loca = 'best'):
	plt.legend(frameon= False, fontsize = s1,ncols=n,loc=loca )
	plt.tick_params(labelsize=s2)

fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(LogLocator(subs='all'))


plt.plot(X,[d(x)-db(x) for x in X],label=r'$xd_{\text{val}}$')
# plt.plot(X,[db(x) for x in X],label=r'$x\bar{d}_p(x)$')
plt.plot(X,[u(x)-ub(x) for x in X],label=r'$xu_{\text{val}}$')
# plt.plot(X,[ub(x) for x in X],label=r'$x\bar{u}_p(x)$')
# plt.plot(X,[s(x) for x in X],label=r'$xs_p(x)$')
# plt.plot(X,[sb(x) for x in X],label=r'$x\bar{s}_p(x)$')
plt.plot(X,[g(x) for x in X],label=r'$xg(\times 1/4)$')
plt.plot(X,[F_2p(x) for x in X], label = r'$F_{2,p}$')
plt.plot(X,[2*(ub(x)+db(x))+s(x)+sb(x) for x in X],label=r'$xq_{\text{sea}}$')
# plt.plot(X,[F_2p_v(x) for x in X], label = r'$F_2(x)$ valence')
plt.xscale('log')
plt.yscale('log')
plt.axvline(x=1/3,color='grey',alpha=0.3)
plt.ylim(10**(-5),50)
# plt.xlim(top=1)
plot_usuals(n=2)
plt.xlabel(r'$x$',fontsize=f_size)
plt.ylabel(r'$xf_p(x,Q)$',fontsize=f_size)
plt.text(0.1, 0.75,r'$Q =$ '+str(Q)+'GeV ',horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
plt.text(0.1,0.65,p,horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
n_fig = 'proton_pdf_Q.pdf'
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
plt.show()
plt.close()

fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(LogLocator(subs='all'))
plt.plot(X,[F_2pA(x)for x in X],label=r'$F_{2,pA}$')
plt.axhline(y=1,color='grey',alpha=0.3)
plt.xscale('log')
# plt.yscale('log')
plt.ylim(0.5,1.5)
# plt.xlim(top=1)
plot_usuals(n=1)
plt.xlabel(r'$x$',fontsize=f_size)
# plt.ylabel(r'$xf_p(x,Q)$',fontsize=f_size)
plt.text(0.1, 0.85,r'$Q =$ '+str(Q)+'GeV ',horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
# plt.text(0.1,0.65,p,horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
n_fig = 'F_2_pA_pdf_Q.pdf'
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
plt.show()
plt.close()

# try with the errors

fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_minor_locator(LogLocator(subs='all'))
plt.plot(X,[g_pA(x)for x in X],label=r'$g_{pA}$')
plt.plot(X,[u_pA(x)for x in X],label=r'$u_{pA}$')
# plt.plot(X,[d_pA(x)for x in X],label=r'$d_{pA}$')
plt.axhline(y=1,color='grey',alpha=0.3)
plt.xscale('log')
# plt.yscale('log')
plt.ylim(0,3)
# plt.xlim(top=1)
plot_usuals(n=1)
plt.xlabel(r'$x$',fontsize=f_size)
# plt.ylabel(r'$xf_p(x,Q)$',fontsize=f_size)
plt.text(0.1, 0.85,r'$Q =$ '+str(Q)+'GeV ',horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
# plt.text(0.1,0.65,p,horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,fontsize = f_size)
n_fig = 'g_pA_pdf_Q.pdf'
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.savefig(os.path.join(plots_dir,n_fig),bbox_inches="tight")# bbox_inches="tight"
plt.show()
plt.close()