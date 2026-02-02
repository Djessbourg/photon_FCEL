import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from src import sigma as sig
from src import Collision
from scipy.integrate import quad

### Initialisation of a sigma object ###
rs = 8800.
s = (rs)**2 															# CM energy in Gev2
p = "NNPDF40_lo_as_01180"
Pb = "nNNPDF30_nlo_as_0118_A208_Z82"
Z = 82.
A = 208.

pPb_cross_section = sig.Sigma(p,Pb,s,Z,A)

n_p = len(pPb_cross_section.p_list)
n_A = len(pPb_cross_section.A_list)

# By taking G_A(num=0)/G_p(num=0) as the central value.
Q2 = 10 #GeV^2 
X = np.logspace(-3,0,1000)
G_p = []
G_A = []

for i in range(n_p):
    Gp = []
    for x in X:
        Gp.append(pPb_cross_section.Gluon_p(x,Q2,i))
    G_p.append(Gp)

for i in range(n_A):
    GA = []
    for x in X:
        GA.append(pPb_cross_section.Gluon_A(x,Q2,i))
    G_A.append(GA)