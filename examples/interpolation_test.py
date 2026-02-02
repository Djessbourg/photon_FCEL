# =============================================================================
# File for testing interpolation methods 
# Last modified: 4/11/2025
# =============================================================================
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'plots'))
examples_dir =  os.path.abspath(os.path.dirname(__file__))						# current directory
pawres_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'copy_pawres')) # put here your jetphox data
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'RpA_tot_dir_onef')) # the direcory to save data from this file

import numpy as np
import matplotlib.pyplot as plt
from src import sigma as sig
from src import Probability as p
from scipy.integrate import quad

import uproot
from scipy.interpolate import CubicSpline
from scipy.interpolate import make_smoothing_spline

lam = None

process ={
	"1":{"rho":{"3bar":2./3.,"6":1./3.},"Fc":{"3bar":4./3.,"6":2.}},
	"3":{"rho":{"1":8./9.,"8":1./9.},"Fc":{"1":0.,"8":3.}},
	"5":{"rho":{"3bar":4./22.,"6":18./22.},"Fc":{"3bar":4./3.,"6":2.}},
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

file = uproot.open(os.path.join(pawres_dir,"LO_frag_ratios.root"))

# retrieve all fragmentaion histograms
histo = {k: file[k] for k in file.keys() if file[k].classname.startswith("TH1")} #labels are the same as the process dict but with a ';1' in addition
f_alpha_cs = {}
f_alpha_smooth = {}

def prefix(string):																# useful for normalizing root histograms
	'''return the prefix of a string with a ; separator'''
	s = ''
	for i in string:
		if i ==';':
			break
		s = s + i
	return s

# Create a dictionary of cubic sline fits of the f_alpha
for k in histo.keys():
	values = histo[k].values()
	edges = histo[k].axes[0].edges()
	centers = (edges[:-1] + edges[1:]) / 2.
	cs = CubicSpline(centers, values)
	smooth = make_smoothing_spline(centers, values,lam=lam)
	f_alpha_cs[prefix(k)]=cs 																# cs gives an array for a float or an array of float given
	f_alpha_smooth[prefix(k)]=smooth

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

Processus_cs = {} 																	# dict for all cross sections interpolated from jetphox data
Processus_smooth = {}
Processus_smooth2 = {}

for p_num in list(process.keys()):
	fn = 'ggoLOdiag_'+p_num+'.root'
	file = uproot.open(os.path.join(pawres_dir,fn))
	histo = file[fname]
	Values = histo.values()
	edges = histo.axes[0].edges()
	rapidity = (edges[:-1] + edges[1:]) / 2.
	cs = CubicSpline(rapidity,Values)
	smooth = make_smoothing_spline(rapidity, Values,lam=lam)
	lil_Y = np.linspace(-6, 6,31)
	smooth2 = make_smoothing_spline(lil_Y, smooth(lil_Y),lam=lam)
	Processus_cs[p_num] = cs
	Processus_smooth[p_num]= smooth
	Processus_smooth2[p_num] = smooth2

p_num = '14a'

def plot_usuals():
	plt.legend(frameon= False, fontsize = 15)
	plt.tick_params(labelsize=15)


Y = np.linspace(-6, 6,100)
# fig, ax = plt.subplots(constrained_layout=True)
# plt.plot(Y,f_alpha_cs[p_num](Y),label='cs')
# plt.plot(Y,f_alpha_smooth[p_num](Y),label='smooth')
# plot_usuals()
# plt.title(p_num,fontsize=15)
# plt.xlabel('y',fontsize=15)
# plt.ylabel(r'$R_{pA}^{frag,R_\alpha}$ for $\alpha=$ '+p_num,fontsize=15)
# plt.show()
# plt.close()

fig, ax = plt.subplots(constrained_layout=True)
plt.plot(Y,Processus_cs[p_num](Y),label='cs')
plt.plot(Y,Processus_smooth[p_num](Y),label='smooth')
plt.plot(Y,Processus_smooth2[p_num](Y),label='smooth2')
plot_usuals()
plt.title(p_num)
plt.show()