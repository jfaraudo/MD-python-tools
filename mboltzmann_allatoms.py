import numpy
from numpy import sqrt, pi, exp, linspace, loadtxt, power
import matplotlib.pyplot as plt 
import pylab as pl
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit
k=0.00198657 #Boltzmann constant in kcal/mol

#Data loading (kinetic energies from the MD simulation)
rawdata = loadtxt('kineticall_cada25frames_sel0_1_2.dat')
r1 = loadtxt('kineticall_cada25frames_sel3_4_5.dat')
rawdata = numpy.concatenate((rawdata,r1))
del r1
r1 = loadtxt('kineticall_cada25frames_sel6_7_8.dat')
rawdata = numpy.concatenate((rawdata,r1))
del r1
r1 = loadtxt('kineticall_cada25frames_sel9_10_11.dat')
rawdata = numpy.concatenate((rawdata,r1))
del r1
r1 = loadtxt('kineticall_cada25frames_sel12_13_14.dat')
rawdata = numpy.concatenate((rawdata,r1))
del r1
r1 = loadtxt('kineticall_cada25frames_sel15_16_17.dat')
rawdata = numpy.concatenate((rawdata,r1))
del r1
r1 = loadtxt('kineticall_cada25frames_sel18_19_n.dat')
rawdata = numpy.concatenate((rawdata,r1))
del r1
r1 = loadtxt('kineticall_cada25frames_sel20_21_n.dat')
rawdata = numpy.concatenate((rawdata,r1))
del r1

#Create a normalised histogram (normalised area) with 300 bins
pl.figure(figsize=(3,2.8))
freq, energies, patches = pl.hist(rawdata,bins=300,range=(0,4),normed=True,histtype='bar',align='mid',color='#0072B2',linewidth=0.2,alpha=0.8) 

freq=numpy.append(0,freq)
s=numpy.size(freq)

#Plot the theoretical Maxwell-Boltzmann distribution at 300 K
def maxwell(e,T):
	return 2/sqrt(pi)*power(k*T,-3./2.)*sqrt(e)*exp(-e/k/T)
pl.plot(energies,maxwell(energies,300),'-k',linewidth=1.3,label='Theoretical distribution \nat 300K')

#Fit the kinetic energies histogram to the M-B distribution at the temperature which best reproduces the M-B distribution
#popt, pcov = curve_fit(maxwell, energies, freq ,p0=200)
# The popt argument are the best-fit paramters for T
#The pcov variable contains the covariance matrix, which indicates the uncertainties and correlations between parameters. This is mostly useful when the data has uncertainties.


pl.legend(fontsize=8,frameon=False), pl.axis('tight')
pl.xlabel("Kinetic Energy (kcal/mol)",fontsize=8), pl.ylabel("Normalised frequency",fontsize=8), pl.xlim(0,4), pl.ylim(0)
ax = plt.gca() #get current axis
ax.xaxis.set_major_locator( MaxNLocator(nbins=5) ) #place 5 bins in the x axis
pl.tick_params(axis='both',labelsize=7) #change x and y ticks size
plt.tight_layout() 
pl.savefig('histogramboltzmann_allatoms_all_s4',dpi=1200)
