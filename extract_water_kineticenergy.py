import numpy
from numpy import sqrt, pi, exp, linspace, loadtxt, power, size, sqrt
import matplotlib.pyplot as plt 
import pylab as pl
#from lmfit import  Model
from scipy.optimize import curve_fit
k=0.00198657 #Boltzmann constant in kcal/mol

f = open ("kinetic_energy_water_s4.dat","w") #output file

mass = loadtxt('watermass_cada25frame_sel0_1_2_n.dat') #open mass file
velo = loadtxt('watervelo_cada25frame_sel0_1_2_n.dat') #open velocities file

max = size(mass)

#While loop, for computing the kinetic energy of the water center of mass
a=0
while a < max : 
	VCMx = mass[a]*velo[a,0] + mass[a+1]*velo[a+1,0] + mass[a+2]*velo[a+2,0]  
	VCMy = mass[a]*velo[a,1] + mass[a+1]*velo[a+1,1] + mass[a+2]*velo[a+2,1]  
	VCMz = mass[a]*velo[a,2] + mass[a+1]*velo[a+1,2] + mass[a+2]*velo[a+2,2]  
	MCM = mass[a] + mass[a+1] + mass[a+2]
	VCM = sqrt( VCMx*VCMx/MCM/MCM + VCMy*VCMy/MCM/MCM + VCMz*VCMz/MCM/MCM ) 
	kinetic = 0.5*MCM*VCM*VCM
	f.write("%s\n"%(str(kinetic)))
	a = a + 3
f.close()

