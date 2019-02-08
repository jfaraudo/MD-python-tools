# -------------------------------------------------------
# INTRODUCTION TO MOLECULAR DYNAMICS
# May 2018 
# 
# Very simple example of MD molecular dynamics
# Harmonic oscillator
# By Jordi Faraudo 2018
# -------------------------------------------------------

# Here we import the mathematical library and the plots library
import numpy as np
import matplotlib.pyplot as plt

# Particle submitted to harmonic force
# mass
m = 1
# Period in ns
T = 1.0
#frequency
w = 2.0*np.pi
#Force constant
k=m*w*w

#Initial condition (position and velocity)
x0=0.1
v0=0.0

#Show data of the program

print '\n------------------------------------------------------'
print 'SIMPLE MD SIMULATION OF A SINGLE ATOM IN HARMONIC TRAP'
print '\n------------------------------------------------------'
print 'Atom of mass:',m,' ng',' '
print 'Period:',T,' ns'

# input time step
dt = float(raw_input("\n Time step dt (in ns):\n>"))
# Final time
ntot = int(raw_input("\n Number of time steps:\n>"))
print 'Simulation time will be',dt*ntot,' ns'

# create empty time, position and velocity arrays (x,v,t) arrays
x = np.empty(shape=(0,ntot))
v = np.empty(shape=(0,ntot))
t = np.empty(shape=(0,ntot))

print x[1]

#Initial conditions
x[0] = x0
v[0] = v0
   
# Time evolution
print '\n Calculating time evolution...'
for i in range(0, ntot):
    print i
    #Calculate Force over the particle
    f = -k*x[i]
    #Calculate acceleration from 2nd Law
    a = f/m 
    # New velocity after time dt
    v[i] = v[i-1]+a*dt
    #Average velocity from t to t+dt
    v_av= (v[i]+v[i+1])/2.0
    # New position
    x[i+1] = x[i]+v_av*dt
    #Update time
    t[i+1] = t[i]+dt
    
print 'Calculation finished. Showing plot with results'

# plot output
plt.plot(t,x, 'ro', t, v, 'bv')
#create axis
plt.axhline(0, color='black')
plt.axvline(0, color='black')
#Show plot in screen
plt.show()
