#import python modules
import numpy as np
import time

#import compiled fortran library
import ljlib

#import pdb writing
import atomwrite


#whether or not to save the trajectory as a PdbFile
SaveTrajectory = True

#distance cutoff for pairwise interactions
Cut = 3.0


#NOTE:
#everything below assumes unit atomic masses,
#such that forces = accelerations.


def InitPositions(N, L):
    """Returns an array of initial positions of each atom,
placed on a cubic lattice for convenience.
Input:
    N: number of atoms
    L: box length
Output:
    Pos: (N,3) array of positions
"""
    #make the position array
    Pos = np.zeros((N,3), float)
    #compute integer grid # of locations for cubic lattice
    NLat = int(N**(1./3.) + 1.)
    #make an array of lattice sites
    r = L * (np.arange(NLat, dtype=float)/NLat - 0.5)
    #loop through x, y, z positions in lattice until done
    #for every atom in the system
    i = 0
    for x in r:
        for y in r:
            for z in r:
                Pos[i] = np.array([x,y,z], float)
                i += 1
                #if done placing atoms, return
                if i >= N:
                    return Pos
    return Pos


def RescaleVelocities(Vel, T):
    """Rescales velocities in the system to the target temperature.
Input:
    Vel: (N,3) array of atomic velocities
    T: target temperature
Output:
    Vel: same as above
"""
    #recenter to zero net momentum (assuming all masses same)
    Vel = Vel - Vel.mean(axis=0)
    #find the total kinetic energy
    KE = 0.5 * np.sum(Vel * Vel)
    #find velocity scale factor from ratios of kinetic energy
    VScale = np.sqrt(1.5 * len(Vel) * T / KE)
    Vel = Vel * VScale
    return Vel  


def InitVelocities(N, T):
    """Returns an initial random velocity set.
Input:
    N: number of atoms
    T: target temperature
Output:
    Vel: (N,3) array of atomic velocities
"""
    Vel = np.random.rand(N, 3)
    Vel = RescaleVelocities(Vel, T)
    return Vel


def InitAccel(Pos, L):
    """Returns the initial acceleration array.
Input:
    Pos: (N,3) array of atomic positions
    L: simulation box length
Output:
    Accel: (N,3) array of acceleration vectors
"""
    Accel = np.zeros_like(Pos)
    #get the acceleration from the forces
    PEnergy, Accel = ljlib.energyforces(Pos, L, Cut, Accel)
    return Accel  


def RunTest():
    #set the init box width, number of particles, temperature, and timestep
    N = 108
    rho = 0.80
    L = (N / rho)**(1./3.)
    Temp = 1.0
    dt = 0.001

    #set the frequency in seconds to update the display    
    DisplayFreq = 0.1
    #set the frequency in md steps to rescale velocities
    RescaleSteps = 1000
    #set the frequency in md steps to write coordinates
    WriteSteps = 100

    #set the max number of md steps; 0 for infinite loop
    MaxSteps = WriteSteps*1000

    #set the random number seed; useful for debugging
    np.random.seed = 342324

    #get the initial positions, velocities, and acceleration (forces)    
    Pos = InitPositions(N, L)
    Vel = InitVelocities(N, Temp)
    Accel = InitAccel(Pos, L)

    #MD steps
    StartTime = time.time()
    LastTime = time.time()
    i = 0
    if SaveTrajectory:
        Pdb = atomwrite.pdbfile("anim.pdb", L)
        
    while i < MaxSteps or MaxSteps <= 0:
        #do one step of the integration by calling the Fortran libraries
        Pos, Vel, Accel, KEnergy, PEnergy = ljlib.vvintegrate(Pos, Vel, Accel, L, Cut, dt)
        i += 1
        
        #check if we need to rescale the velocities 
        if i % RescaleSteps == 0:
            Vel = RescaleVelocities(Vel, Temp)

        #check if we need to output the positions 
        if SaveTrajectory and WriteSteps > 0 and i % WriteSteps == 0:
            Pdb.write(Pos)            

        #check if we need to update the display            
        if time.time() - LastTime > DisplayFreq:
            print "%d  %11.4f  %11.4f  %11.4f" % (i, PEnergy + KEnergy, PEnergy, KEnergy)
            LastTime = time.time()

    if SaveTrajectory:
        Pdb.close()            

    #do one last update of the display         
    print "%d  %11.4f  %11.4f  %11.4f" % (i, PEnergy + KEnergy, PEnergy, KEnergy)

    StopTime = time.time()
    print "Total time: %.1f s" % (StopTime - StartTime)


#check to see if we were run at the command line
if __name__ == '__main__':
    #run the test simulation
    RunTest()

    