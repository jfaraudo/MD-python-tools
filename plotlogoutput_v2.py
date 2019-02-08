import numpy
import sys #library for making script with input parameters
from numpy import sqrt, pi, exp, linspace, loadtxt, power
import matplotlib.pyplot as plt 
import pylab as pl
k = 0.00198657 #Boltzmann constant in kcal/mol

nu = len(sys.argv) #number of arguments
lista = sys.argv

f = open('averageproperties.txt',"w")

def plotwrite(prop,units):
	array = loadtxt('%s.dat'%(str(prop)))
	array = array[:,1] #First column only contains TS
	avg=numpy.median(array)
	dev=numpy.std(array)
	print('%s %s+-%s %s\n'%(str(prop),str(avg),str(dev),str(units)))
	f.write('%s %s+-%s %s\n'%(str(prop),str(avg),str(dev),str(units)))
	pl.figure()
	pl.ylabel('%s (%s)'%(str(prop),str(units))), pl.xlabel('TS')
	pl.plot(array,label='Average: %s+-%s%s'%(str(avg),str(dev),str(units)))
	pl.legend(), pl.savefig('%s.png'%(str(prop)))

#"-h parameter: HELP (information for callable properties)
for a in range(nu):
	if lista[a] == '-h':
		print('########NAMD LOG FILE Output generator########')
		print('Script prints average value, and plots the property\n\n')	
		print('Callable variables:')
		print('p: PRESSURE\n')
		print('to: TOTAL\n')
		print('gp: GPRESSURE\n')
		print('pa: PRESSAVG\n')
		print('gpa: GPRESSAVG\n')
		print('k: KINETIC\n')
		print('t: TEMP\n')


print('##Average properties and STD will be printed:##\n')

for a in range(nu):
	##PRESSURE
	if lista[a] == 'p':
		plotwrite('pressure','Bar')
	##TOTAL
	if lista[a] == 'to':
		plotwrite('total','kcal/mol')
	##GPRESSURE
	if lista[a] == 'gp':
		plotwrite('gpressure','Bar')
	##VOLUME
	if lista[a] == 'v':
		plotwrite('volume','A^3')
	##PRESSAVG
	if lista[a] == 'pa':
		plotwrite('pressavg','Bar')
	##GPRESSAVG
	if lista[a] == 'gpa':
		plotwrite('gpressavg','Bar')
	##KINETIC
	if lista[a] == 'k':
		plotwrite('kinetic','kcal/mol')
	##TEMP
	if lista[a] == 't':
		plotwrite('temp','K')
f.close()
