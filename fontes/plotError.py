#!/usr/bin/python
# -*- coding: utf-8 -*-

# Post-processing program to plot results from linear elasticity problem
# Author: Diego Volpatto

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.mlab import griddata
from matplotlib.collections import EllipseCollection
import sys

def lamb(E, ni):
  return (ni*E)/((1.0 + ni)*(1.0-2.0*ni))

def mi(E, ni):
  return (E)/(2.0*(1.0 + ni))

def solU_grav(z,L,rho,cgrav,lamb,mi):
  return ((rho*cgrav)/(2.0*(lamb + 2.0*mi)))*(L**2.0 - z**2.0)

def solU_strain(z,t,lamb,mi):
  return (t/(lamb + 2.0*mi))*z

E = 40.0e9
ni = 0.23
t = -3.40e6/(10.0/200.0)
lambU = lamb(E,ni)
miU = mi(E,ni)

filepath = sys.argv[1]

filename = './' + filepath + 'passosPressaoBlocoMacro.dat'

inDataLegendP = np.loadtxt(filename,unpack=True)
dataLegendP = inDataLegendP[1:]
pnum = len(dataLegendP)

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameErrorL = './' + filepath + ('eqdivL.%d' % i)
	Y, Ex, Ey = np.loadtxt(inDataNameErrorL,unpack=True,usecols=[2,3,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_{E_y}}$',fontsize=18)
	#ax.set_ylabel(r'$\sigma_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Ey,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'eqdivLy.png')

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1
	
	inDataNameErrorC = './' + filepath + ('eqdivC.%d' % i)
	Y, Ex, Ey = np.loadtxt(inDataNameErrorC,unpack=True,usecols=[2,3,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	#ax.set_ylabel(r'$U_y$',fontsize=18)
	ax.set_ylabel(r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_{E_y}}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Ey,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'eqdivCy.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1
	
	inDataNameErrorR = './' + filepath + ('eqdivR.%d' % i)
	Y, Ex, Ey = np.loadtxt(inDataNameErrorR,unpack=True,usecols=[2,3,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	#ax.set_ylabel(r'$U_y$',fontsize=18)
	ax.set_ylabel(r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_{E_y}}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Ey,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'eqdivRy.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameErrorL = './' + filepath + ('eqdivL.%d' % i)
	Y, Ex, Ey = np.loadtxt(inDataNameErrorL,unpack=True,usecols=[2,3,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\frac{\partial_x \sigma_E - \alpha \nabla p}{\sigma_{E_x}}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Ex,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'eqdivLx.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1
	
	inDataNameErrorC = './' + filepath + ('eqdivC.%d' % i)
	Y, Ex, Ey = np.loadtxt(inDataNameErrorC,unpack=True,usecols=[2,3,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	#ax.set_ylabel(r'$U_x$',fontsize=18)
	ax.set_ylabel(r'$\frac{\partial_x \sigma_E - \alpha \nabla p}{\sigma_{E_x}}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Ex,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'eqdivCx.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1
	
	inDataNameErrorR = './' + filepath + ('eqdivR.%d' % i)
	Y, Ex, Ey = np.loadtxt(inDataNameErrorR,unpack=True,usecols=[2,3,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	#ax.set_ylabel(r'$U_x$',fontsize=18)
	ax.set_ylabel(r'$\frac{\partial_x \sigma_E - \alpha \nabla p}{\sigma_{E_x}}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Ex,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'eqdivRx.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameUL = './' + filepath + ('UL.%d' % i)
	Y, Uy = np.loadtxt(inDataNameUL,unpack=True,usecols=[2,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$U_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Uy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
solU_disp = solU_strain(Y,t,lambU,miU)
#plt.plot(Y,solU_disp,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'ULy.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameUC = './' + filepath + ('UC.%d' % i)
	Y, Uy = np.loadtxt(inDataNameUC,unpack=True,usecols=[2,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$U_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Uy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')

#solU_disp = solU_strain(Y,t,lambU,miU)
#plt.plot(Y,solU_disp,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'UCy.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameUR = './' + filepath + ('UR.%d' % i)
	Y, Uy = np.loadtxt(inDataNameUR,unpack=True,usecols=[2,4])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$U_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(Y,Uy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

#ax.legend(loc='best', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
solU_disp = solU_strain(Y,t,lambU,miU)
#plt.plot(Y,solU_disp,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'URy.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameUCx = './' + filepath + ('UCx.%d' % i)
	X, Uy = np.loadtxt(inDataNameUCx,unpack=True,usecols=[1,4])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$U_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(X,Uy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.70, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#solU_disp = solU_strain(Y,t,lambU,miU)
#plt.plot(Y,solU_disp,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'UCx.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameUT = './' + filepath + ('UTx.%d' % i)
	X, Uy = np.loadtxt(inDataNameUT,unpack=True,usecols=[1,4])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$U_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(X,Uy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.70, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
solU_disp = solU_strain(Y,t,lambU,miU)
#plt.plot(Y,solU_disp,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'UT.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

i = 1
inDataNameSigyL = './' + filepath + ('sigmasL.%d' % i)
ys = np.loadtxt(inDataNameSigyL,unpack=True,usecols=[0])
numElem = len(ys)
Ly = 50.0
Lx = 10.0
y = np.linspace(0.0, Ly, numElem)
x = np.linspace(0.0, Lx, numElem)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameSigyL = './' + filepath + ('sigmasL.%d' % i)
	sigy = np.loadtxt(inDataNameSigyL,unpack=True,usecols=[2])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\sigma_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(y,sigy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'sigyL.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameSigyC = './' + filepath + ('sigmasC.%d' % i)
	sigy = np.loadtxt(inDataNameSigyC,unpack=True,usecols=[2])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\sigma_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(y,sigy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
#plt.plot(y,ty,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'sigyC.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameSigyR = './' + filepath + ('sigmasR.%d' % i)
	sigy = np.loadtxt(inDataNameSigyR,unpack=True,usecols=[2])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\sigma_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(y,sigy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
#plt.plot(y,ty,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'sigyR.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameTauyL = './' + filepath + ('sigmasL.%d' % i)
	tauy = np.loadtxt(inDataNameTauyL,unpack=True,usecols=[3])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\tau_{xy}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(y,tauy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'tauyL.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameTauyC = './' + filepath + ('sigmasC.%d' % i)
	tauy = np.loadtxt(inDataNameTauyC,unpack=True,usecols=[3])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\tau_{xy}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(y,tauy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
#plt.plot(y,ty,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'tauyC.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameTauyR = './' + filepath + ('sigmasR.%d' % i)
	tauy = np.loadtxt(inDataNameTauyR,unpack=True,usecols=[3])
	
	ax.set_xlabel(r'$y\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\tau_{xy}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(y,tauy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.3*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
#plt.plot(y,ty,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'tauyR.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

i = 1
inDataNameSigyB = './' + filepath + ('sigmasLx.%d' % i)
xs = np.loadtxt(inDataNameSigyB,unpack=True,usecols=[0])
numElem = len(xs)
Lx = 10.0
x = np.linspace(0.0, Lx, numElem)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameSigyL = './' + filepath + ('sigmasLx.%d' % i)
	sigyB = np.loadtxt(inDataNameSigyL,unpack=True,usecols=[2])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\sigma_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(x,sigyB,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'sigyB.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

#i = 1
#inDataNameSigxC = './' + filepath + ('sigmasCx.%d' % i)
#xs = np.loadtxt(inDataNameSigxC,unpack=True,usecols=[0])
#numElem = len(xs)
#Lx = 10.0
#x = np.linspace(0.0, Lx, numElem)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameSigyC = './' + filepath + ('sigmasCx.%d' % i)
	sigyC = np.loadtxt(inDataNameSigyC,unpack=True,usecols=[2])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\sigma_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(x,sigyC,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'sigyCx.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameSigyT = './' + filepath + ('sigmasRx.%d' % i)
	sigyT = np.loadtxt(inDataNameSigyT,unpack=True,usecols=[2])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\sigma_y$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(x,sigyT,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'sigyT.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameTauxL = './' + filepath + ('sigmasLx.%d' % i)
	tauxB = np.loadtxt(inDataNameTauxL,unpack=True,usecols=[3])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\tau_{xy}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(x,tauxB,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'tauxB.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

#i = 1
#inDataNameSigxC = './' + filepath + ('sigmasCx.%d' % i)
#xs = np.loadtxt(inDataNameSigxC,unpack=True,usecols=[0])
#numElem = len(xs)
#Lx = 10.0
#x = np.linspace(0.0, Lx, numElem)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameTauxC = './' + filepath + ('sigmasCx.%d' % i)
	tauxC = np.loadtxt(inDataNameTauxC,unpack=True,usecols=[3])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\tau_{xy}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(x,tauxC,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'tauxCx.png')
plt.close()

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameTauxT = './' + filepath + ('sigmasRx.%d' % i)
	tauxT = np.loadtxt(inDataNameTauxT,unpack=True,usecols=[3])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\tau_{xy}$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(x,tauxT,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.25*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
plt.grid(True)
#plt.plot(y,ty,color='green')
plt.savefig('./' + filepath + 'tauxT.png')
plt.close()

i = 1
inDataNameSigy = './' + filepath + ('sigy.%d' % i)
ys = np.loadtxt(inDataNameSigy,unpack=True,usecols=[0])
numElem = len(ys)
Lx = 10.0
x = np.linspace(0.0, Lx, numElem)

fig = plt.figure()
ax = plt.subplot(111)

for i in range(1,pnum+1):
	#i = 1	
	
	inDataNameSigy = './' + filepath + ('sigy.%d' % i)
	sigy = np.loadtxt(inDataNameSigy,unpack=True,usecols=[1])
	
	ax.set_xlabel(r'$x\,(m)$',fontsize=16)
	ax.set_ylabel(r'$\sigma_y - \alpha p$',fontsize=18)
	if dataLegendP[i-1] == 1.:
		leg = (u"%d mês" % dataLegendP[i-1])
	else:
		leg = ("%d meses" % dataLegendP[i-1])
	ax.plot(x,sigy,'o-',label=leg)#, label=r'$\frac{\partial_y \sigma_E - \alpha \nabla p}{\sigma_E_y}$')

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ty = t*np.ones(numElem)
#plt.plot(x,ty,color='green')
plt.grid(True)
plt.savefig('./' + filepath + 'sigyBC.png')
plt.close()
print "Error's plots OK"