#!/usr/bin/python
# -*- coding: utf-8 -*-

# Post-processing program to plot results from linear elasticity problem
# Author: Diego Volpatto

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.mlab import griddata
import sys

filepath = sys.argv[1]

def lamb(E, ni):
  return (ni*E)/((1.0 + ni)*(1.0-2.0*ni))

def mi(E, ni):
  return (E)/(2.0*(1.0 + ni))

def solU_grav(z,L,rho,cgrav,lamb,mi):
  return ((rho*cgrav)/(2.0*(lamb + 2.0*mi)))*(L**2.0 - z**2.0)

def solU_strain(z,t,lamb,mi):
  return (t/(lamb + 2.0*mi))*z

def solSigma_grav(z,rho,cgrav):
  return -(rho*cgrav*z)

def solSigma_strain(z,t):
  return t + 0.0*z

#inDataYoung = 'youngdata.dat'
#E = np.loadtxt(inDataYoung)
E = 1.0e4
#inDataPoisson = 'poissondata.dat'
#ni = np.loadtxt(inDataPoisson)
ni = 0.0
#inDataRho = 'rhodata.dat'
#rho = np.loadtxt(inDataRho)
rho = 1.0
clamb = lamb(E,ni)
cmi = mi(E,ni)
cgrav = -9.8

fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)
#fig, ax = plt.subplots()

inDataNameDisp = './' + filepath + 'disp.1'
X, Y, U, V = np.loadtxt(inDataNameDisp,unpack=True,usecols=[1,2,3,4])
#inDataNameNodes = './dxgeo01/fnodes.stoc'
#U, V = np.loadtxt(inDataNameNodes,unpack=True,usecols=[3,4])
#skip = (slice(None, None, 3), slice(None, None, 3))
UU, VV = np.meshgrid(U,V)
XX, YY = np.meshgrid(X,Y)
#disp = np.sqrt(UU**2. + VV**2.)
disp = np.sqrt(U**2. + V**2.)
#print len(disp)
#im = ax.imshow(disp, extent=[X.min(), X.max(), Y.min(), Y.max()])
#ax.set_xlabel(r'$x\,(m)$',fontsize=18)
plt.xlabel(r'$x$',fontsize=18)
#ax.set_ylabel(r'$y\,(m)$',fontsize=18)
plt.ylabel(r'$y$',fontsize=18)
#ax.quiver(X,Y,U,V, disp, cmap='Blues' ,angles='xy', scale_units='xy', scale=1.2)
#ax.quiver(X,Y,U,V, disp, cmap='Blues' ,angles='xy', scale_units='xy')
plt.quiver(X,Y,U,V, disp,cmap=cm.afmhot_r)
#plt.imshow(disp)
#plt.plot(XX,YY)
plt.colorbar()
#plt.streamplot(X,Y,U,V)

#print X[:51]

#print U[0:51]

#Y, X = np.mgrid[-3:3:100j, -3:3:100j]
#U = -1 - X**2 + Y
#V = 1 + X - Y**2
#speed = np.sqrt(U*U + V*V)

##speed = np.sqrt(U*U + V*V)

#plt.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.autumn)
#plt.colorbar()

#f, (ax1, ax2) = plt.subplots(ncols=2)
#ax1.streamplot(X, Y, U, V, density=[0.5, 1])

#lw = 5*speed/speed.max()
#ax2.streamplot(X, Y, U, V, density=0.6, color='k', linewidth=lw)

#plt.show()




#plt.colorbar()
#plt.pcolor(X, Y, disp, cmap='RdBu', vmin=disp_min, vmax=disp_max)
#plt.pcolor(X,Y,disp)
#ax.set(aspect=1, title='Quiver Plot')
#ax.streamplot(X, Y, U, V, disp)
#fig.colorbar(im)

#box = ax.get_position()

#plt.colorbar(disp)

#ax.set_position([0.1*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 1.0, box.height])

#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('./'+ filepath + 'quiverUgrav.pdf')
#plt.show()

# Deslocamentos em y

fig = plt.figure()
ax = plt.subplot(111)
ax.set_xlabel(r'$y\,(m)$',fontsize=18)
ax.set_ylabel(r'$U_y$',fontsize=18)
L = Y.max()
t = 1.0e3
#solUg = solU_grav(Y,L,rho,cgrav,clamb,cmi)
#solUt = -1.*solU_strain(Y,t,clamb,cmi)
ax.plot(Y,V,'o', label=u'Sol Numérica')
#ax.plot(Y,solUt,'--', label=u'Sol Analítica',color='r')
box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 1.0, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.savefig('./' + filepath + 'Uy.pdf')
#plt.show()

XX, YY = np.meshgrid(X,Y)
UU,VV = np.meshgrid(U,V)

#print XX

#plt.streamplot(XX, YY, U, V, color=U, linewidth=2, cmap=plt.cm.autumn)
#plt.colorbar()
#plt.plot(XX,YY,'.',zorder=1)
#plt.scatter(XX,YY,zorder=2)

#f, (ax1, ax2) = plt.subplots(ncols=2)
#ax1.streamplot(X, Y, U, V, density=[0.5, 1])

#lw = 5*speed/speed.max()
#ax2.streamplot(X, Y, U, V, density=0.6, color='k', linewidth=lw)

#plt.show()

# Tensões

#fig = plt.figure()
#ax = plt.subplot(111)
#inDataSigy = './dxgeo01/sigy000.stoc'
#sigy = np.loadtxt(inDataSigy,unpack=True)
#it = len(sigy)
#zSig = np.linspace(Y.min(),Y.max(),it)
##print "len zSig = ", len(zSig)
##h = np.array([])
##for i in range(len(sigy)):
  ##h = np.append(h,())
#ax.set_xlabel(r'$y\,(m)$',fontsize=18)
#ax.set_ylabel(r'$\sigma_y$',fontsize=18)
##L = Y.max()
#solSigY = solSigma_grav(zSig,rho,cgrav)
#ax.plot(zSig,sigy,'-o', label=u'Sol Numérica')
#ax.plot(zSig,solSigY,'--', label=u'Sol Analítica',color='r')
#box = ax.get_position()
#ax.set_position([0.2*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 1.0, box.height])
##ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#ax.legend(loc='best')
#plt.savefig('Sigmay.pdf')
#plt.show()


fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)
#fig, ax = plt.subplots()

inDataNameSigma = './' + filepath + 'sigma.3'
X, Y, sigmaX, sigmaY = np.loadtxt(inDataNameSigma,unpack=True,usecols=[1,2,3,4])

#UU, VV = np.meshgrid(U,V)
#XX, YY = np.meshgrid(X,Y)

sig = np.sqrt(sigmaX**2. + sigmaY**2.)

plt.xlabel(r'$x$',fontsize=18)

plt.ylabel(r'$y$',fontsize=18)

plt.quiver(X,Y,sigmaX,sigmaY, sig,cmap=cm.afmhot_r)

plt.colorbar()

plt.savefig('./'+ filepath + 'quiverSig.pdf')
#plt.show()

