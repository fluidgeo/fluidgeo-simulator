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

#fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

#ax = plt.subplot(111)
#fig, ax = plt.subplots()

filename = './' + filepath + 'passosPressaoBlocoMacro.dat'

inDataLegendP = np.loadtxt(filename,unpack=True)
dataLegendP = inDataLegendP[1:]
pnum = len(dataLegendP)

inDataNameDisp = './' + filepath + 'disp.1'
X, Y = np.loadtxt(inDataNameDisp,unpack=True,usecols=[1,2])

itU = 51
xU = np.linspace(X.min(),X.max(),itU)
yU = np.linspace(Y.min(),Y.max(),itU)
xxU, yyU = np.meshgrid(xU, yU)

sxU = np.zeros((itU,itU))
syU = np.zeros((itU,itU))
#i = 5
#fig = plt.figure()
	#ax = fig.gca(projection='3d')
	
	#ax = fig.add_subplot(111, projection='3d')
	
#ax = plt.subplot(111)
#inDataNameDisp = './' + filepath + ('disp.%d' % i)
#X, Y, U, V = np.loadtxt(inDataNameDisp,unpack=True,usecols=[1,2,3,4])
	#inDataNameNodes = './dxgeo01/fnodes.stoc'
	#U, V = np.loadtxt(inDataNameNodes,unpack=True,usecols=[3,4])
	#skip = (slice(None, None, 3), slice(None, None, 3))
	#UU, VV = np.meshgrid(U,V)
	#XX, YY = np.meshgrid(X,Y)
	#disp = np.sqrt(UU**2. + VV**2.)
#disp = np.sqrt(U**2. + V**2.)
	#print len(disp)
	#im = ax.imshow(disp, extent=[X.min(), X.max(), Y.min(), Y.max()])
	#ax.set_xlabel(r'$x\,(m)$',fontsize=18)
#plt.xlabel(r'$x$',fontsize=18)
	#ax.set_ylabel(r'$y\,(m)$',fontsize=18)
#plt.ylabel(r'$y$',fontsize=18)
	#ax.quiver(X,Y,U,V, disp, cmap='Blues' ,angles='xy', scale_units='xy', scale=1.2)
	#ax.quiver(X,Y,U,V, disp, cmap='Blues' ,angles='xy', scale_units='xy')
#plt.quiver(X,Y,U,V, disp,cmap=cm.afmhot_r)
	#plt.imshow(disp)
	#plt.plot(XX,YY)
#plt.colorbar()
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
#plt.savefig('./'+ filepath + ('quiverUgrav%d.pdf' % i))
	#plt.show()

for i in range(1,pnum+1):
	fig = plt.figure()
	#ax = fig.gca(projection='3d')

	#ax = fig.add_subplot(111, projection='3d')

	ax = plt.subplot(111)
	inDataNameDisp = './' + filepath + ('disp.%d' % i)
	U, V = np.loadtxt(inDataNameDisp,unpack=True,usecols=[3,4])
	k = 0
	for ii in range(itU):
		for jj in range(itU):
			#print "k = ", k
			#print ("phi[%d] = " % k), phin[k]
			sxU[ii,jj] = U[k]
			syU[ii,jj] = V[k]
			#print ("sX[%d,%d] = " % (i,j)), sX[i,j]
			k = k + 1 
	#print sxU
	#sys.exit("sxU")
	disp = np.sqrt(sxU*sxU + syU*syU)

	plt.streamplot(xxU, yyU, sxU, syU, color=disp, linewidth=2, cmap=plt.cm.autumn)
	plt.colorbar()
	#inDataNameNodes = './dxgeo01/fnodes.stoc'
	#U, V = np.loadtxt(inDataNameNodes,unpack=True,usecols=[3,4])
	#skip = (slice(None, None, 3), slice(None, None, 3))
	#UU, VV = np.meshgrid(U,V)
	#XX, YY = np.meshgrid(X,Y)
	#disp = np.sqrt(UU**2. + VV**2.)
	#disp = np.sqrt(U**2. + V**2.)
	#print len(disp)
	#im = ax.imshow(disp, extent=[X.min(), X.max(), Y.min(), Y.max()])
	#ax.set_xlabel(r'$x\,(m)$',fontsize=18)
	plt.xlabel(r'$x$',fontsize=18)
	#ax.set_ylabel(r'$y\,(m)$',fontsize=18)
	plt.ylabel(r'$y$',fontsize=18)
	#ax.quiver(X,Y,U,V, disp, cmap='Blues' ,angles='xy', scale_units='xy', scale=1.2)
	#ax.quiver(X,Y,U,V, disp, cmap='Blues' ,angles='xy', scale_units='xy')
	
	#plt.quiver(X,Y,U,V, disp,cmap=cm.afmhot_r)
	
	#plt.imshow(disp)
	#plt.plot(XX,YY)
	#plt.colorbar()
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
	plt.savefig('./'+ filepath + ('streamUgrav%d.png' % i))
	#sys.exit("streamlines")
	#plt.show()
#sys.exit("streamlines")
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

#XX, YY = np.meshgrid(X,Y)
#UU,VV = np.meshgrid(U,V)

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


#fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

#ax = plt.subplot(111)
#fig, ax = plt.subplots()

inDataNameSigma = './' + filepath + 'porosity.5'
phin = np.loadtxt(inDataNameSigma,unpack=True,usecols=[1])

it = 50
xSig = np.linspace(X.min(),X.max(),it)
ySig = np.linspace(Y.min(),Y.max(),it)

XX, YY = np.meshgrid(xSig, ySig)
XY = np.hstack((XX.ravel()[:,np.newaxis], YY.ravel()[:,np.newaxis]))

sX = np.zeros((it,it))
k = 0
for i in range(it):
	for j in range(it):
		#print "k = ", k
		#print ("phi[%d] = " % k), phin[k]
		sX[i,j] = phin[k]
		#print ("sX[%d,%d] = " % (i,j)), sX[i,j]
		k = k + 1 

#print sX
#sys.exit("Aqui")



#plt.figure()
#CS = plt.contour(XX, YY, sX, 4)
#plt.clabel(CS, inline=1, fontsize=10)
#plt.title('Porosidade')

#plt.figure()

#CS = plt.contour(XX, YY, sX, 3,
                 #colors='k', # negative contours will be dashed by default
                 #)
#plt.clabel(CS, fontsize=9, inline=1)
#plt.title('Porosidade')
#plt.savefig('./'+ filepath + 'phi_n.pdf')
#sys.exit("Porosity")

#plt.figure()

#z_min, z_max = np.abs(phin).min(), np.abs(phin).max()
#plt.pcolor(XX, YY, sX, cmap='RdBu', vmin=z_min, vmax=z_max)
#plt.title('pcolor')
##set the limits of the plot to the limits of the data
#plt.axis([XX.min(), XX.max(), YY.min(), YY.max()])
#plt.colorbar()

plt.figure()

plt.imshow(sX, interpolation='lanczos', cmap=cm.RdYlGn,
                origin='lower', extent=[XX.min(),XX.max(),YY.min(),YY.max()],
                vmax=abs(sX).max(), vmin=abs(sX).min())
plt.colorbar()
#plt.set_label(r"$\phi_n$")

plt.savefig('./'+ filepath + 'phi_n5.png')
plt.show()
sys.exit("Porosity")


fig, ax = plt.subplots()

ec = EllipseCollection(
                        XX,
                        YY,
                        sX,
                        units='x',
                        offsets=XY,
                        transOffset=ax.transData)
ec.set_array((sX).ravel())
ax.add_collection(ec)
ax.autoscale_view()
ax.set_xlabel('X')
ax.set_ylabel('y')
cbar = plt.colorbar(ec)
cbar.set_label('X+Y')
#plt.show()
plt.savefig('./'+ filepath + 'phi_n.pdf')
sys.exit("EllipseCollection")
#sys.exit("Parei aqui")
#UU, VV = np.meshgrid(U,V)
#XX, YY = np.meshgrid(xSig,ySig)

#sig = np.sqrt(sigmaX**2. + sigmaY**2.)
sig = np.sqrt(phin)

plt.xlabel(r'$x$',fontsize=18)

plt.ylabel(r'$y$',fontsize=18)

plt.quiver(xSig,ySig,sigmaX,sigmaY, sig,cmap=cm.afmhot_r)

plt.colorbar()

plt.savefig('./'+ filepath + 'phi_n.pdf')
#plt.show()

print "Displacement's plots OK"