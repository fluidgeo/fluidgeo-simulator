#!/usr/bin/python
# -*- coding: utf-8 -*-

# Post-processing program to plot results from fluid dynamics
# reservoir problem.
# Author: Diego Volpatto

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import sys

filepath = sys.argv[1]
filepath2 = sys.argv[2]
filepath3 = sys.argv[3]
filepath4 = sys.argv[4]
filepath5 = sys.argv[5]
filepath6 = sys.argv[6]

#plt.rcParams['legend.loc'] = 'best'
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

ano = 30.*12.*86400.
mes = 30.*86400.

# Perfil de pressão

filename = './' + filepath + 'passosPressaoBlocoMacro.dat'
filename2 = './' + filepath2 + 'passosPressaoBlocoMacro.dat'
filename3 = './' + filepath3 + 'passosPressaoBlocoMacro.dat'
filename4 = './' + filepath4 + 'passosPressaoBlocoMacro.dat'
filename5 = './' + filepath5 + 'passosPressaoBlocoMacro.dat'
filename6 = './' + filepath6 + 'passosPressaoBlocoMacro.dat'

inDataLegendP = np.loadtxt(filename,unpack=True)
dataLegendP = inDataLegendP[1:]
inDataLegendP2 = np.loadtxt(filename2,unpack=True)
dataLegendP2 = inDataLegendP2[1:]
inDataLegendP3 = np.loadtxt(filename3,unpack=True)
dataLegendP3 = inDataLegendP3[1:]
inDataLegendP4 = np.loadtxt(filename,unpack=True)
dataLegendP4 = inDataLegendP[1:]
inDataLegendP5 = np.loadtxt(filename2,unpack=True)
dataLegendP5 = inDataLegendP2[1:]
inDataLegendP6 = np.loadtxt(filename3,unpack=True)
dataLegendP6 = inDataLegendP3[1:]
#if (dataLegendP.any() != dataLegendP2.any() or dataLegendP.any() != dataLegendP3.any() or dataLegendP2.any() != dataLegendP3.any()):
#	sys.exit("Erro: Passos de tempo diferentes")
pnum = len(dataLegendP)
Px0 = np.zeros(pnum)
#pData = np.loadtxt('./expTeste/fort.1111',unpack=True)
#Pr = pData[0]
#Pw = pData[1]
#PI = pData[2]

fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.adsd_subplot(111, projection='3d')

ax = plt.subplot(111)

for i in range(1,pnum+1):
    #i_mask = 10*i
    inDataNameP = './' + filepath + ('solP_C.%d' % i)
    inDataNameP2 = './' + filepath2 + ('solP_C.%d' % i)
    inDataNameP3 = './' + filepath3 + ('solP_C.%d' % i)
    inDataNameP4 = './' + filepath4 + ('solP_C.%d' % i)
    inDataNameP5 = './' + filepath5 + ('solP_C.%d' % i)
    inDataNameP6 = './' + filepath6 + ('solP_C.%d' % i)
    xx, yy, P = np.loadtxt(inDataNameP,unpack=True,usecols=[2,3,4])
    xx2, yy2, P2 = np.loadtxt(inDataNameP2,unpack=True,usecols=[2,3,4])
    xx3, yy3, P3 = np.loadtxt(inDataNameP3,unpack=True,usecols=[2,3,4])
    xx4, yy4, P4 = np.loadtxt(inDataNameP4,unpack=True,usecols=[2,3,4])
    xx5, yy5, P5 = np.loadtxt(inDataNameP5,unpack=True,usecols=[2,3,4])
    xx6, yy6, P6 = np.loadtxt(inDataNameP6,unpack=True,usecols=[2,3,4])
    #print xx
    #print yy
    #print P
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    '''
    if dataLegendP[i-1] == 1.:
        leg = (u"%d mês" % dataLegendP[i-1])
    else:
        leg = ("%d meses" % dataLegendP[i-1])
    '''
    Px0[i-1] = P[0]
    #ax.plot(xx,yy,P,'o',label=leg)
    ax.plot(xx,P,label=r'Desacoplado')
    ax.plot(xx2,P2,label=r'Uma-via')
    ax.plot(xx3,P3,label=r'Duas-vias')
    ax.plot(xx4,P4,'-.',label=r'Desacoplado $(raso)$')
    ax.plot(xx5,P5,'-.',label=r'Uma-via $(raso)$')
    ax.plot(xx6,P6,'-.',label=r'Duas-vias $(raso)$')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel(r'$x$',fontsize=18)
#ax.set_ylabel('$y^{*}$',fontsize=18)
ax.set_ylabel(r'$p$',fontsize=18)
ax.grid(True)

plt.savefig('./' + filepath + 'tmpPressaoC_compH.pdf')
#plt.show()

fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.adsd_subplot(111, projection='3d')

ax = plt.subplot(111)

for i in range(1,pnum+1):
    #i_mask = 10*i
    inDataNameP = './' + filepath + ('solP_BC.%d' % i)
    inDataNameP2 = './' + filepath2 + ('solP_BC.%d' % i)
    inDataNameP3 = './' + filepath3 + ('solP_BC.%d' % i)
    inDataNameP4 = './' + filepath4 + ('solP_BC.%d' % i)
    inDataNameP5 = './' + filepath5 + ('solP_BC.%d' % i)
    inDataNameP6 = './' + filepath6 + ('solP_BC.%d' % i)
    xx, yy, P = np.loadtxt(inDataNameP,unpack=True,usecols=[2,3,4])
    xx2, yy2, P2 = np.loadtxt(inDataNameP2,unpack=True,usecols=[2,3,4])
    xx3, yy3, P3 = np.loadtxt(inDataNameP3,unpack=True,usecols=[2,3,4])
    xx4, yy4, P4 = np.loadtxt(inDataNameP4,unpack=True,usecols=[2,3,4])
    xx5, yy5, P5 = np.loadtxt(inDataNameP5,unpack=True,usecols=[2,3,4])
    xx6, yy6, P6 = np.loadtxt(inDataNameP6,unpack=True,usecols=[2,3,4])
    #print xx
    #print yy
    #print P
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    
    if dataLegendP[i-1] == 1.:
        leg = (u"%d mês" % dataLegendP[i-1])
    else:
        leg = ("%d meses" % dataLegendP[i-1])
    Px0[i-1] = P[0]
    #ax.plot(xx,yy,P,'o',label=leg)
    ax.plot(xx,P,label=r'Desacoplado')
    ax.plot(xx2,P2,label=r'Uma-via')
    ax.plot(xx3,P3,label=r'Duas-vias')
    ax.plot(xx4,P4,'-.',label=r'Desacoplado $(raso)$')
    ax.plot(xx5,P5,'-.',label=r'Uma-via $(raso)$')
    ax.plot(xx6,P6,'-.',label=r'Duas-vias $(raso)$')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel(r'$x$',fontsize=18)
#ax.set_ylabel('$y^{*}$',fontsize=18)
ax.set_ylabel(r'$p$',fontsize=18)
ax.grid(True)

plt.savefig('./' + filepath + 'tmpPressaoBC_compH.pdf')
#plt.show()
'''

fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)

for i in range(1,pnum+1):
    #i_mask = 11*i
    inDataNameGradP = './' + filepath + ('gradPx.%d' % i)
    inDataNameGradP2 = './' + filepath2 + ('gradPx.%d' % i)
    inDataNameGradP3 = './' + filepath3 + ('gradPx.%d' % i)
    xx, yy, gradP = np.loadtxt(inDataNameGradP,unpack=True,usecols=[2,3,4])
    xx2, yy2, gradP2 = np.loadtxt(inDataNameGradP2,unpack=True,usecols=[2,3,4])
    xx3, yy3, gradP3 = np.loadtxt(inDataNameGradP3,unpack=True,usecols=[2,3,4])
    #print xx
    #print yy
    #print gradP
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    if dataLegendP[i-1] == 1.:
        leg = (u"%d mês" % dataLegendP[i-1])
    else:
        leg = ("%d meses" % dataLegendP[i-1])
    Px0[i-1] = P[0]
    #ax.plot(xx,yy,gradP,'o',label=leg)
    ax.plot(xx,gradP,'-.',label=leg)
    ax.plot(xx2,gradP2,label=leg)
    ax.plot(xx3,gradP3,'x',label=leg)
    #surf = ax.plot_surface(xx, yy, P, rstride=1, cstride=1, cmap=cm.coolwarm,
        #linewidth=0, antialiased=True)
    #surf = ax.plot_surface(xx, yy, P)
    #plt.show()
    #plt.pcolor(xx,yy,P)
    #ax.contour(xx, yy, P)
    #plt.savefig('./teste/tmpPressao%dH.pdf' % i_mask)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel(r'$x$',fontsize=18)
#ax.set_ylabel('$y^{*}$',fontsize=18)
ax.set_ylabel(r"$\nabla p$",fontsize=18)
#plt.xlabel('$x^{*}$',fontsize=18)
#plt.ylabel('$p^{*}$',fontsize=18,rotation='horizontal')
#plt.grid(b=True, which='major', color='k', linestyle='--')
#ax.legend()
plt.grid(True)
plt.savefig('./' + filepath + 'tmpgradPressao_compH.pdf')
#plt.show()

fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)

for i in range(1,pnum+1):
    #i_mask = 14*i
    inDataNameV = './' + filepath + ('solVelocity_x.%d' % i)
    inDataNameV2 = './' + filepath2 + ('solVelocity_x.%d' % i)
    inDataNameV3 = './' + filepath3 + ('solVelocity_x.%d' % i)
    xx, yy, V = np.loadtxt(inDataNameV,unpack=True,usecols=[2,3,4])
    xx2, yy2, V2 = np.loadtxt(inDataNameV2,unpack=True,usecols=[2,3,4])
    xx3, yy3, V3 = np.loadtxt(inDataNameV3,unpack=True,usecols=[2,3,4])
    #print xx
    #print yy
    #print V
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    if dataLegendP[i-1] == 1.:
        leg = (u"%d mês" % dataLegendP[i-1])
    else:
        leg = ("%d meses" % dataLegendP[i-1])
    Px0[i-1] = P[0]
    #ax.plot(xx,yy,V,'o',label=leg)
    ax.plot(xx,V,'-.',label=leg)
    ax.plot(xx2,V2,label=leg)
    ax.plot(xx3,V3,'x',label=leg)
    #surf = ax.plot_surface(xx, yy, P, rstride=1, cstride=1, cmap=cm.coolwarm,
        #linewidth=0, antialiased=True)
    #surf = ax.plot_surface(xx, yy, P)
    #plt.show()
    #plt.pcolor(xx,yy,P)
    #ax.contour(xx, yy, P)
    #plt.savefig('./teste/tmpPressao%dH.pdf' % i_mask)

box = ax.get_position()
ax.set_position([0.45*box.x0+box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#ax.legend(loc='best')
ax.set_xlabel(r'$x$',fontsize=18)
#ax.set_ylabel('$y^{*}$',fontsize=18)
ax.set_ylabel(r"$u_D\,(m/s)$",fontsize=18)
#plt.xlabel('$x^{*}$',fontsize=18)
#plt.ylabel('$p^{*}$',fontsize=18,rotation='horizontal')
#plt.grid(b=True, which='major', color='k', linestyle='--')
#ax.legend()
plt.grid(True)
plt.savefig('./' + filepath + 'tmpV_compH.pdf')
#plt.show()

fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)


for i in range(1,pnum+1):
    #i_mask = 19*i
    inDataNameJ = './' + filepath + ('nodeFlux_x.%d' % i)
    inDataNameJ2 = './' + filepath2 + ('nodeFlux_x.%d' % i)
    inDataNameJ3 = './' + filepath3 + ('nodeFlux_x.%d' % i)
    xx, yy, J = np.loadtxt(inDataNameJ,unpack=True,usecols=[2,3,4])
    xx2, yy2, J2 = np.loadtxt(inDataNameJ2,unpack=True,usecols=[2,3,4])
    xx3, yy3, J3 = np.loadtxt(inDataNameJ3,unpack=True,usecols=[2,3,4])
    #print xx
    #print yy
    #print J
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    if dataLegendP[i-1] == 1.:
        leg = (u"%d mês" % dataLegendP[i-1])
    else:
        leg = ("%d meses" % dataLegendP[i-1])
    Px0[i-1] = P[0]
    #ax.plot(xx,yy,J,'.',label=leg)
    #ax.plot(xx,yy,J,'o',label=leg)
    ax.plot(xx,J,'-.',label=leg)
    ax.plot(xx2,J2,label=leg)
    ax.plot(xx3,J3,'x',label=leg)
    #surf = ax.plot_surface(xx, yy, P, rstride=1, cstride=1, cmap=cm.coolwarm,
        #linewidth=0, antialiased=True)
    #surf = ax.plot_surface(xx, yy, P)
    #plt.show()
    #plt.pcolor(xx,yy,P)
    #ax.contour(xx, yy, P)
    #plt.savefig('./teste/tmpPressao%dH.pdf' % i_mask)

box = ax.get_position()
ax.set_position([0.45*box.x0+box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
ax.set_xlabel('$x$',fontsize=18)
#ax.set_ylabel(r'$y^{*}$',fontsize=18)
ax.set_ylabel(r'$J \,\left(\frac{kg}{m^2 s}\right)$',fontsize=18)
#plt.xlabel('$x^{*}$',fontsize=18)
#plt.ylabel('$p^{*}$',fontsize=18,rotation='horizontal')
#plt.grid(b=True, which='major', color='k', linestyle='--')
#ax.legend()
plt.grid(True)
plt.savefig('./' + filepath + 'tmpJ_compH.pdf')
#plt.show()

fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)


for i in range(1,pnum+1):
    #i_mask = 23*i
    inDataNameResid = './' + filepath + ('residueFlux_x.%d' % i)
    inDataNameResid2 = './' + filepath2 + ('residueFlux_x.%d' % i)
    inDataNameResid3 = './' + filepath3 + ('residueFlux_x.%d' % i)
    xx, yy, Resid = np.loadtxt(inDataNameResid,unpack=True,usecols=[2,3,4])
    xx2, yy2, Resid2 = np.loadtxt(inDataNameResid2,unpack=True,usecols=[2,3,4])
    xx3, yy3, Resid3 = np.loadtxt(inDataNameResid3,unpack=True,usecols=[2,3,4])
    #print xx
    #print yy
    #print Resid
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    if dataLegendP[i-1] == 1.:
	leg = (u"%d mês" % dataLegendP[i-1])
    else:
	leg = ("%d meses" % dataLegendP[i-1])
    Px0[i-1] = P[0]
    #ax.plot(xx,yy,J,'.',label=leg)
    #ax.plot(xx,yy,J,'o',label=leg)
    ax.plot(xx,Resid,'-.',label=leg)
    ax.plot(xx2,Resid2,label=leg)
    ax.plot(xx3,Resid3,'x',label=leg)
    #surf = ax.plot_surface(xx, yy, P, rstride=1, cstride=1, cmap=cm.coolwarm,
        #linewidth=0, antialiased=True)
    #surf = ax.plot_surface(xx, yy, P)
    #plt.show()
    #plt.pcolor(xx,yy,P)
    #ax.contour(xx, yy, P)
    #plt.savefig('./teste/tmpPressao%dH.pdf' % i_mask)

box = ax.get_position()
ax.set_position([0.2*box.x0+box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel('$x$',fontsize=18)
#ax.set_ylabel(r'$y^{*}$',fontsize=18)
ax.set_ylabel(u'Resíduo',fontsize=18)
#plt.xlabel('$x^{*}$',fontsize=18)
#plt.ylabel('$p^{*}$',fontsize=18,rotation='horizontal')
#plt.grid(b=True, which='major', color='k', linestyle='--')
#ax.legend()
plt.grid(True)
plt.savefig('./' + filepath + 'tmpResid_compH.pdf')
#plt.show()
'''
## Produção do bloco
fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)

inDataNameJprod = './' + filepath + 'echoProducao.dat'
inDataNameJprod2 = './' + filepath2 + 'echoProducao.dat'
inDataNameJprod3 = './' + filepath3 + 'echoProducao.dat'
inDataNameJprod4 = './' + filepath4 + 'echoProducao.dat'
inDataNameJprod5 = './' + filepath5 + 'echoProducao.dat'
inDataNameJprod6 = './' + filepath6 + 'echoProducao.dat'
dt, Jprod = np.loadtxt(inDataNameJprod,unpack=True,usecols=[1,5])
dt2, Jprod2 = np.loadtxt(inDataNameJprod2,unpack=True,usecols=[1,5])
dt3, Jprod3 = np.loadtxt(inDataNameJprod3,unpack=True,usecols=[1,5])
dt4, Jprod4 = np.loadtxt(inDataNameJprod4,unpack=True,usecols=[1,5])
dt5, Jprod5 = np.loadtxt(inDataNameJprod5,unpack=True,usecols=[1,5])
dt6, Jprod6 = np.loadtxt(inDataNameJprod6,unpack=True,usecols=[1,5])
dt = dt/mes
dt2 = dt2/mes
dt3 = dt3/mes
dt4 = dt4/mes
dt5 = dt5/mes
dt6 = dt6/mes
ax.set_xlabel(r'$t\,(meses)$',fontsize=18)
ax.set_ylabel(r'$Produ \c c \~ a o\, \left(kg\right)$',fontsize=16)
ax.plot(dt,Jprod,label=u'Desacoplado')
ax.plot(dt2,Jprod2,label=u'Uma-via')
ax.plot(dt3,Jprod3,label=u'Duas-vias')
ax.plot(dt4,Jprod4,'-.',label=r'Desacoplado $(raso)$')
ax.plot(dt5,Jprod5,'-.',label=r'Uma-via $(raso)$')
ax.plot(dt6,Jprod6,'-.',label=r'Duas-vias $(raso)$')
box = ax.get_position()
ax.set_position([0.1*box.x0+box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'Prod_compH.pdf')
#plt.show()
'''
# RF barra
fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)

inDataNameRF_ = './' + filepath + 'echoProducao.dat'
inDataNameRF_2 = './' + filepath2 + 'echoProducao.dat'
inDataNameRF_3 = './' + filepath3 + 'echoProducao.dat'
dt, RF_ = np.loadtxt(inDataNameRF_,unpack=True,usecols=[1,6])
dt2, RF_2 = np.loadtxt(inDataNameRF_2,unpack=True,usecols=[1,6])
dt3, RF_3 = np.loadtxt(inDataNameRF_3,unpack=True,usecols=[1,6])
dt = dt/mes
dt2 = dt2/mes
dt3 = dt3/mes
ax.set_xlabel(r'$t\,(meses)$',fontsize=18)
ax.set_ylabel(r'$RF_{recuper\'avel}$',fontsize=16)
ax.plot(dt,RF_,'-.',label=u'RF')
ax.plot(dt2,RF_2,label=u'RF')
ax.plot(dt3,RF_3,'x',label=u'RF')
box = ax.get_position()
ax.set_position([0.1*box.x0+box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid(True)
plt.savefig('./' + filepath + 'RF_compH.pdf')
#plt.show()

# RF
fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)

inDataNameRF = './' + filepath + 'echoProducao.dat'
inDataNameRF2 = './' + filepath2 + 'echoProducao.dat'
inDataNameRF3 = './' + filepath3 + 'echoProducao.dat'
dt, RF = np.loadtxt(inDataNameRF,unpack=True,usecols=[1,7])
dt2, RF2 = np.loadtxt(inDataNameRF2,unpack=True,usecols=[1,7])
dt3, RF3 = np.loadtxt(inDataNameRF3,unpack=True,usecols=[1,7])
dt = dt/mes
dt2 = dt2/mes
dt3 = dt3/mes
ax.set_xlabel(r'$t\,(meses)$',fontsize=18)
ax.set_ylabel(r'$RF_{total}$',fontsize=16)
ax.plot(dt,RF,'-.',label=u'RF')
ax.plot(dt2,RF2,label=u'RF')
ax.plot(dt3,RF3,label=u'RF')
box = ax.get_position()
ax.set_position([0.1*box.x0+box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid(True)
plt.savefig('./' + filepath + 'RFcompH.pdf')
#plt.show()
'''
fig = plt.figure()
#ax = fig.gca(projection='3d')

#ax = fig.add_subplot(111, projection='3d')

ax = plt.subplot(111)

inDataNameRF = './' + filepath + 'echoProducao.dat'
inDataNameRF2 = './' + filepath2 + 'echoProducao.dat'
inDataNameRF3 = './' + filepath3 + 'echoProducao.dat'
inDataNameRF4 = './' + filepath4 + 'echoProducao.dat'
inDataNameRF5 = './' + filepath5 + 'echoProducao.dat'
inDataNameRF6 = './' + filepath6 + 'echoProducao.dat'
dt, Flux = np.loadtxt(inDataNameRF,unpack=True,usecols=[1,3])
dt2, Flux2 = np.loadtxt(inDataNameRF2,unpack=True,usecols=[1,3])
dt3, Flux3 = np.loadtxt(inDataNameRF3,unpack=True,usecols=[1,3])
dt4, Flux4 = np.loadtxt(inDataNameRF4,unpack=True,usecols=[1,3])
dt5, Flux5 = np.loadtxt(inDataNameRF5,unpack=True,usecols=[1,3])
dt6, Flux6 = np.loadtxt(inDataNameRF6,unpack=True,usecols=[1,3])
dt = dt/mes
dt2 = dt2/mes
dt3 = dt3/mes
dt4 = dt4/mes
dt5 = dt5/mes
dt6 = dt6/mes
ax.set_xlabel(r'$t\,(meses)$',fontsize=18)
ax.set_ylabel(r'J $\left(\frac{kg}{m^2 s}\right)$',fontsize=16)
ax.plot(dt,Flux,label=r'Desacoplado')
ax.plot(dt2,Flux2,label=r'Uma-via')
ax.plot(dt3,Flux3,label=r'Duas-vias')
ax.plot(dt4,Flux4,'-.',label=r'Desacoplado $(raso)$')
ax.plot(dt5,Flux5,'-.',label=r'Uma-via $(raso)$')
ax.plot(dt6,Flux6,'-.',label=r'Duas-vias $(raso)$')
box = ax.get_position()
ax.set_position([0.45*box.x0+box.x0, box.y0, box.width * 0.95, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.legend(loc='best')
plt.grid(True)
plt.savefig('./' + filepath + 'FluxoH.pdf')
#plt.show()

print "Comparative pressure's plots OK"
