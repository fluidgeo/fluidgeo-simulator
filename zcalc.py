#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# Definindo a precisão das impressões
np.set_printoptions(precision=50)
plt.clf()

# Constante de Boltzmann
kB = 1.381e-23
# Número de Avogrado
Nav = 6.022e23
# Constante dos Gases
R   = Nav*kB

# Propriedades do Metano
T   = 320.0            # Temperatura (Kelvin)
Tb  = 111.63           # Ponto de bolha normal (Kelvin)
Tc  = 190.6            # Temperatura crítica (Kelvin)
Tr  = T/Tc             # Temperatura reduzida (adim)
Tbr = Tb/Tr            # Ponto de bolha reduzido (adim)
Pc  = 4.599e6          # Pressão crítica (Pascal)

for i in ('vdW','RK','SRK','PR', 'mRK', 'mPR', 'mPRnew'):
	print i
	if (i == 'vdW'):
		# Parâmetros da Eq de Estado de vdW

		epsilon = 0.0
		sigma   = 0.0
		psi = 27.0/64.0
		omega = 1.0/8.0
		alpha = 1.0

	if (i == 'RK'):
		# Parâmetros da Eq de Estado de Redlich/Kwong

		epsilon = 0.0
		sigma = 1.0
		psi = 0.42748023354
		omega = 0.086640349965
		alpha = np.sqrt(Tc/T)
		
	if (i == 'SRK'):
		# Parâmetros da Eq de Estado de Soave/Redlich/Kwong
		
		epsilon = 0.0
		sigma = 1.0
		psi = 0.42748023354
		omega = 0.086640349965
		w = 0.0113
		alpha = (1.0+(0.480+1.574*w-0.176*w**2.0)*(1.0-np.sqrt(Tr)))**2.0
	
	if (i == 'mRK'):
		# Parâmetros da Eq de Estado de Redlich/Kwong modificada

		epsilon = 0.0
		sigma = 1.0
		psi = 0.42748023354
		omega = 0.086640349965
		w = 0.0113
		#alpha0 = Tr**(-0.201158)*np.exp(0.141599*(1.0-Tr**(2.29528)))
		alpha0 = Tr**(-1.1)*np.exp(0.441411*(1.0-Tr**(-1.3)))
		#alpha1 = Tr**(-0.660145)*np.exp(0.500315*(1.0-Tr**(2.63165)))
		alpha1 = Tr**(-2.31278)*np.exp(0.032580*(1.0-Tr**(-10.3128)))
		alpha = alpha0 + w*(alpha1 - alpha0)
	
	if (i == 'PR'):
		# Parâmetros da Eq de Estado de Peng-Robinson
		
		epsilon = 1.0 - np.sqrt(2.0)
		sigma = 1.0 + np.sqrt(2.0)
		psi = 0.457235528921
		#omega = 0.07780 # O.0777960739039
		omega = 0.0777960739039
		w = 0.0113
		alpha = (1.0+(0.37464+1.54226*w-0.26992*w**2.0)*(1.0-np.sqrt(Tr)))**2.0
		
	if (i == 'mPRnew'):
		# Parâmetros da Eq de Estado de Peng-Robinson modificada 2015
		
		epsilon = 1.0 - np.sqrt(2.0)
		sigma = 1.0 + np.sqrt(2.0)
		#psi = 0.45724
		psi = 0.457235528921
		#omega = 0.07780
		omega = 0.0777960739039
		w = 0.0113
		a1 = -0.69560
		a2 = 0.786752
		a3 = 0.710102
		a4 = 0.140024
		b1 = -1.07370
		b2 = 1.604866
		b3 = -0.93029
		b4 = -0.44313
		c1 = -0.56996
		c2 = 1.0
		d = 0.594998
		beta = (np.abs(1.0-Tr**a1)**a2)/(np.abs(1.0-Tbr**a3)**a4)
		gamma = Tr**(b1*w**3.0+b2*w**2.0+b3*w+b4)
		alpha = (c1*beta + c2*gamma)**d
		
	if (i == 'mPR'):
		# Parâmetros da Eq de Estado de Peng-Robinson modificada
		
		epsilon = 1.0 - np.sqrt(2.0)
		sigma = 1.0 + np.sqrt(2.0)
		psi = 0.457235528921
		#omega = 0.07780 # O.0777960739039
		omega = 0.0777960739039
		w = 0.0113
		#w = 0.001
		#alpha0 = Tr**(-0.171813)*np.exp(0.125283*(1.0-Tr**(1.77634)))
		alpha0 = Tr**(-0.792615)*np.exp(0.401219*(1.0-Tr**(-0.992615)))
		#alpha1 = Tr**(-0.607352)*np.exp(0.511614*(1.0-Tr**(2.20517)))
		alpha1 = Tr**(-1.98471)*np.exp(0.024955*(1.0-Tr**(-9.98471)))
		alpha = alpha0 + w*(alpha1 - alpha0)
	
	# Definindo a faixa de variação da Concentração
	C = np.arange(1.0,22000.0,1.0)

	# Calculando P e Z para uma mesma faixa de C
	a = psi*(alpha*(R**2.0)*(Tc**2.0))/Pc
	b = omega*(R*Tc)/Pc
	P = (R*T*C)/(1.0-C*b) - (a*C*C)/((1.0 + epsilon*b*C)*(1.0 + sigma*b*C))
	Z = 1.0/(1.0-C*b) - (a*C)/(R*T*(1.0 + epsilon*b*C)*(1.0 + sigma*b*C))
	
	# Aproximação polinomial por mínimos quadrados
	Pb = np.arange(0.05,6.8,0.01)*1e7
	Cb = np.interp(Pb,P,C)
	Zb = 1.0/(1.0-Cb*b) - (a*Cb)/(R*T*(1.0 + epsilon*b*Cb)*(1.0 + sigma*b*Cb))
	Zap, residuals, rank, singular_values, rcond = np.polyfit(Pb,Zb,9,full=True)
	print "Residuos L2 do MMQ = ", residuals
	np.savetxt(("coef_%s.dat" % i),Zap,fmt='%.50e')

	# Plotando Z em função de P
	plt.plot(P,Z,label=i,linewidth=2.0)
	
axes = plt.gca()
axes.set_xlim([5.e5,6.8e7])
axes.set_ylim([Z.min()-0.1*Z.min(),Z.max()+0.2*Z.max()])
plt.grid(True)
plt.legend(loc='best')
plt.title("Z(p) Metano")
plt.xlabel("p (Pa)")
plt.ylabel("Z")
plt.savefig("Zpy.pdf")
#plt.show()
