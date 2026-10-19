"""
Script containing code for matrix element manipulation (reweighting, computing optimal observables)
Based on Fortran code automatically produced by Madgraph5/ALOHA/HELAS, and FeynRules model file stuff
converted to Python by ALOHA.
"""

from __future__ import division 
from __future__ import absolute_import
import math
from math import sqrt, pow
from itertools import product
from six.moves import range

import numpy as np

#------------------------IMPORTED FROM wavefunctions.py------------------------

class WaveFunction(list):


	"""a objet for a WaveFunction"""
	spin_to_size={0:1,
				  1:3,
				  2:6,
				  3:6,
				  4:18,
				  5:18}
    
	def __init__(self, spin= None, size=None):
		"""Init the list with zero value"""

		if spin:
			size = self.spin_to_size[spin]
		list.__init__(self, [0]*size)  

def ixxxxx(p,fmass,nhel,nsf):
	"""Defines an inflow fermion."""

	fi = WaveFunction(2)
    
	fi[0] = complex(-p[0]*nsf,-p[3]*nsf)
	fi[1] = complex(-p[1]*nsf,-p[2]*nsf) 
	nh = nhel*nsf 
	if (fmass != 0.):
		pp = min(p[0],sqrt(p[1]**2 + p[2]**2 + p[3]**2 ))
		if (pp == 0.): 
			sqm = [sqrt(abs(fmass))]
			sqm.append(sign(sqm[0],fmass)) 
			ip = (1+nh)//2 
			im = (1-nh)//2 

			fi[2] = ip*sqm[ip]
			fi[3] = im*nsf*sqm[ip]
			fi[4] = ip*nsf*sqm[im]
			fi[5] = im*sqm[im]

		else:
			sf = [(1+nsf+(1-nsf)*nh)*0.5,(1+nsf-(1-nsf)*nh)*0.5]
			omega = [sqrt(p[0]+pp),fmass/(sqrt(p[0]+pp))]
			ip = int((1+nh)//2)
			im = int((1-nh)//2)
			sfomeg = [sf[0]*omega[ip],sf[1]*omega[im]]
			pp3 = max(pp+p[3],0.)
			if (pp3 == 0.):
				chi1 = complex(-nh,0.) 
			else:
				chi1 = complex(nh*p[1]/sqrt(2.*pp*pp3),\
				p[2]/sqrt(2.*pp*pp3))
			chi = [complex(sqrt(pp3*0.5/pp)),chi1]

			fi[2] = sfomeg[0]*chi[im]
			fi[3] = sfomeg[0]*chi[ip]
			fi[4] = sfomeg[1]*chi[im]
			fi[5] = sfomeg[1]*chi[ip] 
    
	else: 
		sqp0p3 = sqrt(max(p[0]+p[3],0.))*nsf
		if (sqp0p3 == 0.):
			chi1 = complex(-nhel*sqrt(2.*p[0]),0.)
		else:
			chi1 = complex(nh*p[1]/sqp0p3,p[2]/sqp0p3)
		chi = [complex(sqp0p3,0.),chi1]
		if (nh == 1):
			fi[2] = complex(0.,0.)
			fi[3] = complex(0.,0.)
			fi[4] = chi[0]
			fi[5] = chi[1] 
		else:
			fi[2] = chi[1]
			fi[3] = chi[0]
			fi[4] = complex(0.,0.)
			fi[5] = complex(0.,0.) 
    
	return fi 

def oxxxxx(p,fmass,nhel,nsf):
	""" initialize an outgoing fermion"""

	fo = WaveFunction(2)
	fo[0] = complex(p[0]*nsf,p[3]*nsf)
	fo[1] = complex(p[1]*nsf,p[2]*nsf)
	nh = nhel*nsf
	if (fmass != 0.):
		pp = min(p[0],sqrt(p[1]**2 + p[2]**2 + p[3]**2 ))
		if (pp == 0.): 
			sqm = [sqrt(abs(fmass))]
			sqm.append( sign(sqm[0], fmass)) 
			ip = int(-((1-nh)//2) * nhel)
			im = int((1+nh)//2 * nhel)

			fo[2] = im*sqm[abs(im)]
			fo[3] = ip*nsf*sqm[abs(im)]
			fo[4] = im*nsf*sqm[abs(ip)]
			fo[5] = ip*sqm[abs(ip)]

		else:
			sf = [(1+nsf+(1-nsf)*nh)*0.5,(1+nsf-(1-nsf)*nh)*0.5]
			omega = [sqrt(p[0]+pp),fmass/(sqrt(p[0]+pp))]
			ip = int((1+nh)//2)
			im = int((1-nh)//2)
			sfomeg = [sf[0]*omega[ip],sf[1]*omega[im]]
			pp3 = max(pp+p[3],0.)
			if (pp3 == 0.):
				chi1 = complex(-nh,0.) 
			else:
				chi1 = complex(nh*p[1]/sqrt(2.*pp*pp3),\
				-p[2]/sqrt(2.*pp*pp3))
			chi = [complex(sqrt(pp3*0.5/pp)),chi1]

			fo[2] = sfomeg[1]*chi[im]
			fo[3] = sfomeg[1]*chi[ip]
			fo[4] = sfomeg[0]*chi[im]
			fo[5] = sfomeg[0]*chi[ip] 
            
	else: 
		sqp0p3 = sqrt(max(p[0]+p[3],0.))*nsf
		if (sqp0p3 == 0.):
			chi1 = complex(-nhel*sqrt(2.*p[0]),0.)
		else:
			chi1 = complex(nh*p[1]/sqp0p3,-p[2]/sqp0p3)
		chi = [complex(sqp0p3,0.),chi1]
		if (nh == 1):
			fo[2] = chi[0]
			fo[3] = chi[1]
			fo[4] = complex(0.,0.)
			fo[5] = complex(0.,0.) 
		else:
			fo[2] = complex(0.,0.)
			fo[3] = complex(0.,0.)
			fo[4] = chi[1]
			fo[5] = chi[0] 
    
	return fo

def vxxxxx(p,vmass,nhel,nsv):
	""" initialize a vector wavefunction. nhel=4 is for checking BRST"""

	vc = WaveFunction(3)

	sqh = sqrt(0.5)
	nsvahl = nsv*abs(nhel)
	pt2 = p[1]**2 + p[2]**2 
	pp = min(p[0],sqrt(pt2 + p[3]**2))
	pt = min(pp,sqrt(pt2))

	vc[0] = complex(p[0]*nsv,p[3]*nsv)
	vc[1] = complex(p[1]*nsv,p[2]*nsv)

	if (nhel == 4):
		if (vmass == 0.):
			vc[2] = 1.
			vc[3]=p[1]/p[0]
			vc[4]=p[2]/p[0]
			vc[5]=p[3]/p[0]
		else:
			vc[2] = p[0]/vmass
			vc[3] = p[1]/vmass
			vc[4] = p[2]/vmass
			vc[5] = p[3]/vmass
        
		return vc 

	if (vmass != 0.):
		hel0 = 1.-abs(nhel) 

		if (pp == 0.):
			vc[2] = complex(0.,0.)
			vc[3] = complex(-nhel*sqh,0.)
			vc[4] = complex(0.,nsvahl*sqh) 
			vc[5] = complex(hel0,0.)

		else:
			emp = p[0]/(vmass*pp)
			vc[2] = complex(hel0*pp/vmass,0.)
			vc[5] = complex(hel0*p[3]*emp+nhel*pt/pp*sqh)
			if (pt != 0.):
				pzpt = p[3]/(pp*pt)*sqh*nhel
				vc[3] = complex(hel0*p[1]*emp-p[1]*pzpt, -nsvahl*p[2]/pt*sqh)
				vc[4] = complex(hel0*p[2]*emp-p[2]*pzpt, nsvahl*p[1]/pt*sqh) 
			else:
				vc[3] = complex(-nhel*sqh,0.)
				vc[4] = complex(0.,nsvahl*sign(sqh,p[3]))
	else: 
		pp = p[0]
		pt = sqrt(p[1]**2 + p[2]**2)
		vc[2] = complex(0.,0.)
		vc[5] = complex(nhel*pt/pp*sqh)
		if (pt != 0.):
			pzpt = p[3]/(pp*pt)*sqh*nhel
			vc[3] = complex(-p[1]*pzpt,-nsv*p[2]/pt*sqh)
			vc[4] = complex(-p[2]*pzpt,nsv*p[1]/pt*sqh)
		else:
			vc[3] = complex(-nhel*sqh,0.)
			vc[4] = complex(0.,nsv*sign(sqh,p[3]))
    
	return vc

def sign(x,y):
	"""Fortran's sign transfer function"""
	try:
		cmp = (y < 0.)
	except TypeError:
		# y should be complex
		if abs(y.imag) < 1e-6 * abs(y.real):
			y = y.real
		else:
			raise
	finally:
		if (y < 0.):
			return -abs(x) 
		else:
			return abs(x) 

#---------------------------MATRIX ELEMENT FUNCTIONS---------------------------

def FFV1_0(F1,F2,V3,COUP):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#Gamma(3,2,1)
	TMP2 = (F1[2]*(F2[4]*(V3[2]+V3[5])+F2[5]*(V3[3]+1j*(V3[4])))+(F1[3]*(F2[4]*(V3[3]-1j*(V3[4]))+F2[5]*(V3[2]-V3[5]))+(F1[4]*(F2[2]*(V3[2]-V3[5])-F2[3]*(V3[3]+1j*(V3[4])))+F1[5]*(F2[2]*(+1j*(V3[4])-V3[3])+F2[3]*(V3[2]+V3[5])))))
	vertex = COUP*-1j * TMP2
	return vertex

def FFV2_0(F1,F2,V3,COUP):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#-(P(-1,3)*Gamma(-1,2,-2)*Gamma(3,-2,1)) + P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2)

	P3 = [complex(V3[0]).real, complex(V3[1]).real, complex(V3[1]).imag, complex(V3[0]).imag]
	TMP4 = (F1[2]*(F2[2]*(P3[0]*(V3[2]+V3[5])+(P3[1]*(-1)*(V3[3]+1j*(V3[4]))+(P3[2]*(+1j*(V3[3])-V3[4])-P3[3]*(V3[2]+V3[5]))))+F2[3]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4]))))))+(F1[3]*(F2[2]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F2[3]*(P3[0]*(V3[2]-V3[5])+(P3[1]*(-V3[3]+1j*(V3[4]))+(P3[2]*(-1)*(+1j*(V3[3])+V3[4])+P3[3]*(V3[2]-V3[5])))))+(F1[4]*(F2[4]*(P3[0]*(V3[2]-V3[5])+(P3[1]*(-1)*(V3[3]+1j*(V3[4]))+(P3[2]*(+1j*(V3[3])-V3[4])+P3[3]*(V3[2]-V3[5]))))+F2[5]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4]))))))+F1[5]*(F2[4]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F2[5]*(P3[0]*(V3[2]+V3[5])+(P3[1]*(-V3[3]+1j*(V3[4]))+(P3[2]*(-1)*(+1j*(V3[3])+V3[4])-P3[3]*(V3[2]+V3[5]))))))))
	TMP5 = (F1[2]*(F2[2]*(P3[0]*(V3[2]-V3[5])+(P3[1]*(-V3[3]+1j*(V3[4]))+(P3[2]*(-1)*(+1j*(V3[3])+V3[4])+P3[3]*(V3[2]-V3[5]))))+F2[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4]))))))+(F1[3]*(F2[2]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F2[3]*(P3[0]*(V3[2]+V3[5])+(P3[1]*(-1)*(V3[3]+1j*(V3[4]))+(P3[2]*(+1j*(V3[3])-V3[4])-P3[3]*(V3[2]+V3[5])))))+(F1[4]*(F2[4]*(P3[0]*(V3[2]+V3[5])+(P3[1]*(-V3[3]+1j*(V3[4]))+(P3[2]*(-1)*(+1j*(V3[3])+V3[4])-P3[3]*(V3[2]+V3[5]))))+F2[5]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))-P3[3]*(V3[3]+1j*(V3[4]))))))+F1[5]*(F2[4]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F2[5]*(P3[0]*(V3[2]-V3[5])+(P3[1]*(-1)*(V3[3]+1j*(V3[4]))+(P3[2]*(+1j*(V3[3])-V3[4])+P3[3]*(V3[2]-V3[5]))))))))
	vertex = COUP*(+1j*(TMP4)-1j*(TMP5))
	return vertex

def FFV3_0(F1,F2,V3,COUP):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#-(P(-1,3)*Gamma5(-2,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)) + P(-1,3)*Gamma5(-2,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)

	P3 = [complex(V3[0]).real, complex(V3[1]).real, complex(V3[1]).imag, complex(V3[0]).imag]
	TMP7 = (-1)*(F1[2]*(F2[2]*(P3[0]*(V3[2]+V3[5])+(P3[1]*(-1)*(V3[3]+1j*(V3[4]))+(P3[2]*(+1j*(V3[3])-V3[4])-P3[3]*(V3[2]+V3[5]))))+F2[3]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4]))))))+(F1[3]*(F2[2]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F2[3]*(P3[0]*(V3[2]-V3[5])+(P3[1]*(-V3[3]+1j*(V3[4]))+(P3[2]*(-1)*(+1j*(V3[3])+V3[4])+P3[3]*(V3[2]-V3[5])))))+(F1[4]*(F2[4]*(P3[0]*(-V3[2]+V3[5])+(P3[1]*(V3[3]+1j*(V3[4]))+(P3[2]*(-1j*(V3[3])+V3[4])+P3[3]*(-V3[2]+V3[5]))))+F2[5]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))-P3[3]*(V3[3]+1j*(V3[4]))))))+F1[5]*(F2[4]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F2[5]*(P3[0]*(-1)*(V3[2]+V3[5])+(P3[1]*(V3[3]-1j*(V3[4]))+(P3[2]*(+1j*(V3[3])+V3[4])+P3[3]*(V3[2]+V3[5]))))))))
	TMP8 = (-1)*(F1[2]*(F2[2]*(P3[0]*(V3[2]-V3[5])+(P3[1]*(-V3[3]+1j*(V3[4]))+(P3[2]*(-1)*(+1j*(V3[3])+V3[4])+P3[3]*(V3[2]-V3[5]))))+F2[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4]))))))+(F1[3]*(F2[2]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F2[3]*(P3[0]*(V3[2]+V3[5])+(P3[1]*(-1)*(V3[3]+1j*(V3[4]))+(P3[2]*(+1j*(V3[3])-V3[4])-P3[3]*(V3[2]+V3[5])))))+(F1[4]*(F2[4]*(P3[0]*(-1)*(V3[2]+V3[5])+(P3[1]*(V3[3]-1j*(V3[4]))+(P3[2]*(+1j*(V3[3])+V3[4])+P3[3]*(V3[2]+V3[5]))))+F2[5]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4]))))))+F1[5]*(F2[4]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F2[5]*(P3[0]*(-V3[2]+V3[5])+(P3[1]*(V3[3]+1j*(V3[4]))+(P3[2]*(-1j*(V3[3])+V3[4])+P3[3]*(-V3[2]+V3[5]))))))))
	vertex = COUP*(+1j*(TMP7)-1j*(TMP8))
	return vertex

def FFV1_1(F2,V3,COUP,M1,W1):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#Gamma(3,2,1)
	F1 = WaveFunction(size=6)
	F1[0] = +F2[0]+V3[0]
	F1[1] = +F2[1]+V3[1]
	P1 = [-complex(F1[0]).real, -complex(F1[1]).real, -complex(F1[1]).imag, -complex(F1[0]).imag]
	denom = COUP/(P1[0]**2-P1[1]**2-P1[2]**2-P1[3]**2 - M1 * (M1 -1j* W1))
	F1[2]= denom*1j*(F2[2]*(P1[0]*(V3[5]-V3[2])+(P1[1]*(V3[3]-1j*(V3[4]))+(P1[2]*(V3[4]+1j*(V3[3]))+P1[3]*(V3[5]-V3[2]))))+(F2[3]*(P1[0]*(V3[3]+1j*(V3[4]))+(P1[1]*(-1)*(V3[2]+V3[5])+(P1[2]*(-1)*(+1j*(V3[2]+V3[5]))+P1[3]*(V3[3]+1j*(V3[4])))))+M1*(F2[4]*(V3[2]+V3[5])+F2[5]*(V3[3]+1j*(V3[4])))))
	F1[3]= denom*(-1j)*(F2[2]*(P1[0]*(+1j*(V3[4])-V3[3])+(P1[1]*(V3[2]-V3[5])+(P1[2]*(-1j*(V3[2])+1j*(V3[5]))+P1[3]*(V3[3]-1j*(V3[4])))))+(F2[3]*(P1[0]*(V3[2]+V3[5])+(P1[1]*(-1)*(V3[3]+1j*(V3[4]))+(P1[2]*(+1j*(V3[3])-V3[4])-P1[3]*(V3[2]+V3[5]))))+M1*(F2[4]*(+1j*(V3[4])-V3[3])+F2[5]*(V3[5]-V3[2]))))
	F1[4]= denom*(-1j)*(F2[4]*(P1[0]*(V3[2]+V3[5])+(P1[1]*(+1j*(V3[4])-V3[3])+(P1[2]*(-1)*(V3[4]+1j*(V3[3]))-P1[3]*(V3[2]+V3[5]))))+(F2[5]*(P1[0]*(V3[3]+1j*(V3[4]))+(P1[1]*(V3[5]-V3[2])+(P1[2]*(-1j*(V3[2])+1j*(V3[5]))-P1[3]*(V3[3]+1j*(V3[4])))))+M1*(F2[2]*(V3[5]-V3[2])+F2[3]*(V3[3]+1j*(V3[4])))))
	F1[5]= denom*1j*(F2[4]*(P1[0]*(+1j*(V3[4])-V3[3])+(P1[1]*(V3[2]+V3[5])+(P1[2]*(-1)*(+1j*(V3[2]+V3[5]))+P1[3]*(+1j*(V3[4])-V3[3]))))+(F2[5]*(P1[0]*(V3[5]-V3[2])+(P1[1]*(V3[3]+1j*(V3[4]))+(P1[2]*(V3[4]-1j*(V3[3]))+P1[3]*(V3[5]-V3[2]))))+M1*(F2[2]*(+1j*(V3[4])-V3[3])+F2[3]*(V3[2]+V3[5]))))
	return F1

def FFV2_1(F2,V3,COUP,M1,W1):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#-(P(-1,3)*Gamma(-1,2,-2)*Gamma(3,-2,1)) + P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2)

	P3 = [complex(V3[0]).real, complex(V3[1]).real, complex(V3[1]).imag, complex(V3[0]).imag]
	F1 = WaveFunction(size=6)
	F1[0] = +F2[0]+V3[0]
	F1[1] = +F2[1]+V3[1]
	P1 = [-complex(F1[0]).real, -complex(F1[1]).real, -complex(F1[1]).imag, -complex(F1[0]).imag]
	denom = COUP/(P1[0]**2-P1[1]**2-P1[2]**2-P1[3]**2 - M1 * (M1 -1j* W1))
	F1[2]= denom*2j*(F2[4]*(P1[1]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[2]*(P3[0]*(-1)*(+1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2]+V3[5]))+(P3[2]*(V3[2]+V3[5])-P3[3]*(+1j*(V3[3])+V3[4]))))+(P1[0]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+(F2[5]*(P1[0]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[1]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M1*(F2[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+F2[2]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	F1[3]= denom*2j*(F2[4]*(P1[0]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P1[1]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F2[5]*(P1[1]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[2]*(P3[0]*(+1j*(V3[3])-V3[4])+(P3[1]*(-1j*(V3[2])+1j*(V3[5]))+(P3[2]*(V3[2]-V3[5])+P3[3]*(-1j*(V3[3])+V3[4]))))+(P1[0]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+M1*(F2[2]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F2[3]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	F1[4]= denom*2j*(F2[2]*(P1[1]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P1[2]*(P3[0]*(-1)*(+1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2])-1j*(V3[5]))+(P3[2]*(V3[2]-V3[5])+P3[3]*(+1j*(V3[3])+V3[4]))))+(P1[0]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+(F2[3]*(P1[0]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P1[1]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M1*(F2[5]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+F2[4]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	F1[5]= denom*2j*(F2[2]*(P1[0]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[1]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F2[3]*(P1[1]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P1[2]*(P3[0]*(+1j*(V3[3])-V3[4])+(P3[1]*(-1)*(+1j*(V3[2]+V3[5]))+(P3[2]*(V3[2]+V3[5])+P3[3]*(+1j*(V3[3])-V3[4]))))+(P1[0]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+M1*(F2[4]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F2[5]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	return F1

def FFV3_1(F2,V3,COUP,M1,W1):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#-(P(-1,3)*Gamma5(-2,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)) + P(-1,3)*Gamma5(-2,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)

	P3 = [complex(V3[0]).real, complex(V3[1]).real, complex(V3[1]).imag, complex(V3[0]).imag]
	F1 = WaveFunction(size=6)
	F1[0] = +F2[0]+V3[0]
	F1[1] = +F2[1]+V3[1]
	P1 = [-complex(F1[0]).real, -complex(F1[1]).real, -complex(F1[1]).imag, -complex(F1[0]).imag]
	denom = COUP/(P1[0]**2-P1[1]**2-P1[2]**2-P1[3]**2 - M1 * (M1 -1j* W1))
	F1[2]= denom*2j*(F2[4]*(P1[1]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[2]*(P3[0]*(-1)*(+1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2]+V3[5]))+(P3[2]*(V3[2]+V3[5])-P3[3]*(+1j*(V3[3])+V3[4]))))+(P1[0]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+(F2[5]*(P1[0]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[1]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M1*(F2[3]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+F2[2]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	F1[3]= denom*2j*(F2[4]*(P1[0]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P1[1]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F2[5]*(P1[1]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[2]*(P3[0]*(+1j*(V3[3])-V3[4])+(P3[1]*(-1j*(V3[2])+1j*(V3[5]))+(P3[2]*(V3[2]-V3[5])+P3[3]*(-1j*(V3[3])+V3[4]))))+(P1[0]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+M1*(F2[2]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F2[3]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	F1[4]= denom*(-2j)*(F2[2]*(P1[1]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P1[2]*(P3[0]*(-1)*(+1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2])-1j*(V3[5]))+(P3[2]*(V3[2]-V3[5])+P3[3]*(+1j*(V3[3])+V3[4]))))+(P1[0]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+(F2[3]*(P1[0]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P1[1]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M1*(F2[5]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+F2[4]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	F1[5]= denom*(-2j)*(F2[2]*(P1[0]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P1[1]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F2[3]*(P1[1]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P1[2]*(P3[0]*(+1j*(V3[3])-V3[4])+(P3[1]*(-1)*(+1j*(V3[2]+V3[5]))+(P3[2]*(V3[2]+V3[5])+P3[3]*(+1j*(V3[3])-V3[4]))))+(P1[0]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P1[3]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))+M1*(F2[4]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F2[5]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	return F1

def FFV1_2(F1,V3,COUP,M2,W2):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#Gamma(3,2,1)
	F2 = WaveFunction(size=6)
	F2[0] = +F1[0]+V3[0]
	F2[1] = +F1[1]+V3[1]
	P2 = [-complex(F2[0]).real, -complex(F2[1]).real, -complex(F2[1]).imag, -complex(F2[0]).imag]
	denom = COUP/(P2[0]**2-P2[1]**2-P2[2]**2-P2[3]**2 - M2 * (M2 -1j* W2))
	F2[2]= denom*1j*(F1[2]*(P2[0]*(V3[2]+V3[5])+(P2[1]*(-1)*(V3[3]+1j*(V3[4]))+(P2[2]*(+1j*(V3[3])-V3[4])-P2[3]*(V3[2]+V3[5]))))+(F1[3]*(P2[0]*(V3[3]-1j*(V3[4]))+(P2[1]*(V3[5]-V3[2])+(P2[2]*(-1j*(V3[5])+1j*(V3[2]))+P2[3]*(+1j*(V3[4])-V3[3]))))+M2*(F1[4]*(V3[2]-V3[5])+F1[5]*(+1j*(V3[4])-V3[3]))))
	F2[3]= denom*(-1j)*(F1[2]*(P2[0]*(-1)*(V3[3]+1j*(V3[4]))+(P2[1]*(V3[2]+V3[5])+(P2[2]*(+1j*(V3[2]+V3[5]))-P2[3]*(V3[3]+1j*(V3[4])))))+(F1[3]*(P2[0]*(V3[5]-V3[2])+(P2[1]*(V3[3]-1j*(V3[4]))+(P2[2]*(V3[4]+1j*(V3[3]))+P2[3]*(V3[5]-V3[2]))))+M2*(F1[4]*(V3[3]+1j*(V3[4]))-F1[5]*(V3[2]+V3[5]))))
	F2[4]= denom*(-1j)*(F1[4]*(P2[0]*(V3[5]-V3[2])+(P2[1]*(V3[3]+1j*(V3[4]))+(P2[2]*(V3[4]-1j*(V3[3]))+P2[3]*(V3[5]-V3[2]))))+(F1[5]*(P2[0]*(V3[3]-1j*(V3[4]))+(P2[1]*(-1)*(V3[2]+V3[5])+(P2[2]*(+1j*(V3[2]+V3[5]))+P2[3]*(V3[3]-1j*(V3[4])))))+M2*(F1[2]*(-1)*(V3[2]+V3[5])+F1[3]*(+1j*(V3[4])-V3[3]))))
	F2[5]= denom*1j*(F1[4]*(P2[0]*(-1)*(V3[3]+1j*(V3[4]))+(P2[1]*(V3[2]-V3[5])+(P2[2]*(-1j*(V3[5])+1j*(V3[2]))+P2[3]*(V3[3]+1j*(V3[4])))))+(F1[5]*(P2[0]*(V3[2]+V3[5])+(P2[1]*(+1j*(V3[4])-V3[3])+(P2[2]*(-1)*(V3[4]+1j*(V3[3]))-P2[3]*(V3[2]+V3[5]))))+M2*(F1[2]*(V3[3]+1j*(V3[4]))+F1[3]*(V3[2]-V3[5]))))
	return F2

def FFV2_2(F1,V3,COUP,M2,W2):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#-(P(-1,3)*Gamma(-1,2,-2)*Gamma(3,-2,1)) + P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2)

	P3 = [complex(V3[0]).real, complex(V3[1]).real, complex(V3[1]).imag, complex(V3[0]).imag]
	F2 = WaveFunction(size=6)
	F2[0] = +F1[0]+V3[0]
	F2[1] = +F1[1]+V3[1]
	P2 = [-complex(F2[0]).real, -complex(F2[1]).real, -complex(F2[1]).imag, -complex(F2[0]).imag]
	denom = COUP/(P2[0]**2-P2[1]**2-P2[2]**2-P2[3]**2 - M2 * (M2 -1j* W2))
	F2[2]= denom*(-2j)*(F1[4]*(P2[1]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P2[2]*(P3[0]*(-1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2])-1j*(V3[5]))+(P3[2]*(-V3[2]+V3[5])+P3[3]*(+1j*(V3[3])-V3[4]))))+(P2[0]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+(F1[5]*(P2[0]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P2[1]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M2*(F1[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F1[2]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	F2[3]= denom*(-2j)*(F1[4]*(P2[0]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[1]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F1[5]*(P2[1]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P2[2]*(P3[0]*(+1j*(V3[3])+V3[4])+(P3[1]*(-1)*(+1j*(V3[2]+V3[5]))+(P3[2]*(-1)*(V3[2]+V3[5])+P3[3]*(+1j*(V3[3])+V3[4]))))+(P2[0]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+M2*(F1[2]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+F1[3]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	F2[4]= denom*(-2j)*(F1[2]*(P2[1]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[2]*(P3[0]*(-1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2]+V3[5]))+(P3[2]*(-1)*(V3[2]+V3[5])+P3[3]*(-1j*(V3[3])+V3[4]))))+(P2[0]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+(F1[3]*(P2[0]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[1]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M2*(F1[5]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+F1[4]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	F2[5]= denom*(-2j)*(F1[2]*(P2[0]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P2[1]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F1[3]*(P2[1]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[2]*(P3[0]*(+1j*(V3[3])+V3[4])+(P3[1]*(-1j*(V3[2])+1j*(V3[5]))+(P3[2]*(-V3[2]+V3[5])-P3[3]*(+1j*(V3[3])+V3[4]))))+(P2[0]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+M2*(F1[4]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+F1[5]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	return F2

def FFV3_2(F1,V3,COUP,M2,W2):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#-(P(-1,3)*Gamma5(-2,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)) + P(-1,3)*Gamma5(-2,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)

	P3 = [complex(V3[0]).real, complex(V3[1]).real, complex(V3[1]).imag, complex(V3[0]).imag]
	F2 = WaveFunction(size=6)
	F2[0] = +F1[0]+V3[0]
	F2[1] = +F1[1]+V3[1]
	P2 = [-complex(F2[0]).real, -complex(F2[1]).real, -complex(F2[1]).imag, -complex(F2[0]).imag]
	denom = COUP/(P2[0]**2-P2[1]**2-P2[2]**2-P2[3]**2 - M2 * (M2 -1j* W2))
	F2[2]= denom*(-2j)*(F1[4]*(P2[1]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P2[2]*(P3[0]*(-1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2])-1j*(V3[5]))+(P3[2]*(-V3[2]+V3[5])+P3[3]*(+1j*(V3[3])-V3[4]))))+(P2[0]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+(F1[5]*(P2[0]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P2[1]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M2*(F1[3]*(P3[0]*(-V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F1[2]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	F2[3]= denom*(-2j)*(F1[4]*(P2[0]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]-V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[1]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])+V3[3]*P3[2]-V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F1[5]*(P2[1]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+(P2[2]*(P3[0]*(+1j*(V3[3])+V3[4])+(P3[1]*(-1)*(+1j*(V3[2]+V3[5]))+(P3[2]*(-1)*(V3[2]+V3[5])+P3[3]*(+1j*(V3[3])+V3[4]))))+(P2[0]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+M2*(F1[2]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+F1[3]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	F2[4]= denom*2j*(F1[2]*(P2[1]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[2]*(P3[0]*(-1j*(V3[3])+V3[4])+(P3[1]*(+1j*(V3[2]+V3[5]))+(P3[2]*(-1)*(V3[2]+V3[5])+P3[3]*(-1j*(V3[3])+V3[4]))))+(P2[0]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+(F1[3]*(P2[0]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[1]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+M2*(F1[5]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]-1j*(V3[4])))))+F1[4]*(-V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])+V3[5]*P3[0]))))
	F2[5]= denom*2j*(F1[2]*(P2[0]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-1)*(V3[2]+V3[5])+(P3[2]*(-1)*(+1j*(V3[2]+V3[5]))+P3[3]*(V3[3]+1j*(V3[4])))))+(P2[3]*(P3[0]*(-1)*(V3[3]+1j*(V3[4]))+(P3[1]*(V3[2]+V3[5])+(P3[2]*(+1j*(V3[2]+V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+(P2[1]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0])+P2[2]*(-1j*(V3[2]*P3[3])-V3[3]*P3[2]+V3[4]*P3[1]+1j*(V3[5]*P3[0])))))+(F1[3]*(P2[1]*(P3[0]*(V3[3]-1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(+1j*(V3[2])-1j*(V3[5]))+P3[3]*(-V3[3]+1j*(V3[4])))))+(P2[2]*(P3[0]*(+1j*(V3[3])+V3[4])+(P3[1]*(-1j*(V3[2])+1j*(V3[5]))+(P3[2]*(-V3[2]+V3[5])-P3[3]*(+1j*(V3[3])+V3[4]))))+(P2[0]*(V3[2]*P3[3]-1j*(V3[3]*P3[2])+1j*(V3[4]*P3[1])-V3[5]*P3[0])+P2[3]*(-V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])+V3[5]*P3[0]))))+M2*(F1[4]*(P3[0]*(V3[3]+1j*(V3[4]))+(P3[1]*(-V3[2]+V3[5])+(P3[2]*(-1j*(V3[2])+1j*(V3[5]))-P3[3]*(V3[3]+1j*(V3[4])))))+F1[5]*(V3[2]*P3[3]+1j*(V3[3]*P3[2])-1j*(V3[4]*P3[1])-V3[5]*P3[0]))))
	return F2

def FFV1_2_3_0(F1,F2,V3,COUP1,COUP2,COUP3):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#Gamma(3,2,1)

	vertex = FFV1_0(F1,F2,V3,COUP1)
	tmp = FFV2_0(F1,F2,V3,COUP2)
	vertex += tmp
	tmp = FFV3_0(F1,F2,V3,COUP3)
	vertex += tmp
	return vertex

def FFV1_2_3_1(F2,V3,COUP1,COUP2,COUP3,M1,W1):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#Gamma(3,2,1)

	F1 = FFV1_1(F2,V3,COUP1,M1,W1)
	tmp = FFV2_1(F2,V3,COUP2,M1,W1)
	for i in range(2,6):
		F1[i] += tmp[i]
	tmp = FFV3_1(F2,V3,COUP3,M1,W1)
	for i in range(2,6):
		F1[i] += tmp[i]
	return F1

def FFV1_2_3_2(F1,V3,COUP1,COUP2,COUP3,M2,W2):
	#This File is Automatically generated by ALOHA 
	#The process calculated in this file is: 
	#Gamma(3,2,1)

	F2 = FFV1_2(F1,V3,COUP1,M2,W2)
	tmp = FFV2_2(F1,V3,COUP2,M2,W2)
	for i in range(2,6):
		F2[i] += tmp[i]
	tmp = FFV3_2(F1,V3,COUP3,M2,W2)
	for i in range(2,6):
		F2[i] += tmp[i]
	return F2

#-------------------------------HELPER FUNCTIONS-------------------------------

def do_z_boost(fourmom, beta):
	"""
	Do a longitudinal boost of a 4-momentum
	Arguments:
		-- fourmom : 4momentum in the [E,px,py,pz] format
		-- beta : boost factor
	"""
	
	# Compute the gamma factor
	gamma = np.sqrt(1.0/(1-beta*beta))

	# Compute the longitudinally boosted E and pz, with label 'prime' as '_pr'
	E_pr  =       gamma*fourmom[0] + -beta*gamma*fourmom[3];
	pz_pr = -beta*gamma*fourmom[0] +       gamma*fourmom[3];

	# Return the boosted 4momentum
	boosted_fourmom = [E_pr, fourmom[1], fourmom[2], pz_pr]
	return(boosted_fourmom)

def compute_matrix_element_sq(V1, V2, F1, F2, COUP1, COUP2, COUP3, mtau):
	"""
	Compute the matrix element squared of the gamma gamma -> tau+ tau- process
	Arguments:
		-- V1, V2 : incoming photon wavefunctions
		-- F1, F2 : outgoing tau lepton wavefunctions
		-- COUP1, COUP2, COUP3 : coupling values associated with the F1, F2, F3 form factors,
		-- mtau : mass of tau lepton
	"""

	# Complex matrix elements of the two channels (don't really know which is t/u channel but this doesn't matter)
	matr_chan1 = FFV1_2_3_0(FFV1_2_3_2(F1, V1, COUP1, COUP2, COUP3, mtau, 0.0), F2, V2, COUP1, COUP2, COUP3)
	matr_chan2 = FFV1_2_3_0(F1, FFV1_2_3_1(F2, V1, COUP1, COUP2, COUP3, mtau, 0.0), V2, COUP1, COUP2, COUP3)

	# Add the two channels, square and return
	matr = matr_chan1 + matr_chan2
	matr2 = abs(matr)**2
	return(matr2)

#---------------------------------REWEIGHTING----------------------------------

def reweight_event(event, args):
	"""
	Function to reweight an event with given old & new values of vertexF2 and vertexF3
	Arguments:
		-- event : object of class 'Event'
		-- args : ArgumentParser object from wherever this function has been called
	Returns:
		-- [event, reweight_factor] : event with updated weight, and the reweight factor
	"""

	# Get the tau mass
	mtau = event.ta1.mass

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# Build the wavefunctions
	wf_y1  = vxxxxx(fourmom_y1_com, 0,      event.y1.hel, -1)
	wf_y2  = vxxxxx(fourmom_y2_com, 0,      event.y2.hel, -1)
	wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, event.ta1.hel, -1)
	wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, event.ta2.hel, +1)

	# Get the old and new couplings, where:
	# COUP1 = -(ee*complex(0,1))
	# COUP2 = -0.25*(ee*f2tau*complex(0,1))/ymtau
	# COUP3 = -0.25*(ee*f3tau)/ymtau

	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1_old = -1j*ee
	COUP2_old = -0.25*(ee*args.vertexf2old*1j)/mtau
	COUP3_old = -0.25*(ee*args.vertexf3old)/mtau

	COUP1_new = -1j*ee
	COUP2_new = -0.25*(ee*args.vertexf2new*1j)/mtau
	COUP3_new = -0.25*(ee*args.vertexf3new)/mtau

	# Compute the old and new matrix element squared in the COM and update the event weight
	old_matr2 = compute_matrix_element_sq(wf_y1, wf_y2, wf_ta1, wf_ta2, COUP1_old, COUP2_old, COUP3_old, mtau)
	new_matr2 = compute_matrix_element_sq(wf_y1, wf_y2, wf_ta1, wf_ta2, COUP1_new, COUP2_new, COUP3_new, mtau)
	reweight_factor = (new_matr2)/(old_matr2)
	new_weight = reweight_factor*event.weight

	# Return both the reweighted event, and the reweighting factor as a list
	return(new_weight, reweight_factor)

#-----------------------------OPTIMAL OBSERVABLES------------------------------

def compute_optobs_atau_fullpol(event):
	"""
	Compute and return the optimal observable for atau (all 4 orders)
	Matrix elements are computed fully polarized
	Arguments:
		-- event : object of class 'Event' described elsewhere
	Returns:
		-- [oo1, oo2, oo3, oo4] : all 4 orders of optimal observable for atau
	"""
	
	mtau = 1.77686 # GeV

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# Build the wave functions
	wf_y1  = vxxxxx(fourmom_y1_com,   0.0,  event.y1.hel, -1)
	wf_y2  = vxxxxx(fourmom_y2_com,   0.0,  event.y2.hel, -1)
	wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, event.ta1.hel, -1)
	wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, event.ta2.hel, +1)

	# Get the old and new couplings, where:
	# COUP1 = -(ee*complex(0,1))
	# COUP2 = -0.25*(ee*f2tau*complex(0,1))/ymtau
	# COUP3 = -0.25*(ee*f3tau)/ymtau
	# For OO, we take dummy f2, f3 = 1.0, 0.0 since the couplings are PURELY multiplicative

	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1 = -1j*ee
	COUP2 = -0.25*(ee*1.0*1j)/mtau
	COUP3 = 0.0

	# Compute individual contributions to the matrix elements in powers of atau and dtau
	# Notation follows Fortran code automatically produced by Madgraph

	# t-channel diagrams
	amp1 = FFV1_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
	amp2 = FFV2_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP2)
	amp3 = FFV1_0(FFV2_2(wf_ta1, wf_y1, COUP2, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
	amp4 = FFV2_0(FFV2_2(wf_ta1, wf_y1, COUP2, mtau, 0.0) , wf_ta2, wf_y2, COUP2)

	# u-channel diagrams
	amp5 = FFV1_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP1)
	amp6 = FFV2_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP2)
	amp7 = FFV1_0(wf_ta1, FFV2_1(wf_ta2, wf_y1, COUP2, mtau, 0.0), wf_y2, COUP1)
	amp8 = FFV2_0(wf_ta1, FFV2_1(wf_ta2, wf_y1, COUP2, mtau, 0.0), wf_y2, COUP2)

	# combine the matrix element terms according to powers of atau
	mat_0 = amp1 + amp5 				# F1^2
	mat_1 = amp2 + amp6 + amp3 + amp7	# F1*F2 + F2*F1
	mat_2 = amp4 + amp8					# F2^2

	# optimal observable for atau^1, atau^2, atau^3, atau^4
	oo1 = 2*np.real(np.conj(mat_0)*mat_1)/abs(mat_0)**2
	oo2 = (2*np.real(np.conj(mat_0)*mat_2) + abs(mat_1)**2)/abs(mat_0)**2
	oo3 = 2*np.real(np.conj(mat_1)*mat_2)/abs(mat_0)**2
	oo4 = abs(mat_2)**2/abs(mat_0)**2

	# return the optimal observables
	return([oo1, oo2, oo3, oo4])

def compute_optobs_atau_inithelavg(event):
	"""
	Compute and return the optimal observable for atau (all 4 orders)
	Matrix elements are computed with average over photon helicities
	Arguments:
		-- event : object of class 'Event' described elsewhere
	Returns:
		-- [oo1, oo2, oo3, oo4] : all 4 orders of optimal observable for atau
	"""
	
	mtau = 1.77686 # GeV

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# Get the old and new couplings, where:
	# COUP1 = -(ee*complex(0,1))
	# COUP2 = -0.25*(ee*f2tau*complex(0,1))/ymtau
	# COUP3 = -0.25*(ee*f3tau)/ymtau
	# For OO, we take dummy f2, f3 = 1.0, 0.0 since the couplings are PURELY multiplicative

	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1 = -1j*ee
	COUP2 = -0.25*(ee*1.0*1j)/mtau
	COUP3 = 0.0

	# List possible photon helicity combinations
	init_hels = [[+1,+1], [+1,-1], [-1,+1], [-1,-1]]

	# Initialize mat^2 (in powers of atau) for summing over photon helicities
	mat_sq = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

	# loop over initial photon helicities
	for hel in init_hels:

		# Build the wave functions
		wf_y1  = vxxxxx(fourmom_y1_com,   0.0,        hel[0], -1)
		wf_y2  = vxxxxx(fourmom_y2_com,   0.0,        hel[1], -1)
		wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, event.ta1.hel, -1)
		wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, event.ta2.hel, +1)

		# Compute individual contributions to the matrix elements in powers of atau and dtau
		# Notation follows Fortran code automatically produced by Madgraph

		# t-channel diagrams
		amp1 = FFV1_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp2 = FFV2_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP2)
		amp3 = FFV1_0(FFV2_2(wf_ta1, wf_y1, COUP2, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp4 = FFV2_0(FFV2_2(wf_ta1, wf_y1, COUP2, mtau, 0.0) , wf_ta2, wf_y2, COUP2)

		# u-channel diagrams
		amp5 = FFV1_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP1)
		amp6 = FFV2_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP2)
		amp7 = FFV1_0(wf_ta1, FFV2_1(wf_ta2, wf_y1, COUP2, mtau, 0.0), wf_y2, COUP1)
		amp8 = FFV2_0(wf_ta1, FFV2_1(wf_ta2, wf_y1, COUP2, mtau, 0.0), wf_y2, COUP2)

		# combine the matrix element terms according to powers of atau
		mat_0 = amp1 + amp5 				# F1^2
		mat_1 = amp2 + amp6 + amp3 + amp7	# F1*F2 + F2*F1
		mat_2 = amp4 + amp8					# F2^2

		# compute the mat^2 terms in powers of atau, add to the averaged quantity
		mat_sq[0] += abs(mat_0)**2
		mat_sq[1] += 2*np.real(np.conj(mat_0)*mat_1)
		mat_sq[2] += 2*np.real(np.conj(mat_0)*mat_2) + abs(mat_1)**2
		mat_sq[3] += 2*np.real(np.conj(mat_1)*mat_2)
		mat_sq[4] += abs(mat_2)**2

	# compute the optimal observables (no need to divide mat^2 by a constant for averaging)
	oo1 = mat_sq[1]/mat_sq[0]
	oo2 = mat_sq[2]/mat_sq[0]
	oo3 = mat_sq[3]/mat_sq[0]
	oo4 = mat_sq[4]/mat_sq[0]

	# return the optimal observables
	return([oo1, oo2, oo3, oo4])

def compute_optobs_atau_avghel(event):
	"""
	Compute and return the optimal observable for atau (all 4 orders)
	Matrix elements are computed with average over all helicities
	Arguments:
		-- event : object of class 'Event' described elsewhere
	Returns:
		-- [oo1, oo2, oo3, oo4] : all 4 orders of optimal observable for atau
	"""
	
	mtau = 1.77686 # GeV

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# Get the old and new couplings, where:
	# COUP1 = -(ee*complex(0,1))
	# COUP2 = -0.25*(ee*f2tau*complex(0,1))/ymtau
	# COUP3 = -0.25*(ee*f3tau)/ymtau
	# For OO, we take dummy f2, f3 = 1.0, 0.0 since the couplings are PURELY multiplicative

	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1 = -1j*ee
	COUP2 = -0.25*(ee*1.0*1j)/mtau
	COUP3 = 0.0

	# List possible helicity combinations
	hels = [ [+1,+1,+1,+1]
		  ,[+1,+1,+1,-1]
		  ,[+1,+1,-1,+1]
		  ,[+1,-1,+1,+1]
		  ,[-1,+1,+1,+1]
		  ,[+1,+1,-1,-1]
		  ,[+1,-1,+1,-1]
		  ,[-1,+1,+1,-1]
		  ,[+1,-1,-1,+1]
		  ,[-1,+1,-1,+1]
		  ,[-1,-1,+1,+1]
		  ,[+1,-1,-1,-1]
		  ,[-1,+1,-1,-1]
		  ,[-1,-1,+1,-1]
		  ,[-1,-1,-1,+1]
		  ,[-1,-1,-1,-1] ]

	# Initialize mat^2 (in powers of atau) for summing over photon helicities
	mat_sq = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

	# loop over initial photon helicities
	for hel in hels:

		# Build the wave functions
		wf_y1  = vxxxxx(fourmom_y1_com,   0.0, hel[0], -1)
		wf_y2  = vxxxxx(fourmom_y2_com,   0.0, hel[1], -1)
		wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, hel[2], -1)
		wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, hel[3], +1)

		# Compute individual contributions to the matrix elements in powers of atau and dtau
		# Notation follows Fortran code automatically produced by Madgraph

		# t-channel diagrams
		amp1 = FFV1_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp2 = FFV2_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP2)
		amp3 = FFV1_0(FFV2_2(wf_ta1, wf_y1, COUP2, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp4 = FFV2_0(FFV2_2(wf_ta1, wf_y1, COUP2, mtau, 0.0) , wf_ta2, wf_y2, COUP2)

		# u-channel diagrams
		amp5 = FFV1_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP1)
		amp6 = FFV2_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP2)
		amp7 = FFV1_0(wf_ta1, FFV2_1(wf_ta2, wf_y1, COUP2, mtau, 0.0), wf_y2, COUP1)
		amp8 = FFV2_0(wf_ta1, FFV2_1(wf_ta2, wf_y1, COUP2, mtau, 0.0), wf_y2, COUP2)

		# combine the matrix element terms according to powers of atau
		mat_0 = amp1 + amp5 				# F1^2
		mat_1 = amp2 + amp6 + amp3 + amp7	# F1*F2 + F2*F1
		mat_2 = amp4 + amp8					# F2^2

		# compute the mat^2 terms in powers of atau, add to the averaged quantity
		mat_sq[0] += abs(mat_0)**2
		mat_sq[1] += 2*np.real(np.conj(mat_0)*mat_1)
		mat_sq[2] += 2*np.real(np.conj(mat_0)*mat_2) + abs(mat_1)**2
		mat_sq[3] += 2*np.real(np.conj(mat_1)*mat_2)
		mat_sq[4] += abs(mat_2)**2

	# compute the optimal observables (no need to divide mat^2 by a constant for averaging)
	oo1 = mat_sq[1]/mat_sq[0]
	oo2 = mat_sq[2]/mat_sq[0]
	oo3 = mat_sq[3]/mat_sq[0]
	oo4 = mat_sq[4]/mat_sq[0]
	
	# return the optimal observables
	return([oo1, oo2, oo3, oo4])

def compute_optobs_dtau_fullpol(event):
	"""
	Compute and return the optimal observable for dtau (all 4 orders)
	Matrix elements are computed fully polarized
	Arguments:
		-- event : object of class 'Event' described elsewhere
	Returns:
		-- [oo1, oo2, oo3, oo4] : all 4 orders of optimal observable for dtau
	"""

	mtau = 1.77686 # GeV

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# Build the wave functions
	wf_y1  = vxxxxx(fourmom_y1_com,   0.0,  event.y1.hel, -1)
	wf_y2  = vxxxxx(fourmom_y2_com,   0.0,  event.y2.hel, -1)
	wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, event.ta1.hel, -1)
	wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, event.ta2.hel, +1)	

	# Get the old and new couplings, where:
	# COUP1 = -(ee*complex(0,1))
	# COUP2 = -0.25*(ee*f2tau*complex(0,1))/ymtau
	# COUP3 = -0.25*(ee*f3tau)/ymtau
	# For OO, we take dummy f2, f3 = 0.0, 1.0 since the couplings are PURELY multiplicative

	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1 = -1j*ee
	COUP2 = 0.0
	COUP3 = -0.25*(ee*1.0)/mtau

	# Compute individual contributions to the matrix elements in powers of atau and dtau
	# Notation follows Fortran code automatically produced by Madgraph

	# t-channel diagrams
	amp1 = FFV1_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
	amp2 = FFV3_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP3)
	amp3 = FFV1_0(FFV3_2(wf_ta1, wf_y1, COUP3, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
	amp4 = FFV3_0(FFV3_2(wf_ta1, wf_y1, COUP3, mtau, 0.0) , wf_ta2, wf_y2, COUP3)

	# u-channel diagrams
	amp5 = FFV1_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP1)
	amp6 = FFV3_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP3)
	amp7 = FFV1_0(wf_ta1, FFV3_1(wf_ta2, wf_y1, COUP3, mtau, 0.0), wf_y2, COUP1)
	amp8 = FFV3_0(wf_ta1, FFV3_1(wf_ta2, wf_y1, COUP3, mtau, 0.0), wf_y2, COUP3)

	# combine the matrix element terms according to powers of dtau
	mat_0 = amp1 + amp5 				# F1^2
	mat_1 = amp2 + amp6 + amp3 + amp7	# F1*F3 + F3*F1
	mat_2 = amp4 + amp8					# F3^2

	# optimal observable for dtau^1, dtau^2, dtau^3, dtau^4
	oo1 = 2*np.real(np.conj(mat_0)*mat_1)/abs(mat_0)**2
	oo2 = (2*np.real(np.conj(mat_0)*mat_2) + abs(mat_1)**2)/abs(mat_0)**2
	oo3 = 2*np.real(np.conj(mat_1)*mat_2)/abs(mat_0)**2
	oo4 = abs(mat_2)**2/abs(mat_0)

	# return the optimal observables
	return([oo1, oo2, oo3, oo4])

def compute_optobs_dtau_inithelavg(event):
	"""
	Compute and return the optimal observable for dtau (all 4 orders)
	Matrix elements are computed with average over photon helicities
	Arguments:
		-- event : object of class 'Event' described elsewhere
	Returns:
		-- [oo1, oo2, oo3, oo4] : all 4 orders of optimal observable for dtau
	"""
	
	mtau = 1.77686 # GeV

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# Get the old and new couplings, where:
	# COUP1 = -(ee*complex(0,1))
	# COUP2 = -0.25*(ee*f2tau*complex(0,1))/ymtau
	# COUP3 = -0.25*(ee*f3tau)/ymtau
	# For OO, we take dummy f2, f3 = 0.0, 1.0 since the couplings are PURELY multiplicative

	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1 = -1j*ee
	COUP2 = 0.0
	COUP3 = -0.25*(ee*1.0)/mtau

	# List possible photon helicity combinations
	init_hels = [[+1,+1], [+1,-1], [-1,+1], [-1,-1]]

	# Initialize mat^2 (in powers of dtau) for summing over photon helicities
	mat_sq = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

	# loop over initial photon helicities
	for hel in init_hels:

		# Build the wave functions
		wf_y1  = vxxxxx(fourmom_y1_com,   0.0,        hel[0], -1)
		wf_y2  = vxxxxx(fourmom_y2_com,   0.0,        hel[1], -1)
		wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, event.ta1.hel, -1)
		wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, event.ta2.hel, +1)

		# Compute individual contributions to the matrix elements in powers of atau and dtau
		# Notation follows Fortran code automatically produced by Madgraph

		# t-channel diagrams
		amp1 = FFV1_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp2 = FFV3_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP3)
		amp3 = FFV1_0(FFV3_2(wf_ta1, wf_y1, COUP3, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp4 = FFV3_0(FFV3_2(wf_ta1, wf_y1, COUP3, mtau, 0.0) , wf_ta2, wf_y2, COUP3)

		# u-channel diagrams
		amp5 = FFV1_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP1)
		amp6 = FFV3_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP3)
		amp7 = FFV1_0(wf_ta1, FFV3_1(wf_ta2, wf_y1, COUP3, mtau, 0.0), wf_y2, COUP1)
		amp8 = FFV3_0(wf_ta1, FFV3_1(wf_ta2, wf_y1, COUP3, mtau, 0.0), wf_y2, COUP3)

		# combine the matrix element terms according to powers of dtau
		mat_0 = amp1 + amp5 				# F1^2
		mat_1 = amp2 + amp6 + amp3 + amp7	# F1*F3 + F3*F1
		mat_2 = amp4 + amp8					# F3^2

		# compute the mat^2 terms in powers of dtau, add to the averaged quantity
		mat_sq[0] += abs(mat_0)**2
		mat_sq[1] += 2*np.real(np.conj(mat_0)*mat_1)
		mat_sq[2] += 2*np.real(np.conj(mat_0)*mat_2) + abs(mat_1)**2
		mat_sq[3] += 2*np.real(np.conj(mat_1)*mat_2)
		mat_sq[4] += abs(mat_2)**2

	# compute the optimal observables (no need to divide mat^2 by a constant for averaging)
	oo1 = mat_sq[1]/mat_sq[0]
	oo2 = mat_sq[2]/mat_sq[0]
	oo3 = mat_sq[3]/mat_sq[0]
	oo4 = mat_sq[4]/mat_sq[0]

	# return the optimal observables
	return([oo1, oo2, oo3, oo4])

def compute_optobs_dtau_avghel(event):
	"""
	Compute and return the optimal observable for dtau (all 4 orders)
	Matrix elements are computed with average over all helicities
	Arguments:
		-- event : object of class 'Event' described elsewhere
	Returns:
		-- [oo1, oo2, oo3, oo4] : all 4 orders of optimal observable for dtau
	"""
	
	mtau = 1.77686 # GeV

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# Get the old and new couplings, where:
	# COUP1 = -(ee*complex(0,1))
	# COUP2 = -0.25*(ee*f2tau*complex(0,1))/ymtau
	# COUP3 = -0.25*(ee*f3tau)/ymtau
	# For OO, we take dummy f2, f3 = 0.0, 1.0 since the couplings are PURELY multiplicative

	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1 = -1j*ee
	COUP2 = 0.0
	COUP3 = -0.25*(ee*1.0)/mtau

	# List possible helicity combinations
	hels = [ [+1,+1,+1,+1]
		  ,[+1,+1,+1,-1]
		  ,[+1,+1,-1,+1]
		  ,[+1,-1,+1,+1]
		  ,[-1,+1,+1,+1]
		  ,[+1,+1,-1,-1]
		  ,[+1,-1,+1,-1]
		  ,[-1,+1,+1,-1]
		  ,[+1,-1,-1,+1]
		  ,[-1,+1,-1,+1]
		  ,[-1,-1,+1,+1]
		  ,[+1,-1,-1,-1]
		  ,[-1,+1,-1,-1]
		  ,[-1,-1,+1,-1]
		  ,[-1,-1,-1,+1]
		  ,[-1,-1,-1,-1] ]

	# Initialize mat^2 (in powers of dtau) for summing over photon helicities
	mat_sq = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

	# loop over initial photon helicities
	for hel in hels:

		# Build the wave functions
		wf_y1  = vxxxxx(fourmom_y1_com,   0.0, hel[0], -1)
		wf_y2  = vxxxxx(fourmom_y2_com,   0.0, hel[1], -1)
		wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, hel[2], -1)
		wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, hel[3], +1)

		# Compute individual contributions to the matrix elements in powers of atau and dtau
		# Notation follows Fortran code automatically produced by Madgraph

		# t-channel diagrams
		amp1 = FFV1_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp2 = FFV3_0(FFV1_2(wf_ta1, wf_y1, COUP1, mtau, 0.0) , wf_ta2, wf_y2, COUP3)
		amp3 = FFV1_0(FFV3_2(wf_ta1, wf_y1, COUP3, mtau, 0.0) , wf_ta2, wf_y2, COUP1)
		amp4 = FFV3_0(FFV3_2(wf_ta1, wf_y1, COUP3, mtau, 0.0) , wf_ta2, wf_y2, COUP3)

		# u-channel diagrams
		amp5 = FFV1_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP1)
		amp6 = FFV3_0(wf_ta1, FFV1_1(wf_ta2, wf_y1, COUP1, mtau, 0.0), wf_y2, COUP3)
		amp7 = FFV1_0(wf_ta1, FFV3_1(wf_ta2, wf_y1, COUP3, mtau, 0.0), wf_y2, COUP1)
		amp8 = FFV3_0(wf_ta1, FFV3_1(wf_ta2, wf_y1, COUP3, mtau, 0.0), wf_y2, COUP3)

		# combine the matrix element terms according to powers of dtau
		mat_0 = amp1 + amp5 				# F1^2
		mat_1 = amp2 + amp6 + amp3 + amp7	# F1*F3 + F3*F1
		mat_2 = amp4 + amp8					# F3^2

		# compute the mat^2 terms in powers of dtau, add to the averaged quantity
		mat_sq[0] += abs(mat_0)**2
		mat_sq[1] += 2*np.real(np.conj(mat_0)*mat_1)
		mat_sq[2] += 2*np.real(np.conj(mat_0)*mat_2) + abs(mat_1)**2
		mat_sq[3] += 2*np.real(np.conj(mat_1)*mat_2)
		mat_sq[4] += abs(mat_2)**2

	# compute the optimal observables (no need to divide mat^2 by a constant for averaging)
	oo1 = mat_sq[1]/mat_sq[0]
	oo2 = mat_sq[2]/mat_sq[0]
	oo3 = mat_sq[3]/mat_sq[0]
	oo4 = mat_sq[4]/mat_sq[0]
	
	# return the optimal observables
	return([oo1, oo2, oo3, oo4])

#---------------------------------APPROXIMATION--------------------------------

def check_helicities_equality(hels1, hels2):
	"""
	Verify if 2 sets of helicities for [y1,y2,ta1,ta2] are the same either directly or 
	via matrix element symmetries. All symmetries hardcoded based on SM matrix element
	I.E. :  ++++,+-++,-+++,++--,--++,+---,-+--,---- = 0
	+++- = ---+ ,++-+ = --+- ,+-+- = -+-+ ,-++- = +--+      ---> equal and opposite

	"""

	# convert to numpy arrays just to be safe
	hels1, hels2 = np.array(hels1), np.array(hels2)

	# check if they are directly equal
	if (hels1==hels2).all(): return(True)

	# check if they are related by minus sign
	if (hels1==-hels2).all(): return(True)

	# all other cases: false
	else: return(False)

def guess_helicities(event):
	"""
	Guess the helicity for an event, given lab 4-momenta of the event.
	Matrix element sq. are computed for all helicity combinations, the one with max probability is returned.
	Couplings assumed SM (atau, dtau = 0)
	"""

	# Get the tau mass
	mtau = event.ta1.mass

	# get the 4momenta in the "lab" frame
	fourmom_y1_lab  = [event.y1.e,  event.y1.px,  event.y1.py,  event.y1.pz ]
	fourmom_y2_lab  = [event.y2.e,  event.y2.px,  event.y2.py,  event.y2.pz ]
	fourmom_ta1_lab = [event.ta1.e, event.ta1.px, event.ta1.py, event.ta1.pz]
	fourmom_ta2_lab = [event.ta2.e, event.ta2.px, event.ta2.py, event.ta2.pz]

	# Compute the boost parameter boost=(E1-E2)/(E1+E2) required to longituinally boost into COM frame
	e1, e2 = fourmom_y1_lab[0], fourmom_y2_lab[0]
	boost = (e1-e2)/(e1+e2)

	# boost the 4momenta of the particles into the COM frame
	fourmom_y1_com  = do_z_boost(fourmom_y1_lab,  boost)
	fourmom_y2_com  = do_z_boost(fourmom_y2_lab,  boost)
	fourmom_ta1_com = do_z_boost(fourmom_ta1_lab, boost)
	fourmom_ta2_com = do_z_boost(fourmom_ta2_lab, boost)

	# define the couplings
	ee = np.sqrt(4*np.pi*7.297352e-3)

	COUP1 = -1j*ee
	COUP2 = 0.0
	COUP3 = 0.0

	# List possible helicity combinations
	hels = [ [+1,+1,+1,+1]
		  ,[+1,+1,+1,-1]
		  ,[+1,+1,-1,+1]
		  ,[+1,-1,+1,+1]
		  ,[-1,+1,+1,+1]
		  ,[+1,+1,-1,-1]
		  ,[+1,-1,+1,-1]
		  ,[-1,+1,+1,-1]
		  ,[+1,-1,-1,+1]
		  ,[-1,+1,-1,+1]
		  ,[-1,-1,+1,+1]
		  ,[+1,-1,-1,-1]
		  ,[-1,+1,-1,-1]
		  ,[-1,-1,+1,-1]
		  ,[-1,-1,-1,+1]
		  ,[-1,-1,-1,-1] ]

	# compute the matrix element sq. for every helicity combination
	# no need to convert to probability i.e. M^2{hels}/M^2{total} since M^2{total} is constant for every eveny
	matr2_hels = []
	for hel in hels:

		# Build the wavefunctions
		wf_y1  = vxxxxx(fourmom_y1_com, 0,     hel[0], -1)
		wf_y2  = vxxxxx(fourmom_y2_com, 0,     hel[1], -1)
		wf_ta1 = ixxxxx(fourmom_ta1_com, mtau, hel[2], -1)
		wf_ta2 = oxxxxx(fourmom_ta2_com, mtau, hel[3], +1)

		# compute the matrix element squared, get the maximimized helicity combo
		matr2 = compute_matrix_element_sq(wf_y1, wf_y2, wf_ta1, wf_ta2, COUP1, COUP2, COUP3, mtau)
		matr2_hels.append(matr2)
	
	# return the helicity combination that maximizes the matrix element
	return(hels[matr2_hels.index(max(matr2_hels))])