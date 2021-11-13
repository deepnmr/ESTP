#!/usr/bin/env python
import numpy as np
from sys import stderr
from numpy.random import randn
from numpy import sqrt

from estp.est_model_noex import est_model

m2=est_model()

fields = [70.12]
#Ts = [0.4, 0.6]
v1s = [15.0]
Ts = [0.4]
#v1s = [20.0]

errorinv1 = 0.1
v1errs = []
for v1 in v1s:
   v1errs.append(v1*errorinv1)

#kex = 200
#pB = 0.1
#kab = pB*kex
#kba = kex-kab
dGs=[120.0]
#dws=[-5.0]
#dGs = [115.0, 115.0,115.0,115.0,115.0,115.0,115.0,115.0,115.0,115.0]
#dws = [-10.0,-9.0, -8.0, -7.0, -6.0,-5.0, -4.0, -3.0, -2.0, -1.0]
r1 = 0.5
r2a = 20.0
#r2b = 100.0
#dR = r2b

center = 120.0
totalppm = 20.0
numofdata = 50


res = []
res.append({'name' : 'A1', 'dG' : dGs[0]})#, 'dw' : dws[0]})
#res.append({'name' : 'A2', 'dG' : dGs[1], 'dw' : dws[1]})
#res.append({'name' : 'A3', 'dG' : dGs[2], 'dw' : dws[2]})
#res.append({'name' : 'A4', 'dG' : dGs[3], 'dw' : dws[3]})
#res.append({'name' : 'A5', 'dG' : dGs[4], 'dw' : dws[4]})
#res.append({'name' : 'A6', 'dG' : dGs[5], 'dw' : dws[5]})
#res.append({'name' : 'A7', 'dG' : dGs[6], 'dw' : dws[6]})
#res.append({'name' : 'A8', 'dG' : dGs[7], 'dw' : dws[7]})
#res.append({'name' : 'A9', 'dG' : dGs[8], 'dw' : dws[8]})
#res.append({'name' : 'A10', 'dG' : dGs[9], 'dw' : dws[9]})



offset = np.linspace(center-totalppm/2,center+totalppm/2.0,numofdata)
noise = 50.0
signal = 200.0*50.0
#noiseRatio = 0.02

for i in range(0,len(fields)):
   field = fields[i]
   for ii in range(0,len(Ts)):
      T = Ts[ii]
      for iii in range(0,len(v1s)):
         v1 = v1s[iii]
         v1err = v1errs[iii]
         buff = '%5.2f\n' % field
         buff += '%f\n' % T
         buff += '%5.2f %5.2f\n' % (v1, v1err)
         buff += '#offset(ppm)     Intensity     error\n'
         for r in res:
	         buff += '# %s R2a: %2.1f R2b: %2.1f dw: %2.1f\n' % (r['name'],r2a+5.0,0.0,0.0)
	         dG = r['dG'] 
	         #dw = r['dw']
	         for of in offset:
#	            intideal = m2.Baldwin(kab,kba,dG,dw,r1,r2a,dR,of,T,v1,v1err,field)
	            intideal = m2.matrix_calc_noex(dG,r1,r2a,of,T,v1,v1err,field)
	            expsignal = intideal*signal
#	            intstd = intideal*sqrt((noise/signal)**2) 
	            intstd = 0.006 
	            appInt = intideal + intstd*randn()
	            buff += '%8.3f %8.3f %8.3f\n' % (of, appInt, abs(intstd))
         #filex = open('../cest/in/synth-%d-t%d-v%5.2f.txt' % (int(field),int(T*1000),v1), 'w')
         filex = open('synth-%d-t%d-v%5.2f.txt' % (int(field),int(T*1000),v1), 'w')
         filex.write(buff)
         filex.close()
