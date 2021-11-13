#!/usr/bin/env python

######
# 2015 by Donghan Lee
########
# modified by Marta Carneiro to output a file per field and residue;
# modified by Marta Carneiro to read error, R2s and dW from input files and output the same values
# modified by Marta Carneiro to invert the CEST profile (i.e., dips become peaks) before to FT/LP, and revert to CEST profile after FT

from sys import argv, stdout
from estp.est_model import est_model
from estp.lp import lp, lp_1d
from estp.proc_base import fft, ifft, rft, irft, zf
#from scipy.fftpack import fft, ifft, rfft, irfft
from numpy import abs, real, append, max, array, min, linspace, sqrt
import json
import matplotlib.pyplot as plt

m2=est_model()
m2.verbose = True

print '**************'
stdout.write(m2.programName + '\n')
print '**************'

# define datasets

configFile = open(argv[1])
conf = json.load(configFile)
configFile.close()

projectName = conf['Project Name']
datasetsNames = conf['datasets']

for dataset in datasetsNames:
    m2.dataset.addData(dataset)


residues = conf['residues']


for r in residues:
    resid = r['name']
    active = r['flag']
    print 'Residue: ' + resid + ' ' + active
    for r in m2.dataset.res:
        if r.label == resid:
	    if active == 'on':
	        r.active = True
	    elif active == 'off':
	        r.active = False
	    else:
	        print 'Error: wrong flag for the residue' + resid + ' ' + active
		exit()

m2.datapdf('data.pdf')

#for data in m2.dataset:
for r in m2.dataset.res:
    if(r.active):
       for i in range(0, len(r.estSpecs)):
	   ep = r.estSpecs[i]
	   buf = '%8.3f\n' % ep.field
	   buf += '%8.3f\n' % ep.T
	   buf += '%8.3f %8.3f\n' % (ep.v1, ep.v1err)
	   buf += '# offset[ppm] intensity interr\n'
	   
	   buf += '# %s R2a: %8.3f R2b: %8.3f dw: %8.3f\n' % (r.label,ep.initr2a,ep.initr2b,ep.initdw)
	   epint = array(ep.int)
	   
	   #epintmax = max(epint)
	   
	   #epintinv = []
	   #for j in range(0,len(epint)):
	   # epintinv.append(epintmax-epint[j])
	   #epintinv = epintmax-epint
	   
	   #epintinva = array(epintinv)

	   
#	   intmax = max(epint)
	   
#	   epinta = []
#	   for i in range(0,len(epint)):
#	        epinta.append(intmax - epint[i])
	   
	   
#	   epinta = array(epinta)
	   
	   offsetmax = max(ep.offset)
	   offsetmin = min(ep.offset)
	   offsetlen = len(ep.offset)
	   npred = offsetlen
	   
	   numoffset = offsetlen + npred
	   
	   print numoffset
	   offsets = linspace(offsetmin,offsetmax + (ep.offset[1] - ep.offset[0])/2,numoffset)
	   
	   inttime = irft(epint)
	   #inttimexx = irft(epint)
	   #
	   #plt.plot(inttime)
	   #plt.plot(inttimexx)
	   #plt.show()

	   lpinttime = lp(inttime,pred=npred,order=2)
	   #lpinttime = inttime
#	   zflpinttime = zf(lpinttime, pad = (offsetlen+npred))
	   print len(lpinttime)
	   
	   #plt.plot(inttime)
	   #plt.plot(lpinttime)
	   #plt.show()

	   lpint = rft(lpinttime)
	   
	   #lpint=rft(inttime)
	   #lpintmax = max(lpint)
	   
	   #lpintinv = []
	   #for j in range(0,len(lpint)):
	   # lpintinv.append(epintmax-lpint[j])
	   # 
	   #plt.plot(epint)
	   #plt.plot(lpintinv)
	   #
	   #plt.show()
	   
	   ### Alternative error calculation ###
	   
	   #print epint
	   #print lpintinv
	   origoffset = array(ep.offset)
	   #print origoffset
	   comparoffset = offsets[::2]
	   #print offsets
	   #print comparoffset
	   
	   lpinterr = lpint[::2]
	   #print lpintinverr
	   #print epint
	   
	   allerr=[]
	   sumerr = 0

	   for j in range(0,len(epint)):
	    temperr = (epint[j]-lpinterr[j])
	    allerr.append(temperr)
	    sumerr += temperr**2
	   
	   #
	   errlp = sqrt(sumerr/len(epint))
	   
	   print errlp

	   
	   err = sqrt(ep.intstd[i]**2+errlp**2)
	   
	   print err
	   #
	   #

	   ##
	   for ii in range(0,len(lpint)-1):
	    #buf += '%8.4f\t%8.4f\t%8.4f\n' % (offsets[ii], lpintinv[ii],ep.intstd[i])
	    buf += '%8.4f\t%8.4f\t%8.4f\n' % (offsets[ii], lpint[ii],err)
	   
	   filex = open('lp_'+str(int(ep.v1))+'Hz_'+str(r.label)+'.txt', 'w')
	   filex.write(buf)
	   filex.close()

