###
# 2015 by Donghan Lee
# Jan 26, 2015
###

from matplotlib import use
from numpy.lib.function_base import linspace
from numpy.linalg import eig, inv
from numpy import array, dot, append, pi, argmin, max, log, zeros, sort, sqrt, mean, std, transpose, exp, cos, sin, tan, arctan2, nan_to_num, diag, real
from scipy.linalg.matfuncs import expm2
from scipy.optimize import leastsq, minimize
from scipy.stats import norm
from matplotlib.pyplot import plot, errorbar, figure, title, legend, xlabel, ylabel, grid, xlim, ylim, close, subplot, savefig
from matplotlib.backends.backend_pdf import PdfPages
from sys import argv
from os import getlogin, uname
from time import ctime
import json

from est_data import EstDataSet
from leastsqbound import leastsqbound

class est_model:
   def __init__(self):
       self.programName = 'est ver. 1.0'
       self.header = 'Donghan Lee\n' \
                       + ctime() + '\n'
       
       self.method = 'Baldwin'                

       # global parameters
       self.kab = 0.0
       self.kab_std = 0.0
       self.kba = 0.0
       self.kba_std = 0.0
       #
       self.dataset = EstDataSet()
       self.verbose = False


       self.buf = ''
       self.buf += self.programName + '\n\n'
       self.buf += self.header



   # Models
   
   # two-state matrix
   def matrix_calc(self, kab, kba, dG, dw, R1, R2a, R2b, dRF, T, v1, v1err, B0):
    
       R1a = R1
       R1b = R1
       
#       R2b = R2a + dR
 
 
       wRF = dRF * B0 * 2.0  * pi
      

       wG = dG * B0 * 2 * pi
       wa = wG - wRF
    
    
       wE = (dG + dw) * B0 * 2 * pi 
       wb = wE - wRF

       wy = 0.0
        
       kex = kab + kba 
       pB = kab / kex
       pA = kba / kex

       w1 = v1*2*pi
       w1err = v1err*2*pi
       
       wxs = linspace(-2.0*w1err + w1, 2.0*w1err + w1, 10)
       
       w1s = []
       w1errs = []
       for i in range(0,10):
           w1s.append(w1)
           w1errs.append(w1err)

       weights = norm.pdf(wxs,w1s,w1errs)
       refM = pA

       startM = [0,0,0,pA,0,0,pB]
       
       preIcal = 0.0
       for i in range(0,10):
          wx = wxs[i]
          Alist = [[0,0,0,0,0,0,0],
                  [0,-kab-R2a,-wa,wy,kba,0,0],
                  [0,wa,-kab-R2a,-wx,0,kba,0],
                  [2*R1a*pA,-wy,wx,-kab-R1a,0,0,kba],
                  [0,kab,0,0,-kba-R2b,-wb,wy],
                  [0,0,kab,0,wb,-kba-R2b,-wx],
                  [2*R1b*pB,0,0,kab,-wy,wx,-kba-R1b]]
          A = array(Alist)
          eA = expm2(A * T)
          endM = dot(eA,startM)
          preIcal += weights[i]*endM[3]/refM
    
       Icalc = preIcal/sum(weights)
       
       if(Icalc < 0):
          Icalc = 0.0

       return Icalc

   def Baldwin(self, kab, kba, dG, dw, R1, R2a, R2b, dRF, T, v1, v1err, B0):
       R1a = R1
       R1b = R1




       kex = kab + kba 
       pE = kab / kex
       pG = kba / kex

       dR = R2b-R2a
       
       wRF = dRF * B0 * 2.0  * pi
      
       
       wG = dG * B0 * 2.0 * pi
       ddG = wG - wRF
       
       dE = dG + dw
       wE = dE * B0 * 2.0 * pi 
       ddE = wE - wRF
       
       ddw = ddE - ddG
       
       da = pG*ddG + pE*ddE

       w1b = v1*2*pi
       w1err = v1err*2*pi
       
       wxs = linspace(-2.0*w1err + w1b, 2.0*w1err + w1b, 10)
       
       w1s = []
       w1errs = []
       for i in range(0,10):
           w1s.append(w1b)
           w1errs.append(w1err)
       

       weights = norm.pdf(wxs,w1s,w1errs)
       
       
       
       preIcal = 0.0
       for i in range(0,10):
           w1 = wxs[i]
           OG2 = w1**2 + ddG**2
           OE2 = w1**2 + ddE**2
           Oa2 = w1**2 + da**2
           
           tantheta2 = (w1**2)/(da**2)
           sintheta2 = (w1**2)/Oa2
           costheta2 = (da**2)/Oa2
           
           F1p = pG*pE*(ddw**2)
           F2p = (kex**2) + (w1**2) + ((ddG*ddE/da)**2)
           Dp = (kex**2) + OG2*OE2/Oa2
           
           F1 = pE*(OG2 + (kex**2) + dR*pG*kex)
           F2 = 2*kex +  (w1**2)/kex + dR*pG
           F3 = 3*pE*kex + (2*pG*kex + (w1**2)/kex + dR + dR * ((pE*kex)**2)/OG2)*(OG2/(w1**2))
           
           CR1 = (F2p + (F1p + dR*(F3-F2))*tantheta2)/(Dp + dR*F3*sintheta2)
           CR2 = (Dp/sintheta2 - F2p/tantheta2 - F1p + dR*F2)/(Dp + dR*F3*sintheta2)
           Rex = (F1p*kex + dR*F1)/(Dp + dR*F3*sintheta2)
           
           R1r = CR1*R1*costheta2 + (CR2*R2a + Rex)*sintheta2
           
           preIcal += weights[i]*costheta2*exp(-1.0*T*R1r)
       
       Icalc = preIcal/sum(weights)
       
       return Icalc
       
       


   def initGuessAll(self,initConf):
       # grid search
       kex = linspace(initConf['kex']['min'], initConf['kex']['max'], initConf['kex']['nsteps'])
       
       pB = linspace(initConf['pB']['min'], initConf['pB']['max'], initConf['pB']['nsteps'])

#       dws = linspace(initConf['dw']['min'], initConf['dw']['max'], initConf['dw']['nsteps'])
       
       if(initConf['Method'] == 'Baldwin'):
           self.method = 'Baldwin'
       elif(initConf['Method'] == 'Matrix'):
           self.method = 'Matrix'
       else:
           self.method = 'Baldwin'    
       
       minChi2 = 9.0e99
       
       for i_kex in range(0, len(kex)):
           for i_pB in range(0,len(pB)):
#              for i_dw in range(0,len(dws)):
                p0 = []
                            
                kab = pB[i_pB] * kex[i_kex]
                kba = kex[i_kex] - kab
                p0.append(kab)
                p0.append(kba)
                            
                for i in range(0, len(self.dataset.res)):
                    if(self.dataset.res[i].active):
                        dG = self.dataset.res[i].estSpecs[0].offset[argmin(self.dataset.res[i].estSpecs[0].int)]
                        ctT = self.dataset.res[i].estSpecs[0].T 
                        maxint = max(self.dataset.res[i].estSpecs[0].int)
                        r1 = -1.0/ctT*log(maxint)
#                        dw = dws[i_dw] 
                        dw = self.dataset.res[i].estSpecs[0].initdw
                        r2a = self.dataset.res[i].estSpecs[0].initr2a
                        r2b = self.dataset.res[i].estSpecs[0].initr2b
                                    
                                    
                        p0.append(dG)
                        p0.append(dw)
                        p0.append(r1)
                        p0.append(r2a)
                        p0.append(r2b)
                        
#                print p0            
                self.errFunc(p0)
                           
                if(self.chi2 < minChi2):
                    minChi2 = self.chi2
                    p1 = p0
#                else:
#                    p1 = p0


       print 'Initial Guess is finished****************\n'
       
       
       return p1



   def seParam(self,p):
       # parameters
       j = 0
       kab = p[0]
       kba = p[1]
       j = 2

       nres = len(self.dataset.res) # number of residues
       dGs = [] # nRes
       dws = [] # nRes
       r1s = [] # nRes
       r2as = [] # nRes
       r2bs = [] # nRes
                                        
       # get local parameters
       for i in range(0,nres):
            if(self.dataset.res[i].active):
                dGs.append(p[j])
                dws.append(p[j+1])
                r1s.append(p[j+2])
                r2as.append(p[j+3])
                r2bs.append(p[j+4])
                j += 5
            else:
                dGs.append(0.0)
                dws.append(0.0)
                r1s.append(0.0)
                r2as.append(0.0)
                r2bs.append(0.0)

       return [kab,kba,dGs,dws,r1s,r2as,r2bs]

   def comParam(self,par):
       # parameters
       kab = par[0]
       kba = par[1]
       dGs = par[2]
       dws = par[3] # nRes
       r1s = par[4] # nRes
       r2as = par[5] # nRes
       r2bs = par[6] # nRes
       
       
       p = []
       p.append(kab)
       p.append(kba)

       nres = len(self.dataset.res) # number of residues
                                        
       # get local parameters
       for i in range(0,nres):
            if(self.dataset.res[i].active):
                p.append(dGs[i])
                p.append(dws[i])
                p.append(r1s[i])
                p.append(r2as[i])
                p.append(r2bs[i])

       return p





   # chi functions for leastsq
   def errFunc(self, p):
       par = self.seParam(p)
       kab = par[0]
       kba = par[1]
       dGs = par[2]
       dws = par[3]
       r1s = par[4]
       r2as = par[5]
       r2bs = par[6]
       
       
       
       # actual calculations
       residuals = array([])
       nres = len(self.dataset.res)
       for i in range(0,nres):
           res = self.dataset.res[i]
           if(res.active):
               for j in range(0, len(res.estSpecs)):
                  estspec = res.estSpecs[j]
                  for k in range(0,len(estspec.int)):
                      if(self.method == 'Matrix'):
                          est_calc = self.matrix_calc(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                      elif(self.method == 'Baldwin'):
                          est_calc = self.Baldwin(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                      else:
                          est_calc = self.Baldwin(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                      est_calc = nan_to_num(est_calc)
                      est_rsd = (estspec.int[k]-est_calc) / estspec.intstd[k]
                      residuals = append(residuals,est_rsd)

#       print residuals
       self.chi2 = sum(residuals**2)
       self.npar = len(p)
       self.dof = len(residuals) - len(p) - 1.0
       self.nvar = len(residuals)


       if(self.verbose):
          print 'kab=%8.3f kba=%8.3f chi2=%12.3f dof=%d nv=%d np=%d' % (kab, kba, self.chi2, self.dof, self.nvar, self.npar)
       print p

       return residuals
#       return sum(residuals**2)
       

   def fit(self,p0):
       bnds = []
       bnds.append((0,None))
       bnds.append((0,None))
       lenp = len(p0)
       numres = (lenp - 2)/5
       for i in range(0,numres):
           bnds.append((None,None))
           bnds.append((None,None))
           bnds.append((0,None))
           bnds.append((0,None))
           bnds.append((0,None))
       
 #      out = leastsq(self.errFunc, x0=p0, full_output=1,ftol=1e-9, xtol=1e-12,factor=100)
       out = leastsqbound(self.errFunc, x0=p0, bounds = bnds, full_output=1)
       

#       out = minimize(self.errFunc, x0=p0,bounds=bnds,method='SLSQP')
#       if(out[4] != 1 and out[4] != 2 and out[4] != 3 and out[4] != 4 ):
#          exit()       
       
       print p0
       p1 = out[0]
       
       covar= out[1]
       
#       par = self.seParam(p1)
#       kab = par[0]
#       kba = par[1]
#       dGs = par[2]
#       dws = par[3]
#       r1s = par[4]
#       r2as = par[5]
#       r2bs = par[6]
       
#       for i in range(0,len(r2as)):
#          if(r2as[i] < 0):
#             if(r2bs[i] < 0):
#                r2as[i] = abs(r2as[i])
#                r2bs[i] = abs(r2bs[i])
#             else:
#                r2as[i] = abs(r2as[i])
#                r2bs[i] = r2bs[i] + 2*abs(r2as[i])
#          elif(r2bs[i] < 0):
#             r2bs[i] = abs(r2bs[i]) + 2.0*(abs(r2as[i]))
#
#       par = [kab,kba,dGs,dws,r1s,r2as,r2bs]
#       
#       p1 = self.comParam(par)
       

       return p1, covar



   def pdf(self, p1, pdfFileName):
        pdf = PdfPages(pdfFileName)
#        p1 = out[0]

        par = self.seParam(p1)
        kab = par[0]
        kba = par[1]
        dGs = par[2]
        dws = par[3]
        r1s = par[4]
        r2as = par[5]
        r2bs = par[6]
       


        nres = len(self.dataset.res)


        for i in range(0,nres):
            res = self.dataset.res[i]
            if(res.active):
                figure(figsize=(8, 6))
                
                colorsSet = ['b', 'g', 'r', 'c', 'm', 'y']
                colorsSet.reverse()
                
                min_x = []
                max_x = []
                min_y = [] 
                max_y = []
                

                for j in range(0, len(res.estSpecs)):
                    ep = res.estSpecs[j]
                    
                    min_x.append(min(ep.offset))
                    max_x.append(max(ep.offset))
                    
                    
                    # get color
                    c = 'k'
                    try:
                        c = colorsSet.pop()
                    except:
                        pass

                    errorbar(ep.offset, ep.int, yerr=ep.intstd, fmt='%so' % c, markersize=3, label=' exp\n%5.1f Hz\n%5.1f ms' % (ep.v1, ep.T*1000.0))
                    tx = array(ep.offset)
                    tx = sort(tx)
                    ty = zeros(tx.shape)

                    for ii in range(0, len(tx)):
                     if(self.method == 'Matrix'):
                        ty[ii] = self.matrix_calc(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],tx[ii],ep.T,ep.v1,ep.v1err,ep.field)
                     elif(self.method == 'Baldwin'):
                        ty[ii] = self.Baldwin(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],tx[ii],ep.T,ep.v1,ep.v1err,ep.field)
                     else:
                        ty[ii] = self.Baldwin(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],tx[ii],ep.T,ep.v1,ep.v1err,ep.field)
                    plot(tx, ty, '%s-' % c, label='calc\n%5.1f Hz\n%5.1f ms' % (ep.v1, ep.T*1000.0))
                    
                    min_yv = []
                    max_yv = []
                    for k in range(0, len(ep.int)):
                        min_yv.append(ep.int[k] - ep.intstd[k])
                        max_yv.append(ep.int[k] + ep.intstd[k])

                    min_y.append(min(min_yv))
                    max_y.append(max(max_yv))
                    
                title(res.label)
                xlabel('chemical shift [ppm]')
                ylabel('Intensity')
                #
                # Shrink current axis by 20%
                box = subplot(111).get_position()
                subplot(111).set_position([box.x0, box.y0, box.width * 0.8, box.height])
                # Put a legend to the right of the current axis
                legend(loc='center left', bbox_to_anchor=(1, 0.5))
                #
                grid(True)
                #ylim(0, 1.1 * max(max(ty), max(dspCurve.R2exp)))
                ylim(0.95 * min(min_y), 1.05 * max(max_y))
                xlim(min(min_x), max(max_x))
                pdf.savefig()
                close()
                outputFileName = pdfFileName.replace(".pdf", "") + str(res.label)
                savefig(outputFileName+str(i)+".png", dpi=100) # KH 2021

        pdf.close()


   def datapdf(self, pdfFileName):
        pdf = PdfPages(pdfFileName)
	
        for r in self.dataset.res:
            if(r.active):
                figure(figsize=(8, 6))
                
                colorsSet = ['b', 'g', 'r', 'c', 'm', 'y']
                colorsSet.reverse()
                
        		min_x = []
                max_x = []
                min_y = [] 
                max_y = []
                
                for i in range(0, len(r.estSpecs)):
                    ep = r.estSpecs[i]
                    
	        	    min_x.append(min(ep.offset))
                    max_x.append(max(ep.offset))
                    
                    # get color
                    c = 'k'
                    try:
                        c = colorsSet.pop()
                    except:
                        pass

                    errorbar(ep.offset, ep.int, yerr=ep.intstd, fmt='%so' % c, markersize=3, label=' exp\n%5.1f Hz\n%5.1f ms' % (ep.v1,ep.T*1000.0))
                    
                    min_yv = []
                    max_yv = []
                    for i in range(0, len(ep.int)):
                        min_yv.append(ep.int[i] - ep.intstd[i])
                        max_yv.append(ep.int[i] + ep.intstd[i])

                    min_y.append(min(min_yv))
                    max_y.append(max(max_yv))
                    
                title(r.label)
                xlabel('chemical shift [ppm]')
                ylabel('Intensity')
                #
                # Shrink current axis by 20%
                box = subplot(111).get_position()
                subplot(111).set_position([box.x0, box.y0, box.width * 0.8, box.height])
                # Put a legend to the right of the current axis
                legend(loc='center left', bbox_to_anchor=(1, 0.5))
                #
                grid(True)
                #ylim(0, 1.1 * max(max(ty), max(dspCurve.R2exp)))
                ylim(0.95 * min(min_y), 1.05 * max(max_y))
                xlim(min(min_x), max(max_x))
                pdf.savefig()
                outputFileName = pdfFileName.replace(".pdf", "") + str(r.label)
                savefig(outputFileName+"_data.png", dpi=100) # KH 2021
                close()
        
        pdf.close()

   def getLogBuffer(self,out):

       buf  = '***************************************************\n'
       buf += self.programName
       buf += '\n'
       buf += '***************************************************\n\n'
       user = getlogin()
       hostname = uname()[1]
       buf += 'User: %s@%s\n' % (user, hostname)
       buf += '%s\n' % ctime()
       buf += '***************************************************\n'
#       buf += 'Input file:\n'
#       file1 = open(argv[1],'r')
#       for line in file1.readlines():
#           buf += line
#       file1.close()
#       buf += '***************\n\n'

#       buf += 'Datatset(s):\n'
       p = out[0]
       covar = out[1]
       parstd = []
       for i in range(0,len(p)):
           std = sqrt(covar[i][i])
           parstd.append(std)

       par = self.seParam(p)
       kab = par[0]
       kabstd = parstd[0]
       kba = par[1]
       kbastd = parstd[1]
       dGs = par[2]
       dws = par[3]
       r1s = par[4]
       r2as = par[5]
       r2bs = par[6]


       nres = len(self.dataset.res)
                                
       for i in range(0,nres):
            res = self.dataset.res[i]
            if(res.active):
                buf += '# %s\n' % (res.label)
                for j in range(0, len(res.estSpecs)):
                    estspec = res.estSpecs[j]
                    buf += 'B0 Field [MHz]: %8.3f\n' % (estspec.field)
                    buf += 'T [ms]:         %8.3f\n' % (estspec.T*1000.0)
                    buf += 'v1 [Hz]:        %8.3f\n' % (estspec.v1)
                    for k in range(0,len(estspec.int)):
                     if(self.method == 'Matrix'):
                        est_calc = self.matrix_calc(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                     elif(self.method == 'Baldwin'):
                        est_calc = self.Baldwin(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                     else:
                        est_calc = self.Baldwin(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                     buf += '%8.3f %12.3f %12.3f %12.3f\n' % (estspec.offset[k], estspec.int[k], estspec.intstd[k], est_calc)
                    buf += '***************************************************\n'

       xxx = self.errFunc(p)

       buf += '\n\n\n*****************************************\n'
       buf += 'Results using %s method\n\n' % (self.method)
       buf += 'Chi2:     %8.6f\n' % (self.chi2)
       buf += 'red_chi2: %8.6f\n' % (self.chi2/self.dof)         
       buf += 'dof:      %d\n' % (self.dof)
       buf += 'nv:       %d\n' % (self.nvar)
       buf += 'np:       %d\n' % (self.npar)
       buf += '*******************************************\n'

       buf += '*******************************************\n'

       buf += 'fitted Parameters\n'
       buf += 'kab [s-1]: %8.3f +/- %8.3f\n' % (kab, kabstd)
       buf += 'kba [s-1]: %8.3f +/- %8.3f\n' % (kba, kbastd)


       j = 2

       for r in self.dataset.res:
           if(r.active):
               resid = r.label
               dG = p[j]
               dGstd = sqrt(covar[j][j])
               dw = p[j+1]
               dwstd = sqrt(covar[j+1][j+1])
               r1 = p[j+2]
               r1std = sqrt(covar[j+2][j+2])
               r2a = p[j+3]
               r2astd = sqrt(covar[j+3][j+3])
               dR = p[j+4]
               dRstd = sqrt(covar[j+4][j+4])
#               if(p[j+4] < 0):
#                   dR = abs(dR) + 2*r2a
#                   dRstd = dRstd + 2*r2astd

             
               buf += '*******************************************\n'
               buf += 'residue: %s\n' % (resid)
               buf += 'peak position [ppm]: %8.4f +/- %8.4f\n' % (dG, dGstd)
               buf += 'cs difference [ppm]: %8.4f +/- %8.4f\n' % (dw, dwstd)
               buf += 'R1 [s-1]:            %8.4f +/- %8.4f\n' % (r1, r1std)
               buf += 'R2a [s-1]:           %8.4f +/- %8.4f\n' % (r2a, r2astd)
               buf += 'R2b [s-1]:           %8.4f +/- %8.4f\n' % (dR, dRstd)
               buf += '*******************************************\n'

               j += 5

       return buf

   def getLogBufferMC(self,p, allps):

       buf  = '***************************************************\n'
       buf += self.programName
       buf += '\n'
       buf += '***************************************************\n\n'
       user = getlogin()
       hostname = uname()[1]
       buf += 'User: %s@%s\n' % (user, hostname)
       buf += '%s\n' % ctime()
       buf += '***************************************************\n'
#       buf += 'Input file:\n'
#       file1 = open(argv[1],'r')
#       for line in file1.readlines():
#           buf += line
#       file1.close()
#       buf += '***************\n\n'

#       buf += 'Datatset(s):\n'

       
       meanps = []
       parstd = []
       mcmat = transpose(allps)
       for i in range(0,len(mcmat)):
           parstd.append(std(mcmat[i]))
           meanps.append(mean(mcmat[i]))

       par = self.seParam(p)
       kab = par[0]
       kba = par[1]
       dGs = par[2]
       dws = par[3]
       r1s = par[4]
       r2as = par[5]
       r2bs = par[6]


       nres = len(self.dataset.res)
                                
       for i in range(0,nres):
            res = self.dataset.res[i]
            if(res.active):
                buf += '# %s\n' % (res.label)
                for j in range(0, len(res.estSpecs)):
                    estspec = res.estSpecs[j]
                    buf += 'B0 Field [MHz]: %8.3f\n' % (estspec.field)
                    buf += 'T [ms]:         %8.3f\n' % (estspec.T*1000.0)
                    buf += 'v1 [Hz]:        %8.3f\n' % (estspec.v1)
                    for k in range(0,len(estspec.int)):
                        #est_calc = self.matrix_calc(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                        est_calc = self.Baldwin(kab,kba,dGs[i],dws[i],r1s[i],r2as[i],r2bs[i],estspec.offset[k],estspec.T,estspec.v1,estspec.v1err,estspec.field)
                        buf += '%8.3f %12.3f %12.3f %12.3f\n' % (estspec.offset[k], estspec.int[k], estspec.intstd[k], est_calc)
                    buf += '***************************************************\n'


       buf += '\n\n\n*****************************************\n'
       buf += 'Results\n\n'
       buf += 'Chi2:     %8.6f\n' % (self.chi2)
       buf += 'red_chi2: %8.6f\n' % (self.chi2/self.dof)         
       buf += 'dof:      %d\n' % (self.dof)
       buf += 'nv:       %d\n' % (self.nvar)
       buf += 'np:       %d\n' % (self.npar)
       buf += '*******************************************\n'

       buf += '*******************************************\n'

       buf += 'fitted Parameters\n'
       buf += 'kab [s-1]: %8.3f +/- %8.3f\n' % (kab, parstd[0])
       buf += 'kba [s-1]: %8.3f +/- %8.3f\n' % (kba, parstd[1])


       j = 2

       for r in self.dataset.res:
           if(r.active):
               resid = r.label
               dG = p[j]
               dGstd = parstd[j]
               dw = p[j+1]
               dwstd = parstd[j+1]
               r1 = p[j+2]
               r1std = parstd[j+2]
               r2a = p[j+3]
               r2astd = parstd[j+3]
               r2b = p[j+4]
               r2bstd = parstd[j+4]
             
               buf += '*******************************************\n'
               buf += 'residue: %s\n' % (resid)
               buf += 'peak position [ppm]: %8.4f +/- %8.4f\n' % (dG, dGstd)
               buf += 'cs difference [ppm]: %8.4f +/- %8.4f\n' % (dw, dwstd)
               buf += 'R1 [s-1]:            %8.4f +/- %8.4f\n' % (r1, r1std)
               buf += 'R2a [s-1]:           %8.4f +/- %8.4f\n' % (r2a, r2astd)
               buf += 'R2b [s-1]:           %8.4f +/- %8.4f\n' % (r2b, r2bstd)
               buf += '*******************************************\n'

               j += 5

       print len(allps)
       return [buf, meanps]


