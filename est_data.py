###
# 2015 by Donghan Lee
###

from re import compile
from numpy import argmin, max, log, std, mean
from numpy.random import randn


# EstSpec class is for an est spectrum which depends on B0 (field), T, and v1.
# Thus, class has field, T, v1, offset, intensity (int), and stadard deviation of intensity (intstd). 
class EstSpec:
   def __init__(self):
       self.field = 0.0   # spectrometer MHz
       self.T = 0.0       # spin-lock field duration in s
       self.centerppm = 0.0       # center of spectrum
       self.v1 = 0.0      # spin-lock field strengths in Hz
       self.v1err = 0.0   # error of spin-lock field strengths in Hz
       self.offset = []   # offset in ppm
       self.int = []      # normalized intensity
       self.intstd = []   # error in the intensity
       self.initdw = 0.0
       self.initr2a = 1.0
       self.initr2b = 200.0

   def info(self):
       print ('--- B0 field: %8.3f [MHz]' % self.field)
       print ('-          T: %8.3f [ms]' %(1000.0*self.T))
       print ('-         v1: %8.3f %8.3f [Hz]' % (self.v1, self.v1err))
       for i in range(0,len(self.offset)):
           print( '%8.3f  %8.3f   %8.3f' % (self.offset[i], self.int[i], self.intstd[i]))
       print ('---')



# Residue class is for residue.
# It has label and active.
class Residue:
   def __init__(self):
      self.label = ''
      self.estSpecs = []
      self.active = True
      
      
      
# EstDataSet class is for all dataset.
# It has array of residue, field, T, v1, and EstSpecs.
class EstDataSet:
   def __init__(self):
      self.res = []    # residue
      self.fields = [] # fields
      self.Ts = []     # Ts
      self.centerppms = []     # Ts
      self.v1s = []    # v1s
      self.v1errs = [] # v1errs
      self.initR2 = False


   def addData(self, fileName):
      p1 = compile('^\s*(\S+)\s+')
      p2 = compile('^\s*(\S+)\s+')
      p3 = compile('^\s*(\S+)\s+(\S+)\s*$')

      inFile = open(fileName, 'r')


      # Read B0 field
      line = inFile.readline()
      try:
          result = p1.match(line)
          currentField = float(result.group(1))
          self.fields.append(currentField)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # Read T
      line = inFile.readline()
      try:
          result = p2.match(line)
          currentT = float(result.group(1))
          self.Ts.append(currentT)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()


      # Read v1
      line = inFile.readline()
      try:
          result = p3.match(line)
          currentV1 = float(result.group(1))
          currentV1err = float(result.group(2))
          self.v1s.append(currentV1)
          self.v1errs.append(currentV1err)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # a line is ignored
      line = inFile.readline()

      # read the data
      p4 = compile('^#\s+(\w+\d+)\s*$')
      p4a = compile('^#\s+(\w+\d+)\s+R2a:\s+(\S+)\s+R2b:\s+(\S+)\s+dw:\s+(\W*\S+)\s*$')

      line = inFile.readline()
      # initial R2 check
      if(p4a.match(line)):
          self.initR2 = True
          p4 = p4a

      p5 = compile('^\s*(\S+)\s+(\S+)\s+(\S+)\s*$')
      # residues

      while(line):
          try:
               result = p4.match(line)
               resLabel = result.group(1)
               if(self.initR2):
                   initR2a = float(result.group(2))
                   initR2b = float(result.group(3))
                   initDW = float(result.group(4))
          except:
               print ("\n Wrong input line:\n")
               print (line)
               exit()

          # create new residue
          resid = None
          for r in self.res:
           if(r.label == resLabel):
               resid = r
               break
          if(resid == None):
           resid = Residue()
           resid.label = resLabel
           self.res.append(resid)

          # read est profile
          line = inFile.readline()
          if(line):
           try:
               result = p5.match(line)
           except:
               print ("\nWrong input data line:\n")
               print (line)
               exit()

          # an est spectrum
          ep = EstSpec()
          ep.field = currentField
          ep.T = currentT
          ep.v1 = currentV1
          ep.v1err = currentV1err
          if(self.initR2):
            ep.initr2a = initR2a
            ep.initr2b = initR2b
            ep.initdw = initDW

          while(result):
           ep.offset.append(float(result.group(1)))
           ep.int.append(float(result.group(2)))
           ep.intstd.append(float(result.group(3)))
           line = inFile.readline()
           if(line):
               try:
                   result = p5.match(line)
               except:
                   print ("\nWrong input data line:\n")
                   print (line)
                   exit()
           else:
               result = False
          	  
	  
          # store to residue.estSpecs
          resid.estSpecs.append(ep)
	  

      #remove empty est profile
      for r in self.res:
           for i in range(len(r.estSpecs)):
               if(len(r.estSpecs[i].offset) == 0):
                   del r.estSpecs[i]

      inFile.close()


   def addDataV1(self, fileName):
      p1 = compile('^\s*(\S+)\s+')
      p2 = compile('^\s*(\S+)\s+')
      p3 = compile('^\s*(\S+)\s+(\S+)\s*$')

      inFile = open(fileName, 'r')

      # Read B0 field
      line = inFile.readline()
      try:
          result = p1.match(line)
          currentField = float(result.group(1))
          self.fields.append(currentField)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # Read T
      line = inFile.readline()
      try:
          result = p2.match(line)
          currentT = float(result.group(1))
          self.Ts.append(currentT)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # Read v1 and add noise for MC runs
      line = inFile.readline()
      try:
          result = p3.match(line)
          currentV1 = float(result.group(1))
          currentV1err = float(result.group(2))
          currentV1 = currentV1 + currentV1err*randn()
          self.v1s.append(currentV1)
          self.v1errs.append(currentV1err)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # a line is ignored
      line = inFile.readline()

      # read the data
      p4 = compile('^#\s+(\w+\d+)\s*$')
      p4a = compile('^#\s+(\w+\d+)\s+R2a:\s+(\S+)\s+R2b:\s+(\S+)\s+dw:\s+(\W*\S+)\s*$')

      line = inFile.readline()
      # initial R2 check
      if(p4a.match(line)):
          self.initR2 = True
          p4 = p4a

      p5 = compile('^\s*(\S+)\s+(\S+)\s+(\S+)\s*$')
      # residues

      while(line):
          try:
               result = p4.match(line)
               resLabel = result.group(1)
               if(self.initR2):
                   initR2a = float(result.group(2))
                   initR2b = float(result.group(3))
                   initDW = float(result.group(4))
          except:
               print ("\n Wrong input line:\n")
               print (line)
               exit()

          # create new residue
          resid = None
          for r in self.res:
           if(r.label == resLabel):
               resid = r
               break
          if(resid == None):
           resid = Residue()
           resid.label = resLabel
           self.res.append(resid)

          # read est profile
          line = inFile.readline()
          if(line):
           try:
               result = p5.match(line)
           except:
               print ("\nWrong input data line:\n")
               print (line)
               exit()

          # an est spectrum
          ep = EstSpec()
          ep.field = currentField
          ep.T = currentT
          ep.v1 = currentV1
          ep.v1err = currentV1err
          if(self.initR2):
            ep.initr2a = initR2a
            ep.initr2b = initR2b
            ep.initdw = initDW

          while(result):
           ep.offset.append(float(result.group(1)))
           ep.int.append(float(result.group(2)))
           ep.intstd.append(float(result.group(3)))
           line = inFile.readline()
           if(line):
               try:
                   result = p5.match(line)
               except:
                   print ("\nWrong input data line:\n")
                   print (line)
                   exit()
           else:
               result = False
          	  
	  
          # store to residue.estSpecs
          resid.estSpecs.append(ep)
	  

      #remove empty est profile
      for r in self.res:
           for i in range(len(r.estSpecs)):
               if(len(r.estSpecs[i].offset) == 0):
                   del r.estSpecs[i]

      inFile.close()

   def addDataWithError(self, fileName):
      p1 = compile('^\s*(\S+)\s+')
      p2 = compile('^\s*(\S+)\s+')
      p3 = compile('^\s*(\S+)\s+(\S+)\s*$')

      inFile = open(fileName, 'r')

      # Read B0 field
      line = inFile.readline()
      try:
          result = p1.match(line)
          currentField = float(result.group(1))
          self.fields.append(currentField)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # Read T
      line = inFile.readline()
      try:
          result = p2.match(line)
          currentT = float(result.group(1))
          self.Ts.append(currentT)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # Read v1 and add noise for MC runs
      line = inFile.readline()
      try:
          result = p3.match(line)
          currentV1 = float(result.group(1))
          currentV1err = float(result.group(2))
          self.v1s.append(currentV1)
          self.v1errs.append(currentV1err)
      except:
          print ("\nWrong input line:\n")
          print (line)
          exit()

      # a line is ignored
      line = inFile.readline()

      # read the data
      p4 = compile('^#\s+(\w+\d+)\s*$')
      p4a = compile('^#\s+(\w+\d+)\s+R2a:\s+(\S+)\s+R2b:\s+(\S+)\s+dw:\s+(\W*\S+)\s*$')

      line = inFile.readline()
      # initial R2 check
      if(p4a.match(line)):
          self.initR2 = True
          p4 = p4a

      p5 = compile('^\s*(\S+)\s+(\S+)\s+(\S+)\s*$')
      # residues

      while(line):
          try:
               result = p4.match(line)
               resLabel = result.group(1)
               if(self.initR2):
                   initR2a = float(result.group(2))
                   initR2b = float(result.group(3))
                   initDW = float(result.group(4))
          except:
               print ("\n Wrong input line:\n")
               print (line)
               exit()

          # create new residue
          resid = None
          for r in self.res:
           if(r.label == resLabel):
               resid = r
               break
          if(resid == None):
           resid = Residue()
           resid.label = resLabel
           self.res.append(resid)

          # read est profile
          line = inFile.readline()
          if(line):
           try:
               result = p5.match(line)
           except:
               print ("\nWrong input data line:\n")
               print (line)
               exit()

          # an est spectrum
          ep = EstSpec()
          ep.field = currentField
          ep.T = currentT
          ep.v1 = currentV1
          ep.v1err = currentV1err
          if(self.initR2):
            ep.initr2a = initR2a
            ep.initr2b = initR2b
            ep.initdw = initDW

          while(result):
           ep.offset.append(float(result.group(1)))
           intind = float(result.group(2))
           intstdind = float(result.group(3))
           intindx = intind + intstdind*randn()
           ep.int.append(intindx)
           ep.intstd.append(float(result.group(3)))
           line = inFile.readline()
           if(line):
               try:
                   result = p5.match(line)
               except:
                   print ("\nWrong input data line:\n")
                   print (line)
                   exit()
           else:
               result = False
          	  
	  
          # store to residue.estSpecs
          resid.estSpecs.append(ep)
	  

      #remove empty est profile
      for r in self.res:
           for i in range(len(r.estSpecs)):
               if(len(r.estSpecs[i].offset) == 0):
                   del r.estSpecs[i]

      inFile.close()



   def getResidues(self):
        """
        Generates residues dictionary
        """

        rdlist = []
        for r in self.res:
            r.active = True
            rd = {'name':r.label}
            if(r.active):
                rd['flag'] = 'on'
            else:
                rd['flag'] = 'off'
            rdlist.append(rd)
        return rdlist
        




   def info(self):
       for i in range(0, len(self.fields)):                     
           print (' % 8.3f' % self.fields[i])
           print (' % 8.3f' % self.Ts[i])
           print (' % 8.3f' % self.v1s[i])
           print ('#csd(ppm)        Intensity      Error')	
           for r in self.res:
               for ep in r.estSpecs:
                    if (ep.field == self.fields[i] and ep.T == self.Ts[i] and ep.v1 == self.v1s[i]):
                        print ('# ' + r.label)
                        for j in range(0, len(ep.offset)):
                            print ('%8.3f %8.3f %8.3f' % (ep.offset[j], ep.int[j], ep.intstd[j]))
        
