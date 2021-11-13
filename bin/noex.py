# -*- coding: utf-8 -*-
######
# 2015 by Donghan Lee
# edited 2021 by Kyungdoe Han
########
import os
from sys import argv, stderr
from json import dumps, encoder
from estp.est_model_noex import est_model
from time import ctime

encoder.FLOAT_REPR = lambda o: format(o, '.4f')

m2 = est_model()
m2.verbose = True

newDict = {'Header' : [],
           'Project Name' : 'default',
           'init': {
                    'kex':{'min' : 10.0, 'max' : 400.0, 'nsteps' : 6},
                     'pB':{'min' : 0.01, 'max' : 0.1, 'nsteps' : 6},
                     'Method' : 'NoEx' 
                    }
 }
newDict['Header'].append(m2.programName)
newDict['Header'].append('Time: ' + ctime())

expList = []

#
if len(argv) < 2:
    stderr.write('Please specify the input file(s):\nprepare.py inputfile1 inputfile2 ...\n')
    exit()

# define datasets

for i in range(1, len(argv)):
    m2.dataset.addData(argv[i])
    expList.append(argv[i])

newDict['datasets'] = expList

rdlist = m2.dataset.getResidues()

newDict['residues'] = rdlist

buf = dumps(newDict, sort_keys=True, indent=4)
print(buf)

stderr.write('Done.\n')
