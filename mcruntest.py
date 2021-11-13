#!/usr/bin/env python

######
# 2015 by Donghan Lee
########

from sys import argv, stdout
from numpy import array
from estp.est_model_noex_test import est_model
import json


nrun = int(argv[2])

m2=est_model()
m2.verbose = True

print ('**************')
stdout.write(m2.programName + '\n')
print ('**************')

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
    print ('Residue: ' + resid + ' ' + active)
    for r in m2.dataset.res:
        if r.label == resid:
            if active == 'on':
                r.active = True
            elif active == 'off':
	            r.active = False
            else:
                print ('Error: wrong flag for the residue' + resid + ' ' + active)
                exit()

m2.datapdf(projectName + '_data.pdf')

p0=m2.initGuessAll(conf['init'])

out=m2.fit(p0)

pp = out[0]
print ('**************')
print (out[0])
print ('**************')
m2.pdf(out[0],projectName + '_mc.pdf')




allps = []
for i in range(0, nrun):
    mmc =est_model()
    mmc.verbose = True
    
    projectName = conf['Project Name']
    datasetsNames = conf['datasets']
    for dataset in datasetsNames:
        mmc.dataset.addDataWithError(dataset)


    residues = conf['residues']

    for r in residues:
        resid = r['name']
        active = r['flag']
        print ('Residue: ' + resid + ' ' + active)
        for r in mmc.dataset.res:
            if r.label == resid:
                if active == 'on':
	                r.active = True
                elif active == 'off':
	                r.active = False
                else:
                    print ('Error: wrong flag for the residue' + resid + ' ' + active)
                    exit()

    mmc.selMethod(conf['init'])
    outmc=mmc.fit(pp)
    
    allps.append(outmc[0])

logBuf = m2.getLogBufferMC(pp,allps)

file1 = open(projectName + '_mc.txt', 'w')
print (logBuf[0])
file1.write(logBuf[0])
file1.close()

m2.pdf(logBuf[1],projectName + '_mcmean.pdf')

print (logBuf[1])

print ('########\n')




