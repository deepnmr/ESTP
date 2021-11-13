#!/usr/bin/env python

######
# 2015 by Donghan Lee
########

from sys import argv, stdout
from est_dl.est_model_noex import est_model

import json

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

m2.datapdf(projectName + '_data.pdf')

p0=m2.initGuessAll(conf['init'])

out=m2.fit(p0)

logBuf = m2.getLogBuffer(out)

file1 = open(projectName + '_result.txt', 'w')
print logBuf
file1.write(logBuf)
file1.close()

m2.pdf(out[0],projectName + '.pdf')

print '########\n'


