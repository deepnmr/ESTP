# ESTP

The ESTP program is for fitting CEST and DEST NMR experiment.

There are two possible ways to use the program.

Following is a quick tutorial using an example in the directory of examples.

* Using terminal

1. go to the directory.

$ cd examples

2. make a configure file

% python ../bin/prepare.py data/* > configure.txt

3. Modifiy configure file, which is a dictionary. By choosing on, fitting will be done on those resigues. If off is chosen, the residue will be skipped for the fitting.

4. actual fitting. It will produce "project_name"_result.txt, "project_name".pdf, and "project_name"_data.pdf.

% python ../bin/run.py configure.txt

5. for MC run, give an additional argument representing how many MC run will be performed. Follwoing is 100 MC run.

% python ../bin/mcrun.py configure.txt 100


* Using gui

% python ESTP.py

------------------------------------
**guiestp1.0.zip**
This contains all.
------------------------------------

Note that proc_base.py and leastbound.py are copied from nmrglue (https://www.nmrglue.com/) written by Dr.Jonathan J. Helmus.
Note that the program was inspired by Dr. Adam Mazur's code ShereKhan
Note that Dr. Marta G Carneiro wrote no exchange portion and output part of linear prediction of cest profile.

Please cite 
1. Marta G Carneiro, Jithender G Reddy, Christian Griesinger, Donghan Lee

Speeding-up exchange-mediated saturation transfer experiments by Fourier transform
J Biomol NMR 2015 Nov;63(3):237-44. doi: 10.1007/s10858-015-9985-9.

2. To be published.



