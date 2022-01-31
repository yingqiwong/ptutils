# ptutils
Companion functions for pantarhei. 

perm/ contains scripts to test and plot permission functions and 
flux and transfer coefficients, as defined in Keller and Suckale
2019. 

makeplots/ contains a script to read pantarhei mat file outputs 
and make various plots to examine the solutions.

utils/ contains anly/ and plot/ functions to analyze and plot
pantarhei outputs respectively. These funtions are called in the 
scripts in perm/ and makeplots/. 



## IMPORTANT NOTE: 
This repo requires dependencies on the pantarhei and ternplot 
git repos. 
Before running codes, PLEASE CHECK Addpaths.m in perm/ and 
makeplots/ to make sure that the path to pantarhei and ternplot
are correctly specified, depending on the file organisation on 
your computer.
The code will break if you don't specify the correct 
paths to pantarhei and ternplot
