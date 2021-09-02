# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
from scipy.optimize import curve_fit
from struct import *

#Variables___________________________________________________________________
stepTot = 2000000
stepInc = 2000
int_data = 4
double_data = 22
int_increment = 4*(int_data - 1)
double_increment = 8*(double_data - 1)

structNum = 5
runNum = 8

#Functions___________________________________________________________________
def BintoArray(datafile, init, datatype, stetot, steInc, dataNum, inc):
    datafile.seek(init, 0)  
    arr = np.fromfile(datafile, dtype=datatype, count=1, sep="")  
    for i in np.arange (1,stetot + 1):     
        if (i%stepInc == 0):
            datafile.seek(np.dtype(datatype).itemsize*dataNum*(steInc-1) + inc, 1) 
            #datafile.seek(inc, 1) 
            arr = np.append(arr, np.fromfile(datafile, datatype, count=1, sep=""))  
    return arr
                         
#Open binary data files with data____________________________________________
data1 = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_double.bin", "rb") 
data2 = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_centre.bin", "rb")

#Fill data arrays________________________________________

#Speed of micelle____________________________________________
speed = np.fromfile(data2, dtype='d', count=-1, sep="")
speed = speed[1:]  
print speed
print speed.size

#Time________________________________________________________
TIME = BintoArray(data1, 0, 'd', 2000000, 2000, 22, 8*(22-1))
TIME = TIME[651:]
print TIME
print TIME.size

#Radius of micelle__________________________________________
rad = BintoArray(data1,112,'d', 2000000, 2000, 22, 8*(22-1))
rad = rad[651:]
print rad
print rad.size

#Average number of solvent particles in mpc box______________
data1.seek(22*8*(2000000 + 1), 0)  
M = np.fromfile(data1, dtype='d', count=1, sep="") 
for i in np.arange (2, 50000 + 1): 
    if (i%50==0):
        if (i==3):
            data1.seek((50-2)*16 + 8, 1)     
        else:
            data1.seek((50-1)*16 + 8, 1) 
        M = np.append(M, np.fromfile(data1, 'd', count=1, sep=""))    
M = M[650:]
print M
print M.size

visc = 0.1*M*(M + np.ones(M.size) - np.exp(-M))*(1.0/(M - np.ones(M.size) + np.exp(-M))) + (5.0/24.0)*(M - np.ones(M.size) + np.exp(-M))

#Scmidt number________________________________________________
sc = M - np.ones(M.size) + np.exp(-M)
sc = sc*sc
sc = (1.0/M)*sc
sc = (1.0/(M + np.ones(M.size) - np.exp(-M)))*sc
sc = sc*25.0/12.0
sc = sc + np.ones(M.size)

#Reynolds number________________________________________________
re = speed*rad
re = (1.0/(0.4*(-0.5*np.ones(M.size) + M*(1.0/(M - np.ones(M.size) + np.exp(-M)))) + (5.0/12.0)*(1.0/M)*(M - np.ones(M.size) + np.exp(-M))))*re

#Peclet number__________________________________________________
pe = 0.2*M*(-0.5*np.ones(M.size) + M/(M - np.ones(M.size) + np.exp(-M))) + (5.0/24)*(M - np.ones(M.size) + np.exp(-M))
pe = 6.0*np.pi*speed*rad*rad*pe



'''
f = py.figure(figsize=(6.503937,8), dpi=1000)

subplots = np.array([f.add_subplot(3,1,1), f.add_subplot(3,1,2), f.add_subplot(3,1,3)])

subplots[0].plot(TIME, sc, label = "sc")   
subplots[1].plot(TIME, re, label = "re") 
subplots[2].plot(TIME, pe, label = "pe")

subplots[0].legend(loc = 'upper right', fontsize = 7)

py.savefig("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/centres.png", dpi = 1000)

py.show()

'''
AVGSIZE=np.average(M)
SPEED=np.average(speed)
visc = np.average(visc)
SC = np.average(sc)
RE = np.average(re)
PE = np.average(pe)

print AVGSIZE
print SPEED
print visc
print SC
print RE
print PE 

'''


plotSize = TIME.size
print TIME
print plotSize

speed = np.zeros(999-651)

sc

re

pe
  
    
#Figure_________________________________________________________________________
f = py.figure(figsize=(6.503937,8), dpi=1000)

subplots = np.array([f.add_subplot(1,3,1), f.add_subplot(1,3,2), f.add_subplot(1,3,3)])

for i in xrange(0,3):
    subplots[i].set_xlabel('Time')
    subplots[i].set_xlim((0.0,10000))
    subplots[i].set_xticks(np.array([0, 2500, 5000, 7500, 10000]))
    subplots[i].set_xticklabels([0,'',5000,'',10000])
subplots[3].set_xlim((0.0,3000))
subplots[3].set_xticks(np.array([0, 1000, 2000, 3000]))
subplots[3].set_xticklabels([0,'','',3000])

subplots[2].set_ylim(0.0,110)
subplots[2].set_yticks(np.arange(0,101,20))

subplots[3].set_ylim(0.0,110)
subplots[3].set_yticks(np.arange(0,101,20))

subplots[4].set_ylim(0.0,110)
subplots[4].set_yticks(np.arange(0,101,20))

subplots[0].set_ylabel('Sc')
subplots[1].set_ylabel('Re')
subplots[2].set_ylabel('Pe')


#Plot________________________________________________________________________

subplots[0].plot(TIME, Sc, label = "2H-1P-3T")   
subplots[1].plot(TIME, Re, label = "1H-1P-4T") 
subplots[2].plot(TIME, Pe, label = "(3H)1H-1P-3T")

subplots[0].legend(loc = 'upper right', fontsize = 7)
        
#f.subplots_adjust(wspace=0.2, hspace=0.6, left=0.05, right=0.96, top=0.96, bottom=0.08)
f.tight_layout()

py.savefig("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/centres.png", dpi = 1000)

py.show()
'''