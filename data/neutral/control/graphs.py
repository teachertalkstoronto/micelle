# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
from struct import *

#Variables___________________________________________________________________
stepTot = 2000000
stepInc = 20
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
    for i in np.arange (0,stetot):
        if (i%stepInc == 0 and i != 0):
            datafile.seek(np.dtype(datatype).itemsize*dataNum*i + inc, 1) 
            #datafile.seek(inc, 1) 
            arr = np.append(arr, np.fromfile(datafile, datatype, count=1, sep=""))  
    return arr
                         
#Open binary data files with data____________________________________________
data = np.empty([structNum,runNum,2])

data[0,0,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_int.bin", "rb")
data[0,1,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/2_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_int.bin", "rb")
data[0,2,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/3_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_int.bin", "rb")
data[0,3,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/4_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_int.bin", "rb")
data[0,4,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_int.bin", "rb")
data[0,5,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/2_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_int.bin", "rb")
data[0,6,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/3_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_int.bin", "rb")
data[0,7,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/4_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_int.bin", "rb")
data[0,0,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_double.bin", "rb")
data[0,1,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/2_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_double.bin", "rb")
data[0,2,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/3_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_double.bin", "rb")
data[0,3,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/4_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000000steps_double.bin", "rb")
data[0,4,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_double.bin", "rb")
data[0,5,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/2_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_double.bin", "rb")
data[0,6,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/3_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_double.bin", "rb")
data[0,7,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/4_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent2000001steps_double.bin", "rb")

data[1,0,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/5_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[1,1,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/6_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[1,2,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/7_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[1,3,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/8_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[1,4,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/5_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[1,5,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/6_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[1,6,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/7_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[1,7,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/8_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[1,0,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/5_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[1,1,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/6_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[1,2,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/7_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[1,3,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/8_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[1,4,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/5_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")
data[1,5,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/6_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")
data[1,6,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/7_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")
data[1,7,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/8_100lipids2headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")

data[2,0,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/9_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[2,1,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/10_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[2,2,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/11_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[2,3,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/12_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_int.bin", "rb")
data[2,4,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/9_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[2,5,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/10_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[2,6,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/11_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[2,7,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/12_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_int.bin", "rb")
data[2,0,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/9_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[2,1,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/10_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[2,2,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/11_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[2,3,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/12_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000000steps_double.bin", "rb")
data[2,4,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/9_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")
data[2,5,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/10_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")
data[2,6,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/11_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")
data[2,7,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/12_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67494solvent2000001steps_double.bin", "rb")

data[3,0,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/13_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[3,1,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/14_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[3,2,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/15_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[3,3,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/16_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[3,4,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/13_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[3,5,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/14_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[3,6,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/15_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[3,7,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/16_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[3,0,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/13_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[3,1,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/14_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[3,2,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/15_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[3,3,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/16_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[3,4,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/13_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")
data[3,5,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/14_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")
data[3,6,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/15_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")
data[3,7,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/16_100lipids2headmonNum1pHsemonNum3tailmonNum2branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")

data[4,0,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/17_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[4,1,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/18_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[4,2,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/19_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[4,3,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/20_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_int.bin", "rb")
data[4,4,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/17_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[4,5,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/18_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[4,6,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/19_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[4,7,0] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/20_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_int.bin", "rb")
data[4,0,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/17_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[4,1,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/18_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[4,2,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/19_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[4,3,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/20_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000000steps_double.bin", "rb")
data[4,4,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/17_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")
data[4,5,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/18_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")
data[4,6,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/19_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")
data[4,7,1] = open("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control/20_100lipids1headmonNum1pHsemonNum4tailmonNum0branchheadmonTot0branchpHsemonTot2branchtailmonTot100drugs67294solvent2000001steps_double.bin", "rb")

#Fill data arrays________________________________________
TIME = BintoArray(data[0,0,1], 0, 'd', stepTot, stepInc, double_data, double_increment)

plotSize = TIME.size

clustNum_avg = np.zeros((5,plotSize))
drugIn_avg = np.zeros((5,plotSize))
uTot_avg = np.zeros((5,plotSize))
clustLipNum_avg_avg = np.zeros((5,plotSize))
clustRad_avg_avg = np.zeros((5,plotSize))
micRad_avg_avg = np.zeros((5,plotSize))
clustDrugNum_avg_avg = np.zeros((5,plotSize))
drugInRad_div_clustRad_avg_avg = np.zeros((5,plotSize))

for i in xrange(0,structNum):
    clustNum_avg[i] = BintoArray(data[i,0,0], 4, 'i4', stepTot, stepInc, int_data, int_increment)
    for j in xrange(1,runNum):
        clustNum_avg[i] += BintoArray(data[i,j,0], 4, 'i4', stepTot, stepInc, int_data, int_increment)
    clustNum_avg[i] *= 1.0/(1.0*runNum)

for i in xrange(0,structNum):
    drugIn_avg[i] = BintoArray(data[i,0,0], 8, 'i4', stepTot, stepInc, int_data, int_increment)
    for j in xrange(1,runNum):
        drugIn_avg[i] += BintoArray(data[i,j,0], 8, 'i4', stepTot, stepInc, int_data, int_increment)
    drugIn_avg[i] *= 1.0/(1.0*runNum)
 
   
for i in xrange(0,structNum):
    uTot_avg[i] = BintoArray(data[i,0,1], 48, 'd', stepTot, stepInc, double_data, double_increment)
    for j in xrange(1,8):
        uTot[i] += BintoArray(data[i,j,1], 48, 'd', stepTot, stepInc, double_data, double_increment)
    uTot[i] *= 1.0/(1.0*runNum)    

for i in xrange(0,structNum):
    clustLipNum_avg_avg[i] = BintoArray(data[i,0,1], 96, 'd', stepTot, stepInc, double_data, double_increment)
    for j in xrange(1,8):
        clustLipNum_avg_avg[i] += BintoArray(data[i,j,1], 96, 'd', stepTot, stepInc, double_data, double_increment)
    clustLipNum_avg_avg[i] *= 1.0/(1.0*runNum)
    
for i in xrange(0,structNum):
    clustRad_avg_avg[i] = BintoArray(data[i,0,1], 112, 'd', stepTot, stepInc, double_data, double_increment)
    for j in xrange(1,8):
        clustRad_avg_avg[i] += BintoArray(data[i,j,1], 112, 'd', stepTot, stepInc, double_data, double_increment)
    clustRad_avg_avg[i] *= 1.0/(1.0*runNum)

for i in xrange(0,structNum):
    micRad_avg_avg[i] = BintoArray(data[i,0,1],128,'d', stepTot, stepInc, double_data, double_increment)
    for j in xrange(1,runNum):
        micRad_avg_avg[i] += BintoArray(data[i,j,1],128,'d', stepTot, stepInc, double_data, double_increment)  
    micRad_avg_avg[i] *= 1.0/(1.0*runNum)

for i in xrange(0,structNum):
    clustDrugNum_avg_avg[i] = BintoArray(data[i,0,1], 144, 'd', stepTot, stepInc, double_data, double_increment)
    for j in xrange(1,8):
        clustDrugNum_avg_avg[i] += BintoArray(data[i,j,1], 144, 'd', stepTot, stepInc, double_data, double_increment)
    clustDrugNum_avg_avg[i] *= 1.0/(1.0*runNum)

for i in xrange(0,structNum):                      
    drugInRad_div_clustRad_avg_avg[i] = BintoArray(data[i,0,1],160,'d', stepTot, stepInc, double_data, double_increment)
    for j in xrange(1,runNum):
        drugInRad_div_clustRad_avg_avg[i] += BintoArray(data[i,j,1],160,'d', stepTot, stepInc, double_data, double_increment)
    drugInRad_div_clustRad_avg_avg[i] *= 1.0/(1.0*runNum)

#Figure_________________________________________________________________________
f = py.figure(figsize=(7.244094,9.2559055), dpi=1000)

f.suptitle('Micelle formation (Neutral): Control, 2H1P3T, 1H1P4T, (3H)1H1P3T, 1H1P3T(3T)')

subplots = np.array([f.add_subplot(4,2,1), f.add_subplot(4,2,1),
                     f.add_subplot(4,2,1), f.add_subplot(4,2,1),
                     f.add_subplot(4,2,1), f.add_subplot(4,2,1),
                     f.add_subplot(4,2,1), f.add_subplot(4,2,1)])
                     
subplots[0].set_title('uTot_avg', y=1.0, fontsize=6, fontweight='bold')
subplots[1].set_title('clustNum_avg', y=1.0, fontsize=6, fontweight='bold')
subplots[2].set_title('clustRad_avg_avg', y=1.0, fontsize=6, fontweight='bold')
subplots[3].set_title('micRad_avg_avg', y=1.0, fontsize=6, fontweight='bold')
subplots[4].set_title('clustLipNum_avg_avg', y=1.0, fontsize=6, fontweight='bold')
subplots[5].set_title('drugIn_avg', y=1.0, fontsize = 6, fontweight='bold')
subplots[6].set_title('clustDrugNum_avg_avg', y=1.0, fontsize=6, fontweight='bold')
subplots[7].set_title('drugInRad_div_clustRad_avg_avg', y=1.0, fontsize=6, fontweight='bold')

for i in xrange(0,8):
    subplots[i].set_xlabel('Time', labelpad = 2, fontsize = 7)

subplots[0].set_ylabel('uTot_avg', labelpad = 0, fontsize = 7)
subplots[1].set_ylabel('clustNum_avg', labelpad = 0, fontsize = 7)
subplots[3].set_ylabel('clustRad_avg_avg', labelpad = 0, fontsize = 7)
subplots[4].set_ylabel('micRad_avg_avg', labelpad = 0, fontsize = 7)
subplots[5].set_ylabel('clustLipNum_avg_avg', labelpad = 0, fontsize = 7)
subplots[6].set_ylabel('drugIn_avg', labelpad = 0, fontsize = 7)
subplots[7].set_ylabel('clustDrugNum_avg_avg', labelpad = 0, fontsize = 7)
subplots[8].set_ylabel('drugInRad_div_clustRad_avg_avg', labelpad = 0, fontsize = 7)


#Plot
subplots[0].plot(TIME, uTot_avg[0], label = "Control")
subplots[0].plot(TIME, uTot_avg[1], label = "2H1P3T")   
subplots[0].plot(TIME, uTot_avg[2], label = "1H1P4T") 
subplots[0].plot(TIME, uTot_avg[3], label = "(3H)1H1P3T")
subplots[0].plot(TIME, uTot_avg[4], label = "1H1P3T(3T)")

subplots[1].plot(TIME, clustNum_avg[0], label = "Control")
subplots[1].plot(TIME, clustNum_avg[1], label = "2H1P3T")   
subplots[1].plot(TIME, clustNum_avg[2], label = "1H1P4T") 
subplots[1].plot(TIME, clustNum_avg[3], label = "(3H)1H1P3T")
subplots[1].plot(TIME, clustNum_avg[4], label = "1H1P3T(3T)")

subplots[2].plot(TIME, clustRad_avg_avg[0], label = "Control")
subplots[2].plot(TIME, clustRad_avg_avg[1], label = "2H1P3T")   
subplots[2].plot(TIME, clustRad_avg_avg[2], label = "1H1P4T") 
subplots[2].plot(TIME, clustRad_avg_avg[3], label = "(3H)1H1P3T")
subplots[2].plot(TIME, clustRad_avg_avg[4], label = "1H1P3T(3T)")

subplots[3].plot(TIME, micRad_avg_avg[0], label = "Control")
subplots[3].plot(TIME, micRad_avg_avg[1], label = "2H1P3T")   
subplots[3].plot(TIME, micRad_avg_avg[2], label = "1H1P4T") 
subplots[3].plot(TIME, micRad_avg_avg[3], label = "(3H)1H1P3T")
subplots[3].plot(TIME, micRad_avg_avg[4], label = "1H1P3T(3T)")

subplots[4].plot(TIME, clustLipNum_avg_avg[0], label = "Control")
subplots[4].plot(TIME, clustLipNum_avg_avg[1], label = "2H1P3T")   
subplots[4].plot(TIME, clustLipNum_avg_avg[2], label = "1H1P4T") 
subplots[4].plot(TIME, clustLipNum_avg_avg[3], label = "(3H)1H1P3T")
subplots[4].plot(TIME, clustLipNum_avg_avg[4], label = "1H1P3T(3T)")

subplots[5].plot(TIME, drugIn_avg[0], label = "Control")
subplots[5].plot(TIME, drugIn_avg[1], label = "2H1P3T")   
subplots[5].plot(TIME, drugIn_avg[2], label = "1H1P4T") 
subplots[5].plot(TIME, drugIn_avg[3], label = "(3H)1H1P3T")
subplots[5].plot(TIME, drugIn_avg[4], label = "1H1P3T(3T)")

subplots[6].plot(TIME, clustDrugNum_avg_avg[0], label = "Control")
subplots[6].plot(TIME, clustDrugNum_avg_avg[1], label = "2H1P3T")   
subplots[6].plot(TIME, clustDrugNum_avg_avg[2], label = "1H1P4T") 
subplots[6].plot(TIME, clustDrugNum_avg_avg[3], label = "(3H)1H1P3T")
subplots[6].plot(TIME, clustDrugNum_avg_avg[4], label = "1H1P3T(3T)")

subplots[7].plot(TIME, drugInRad_div_clustRad_avg_avg[0], label = "Control")
subplots[7].plot(TIME, drugInRad_div_clustRad_avg_avg[1], label = "2H1P3T")   
subplots[7].plot(TIME, drugInRad_div_clustRad_avg_avg[2], label = "1H1P4T") 
subplots[7].plot(TIME, drugInRad_div_clustRad_avg_avg[3], label = "(3H)1H1P3T")
subplots[7].plot(TIME, drugInRad_div_clustRad_avg_avg[4], label = "1H1P3T(3T)")

for i in xrange (0, 8): 
    subplots[i].legend(loc = 2, fontsize=5)
        
#f.subplots_adjust(wspace=0.2, hspace=0.6, left=0.05, right=0.96, top=0.96, bottom=0.08)
#f.tight_layout()

py.savefig("C:/Users/#ballin/Desktop/thesis_kaust/data/neutral/control.png", dpi = 1000)

py.show()

    
