# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
from struct import *

#Variables______________________________________________________________________

stepTot = 100000
stepInc = 2
int_data = 4
double_data = 22
int_increment = 4*(int_data - 1)
double_increment = 8*(double_data - 1)

def BintoArray(datafile, init, datatype, stetot, steInc, dataNum, inc):
    datafile.seek(init, 0)  
    arr = np.fromfile(datafile, dtype=datatype, count=1, sep="")  
    for i in np.arange (0,stetot):
        if (i%stepInc == 0 and i != 0):
            datafile.seek(np.dtype(datatype).itemsize*dataNum*i + inc, 1) 
            #datafile.seek(inc, 1) 
            arr = np.append(arr, np.fromfile(datafile, datatype, count=1, sep=""))  
    return arr
                         
#Open binary data files with similarities.___________________________________________________________________________________________________________________________
data = np.array([open("C:/Users/#ballin/Desktop/thesis_kaust/data/acidic/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent100000steps_int.bin", "rb"), 
                 open("C:/Users/#ballin/Desktop/thesis_kaust/data/acidic/control/1_100lipids1headmonNum1pHsemonNum3tailmonNum0branchheadmonTot0branchpHsemonTot0branchtailmonTot100drugs67594solvent100000steps_double.bin", "rb")])   

#Fill similarities array________________________________________
uTot = BintoArray(data[1], 48, 'd', stepTot, stepInc, double_data, double_increment)
micRad_avg = BintoArray(data[1],128,'d', stepTot, stepInc, double_data, double_increment)
drugInRad_div_clustRad_avg = BintoArray(data[1],160,'d', stepTot, stepInc, double_data, double_increment)
Time = BintoArray(data[1], 0, 'd', stepTot, stepInc, double_data, double_increment)


#Figure_________________________________________________________________________
f = py.figure(figsize=(16.51,22.86), dpi=80)

f.suptitle('Control Acidic')

subplots = np.array([f.add_subplot(1,1,1)])
#subplots[0].set_title('uTot', y=1.08, fontsize=7)
subplots[0].set_xlabel('Time', labelpad = 7, fontsize = 7)
#subplots[0].set_ylabel('uTot', labelpad = 7, fontsize = 7)

#subplots[0].set_ylim(0,1.0)
#subplots[0].set_yticklabels([0.0,'','','','',1.0], fontsize=6)

'''  
ax.set_xticks(ind+width)
xtickNames = ax.set_xticklabels( ['Group'+str(i) for i in range(1,6)])
plt.setp(xtickNames, rotation=45, fontsize=10)
'''


#subplots[0].plot(Time, uTot, label = "uTot")
subplots[0].plot(Time, micRad_avg, label = "micRad_avg")
subplots[0].plot(Time, drugInRad_div_clustRad_avg, label = "drugInRad_div_clustRad_avg")
subplots[0].legend(loc = 2, fontsize=6)
        
#f.subplots_adjust(wspace=0.5, hspace=0.5)


py.savefig("C:/Users/#ballin/Desktop/thesis_kaust/data/acidic/control/graphs.png")


py.show()



    
