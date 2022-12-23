import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv("hpcscan.perf.Grid.log", sep=" ", header=None)
#data = pd.read_csv("hpcscanPerfGridProto27Baseline.log", sep=" ", header=None)
#data = pd.read_csv("hpcscanPerfGridProto27NEC.log", sep=" ", header=None)

print(data.shape)

y_names = []
y_names.append('Fill')
y_names.append('MaxErr')
y_names.append('L1Err')
y_names.append('SumAbs')
y_names.append('SumAbsDiff')
y_names.append('Max')
y_names.append('Min')
y_names.append('UptPressure')

fig, axis = plt.subplots(1,2,figsize=(14,6))

all_y = []
j=0
for i in range(13,28+1,2):  
    all_y.append(((data[i].tolist()),y_names[j]))
    j+=1

for i in range(8):
        line1 = axis[0].bar(all_y[i][1],height=all_y[i][0][0], color='steelblue', width = -0.3,align='edge')
        line2 = axis[0].bar(all_y[i][1],height=all_y[i][0][1], color='firebrick', width = 0.3, align='edge')

line1.set_label('500x500x500') #legend contains only these 2
line2.set_label('1000x1000x1000') #legend contains only these 2
axis[0].set_ylabel('GByte/s')
axis[0].set_title('Test Case Grid')
axis[0].legend()
axis[0].grid()
axis[0].set_axisbelow(True) # to make the grid lines go behind


## Next image

all_y = []
j=0
for i in range(14,29+1,2):
    all_y.append(((data[i].tolist()),y_names[j]))
    j+=1

for i in range(8):
    line1 =axis[1].bar(all_y[i][1],height=all_y[i][0][0], color='steelblue', width = -0.3,align='edge')
    line2 =axis[1].bar(all_y[i][1],height=all_y[i][0][1], color='firebrick', width = 0.3,align='edge')  #linewidth=1, edgecolor='black'
    
line1.set_label('500x500x500') #legend contains only these 2
line2.set_label('1000x1000x1000') #legend contains only these 2
axis[1].set_ylabel('GPoint/s')
axis[1].set_title('Test Case Grid')
axis[1].legend()
axis[1].grid()
axis[1].set_axisbelow(True) # to make the grid lines go behind


plt.tight_layout()
plt.savefig("hpcscanGridPlotData.png")
plt.show()
plt.close()
