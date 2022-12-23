import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv("hpcscan.perf.Memory.log", sep=" ", header=None)
#data = pd.read_csv("hpcscanMemoryShaheen.log", sep=" ", header=None)
#data = pd.read_csv("hpcscanMemoryProto27.log", sep=" ", header=None)

print(data.shape)

x = data[8].tolist()

y_names = []
y_names.append('Fill')
y_names.append('Copy')
y_names.append('Add')
y_names.append('Mult')
y_names.append('AddUpd')

fig, axis = plt.subplots(1,2,figsize=(14,6))

all_y = []
for i in range(13,21+1,2):
    all_y.append(((data[i].tolist()),y_names[int(i%13/2)]))

for i in range(5):
    if(i == 3):
        axis[0].plot(x,all_y[i][0],label=all_y[i][1],marker='o',linestyle='--')
    else:
        axis[0].plot(x,all_y[i][0],label=all_y[i][1],marker='o')


axis[0].set_xlabel('# of threads')
axis[0].set_ylabel('GByte/s')
axis[0].set_title('Test Case Memory / scalability')
axis[0].legend()
axis[0].grid()

## Next image

all_y = []
for i in range(14,22+1,2):
    all_y.append(((data[i].tolist()),y_names[int(i%14/2)]))


for i in range(5):
    if(i == 3):
        axis[1].plot(x,all_y[i][0],label=all_y[i][1],marker='o',linestyle='--')
    else:
        axis[1].plot(x,all_y[i][0],label=all_y[i][1],marker='o')

axis[1].set_xlabel('# of threads')
axis[1].set_ylabel('GPoint/s')
axis[1].set_title('Test Case Memory / scalability')
axis[1].legend()
axis[1].grid()

plt.tight_layout()
plt.savefig("hpcscanMemory.png")
plt.show()
plt.close()
