import matplotlib.pyplot as plt  
import numpy as np

range1 = 2
range2 = 12
prsc = 4

with open('cherndata1d.dat') as f:
    lines = f.readlines()

chern = []
for line in lines:
    chern.append(float(line))  

chern = np.array(chern)
chern = chern.reshape(2*(range1*prsc)+1,range2*prsc+1)

plt.imshow(chern,cmap='turbo',vmin=0,vmax=1)
locs = np.arange(0,range2*prsc,step=prsc)
labels = [i for i in range(0,range2)]
plt.xticks(locs, labels)
locs = np.arange(0,(range1*prsc+1)*2,step=prsc)
labels = [i for i in range(-range1,range1+1)]
plt.yticks(locs, labels)
plt.colorbar()
plt.xlabel("Disorder")
plt.ylabel("Mass")
plt.show()
