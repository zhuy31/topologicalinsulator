import matplotlib.pyplot as plt  
import numpy as np

range1 = 4
range2 = 24
prsc = 4

with open('/home/yuming/dev/research/stronglydisordered/data1.dat') as f:
    lines = f.readlines()

chern = []
for line in lines:
    chern.append(float(line))

chern = np.array(chern)
n = 0
chern = chern[n*(2*(range1*prsc)+1)*(range2*prsc+1):(n+1)*(2*(range1*prsc)+1)*(range2*prsc+1)]
chern = chern.reshape(2*(range1*prsc)+1,range2*prsc+1)

plt.imshow(chern,cmap='turbo')
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
