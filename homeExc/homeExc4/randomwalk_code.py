import numpy as np
import matplotlib.pyplot as plt
import random


def dist(x,y):
    return np.sqrt(x**2+y**2)

def walk(tf):
    T=np.arange(0,tf,1)
    X=np.zeros(T.shape)
    Y=np.zeros(T.shape)
    fwd_bwd=[-1,1]
    maxstep=0
    for i in range(1,len(T)):
        x_or_y=random.randint(0,1)
        step=random.choice(fwd_bwd)
        if x_or_y==0:
            X[i]=X[i-1]+step
            Y[i]=Y[i-1]
        else:
            Y[i]=Y[i-1]+step
            X[i]=X[i-1]
            
        if abs(X[i])>maxstep or abs(Y[i])>maxstep:
            maxstep=max(abs(X[i]),abs(Y[i]))
            
    return X,Y,maxstep

tf=10000
X,Y,scale=walk(tf)

plt.figure(1)
plt.plot(X,Y)
plt.ylim(-scale, scale/3)
plt.xlim(-scale/3, scale/3)
plt.savefig("onewalk.pdf")
plt.show()

Steps=np.arange(10,1000,10)
avdist2=np.empty(Steps.shape)

for i in range(len(Steps)):
    print(Steps[i])
    for j in range(100):
        X,Y,scale=walk(Steps[i])
        R2=dist(X[-1],Y[-1])**2
        avdist2[i]+=R2
    avdist2[i]/=100
    
plt.figure(2)
plt.plot(Steps,np.sqrt(avdist2),label='√<R²>')
plt.plot(Steps,np.sqrt(Steps),label='√N')
plt.xlabel("N")
plt.savefig("comparison_walk.pdf")
plt.legend()
