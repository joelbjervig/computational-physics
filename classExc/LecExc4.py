import matplotlib.pyplot as plt
import numpy as np
from math import *

def euler(f,y0,x):
    y=np.zeros(x.shape)
    y[0]=y0
    for i in range(1,len(y)):
        h=x[i]-x[i-1]
        y[i]=y[i-1]+h*f(x[i-1],y[i-1])
    return y

#dfx and dfy are the derivative of f wrt x and y
def taylorexp(f,dfx,dfy,y0,x):
    y=np.zeros(x.shape)
    y[0]=y0
    for i in range(1,len(y)):
        h=x[i]-x[i-1]
        y[i]=y[i-1]+h*f(x[i-1],y[i-1])
        y[i]+=h*h/2*(dfx(x[i-1],y[i-1])+dfy(x[i-1],y[i-1])*f(x[i-1],y[i-1]))
    return y

#f(x,y)=g(x)*y
def implicit(g,y0,x):
    y=np.zeros(x.shape)
    y[0]=y0
    for i in range(1,len(y)):
        h=x[i]-x[i-1]
        y[i]=y[i-1]*(1+0.5*g(x[i-1])*h)/(1-0.5*g(x[i])*h)
    return y
    
    N=np.arange(1,1000,3)
eul=np.zeros(N.shape)
tayl=np.zeros(N.shape)
imp=np.zeros(N.shape)

a=0
b=1
f= lambda x,y:-x*y
y= lambda x:exp(-x*x/2)
y0=1

dfx = lambda x,y:-y
dfy = lambda x,y:-x

g = lambda x:-x
for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    eul[i]=euler(f,y0,x)[-1]
    tayl[i]=taylorexp(f,dfx,dfy,y0,x)[-1]
    imp[i]= implicit(g,y0,x)[-1]

H=(b-a)/N

plt.figure(1)
plt.xlabel("h")
plt.ylabel("error")
plt.plot(H,abs(eul-y(b)),label="euler")
plt.plot(H,abs(tayl-y(b)),label="taylor")
plt.plot(H,abs(imp-y(b)),label="implicit")
plt.legend()

b=3

for i in range(len(N)):
    x=np.linspace(a,b,N[i])
    eul[i]=euler(f,y0,x)[-1]
    tayl[i]=taylorexp(f,dfx,dfy,y0,x)[-1]
    imp[i]= implicit(g,y0,x)[-1]

H=(b-a)/N

plt.figure(2)
plt.xlabel("h")
plt.ylabel("error")
plt.plot(H,abs(eul-y(b)),label="euler")
plt.plot(H,abs(tayl-y(b)),label="taylor")
plt.plot(H,abs(imp-y(b)),label="implicit")
plt.legend()
