import numpy as np
import matplotlib.pyplot as plt
import time

#bissection function finding one root in one specific interval
def bissection(a,b,f,maxit):
    if f(a)*f(b) >= 0:
        print("Doesn't change sign on interval")
        return None
    c=(a+b)/2
    for i in range(maxit):
        c=(a+b)/2
        if f(c)*f(a)<0:
            b=c
        elif f(c)*f(b)<0:
            a=c
        elif f(c)==0.0:
            print("Exact solution found.")
            return c
    return (a+b)/2

def multipleroots(a,b,f,numroots):
    for scale in range(4,100,2):
        X=np.linspace(a,b,numroots*scale)
        x0=np.zeros((numroots,))
        j=0
        for i in range(1,len(X)):
            if f(X[i-1])*f(X[i])<0:
                x0[j]=bissection(X[i-1],X[i],f,20)
                j+=1
                print(j)
            if f(X[i-1])*f(X[i]) ==0:
                if f(X[i-1]) == 0:
                    x0[j]=X[i-1]
                    j+=1
                    
                else:
                    x0[j]=X[i]
                    j+=1
        if j==numroots:
            return x0
    print("Number of roots not consistent with function provided.")
    return None

# Legendre polynomial 
def P(x,n):  
    if(n == 0): 
        return 1 # P0 = 1 
    elif(n == 1): 
        return x # P1 = x 
    else: 
        return ((2*n-1)*x*P(x, n-1)-(n-1)*P(x, n-2))/(n)

# Derivative of Legendre Polynomial
def P_prime(x,n):
    if(n==0):
        return 0
    elif(n==1):
        return 1
    else:
        return (-n*x*P(x,n)+n*P(x,n-1))/(1-x**2)
        
def simpsons2(f, a, b, n):
    h=(b-a)/n
    k=0.0
    x=a + h
    for i in range(1,int(n/2) + 1):
        k += 4*f(x)
        x += 2*h

    x = a + 2*h
    for i in range(1,int(n/2)):
        k += 2*f(x)
        x += 2*h
    return (h/3)*(f(a)+f(b)+k)



def bode(a,b,N,f):
    if (N-1)%4 != 0:
        print("N is not a multiple of the form 4p+1")
        return None
    s=0
    h=(b-a)/N
    for i in range(0,N-3,4):
        integ=(7*f(a+i*h)+32*f(a+(i+1)*h)+12*f(a+(i+2)*h)+32*f(a+(i+3)*h)+7*f(a+(i+4)*h))
        integ*=(2*h)/45
        s+=integ
    return s

f = lambda x:np.sqrt(1-x**2)
a=-1
b=1

def gausslegendre(f,N):
    a=-1
    b=1
    x=np.linspace(-1,1,1000)
    PN= lambda x:P(x,N)
    weight = lambda x,N:2*(1-x**2)/((N+1)**2*(P(x,N+1))**2)
    
    roots=multipleroots(a,b,PN,N)
    
    I=0
    for xn in roots:
        I+=weight(xn,N)*f(xn)
    return I

N=np.arange(5,25,4)
integ=np.empty(N.shape)
integs=np.empty(N.shape)
integb=np.empty(N.shape)

t=np.empty(N.shape)
ts=np.empty(N.shape)
tb=np.empty(N.shape)

for i in range(len(N)):
    start = time.process_time()
    integ[i]=gausslegendre(f,N[i])
    t[i]=time.process_time() - start
    
    start = time.process_time()
    integs[i]=simpsons2(f,a,b,N[i])
    ts[i]=time.process_time() - start
    
    start = time.process_time()
    integb[i]=bode(a,b,N[i],f)
    tb[i]=time.process_time() - start

plt.figure(1)
plt.semilogy(N,abs(integ-np.pi/2),label="Gauss-Legendre")
plt.semilogy(N,abs(integs-np.pi/2),label="Simpsons")
plt.semilogy(N,abs(integb-np.pi/2),label="Bode")
plt.xlabel("N")
plt.ylabel("Error")
plt.legend()
plt.savefig("compare_integrals.pdf")

plt.figure(2)
plt.semilogy(N,t,label="Gauss-Legendre")
plt.semilogy(N,ts,label="Simpsons")
plt.semilogy(N,tb,label="Bode")
plt.xlabel("N")
plt.ylabel("Runtime")
plt.legend()
plt.savefig("runtime_integrals.pdf")
