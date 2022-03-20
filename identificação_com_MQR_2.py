# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 22:55:55 2019

@author: bruno
"""

import control 
from control.matlab import *
import matplotlib.pyplot as plt 
import numpy as np
import Ident as id

num=[16]
den=[1,2*4*0.7,16]

Gs= tf(num,den)
y,t =step(Gs)


u= np.ones([len(y)])
u[0]=0

plt.plot(t,y,t,u)

Y=y; U=u;
nit=len(y)
#theta=id.MQRE(y,u)
theta=id.MQR(y,u)

a1=theta[0];
a2=theta[1];
b0=theta[2];
b1=theta[3];
#c1=theta[4];
#c2=theta[5];
 


ys=np.zeros(nit);
es=np.zeros(nit);
for k in range(3,nit):
    ys[k]=-a1*Y[k-1]-a2*Y[k-2]+ b0*U[k-1] + b1*U[k-2];
    #+es[k]+ c1*es[k-1] + c2*es[k-2];  
    es[k]=Y[k]-ys[k];
    
plt.Figure()
plt.suptitle('Identificação MQR')
plt.plot(t,y,'r',linewidth=3)
plt.plot(t,ys,'b',linewidth=2)
plt.ylabel('tensão(v)')
plt.xlabel('tempo(s)')
plt.legend(['sinal de entrada','real','Identificado'])