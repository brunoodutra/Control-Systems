# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 15:19:13 2019
https://web.math.princeton.edu/~cwrowley/python-control/intro.html
@author: bruno
"""
import control 
from control.matlab import *
import matplotlib.pyplot as plt 
import numpy as np


num=[16]
den=[1,2*4*0.7,16]

Gs= tf(num,den)
y,t =step(Gs)


u= np.ones([len(y)])
u[0]=0

plt.plot(t,y,t,u)

Y=y; U=u;
#def MQR(Y,U,Ts):
nu=2; ny=2;
K=np.max([ny ,nu]);
nit=len(Y);
##  Mínimos quadrados  recursivo (MQR)
#lam = input('Entre com o fator de esquecimento 0.9<= lam <=1:  ');
lam=1
ordem=ny+nu;
P=1000*np.eye(ordem);
theta=np.zeros(ordem);
I=np.eye(ordem,ordem);
phi=np.zeros(ordem);
yhat=np.zeros(nit);

for k in range(2,nit):
    phi = np.array([-Y[k-1], -Y[k-2], U[k-1], U[k-2]]);
    
    yhat[k]= phi.dot(theta); 
    ek=y[k]-yhat[k];
    
    denominador= lam + phi.dot(P).dot(phi);
    K=(1/denominador)* (P.dot(phi));
    theta= theta + K*ek;  
    
    phi.shape=(ordem,1);
    K.shape=(ordem,1);
     
    P = (I - K.dot(phi.T) ).dot( P/lam)
    P_prace=P.trace
    
a1=theta[0];
a2=theta[1];
b0=theta[2];
b1=theta[3];
 


ys=np.zeros(nit);
es=np.zeros(nit);
for k in range(3,nit):
    ys[k]=-a1*y[k-1]-a2*y[k-2]+ b0*u[k-1] + b1*u[k-2];  
    es[k]=y[k]-ys[k];
    
plt.Figure()
plt.suptitle('Identificação MQR')
plt.plot(t,y,'r',linewidth=3)
plt.plot(t,ys,'b',linewidth=2)
plt.ylabel('tensão(v)')
plt.xlabel('tempo(s)')
plt.grid()
plt.legend(['sinal de entrada','real','Identificado'])
plt.show()
