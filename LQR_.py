# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:14:28 2020
obs falta fazer  a matriz aumentada 
@author: bruno
"""

#!/usr/bin/env python3
# -- coding: utf-8 --
#%% Importing required libraries
import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *  # Load the controls systems library
from control import *
import dlqr as d
#%% Process description and model
Ks = 0.9; # stead-state gain in V/V
zeta = 0.3; # damping factor
wn = 2.0; # natural frequency in rad/s
Gs =tf(Ks*wn*2,[1, 2*zeta*wn, wn*2]) ;
print('G(s) = \n', Gs);

# Discrete model
Ts = 0.1;
Gz = c2d(Gs,Ts,'zoh');
print('G(z) = \n',Gz);

Az = Gz.den[0][0];
print('A(z) = \n',Az);
a1 = Az[1]; a2 = Az[2];

Bz = Gz.num[0][0];
print('B(z) = \n',Bz);
b0 = Bz[0]; b1 = Bz[1];

#%% incremental SS
G_SS=tf2ss(Gz);
A=G_SS.A; B=G_SS.B; C=G_SS.C; D=G_SS.D;
Aa=np.matrix(np.identity(3));
Aa[0:2,0:2]=A;
Aa[2,0:2]=C*A;

Ba=np.matrix(np.zeros([3,1]));
Ba[0:2]=B;
Ba[2:3]=C*B;

Ca=np.matrix(np.zeros([3]));
Ca[0,2]=1;
#%% LQR State Space model 
Q=np.identity(3);
R=np.ones([1,1]);
K, X, eigVals=d.dlqr(Aa,Ba,Q,R);
#K=np.matrix([1.3103  , -0.7264  ,  0.4254])
Kx=K[0,0:2]; Ky=K[0,2];
print('K_LQR = \n',K);
#%% Generating a reference sequence
tfinal = 30; # simulation final time in seconds
N = round(tfinal/Ts);
yr = np.ones([N],float);
yr[0:10] = np.zeros([10],float);
#N = np.shape(yr)[1];

#%% Generating a White Gaussian Noise sequence
# Example: xi = np.random.normal(mean,std deviation,N_samples)
variance = 1e-4;
xi = np.random.normal(0.0,np.sqrt(variance),[N]);

# Defining I/O arrays, error vector, control increment vector
y = np.zeros([N],float);
u = np.zeros([N],float); du = np.zeros([N],float);
e = np.zeros([N],float);
x = np.zeros([N,2],float);

#%% Control-loop
for k in range(N):
 x[k]=A.dot(x[k-1]) + B.dot(u[k-1]).T;
 y[k]=C.dot(x[k]);
 e[k]=y[k]-yr[k];
 dx=x[k]-x[k-1];
 u[k]=u[k-1]  -Kx.dot(dx)  +Ky*-e[k-1];
 # Generating a time vector
t = np.arange(0,N*Ts,Ts);

#%% Ploting results
plt.figure();
plt.subplot(211);
plt.plot(t.T,yr.T,':k',t.T,y.T,'b');
plt.legend(['$y_r(t)$', '$y(t)$']);
plt.ylabel('Output [V]');
plt.title('Feedback control results');

plt.subplot(212);
plt.plot(t.T,u.T,'b');
plt.ylabel('Control [V]');
plt.xlabel('Time [s]');