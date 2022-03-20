# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 14:59:34 2020

@author: bruno
"""

#!/usr/bin/env python3
# -- coding: utf-8 --
#%% Importing required libraries
import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *  # Load the controls systems library

#%% Process description and model
Ks = 0.9; # stead-state gain in V/V
zeta = 0.3; # damping factor
wn = 2.0; # natural frequency in rad/s
Gs = tf(Ks*wn*2,[1, 2*zeta*wn, wn*2]);
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

#%% Digital I-PD control (using Backward approx.)
Kp = 2.0; Ki = 2.0; Kd = 1;
s0 = Kp+Ki*Ts+Kd/Ts;
s1 = -Kp-2*Kd/Ts;
s2 = Kd/Ts;
t0 = Ki*Ts;

#%% Generating a reference sequence
tfinal = 10; # simulation final time in seconds
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

# Control-loop
for k in range(N):
    y[k] = -a1*y[k-1] -a2*y[k-2] +b0*u[k-1] +b1*u[k-2] +xi[k];
    e[k] = yr[k] -y[k];
    du[k] = t0*yr[k] -s0*y[k] -s1*y[k-1] -s2*y[k-2];
    u[k] = u[k-1] +du[k];

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