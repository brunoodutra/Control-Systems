# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 21:32:55 2021

@author: Bruno Dutra
"""

import numpy as np
from sklearn import svm, metrics
from sklearn.metrics import confusion_matrix
from numpy.linalg import inv
import matplotlib.pyplot as plt 
import functions as func
import scipy.io
import scipy.optimize as opt
import daqduino


#%% Importa os dados do controlador sintonizado e as matrizes do sistema
mat = scipy.io.loadmat('controle_protese_MISO_fk.mat')     
A=mat['A'];  B=mat['B']; C=mat['C']; Cfk=mat['Cfk'];
Kx=mat['Kx']; Ky=mat['Ky'];
Ts=mat['Ts']; Ts=float(Ts);
L=mat['Lfk'];

N_inputs=np.shape(B)[1]
N_outputs=C.shape[0]

#%% Gera a sequencia de referencia 
tfinal = 30; # tempo de duração da simulação
N = int(np.round(tfinal/Ts));
r = np.ones([N,N_outputs],float);
r[11:N,0]= 10;
r[0:10] = 0
#N = np.shape(yr)[1];
#%% Newtonw para volts

yr=-(0.9461237455400879 - (2251799813685248*(2.0809623073486785*r/10 + 0.9690904666290612)**(1/2))/2342955267986889);

#%% Gera o ruído gaussiano
# Example: xi = np.random.normal(mean,std deviation,N_samples)
variance = 1e-2;
xi = np.random.normal(0.0,np.sqrt(variance),[N,N_outputs]);

# Definição inicial das variáveis de controle 
y = np.zeros([N,N_outputs],float);
yhat = np.zeros([N,N_outputs],float);
u = np.zeros([N,N_inputs],float); du = np.zeros([N],float);
e = np.zeros([N,N_outputs],float);
x = np.zeros([N,N_inputs],float);
Xf= np.zeros([N,x.shape[1]+y.shape[1]],float);
xa= np.zeros([N,N_inputs],float);
#%% Control-loop
daqduino.start('COM14', 250000) # Open a Serial Port Object

for k in range(N):
    
    
 x[k]=A.dot(x[k-1]) + B.dot(u[k-1]).T;
 #y[k]=C.dot(x[k]) +xi[k];
 Y= daqduino.read() # Measured Output
 y[k]=Y[-1];
 Y_fusion=[Y[7] ,Y[5] ,Y[4], Y[6]]
 e[k]=y[k]-yr[k];
 
 #observador
 xa[k]=((A- L.dot(Cfk)).dot(xa[k-1])) +B.dot(u[k-1]) +L.dot(Y_fusion);
 yhat[k]=C.dot(xa[k]);
 dx=xa[k]-xa[k-1];
 Xf[k]=np.append(dx,y[k])
 
 u[k]=u[k-1]  -Kx.dot(Xf[k])  +Ky.dot(yr[k]);
 daqduino.write(u[k], Ts) # Control Signal
 
daqduino.end() # Close a Serial Port Object
 # Gera o vetor de tempo
t = np.arange(0,N*Ts,Ts);
#%%Volts para newton 
phi_y = np.array([y[:,0]**2, y[:,0]]);
phi_yhat = np.array([y[:,0]**2, y[:,0]]);
theta =np.array([ 0.520240576837170, 0.984423926278238]);

y[:,0] = phi_y.T.dot(theta)*10;
yhat[:,0] = phi_yhat.T.dot(theta)*10;
#%% Ploting results
plt.figure();
plt.subplot(211);
plt.plot(t.T,r[:,0].T,':k',t.T,y[:,0].T,'r');
plt.plot(t.T,r[:,0].T,':k',t.T,yhat[:,0].T,'k');
#plt.plot(t.T,yr[:,2].T,':k',t.T,y[:,2].T,'g');
#plt.plot(t.T,yr[:,3].T,':k',t.T,y[:,3].T,'k');
plt.legend(['$y_r(t)$', '$y(t)$']);
plt.ylabel('Output [V]');
plt.title('Feedback control results');

plt.subplot(212);
plt.plot(t.T,u[:,0].T,'r');
plt.plot(t.T,u[:,1].T,'b');
plt.plot(t.T,u[:,2].T,'g');
plt.plot(t.T,u[:,3].T,'k');
plt.ylabel('Control [V]');
plt.xlabel('Time [s]');
#%%
