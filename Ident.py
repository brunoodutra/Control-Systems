# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 22:48:12 2019

@author: bruno
"""
import matplotlib.pyplot as plt 
import numpy as np


def MQR(Y,U):
    #def MQR(Y,U,Ts):
    nu=2; ny=2;
    K=np.max([ny ,nu]);
    nit=len(Y);
    ##  Mínimos quadrados  recursivo (MQR)
    lam = input('Entre com o fator de esquecimento 0.9<= lam <=1:  ');
    lam=float(lam)
    ordem=ny+nu;
    P=1000*np.eye(ordem);
    theta=np.zeros(ordem);
    I=np.eye(ordem,ordem);
    phi=np.zeros(ordem);
    yhat=np.zeros(nit);
    
    for k in range(2,nit):
        phi = np.array([-Y[k-1], -Y[k-2], U[k-1], U[k-2]]);
        
        yhat[k]= phi.dot(theta); 
        ek=Y[k]-yhat[k];
        
        denominador= lam + phi.dot(P).dot(phi);
        K=(1/denominador)* (P.dot(phi));
        theta= theta + K*ek;  
        
        phi.shape=(ordem,1);
        K.shape=(ordem,1);
         
        P = (I - K.dot(phi.T) ).dot( P/lam)
        P_trace=P.trace
        
    a1=theta[0];
    a2=theta[1];
    b0=theta[2];
    b1=theta[3];
     
    
    
    ys=np.zeros(nit);
    es=np.zeros(nit);
    for k in range(3,nit):
        ys[k]=-a1*Y[k-1]-a2*Y[k-2]+ b0*U[k-1] + b1*U[k-2];  
        es[k]=Y[k]-ys[k];
        
    return theta

def MQRE(Y,U):
    #def MQR(Y,U,Ts):
    nu=2; ny=2; nw=2;
    K=np.max([ny ,nu, nw]);
    nit=len(Y);
    ##  Mínimos quadrados  recursivo (MQR)
    lam = input('Entre com o fator de esquecimento 0.9<= lam <=1:  ');
    lam=float(lam)
    ordem=ny+nu+nw;
    P=1000*np.eye(ordem);
    theta=np.zeros(ordem);
    I=np.eye(ordem,ordem);
    phi=np.zeros(ordem);
    yhat=np.zeros(nit);
    W=np.zeros(nit);
    for k in range(2,nit):
        phi = np.array([-Y[k-1], -Y[k-2], U[k-1], U[k-2],W[k-1], W[k-2]]);
        
        yhat[k]= phi.dot(theta); 
        ek=Y[k]-yhat[k];
        
        denominador= lam + phi.dot(P).dot(phi);
        K=(1/denominador)* (P.dot(phi));
        theta= theta + K*ek;  
        
        phi.shape=(ordem,1);
        K.shape=(ordem,1);
         
        # ruído brando, processo estocástico 
        W[k]=Y[k]-(phi.T.dot(theta));
        
        P = (I - K.dot(phi.T) ).dot( P/lam)
        P_trace=P.trace
        
    a1=theta[0];
    a2=theta[1];
    b0=theta[2];
    b1=theta[3];
    c1=theta[4];
    c2=theta[5];
    
    
    ys=np.zeros(nit);
    es=np.zeros(nit);
    for k in range(3,nit):
        ys[k]=-a1*Y[k-1]-a2*Y[k-2]+ b0*U[k-1] + b1*U[k-2]+es[k]+ c1*es[k-1] + c2*es[k-2];  
        es[k]=Y[k]-ys[k];
        
    return theta

def MQR_ident(Y,U):
    #def MQR(Y,U,Ts):
    ny = input('Qual a ordem de y ? ');
    ny=int(ny)
    nu = input('Qual a ordem de u ? ')
    nu=int(nu)
    #nu2 = input('Qual a ordem de u2 ?')
    
    nit=len(Y);
    ##  Mínimos quadrados  recursivo (MQR)
    lam = input('Entre com o fator de esquecimento 0.9<= lam <=1:  ');
    lam=float(lam)
    ordem=ny+nu;
    P=1000*np.eye(ordem);
    theta=np.zeros([ordem,1]);
    I=np.eye(ordem,ordem);
    phi=np.zeros([ordem,1]);
    K=np.zeros([ordem,1]);
    yhat=np.zeros(nit);
    
    for k in range(2,nit):
        
        n=-1;
        if ny>0 :
            for i in range(1,ny+1):
                n=n+1
                phi[n]=-Y[k-i];
                
        n=ny-1    
        if nu>0 :
            for j in range(1,nu+1):
                n=n+1
                phi[n]=U[k-j];
    
        #phi = np.array([-Y[k-1], -Y[k-2], U[k-1], U[k-2]]);
        
        yhat[k]= phi.T.dot(theta); 
        ek=Y[k]-yhat[k];
        
        denominador= lam + phi.T.dot(P).dot(phi);
        K=(1/denominador)* (P.dot(phi));
        
        theta= theta + K*ek;  
             
        P = (I - K.dot(phi.T) ).dot( P/lam)
        P_trace=P.trace
        
    ys=np.zeros(nit);
    es=np.zeros(nit);
    PHI=np.zeros([ordem,1]);
    for k in range(3,nit):
        n=-1;
        if ny>0 :
            for i in range(1,ny+1):
                n=n+1
                PHI[n]=-Y[k-i];
                
        n=ny-1    
        if nu>0 :
            for j in range(1,nu+1):
                n=n+1
                PHI[n]=U[k-j];
                
        ys[k]=PHI.T.dot(theta);  
        es[k]=Y[k]-ys[k];
        

        
    return theta