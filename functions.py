# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 00:36:38 2019

@author: bruno
"""
import numpy as np
import matplotlib.pyplot as plt 
import Ident as id

def Norm(x,amplitude_MAX, amplitude_MIN):     
# Normalização de dados
        M=np.max(x);
        M1=np.min(x);
        x=(x -M1)*(amplitude_MAX-amplitude_MIN)/(M-M1) + amplitude_MIN;
    
        return x
def Norms(x,amplitude_MAX, amplitude_MIN,MAX_MIN):     
# Normalização de dados
        M=MAX_MIN[0];
        M1=MAX_MIN[1];
        x=(x -M1)*(amplitude_MAX-amplitude_MIN)/(M-M1) + amplitude_MIN;
    
        return x

def Mav(x, Win , incre):
    vec = np.arange(Win,len(x),incre);
    tempsimu=vec[len(vec)-1];
    j=0; k=0; X=0;  tam=1+ np.round((tempsimu-Win)/incre);
    i=0;
    #MAV=np.zeros((tam.astype(np.int64),1));
    MAV=np.zeros(1);
    while i < tempsimu:
        j=j+1;
        #print("i:"+ str(i))
        X=x[i]+X
       
        
        i=i+1;
        
        if j==Win:
            #print("x antes:"+ str(abs(X)/Win))
            #print("MAV"+ str(abs(X)/Win))
            MAV=np.append(MAV, abs(X)/Win)
            #print("MAV"+ str(abs(X)/Win))
            k=k+1
            #print("k:"+ str(k))
            #print("len mav:"+ str(np.size(MAV)))
            X=0
            j=0
            i=incre*k
            
           # print("i depois:"+ str(i))
            
    return MAV

def Mav_2(x, Win):    
    Mav=np.sum(np.abs(x))/Win;
        
    return Mav

def Rms(x, Win , incre):
    vec = np.arange(Win,len(x),incre);
    tempsimu=vec[len(vec)-1];
    j=0; k=0; X=0;  tam=1+ np.round((tempsimu-Win)/incre);
    i=0;
    #MAV=np.zeros((tam.astype(np.int64),1));
    RMS=np.zeros(1);
    #print(tam)
    while i < tempsimu:
        #print(i)
       # print(j)
        j=j+1;
        X=x[i]**2+X;
        i=i+1;
    
      
        if j>=Win:
            k=k+1
            RMS=np.append(RMS, np.sqrt(X/Win))
            X=0; j=0;
            i=incre*k

    return RMS
def Rms_2(x, Win):
    Rms=np.sqrt(np.sum(np.array(x)**2)/Win);
    return Rms

def Wl(x, Win , incre):
    vec = np.arange(Win,len(x),incre);
    tempsimu=vec[len(vec)-1];
    j=0; k=0; X=0;  tam=1+ np.round((tempsimu-Win)/incre);
    i=1;

    WL=np.zeros(1);
    while i < tempsimu:
        j=j+1;
        X=np.append(X,(x[i]+x[i-1]))
        i=i+1;
        
        if j==Win:
            WL=np.append(WL, sum(X))
            k=k+1
            X=0
            j=0
            i=incre*k         
    return WL
def Wl_2(x, Win):
        Wl=np.sum(np.diff(x));
        return Wl
def Var(x, Win , incre):
    vec = np.arange(Win,len(x),incre);
    tempsimu=vec[len(vec)-1];
    j=0; k=0; X=0;  tam=1+ np.round((tempsimu-Win)/incre);
    i=0;
    VAR=np.zeros(1);
    while i < tempsimu:
        j=j+1;
        X=np.append(X,x[i])
       # print("i:"+ str(i))            
        i=i+1;
        if j==Win:
            k=k+1
            VAR=np.append(VAR,np.var(X))
            j=0;
            i=incre*k
            X=np.zeros(1)
           # print("i depois:"+ str(i))
            
    return VAR

def Mov(x, Win , incre):
    vec = np.arange(Win,len(x),incre);
    tempsimu=vec[len(vec)-1];
    j=0; k=0; X=0;  tam=1+ np.round((tempsimu-Win)/incre);
    i=0;
    #MOV=np.zeros((tam.astype(np.int64),1));
    MOV=np.zeros(1);
    while i < tempsimu:
        j=j+1;
        X=x[i]+X
        i=i+1;
        if j==Win:
           # MOV[k]=x[i]
            MOV=np.append(MOV, x[i])
            k=k+1
            X=0
            j=0
            i=incre*k
            
    return MOV
#%%
def ARs(y, ny, Win, incre):
    vec = np.arange(Win,len(y),incre);
    model_AR=np.zeros(ny*np.shape(y)[1])
    THETA=np.zeros(ny*np.shape(y)[1])

    for i in vec:
        sampleEMG=y[i-Win:i,:]
        for j in range(np.shape(sampleEMG)[1]):
            THETA[(j)*ny :(j+1)*ny]=  id.MQ_AR(sampleEMG[:,j],ny)
        model_AR=np.c_[model_AR,THETA];
        
    return model_AR
#%%    
def MSE(y,ys): # Erro Quadrático 
    e1=y-ys;
    emqv = np.sum(e1**2)/np.size(e1); # EMQ 
    return emqv

def R_2(y,ys): # Erro Quadrático 
    e_AR = y-ys; # Erro da saída real e a saída estimada
    e2_AR = e_AR**2; # Erro quadrático
    SEQ_AR = np.sum(e2_AR); # Somatório dos erros quadráticos

## Coeficiente de Correlação Múltipla RMS
    m1_AR = np.sum((y-np.mean(ys))**2); # Soma a diferença entre cada amostra e a média das amostras elevadas ao quadrado armazenadas
    R2_1 = (1-(SEQ_AR/m1_AR)) *100; # Coeficiente de correlação múltipla
    return R2_1

def decimalToBinary(n,byte):
    n=bin(n).replace("0b", "")
    #print(n)
    code=np.zeros([byte])
    for i in range(len(n)):
        code[(byte-len(n))+i]=n[i]
    return code

def to_bin_hot(labels, dimension):
    Bin_labels=np.zeros([len(labels),dimension])
    for i in range(len(labels)):
        Bin_labels[i]= decimalToBinary(labels[i],dimension)
    return Bin_labels

def to_one_hot(labels, dimension):
    results = np.zeros((len(labels), dimension))
    for i, label in enumerate(labels):
        results[i, int(label)] = 1.
    return results
#%%
def removelag(y,lag):
    new_y=y[0:lag];
    new_y=np.append(new_y,y[0:len(y)-lag])
    
    return new_y