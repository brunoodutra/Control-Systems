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
num=[1 ,-0.1]
den=[1,-1.7, 0.72]

Gs= tf(num,den)
t ,y =step(Gs)
u= np.ones([len(y)])
plt.plot(t,y,t,u)

