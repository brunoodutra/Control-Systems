# Nelson Yamaguti

import control
from numpy.linalg import matrix_rank
import matplotlib.pyplot as plt
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
import pyqtgraph.exporters
import numpy as np

#import daqduino

# Identified model

t1 = 2; t2 = 4; t3 = 6; Tfinal = 8;
Ts = 0.0135; nit = round(Tfinal/Ts)

A = np.array([[0.945276502034681 , 0.001671177993333],
              [-4.053592441875310, 0.123790962469117]])

B = np.array([[-0.837055492409139 , 0.480539072489219],
              [-62.004110548821536, 35.595486851051128]])
             
C = np.array([1.000000000000000, 0.000000000000000])

# %% Augmented state

# MÉTODO 1 xa[:, k] = np.array([[  y[k]  ],
#                               [dx[:, k]]])

# Aa = np.array([[1, np.dot(C, A)[0], np.dot(C, A)[1]],
#                [0,     A[0, 0]    ,     A[0, 1]    ],
#                [0,     A[1, 0]    ,     A[1, 1]    ]])

# Ba = np.array([[np.dot(C, B)[0], np.dot(C, B)[1]],
#                [    B[0, 0]    ,     B[0, 1]    ],
#                [    B[1, 0]    ,     B[1, 1]    ]])

# Ca = np.array([1, 0, 0])
# Da = 0

# MÉTODO 2 xa[:, k] = np.array([[dx[:, k]],
                              # [  y[k]  ]])

Aa = np.array([[    A[0, 0]    ,     A[0, 1]    , 0],
                [    A[1, 0]    ,     A[1, 1]    , 0],
                [np.dot(C, A)[0], np.dot(C, A)[1], 1]])

Ba = np.array([[    B[0, 0]    ,     B[0, 1]    ],
                [    B[1, 0]    ,     B[1, 1]    ],
                [np.dot(C, B)[0], np.dot(C, B)[1]]])

Ca = np.array([0, 0, 1])
Da = 0

# %% Controllability and Observability analysis

Co = control.ctrb(Aa, Ba)
if matrix_rank(Co) == len(Aa):
    print('The augmented system is controllable.')

Ob = control.obsv(Aa, Ca)
if matrix_rank(Ob) == len(Aa):
    print('The augmented system is observable.')

# %% LQR controller

#               y | Delta*x1 | Delta*x2
Qlqr = np.diag([1 ,    100   ,     1   ]); Rlqr = 1*np.eye(2);
# K, X, eigVals = dlqr(Aa, Ba, Qlqr, Rlqr)

# Recursive Riccati Difference Equation (RDE)
P = 100*np.eye(3)

for k in range(300):
    P = np.dot(np.dot(Aa.T, P), Aa) - np.dot(np.dot(np.dot(np.dot(Aa.T, P), Ba), np.linalg.inv(np.dot(np.dot(Ba.T, P), Ba) + Rlqr)), np.dot(np.dot(Ba.T, P), Aa)) + Qlqr

K = (np.dot(np.dot(np.dot(Aa.T, P), Ba), np.linalg.inv(np.dot(np.dot(Ba.T, P), Ba) + Rlqr))).T

#        y  | Delta*x1 | Delta*x2
# K = [[k11 ,    k12   ,    k13   ], # Delta*u1
#      [k21 ,    k22   ,    k23   ]] # Delta*u2

# %% Kalman Filter

#              y  | Delta*x1 | Delta*x2
Qkf = np.diag([10 ,    10    ,    10   ]); Rkf = 1000

# Recursive Riccati Difference Equation (RDE)
S = 100*np.eye(3);

for k in range(300):
    #S = np.dot(np.dot(Aa, S), Aa.T) - np.dot(np.dot(np.dot(Aa, S), Ca.T), np.dot(np.linalg(np.dot(np.dot(Ca, S), Ca.T)) + Rkf),np.dot(np.dot(Ca, S), Aa.T)) +Qfk
    aux1=np.dot(np.dot(Aa, S), Aa.T)
    aux2=Aa.dot(S).dot(Ca.T).dot(1/(Ca.dot(S).dot(Ca.T) +Rkf))
    aux2.shape=(3,1)
    aux3=(Ca.T.dot(S).dot(Aa.T))
    aux3.shape=(1,3)
    S =aux1-aux2.dot(aux3)+Qkf
    

L = ((np.dot(np.dot(Aa, S), Ca.T))/(np.dot(np.dot(Ca, S), Ca.T) + Rkf)).T

# %% Reference

yr = np.zeros((nit, 1))
yr[1 : round(t1/Ts)] = 30;
yr[round(t1/Ts) : round(t2/Ts)] = -10;
yr[round(t2/Ts) : round(t3/Ts)] = 20;
yr[round(t3/Ts) : nit] = 0;

# % Initial simulation conditions
y = np.zeros((nit, 1)); x = np.zeros((2, nit)); dx = np.zeros((2, nit))
xa = np.zeros((3, nit)); e = np.zeros((nit, 1)); du = np.zeros((nit, 2))
u = np.zeros((nit, 2)); ya = np.zeros((nit, 1))
y_plot = []; u1_plot = []; u2_plot = []

# %% Initial real-time plot settings

app = QtGui.QApplication([])

win = pg.GraphicsWindow(title = "Linear Quadratic Regulator (LQR)")

pg.setConfigOptions(antialias = True) # Enable antialiasing for prettier plots
p1 = win.addPlot(0, 0, title = "Closed-Loop Response") # Adds the plot in the window 211
p1.showGrid(True, True)
p1.setLabel('left', 'Degree', units = '°')
p1.setRange(xRange = [0, nit]) #yRange = [0,60])
curva_resposta = p1.plot()
curva_referencia = p1.plot()
curva_referencia.setData(yr[:, 0], pen = 'w',name = 'Reference')

p2 = win.addPlot(1, 0, title = "Control Action") # Adds the plot in the window 212
p2.showGrid(True, True)
p2.setLabel('left', 'Amplitude', units = 'V')
p2.setLabel('bottom', 'Iteration')
curve_controle1 = p2.plot()
curve_controle2 = p2.plot()
p2.setRange(xRange = [0, nit], yRange = [0, 5])
# p2.disableAutoRange()

pg.setConfigOptions(antialias = True)
pg.setConfigOption('background', 'k')
pg.setConfigOption('foreground', 'w')

# %% Control loop

daqduino.start('COM4', 250000) # Open a Serial Port Object

for k in range(nit):
    
    # MÉTODO 1 xa[:, k] = np.array([[  y[k]  ],
    #                               [dx[:, k]]])
    y[k] = daqduino.read() # Measured Output
    # x[:, k] = np.dot(A, x[:, k-1]) + np.dot(B, u[k-1, :].T)
    # y[k] = np.dot(C, x[:,k])
    # dx[:, k] = x[:, k] - x[:, k-1]
    # xa[0:4, k] = np.array([float(y[k]), dx[0, k], dx[1, k]]).reshape(3, 1)
    # xa[0, k] = y[k]; xa[1, k] = dx[0, k]; xa[2, k] = dx[1, k]
    
    # MÉTODO 2 xa[:, k] = np.array([[dx[:, k]],
    #                               [  y[k]  ]])
    # y[k] = daqduino.read() # Measured Output
    # x[:, k] = np.dot(A, x[:, k-1]) + np.dot(B, u[k-1, :].T)
    # dx[:, k] = x[:, k] - x[:, k-1]
    # xa[0, k] = dx[0, k]; xa[1, k] = dx[1, k]; xa[2, k] = y[k]

    # xa[:, k] = np.dot(Aa, xa[:, k-1]) + np.dot(Ba, du[k-1, :].T)
    # y[k] = np.dot(Ca, xa[:, k])
    
    e[k] = yr[k] - y[k]
    
    xa[:, k] = np.dot(Aa, xa[:, k-1]) + np.dot(Ba, du[k-1, :].T) + np.dot(L.reshape(3, 1), (y[k-1] - ya[k-1]))
    ya[k] = np.dot(Ca, xa[:, k])
    
    du[k, :] = np.dot(K[:, 0].reshape(2, 1), yr[k]) - np.dot(K, xa[:, k])
    # du[k, :] = K[:, 0].dot(yr[k]) - np.dot(K, xa[:, k])
    u[k, :] = u[k-1, :] + du[k, :]
    
    # Control signal saturation
    if u[k, 0] >= 5:
        u[k, 0] = 5
    elif u[k, 0] <= 0:
        u[k,0] = 0
    
    # Control signal saturation
    if u[k, 1] >= 5:
        u[k, 1] = 5
    elif u[k, 1] <= 0:
        u[k, 1] = 0
    
    # Plotting real-time
    y_plot = np.append(y_plot, y[k])
    curva_resposta.setData(y_plot, pen = 'b',name = 'Reference')
    u1_plot = np.append(u1_plot, u[k, 0])
    curve_controle1.setData(u1_plot, pen='b')
    u2_plot = np.append(u2_plot, u[k, 1])
    curve_controle2.setData(u2_plot, pen = 'r')
    daqduino.write_mimo(u[k, 0], u[k, 1], Ts) # Control Signal
    
    app.processEvents()
    
daqduino.end() # Close a Serial Port Object

# %% Plot results

t = np.arange(0, Tfinal, Ts)
plt.figure(2)
plt.subplot(211)
plt.plot(t, yr, ':k', t, y, 'b'); plt.grid()
plt.title('Closed-Loop Response')
plt.ylabel('Angle (°)'); plt.xlabel('Time (s)')
plt.legend(['Reference (y_r)', 'Angle (°)'])
plt.subplot(212)
plt.plot(t, u[:, 0], 'b', t, u[:, 1], 'r'); plt.grid()
plt.title('Control Action')
plt.ylabel('Amplitude'); plt.xlabel('Time (s)')
plt.legend(['Control (u1)', 'Control (u2)'])
plt.tight_layout()

# Performance indices
J_ISE = float(sum(np.square(e))/nit)
print(f'J (ISE) = {J_ISE:.3f}')

J_ISU = sum(sum(np.square(u)))/nit
print(f'J (ISU) = {J_ISU:.3f}')

J = J_ISE + J_ISU
print(f'J = {J:.3f}')

sigm_e = float(sum(np.square(e - (sum(e)/len(e))))/nit)
print(f'sigma(e)=  {sigm_e:.3f}')

sigm_u = np.zeros((2, 1))
sigm_u[0] = sum(np.square(u[:, 0] - (sum(u[:, 0])/len(u[:, 0]))))/nit
sigm_u[1] = sum(np.square(u[:, 1] - (sum(u[:, 1])/len(u[:, 1]))))/nit
soma = float(sum(sigm_u))
print(f'sigm(u) = {soma:.3f}')