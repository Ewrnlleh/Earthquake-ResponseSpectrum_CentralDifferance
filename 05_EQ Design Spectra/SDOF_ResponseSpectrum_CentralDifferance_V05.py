# -*- coding: utf-8 -*-
" Author: @DrAydinDemir  ;  Date: 21.11.2019 " 
" To create a Response Spectrum for a given EQ record by using Central Differance Method "

import numpy as np
import matplotlib.pyplot as plt
from numpy import array as matrix
from numpy import arange as seri
from numpy import pi as π


#%%
# INPUT MOTION (GROUND MOTION)
'==== Read PEER data ===================='
file = open('PEER_DUZCE_BOLU-0.AT2', 'r')  # Name of the file
HeaderNo = 4                     # Number of Header Lines  not to be read
Δt = 0.01                        # Time interval of the record
SF = 1                        # Scale Factor; if unit is in g, multiply with 9.81 m/s2
lines = file.readlines() 
file.close() 
data = []
for line in lines[HeaderNo:]:
    row = line.split()               # Split the rows
    for i in row:                    # Loop in the each row
        data.append(float(i))
acc = SF * matrix([data]).T               # Acceleration data
t = matrix([seri(0,len(acc)*Δt,Δt)]).T   # Time interval

# Plot Ground Motion
fig, xy = plt.subplots()
xy.plot(t,acc, color="g")
xy.set(xlabel='t [s]', ylabel='ag [g]', title='Time - Ground Acceleration Curve')
xy.grid()
plt.show()


#%%
# Function of Central Differance Method
def SDOFCentDiff(acc,Tn,ξ,m,u0,v0,a0):
  ωn = 2*π/Tn              # Natural angular Frequency
  k = m*ωn**2              # System Stiffness	
  c = 2*m*ξ*ωn             # Damping Coefficient
  ḱ = m/Δt**2 + c/(2*Δt)
  y = m/Δt**2 - c/(2*Δt)
  z = k-2*m/Δt**2
  P = -m*acc
  ## Recursive Relationship
  u = [(u0 - Δt*v0 + Δt**2/2*a0), u0]; v = [v0]; a = [a0]
  i = 1
  for Pi in P[1:-1]:
      Ṕ = Pi-y*u[i-1]-z*u[i] 
      ui = Ṕ/ḱ
      vi = (ui-u[i-1])/(2*Δt)
      ai = (ui-2*u[i]+u[i-1])/Δt**2
      u.append(ui)
      v.append(vi)
      a.append(ai)
      i += 1
  return(u, v, a)


" Input Parameters "
ξ = 0.05    # Damping Ratio  
m = 1.0     # Assumption for Mass 
u0 = 0.0    # Intial displacement
v0 = 0.0    # Intial velocity
a0 = 0.0    # Intial acceleration

# Spectrum Periods
Tmax = 5    # Maximum period on the spectrum
ΔT = 0.1    # Period increment/interval of the spectrum
Tn = matrix([seri(ΔT,Tmax+ΔT,ΔT)]).T   # Period interval


" 1) To obtain Acceleration, Velocity and Displacement response of an SDOF system under an EQ motion"
Tni = 2.5    # for this period
u, v, a = SDOFCentDiff(acc,Tni,ξ,m,u0,v0,a0)
u = matrix([u]).T ; v = matrix([v]).T ; a = matrix([a]).T

# Plot response 
fig, grfk = plt.subplots(nrows=3, ncols=1)
grfk[0].plot(t, u, color='k', label="Umax={} m".format(np.around(np.max(u),decimals=4)))
grfk[0].set(ylabel='Disp. [m]', title='Response of an SDOF system under an EQ motion (Tn:{} s)'.format(Tni))
grfk[0].legend()
grfk[0].grid()
grfk[1].plot(t[:-1], v, color='r', label="Vmax={} m/s".format(np.around(np.max(v),decimals=3)))
grfk[1].set(ylabel='Vel. [m/s]')
grfk[1].legend()
grfk[1].grid()
grfk[2].plot(t[:-1], a, color='b', label="Amax={} m/s2".format(np.around(np.max(a),decimals=2)))
grfk[2].set(xlabel='Time [s]', ylabel='Acc. [m/s2]')
grfk[2].legend()
grfk[2].grid()
fig.tight_layout()
plt.show()


" 2) To create response spectra"
Sd = [] ; Sv = [] ; Sa = []
for Ti in Tn:
    u, v, a = SDOFCentDiff(acc,Ti,ξ,m,u0,v0,a0)     # Recursive Relationship
    Sd.append(np.max(np.abs(u)))
    Sv.append(np.max(np.abs(v)))
    Sa.append(np.max(np.abs(a)))
Sd = matrix(Sd)       # Relative Displacement Response Spectrum
Sv = matrix(Sv)       # Relative Velocity Response Spectrum
Sa = matrix(Sa)       # Relative Acceleration Response Spectrum

# Plot Graphs
fig, grfk = plt.subplots(nrows=3, ncols=1)
grfk[0].plot(Tn, Sd, color='g')
grfk[0].set(ylabel='Sd [m]', title="Displacement Response Spectrum")
grfk[0].grid()
grfk[1].plot(Tn, Sv, color='r')
grfk[1].set(ylabel='Sv [m/s]', title="Velocity Response Spectrum")
grfk[1].grid()
grfk[2].plot(Tn, Sa, color='b')
grfk[2].set(xlabel='Tn [s]', ylabel='Sa [m/s2]', title="Acceleration Response Spectrum")
grfk[2].grid()
fig.tight_layout()
plt.show()


" 3) To create Pseudo response spectra"  
wn = 2*π/Tn
PSv = Sd * wn           # Pseudo Velocity Response Spectrum
PSa = Sd * wn**2        # Pseudo Acceleration Response Spectrum
  
# Plot Graphs
fig, grfk = plt.subplots(nrows=3, ncols=1)
grfk[0].plot(Tn, Sd, color='r', label="Sd")
grfk[0].set(ylabel='Sd [m]', title="Displacement Response Spectrum")
grfk[0].legend()
grfk[0].grid()
grfk[1].plot(Tn, Sv, color='r', label="Sv")
grfk[1].plot(Tn, PSv, color='b', label="PSv")
grfk[1].set(ylabel='Sv, PSv [m/s]', title="Velocity Response Spectrum")
grfk[1].legend()
grfk[1].grid()
grfk[2].plot(Tn, Sa, color='r', label="Sa")
grfk[2].plot(Tn, PSa, color='g', label="PSa")
grfk[2].set(xlabel='Tn [s]', ylabel='Sa, PSa [m/s2]', title="Acceleration Response Spectrum")
grfk[2].legend()
grfk[2].grid()
fig.tight_layout()
plt.show() 
  
  