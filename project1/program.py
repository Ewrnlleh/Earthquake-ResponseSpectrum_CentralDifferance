# -*- coding: utf-8 -*-
"""
Created on Sun May 17 12:41:26 2020

@author: Tarık
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy import array as matrix
from numpy import arange as seri
from numpy import pi as π
from operator import add
import xlrd


#%%
# INPUT MOTION (GROUND MOTION)
'==== Read PEER data ===================='

class Program:
    
    global t
    
    def __init__(self):
         self.Sd = []
         self.Sv = []
         self.Sa = []
         self.acc = []
         # Spectrum Periods
         Tmax = 7    # Maximum period on the spectrum
         ΔT = 0.1    # Period increment/interval of the spectrum
         self.Tn = matrix([seri(ΔT,Tmax+ΔT,ΔT)]).T   # Period interval
         
    def ReadData(self,dosyaIsim):
        file = open(dosyaIsim, 'r')  # Name of the file
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
        self.acc = SF * matrix([data]).T               # Acceleration data
        self.t = matrix([seri(0,len(self.acc)*Δt,Δt)]).T   # Time interval
    
    def Ground(self):
        # Plot Ground Motion
        fig, xy = plt.subplots()
        xy.plot(self.t,self.acc, color="g")
        xy.set(xlabel='t [s]', ylabel='ag [g]', title='Time - Ground Acceleration Curve')
        xy.grid()
        plt.show()
        
    
    #%%
    # Function of Central Differance Method
    def SDOFCentDiff(self,acc,Tn,ξ,m,u0,v0,a0):
      Δt = 0.01  
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
    
    " 1) To obtain Acceleration, Velocity and Displacement response of an SDOF system under an EQ motion"
    def Graph1(self):
        " Input Parameters "
        ξ = 0.05    # Damping Ratio  
        m = 1.0     # Assumption for Mass 
        u0 = 0.0    # Intial displacement
        v0 = 0.0    # Intial velocity
        a0 = 0.0    # Intial acceleration
        Tni = 2.5    # for this period
        u, v, a = self.SDOFCentDiff(self.acc,Tni,ξ,m,u0,v0,a0)
        u = matrix([u]).T ; v = matrix([v]).T ; a = matrix([a]).T
        
        # Plot response 
        fig, grfk = plt.subplots(nrows=3, ncols=1)
        grfk[0].plot(self.t, u, color='k', label="Umax={} m".format(np.around(np.max(u),decimals=4)))
        grfk[0].set(ylabel='Disp. [m]', title='Response of an SDOF system under an EQ motion (Tn:{} s)'.format(Tni))
        grfk[0].legend()
        grfk[0].grid()
        grfk[1].plot(self.t[:-1], v, color='r', label="Vmax={} m/s".format(np.around(np.max(v),decimals=3)))
        grfk[1].set(ylabel='Vel. [m/s]')
        grfk[1].legend()
        grfk[1].grid()
        grfk[2].plot(self.t[:-1], a, color='b', label="Amax={} m/s2".format(np.around(np.max(a),decimals=2)))
        grfk[2].set(xlabel='Time [s]', ylabel='Acc. [m/s2]')
        grfk[2].legend()
        grfk[2].grid()
        fig.tight_layout()
        plt.show()
    
    
    " 2) To create response spectra"
    def Graph2(self):
        " Input Parameters "
        ξ = 0.05    # Damping Ratio  
        m = 1.0     # Assumption for Mass 
        u0 = 0.0    # Intial displacement
        v0 = 0.0    # Intial velocity
        a0 = 0.0    # Intial acceleration
        for Ti in self.Tn:
            u, v, a = self.SDOFCentDiff(self.acc,Ti,ξ,m,u0,v0,a0)     # Recursive Relationship
            self.Sd.append(np.max(np.abs(u)))
            self.Sv.append(np.max(np.abs(v)))
            self.Sa.append(np.max(np.abs(a)))       
        self.Sd = matrix(self.Sd)       # Relative Displacement Response Spectrum
        self.Sv = matrix(self.Sv)       # Relative Velocity Response Spectrum
        self.Sa = matrix(self.Sa)       # Relative Acceleration Response Spectrum
    
        # Plot Graphs
        fig, grfk = plt.subplots(nrows=3, ncols=1)
        grfk[0].plot(self.Tn, self.Sd, color='g')
        grfk[0].set(ylabel='Sd [m]', title="Displacement Response Spectrum")
        grfk[0].grid()
        grfk[1].plot(self.Tn, self.Sv, color='r')
        grfk[1].set(ylabel='Sv [m/s]', title="Velocity Response Spectrum")
        grfk[1].grid()
        grfk[2].plot(self.Tn, self.Sa, color='b')
        grfk[2].set(xlabel='Tn [s]', ylabel='Sa [m/s2]', title="Acceleration Response Spectrum")
        grfk[2].grid()
        fig.tight_layout()
        plt.show()
    
    " 3) To create Pseudo response spectra"
    def Graph3(self):        
        wn = 2*π/self.Tn
        PSv = self.Sd * wn           # Pseudo Velocity Response Spectrum
        PSa = self.Sd * wn**2        # Pseudo Acceleration Response Spectrum 
        # Plot Graphs
        
        fig, grfk = plt.subplots(nrows=2, ncols=1)
        grfk[0].plot(self.Tn, PSv, color='g',label="PSv")
        grfk[0].set(ylabel='PSv [m/s]', title="pseudo velocity response spectrum")
        grfk[0].grid()
        grfk[1].plot(self.Tn, PSa, color='r',label="PSa")
        grfk[1].set(ylabel='PSa [m/s2]', title="pseudo acceleration response spectrum")
        grfk[1].grid()

        fig.tight_layout()
        plt.show()
        
        fig, grfk = plt.subplots(nrows=2, ncols=1)
        
        grfk[0].plot(self.Tn, self.Sv, color='r', label="Sv")
        grfk[0].plot(self.Tn, PSv, color='b', label="PSv")
        grfk[0].set(ylabel='Sv, PSv [m/s]', title="Velocity Response Spectrum")
        grfk[0].legend()
        grfk[0].grid()
        grfk[1].plot(self.Tn, self.Sa, color='r', label="Sa")
        grfk[1].plot(self.Tn, PSa, color='g', label="PSa")
        grfk[1].set(xlabel='Tn [s]', ylabel='Sa, PSa [m/s2]', title="Acceleration Response Spectrum")
        grfk[1].legend()
        grfk[1].grid()
        fig.tight_layout()
        plt.show()
             
dosya1 = Program()
dosya2 = Program()
dosya3 = Program()
dosya4 = Program()

dosya1.ReadData('RSN186_IMPVALL.H_H-NIL090.AT2')
#dosya1.Ground()
#dosya1.Graph1()
dosya1.Graph2()
dosya1.Graph3()

dosya2.ReadData('RSN186_IMPVALL.H_H-NIL360.AT2')
#dosya2.Ground()
#dosya2.Graph1()
dosya2.Graph2()
dosya2.Graph3()

dosya3.ReadData('RSN1165_KOCAELI_IZT090.AT2')
#dosya3.Ground()
#dosya3.Graph1()
dosya3.Graph2()
dosya3.Graph3()

dosya4.ReadData('RSN1165_KOCAELI_IZT180.AT2')
#dosya4.Ground()
#dosya4.Graph1()
dosya4.Graph2()
dosya4.Graph3()

Sar1 = np.sqrt(list(map(add, (i ** 2 for i in dosya1.Sa), (i ** 2 for i in dosya2.Sa)))) 

fig, gr1 = plt.subplots()
gr1.plot(dosya1.Tn,Sar1)
gr1.set(xlabel='t [s]', ylabel='Sar', title='Sar - IMPVALL.H_H-NIL')
gr1.grid()
plt.show() 

Sar2 = np.sqrt(list(map(add, (i ** 2 for i in dosya3.Sa), (i ** 2 for i in dosya4.Sa))))

fig, gr2 = plt.subplots()
gr2.plot(dosya3.Tn,Sar2)
gr2.set(xlabel='t [s]', ylabel='Sar', title='Sar - KOCAELI_IZT')
gr2.grid()
plt.show() 


Sartoplam = list(map(add, Sar1, Sar2))
Saravr = [x / 2 for x in Sartoplam]

fig, gravr = plt.subplots()
gravr.plot(dosya1.Tn,Saravr)
gravr.set(xlabel='t [s]', ylabel='Saravr', title='Saravr - IMPVALL.H_H-NIL,KOCAELI_IZT')
gravr.grid()
plt.show() 


wb = xlrd.open_workbook("YatayElastikTasarimSpektrumu.xlsx")
Sae = []; SaeT = []
sheet = wb.sheet_by_index(0)

for row in range(1,sheet.nrows):
    SaeT.append(sheet.cell_value(row,0))
    Sae.append(sheet.cell_value(row,1))

fig, grsae = plt.subplots()
grsae.plot(SaeT,Sae)
grsae.set(xlabel='t [s]', ylabel='Sae', title='T-Sae')
grsae.grid()
plt.show()

         
fig, grfk = plt.subplots()
grfk.plot(dosya1.Tn,Saravr, color="r", label="Saravr")
grfk.plot(SaeT,Sae, color="k", label="Sae")
grfk.set(xlabel='Time (s)', title='Plot “T – Saravr” and “T – Sae”')
grfk.grid()
grfk.legend(loc="best")
plt.show()
 

  