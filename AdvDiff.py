"""
Astrophysical Fluids - PHYS 643: Problem Set 5

Section 4: Advection-Diffusion Equation

This code uses Lax-Friedrich method for advection and the implicit method for diffusion to solve the Advection-Diffusion Equation using operator splitting.

@author: Samskruthi Ganjam

November 12th, 2020
"""

import numpy as np
import matplotlib.pyplot as plt
  
#Set up the grid

Nx = 50 #Space
Nt = 5000 #Time

dx = 1 #Space step
dt = 1 #Time step

# Define space grid

x = np.arange(Nx)*dx

# Set up the advection parameters

u = - 0.1 # Given bulk velocity u=0.1
alpha = u*dt/(2*dx)

# Set Diffusion parameters

# Choose diffusion coefficients
# Check for stability: dt <= 2dx^2/D. We have xhosen the same timestep dt=1
# For D1 = 0.1: dt=1 <= 2(1)/0.1 and hence stable dt <= 20
# For D2 = 1.0: dt=1 <= 2(1)/1.0 and hence stable dt <= 2
# For the same time step, the solution becomes unstable when D>2 (change D2 to, say 10, to observe this)

D1 = 0.1
D2 = 1.0

beta1 = (D1*dt)/(dx**2)
beta2 = (D2*dt)/(dx**2)

# Define the tri-diagonal matrix A with (1+2*beta) along the diagonals and (-beta) on the upper and lower diagonals

A1 = np.eye(Nx-2)*(1.0 + 2.0*beta1) + np.eye(Nx-2, k=1)*(-beta1) + np.eye(Nx-2, k=-1)*(-beta1)
A2 = np.eye(Nx-2)*(1.0 + 2.0*beta2) + np.eye(Nx-2, k=1)*(-beta2) + np.eye(Nx-2, k=-1)*(-beta2)

# Set boundary conditions for matrix A

A1[0][0] = 1
A1[0][1] = 0
A1[-1][-1] = 1 + beta1

A2[0][0] = 1
A2[0][1] = 0
A2[-1][-1] = 1 + beta2


# Set initial conditions: T(x,t=0) = x over all x  

f0 = x/Nx # First time step
T1 = np.copy(f0)
T2 = np.copy(f0)

# Set up the plot 

plt.ion()
fig, ax = plt.subplots(1,2, figsize=(12,8))
ax[0].set_title('Advection-Diffusion: D = 0.1')
ax[0].set_xlim(0, Nx)
ax[0].set_ylim(0, 2.0)
ax[0].set_xlabel('x')
ax[0].set_ylabel('f(x)')
ax[1].set_title('Advection-Diffusion: D = 1.0')
ax[1].set_xlim(0, Nx)
ax[1].set_ylim(0, 2.0)
ax[1].set_xlabel('x')
ax[1].set_ylabel('f(x)')

#Plot the intial condition for reference in black

ax[0].plot(x, T1, 'k-')
ax[1].plot(x, T2, 'k-')

# Set up plot objects to update 

p1, = ax[0].plot(x, T1, 'ro') # D = 0.1
p2, = ax[1].plot(x, T2, 'ro') # D = 1.0


# Update T in each timestep  

# Fixed boundary conditions are satisfied by keeping first and last cells of f fixed (work with Tj = T[1:Nx-1]) 
# Tj+1 = T[2:] (Shift T by one cell to the right)
# Tj-1 = T[:Nx-2] (Shift T by one cell to the left)

for tstep in range(Nt):
    
    
    # First update Diffusion term - implicit method  
    #T[1:Nx-1] += beta*(T[2:] + T[:Nx-2] - 2 * T[1:Nx-1]) # Explicit     
    T1[1:Nx-1] = np.linalg.solve(A1, T1[1:Nx-1]) #implicit, D = 0.1
    T2[1:Nx-1] = np.linalg.solve(A2, T2[1:Nx-1]) #implicit, D = 1.0

    # Then update Advection term - Lax-Friedrichs method
    T1[1:Nx-1] = 0.5 * (T1[2:] + T1[:Nx-2]) - alpha * (T1[2:] - T1[:Nx-2])# D = 0.1
    T2[1:Nx-1] = 0.5 * (T2[2:] + T2[:Nx-2]) - alpha * (T2[2:] - T2[:Nx-2])# D = 1.0
               
    # Update the plot
    p1.set_ydata(T1) # D = 0.1
    p2.set_ydata(T2) # D = 1.0
    
    fig.canvas.draw()
    plt.pause(0.0001)
        
    

