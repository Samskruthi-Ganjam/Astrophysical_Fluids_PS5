"""
Astrophysical Fluids - PHYS 643: Problem Set 5

Section 3: Advection Equation
Solution to the Advection Equation using FTCS and Lax-Friedrichs methods

@author: Samskruthi Ganjam

November 12th, 2020
"""


import numpy as np
import matplotlib.pyplot as plt


#Set up the grid

Nx = 50 #Space
Nt = 1000 #Time

dx = 1 #Space step
dt = 1 #Time step

# Define space grid
x = np.arange(Nx)*dx

# Set up the parameters for advection

# Given bulk velocity u=0.1
u = - 0.1

alpha = u*dt/(2*dx)

# Set initial conditions: given f(x,t=0) = x over all x

f0 = x/Nx # At first timestep
f1 = np.copy(f0) # Function f1 for FTCS
f2 = np.copy(f0) # Function f2 for Lax-Friedrichs

# Set up the plot

plt.ion()
fig, ax = plt.subplots(1,2, figsize=(14,10))
ax[0].set_title('FTCS')
ax[0].set_xlim(0, Nx)
ax[0].set_ylim(0, 2.0)
ax[0].set_xlabel('x')
ax[0].set_ylabel('f(x)')
ax[1].set_title('Lax-Friedrichs')
ax[1].set_xlim(0, Nx)
ax[1].set_ylim(0, 2.0)
ax[1].set_xlabel('x')
ax[1].set_ylabel('f(x)')

# Plot the intial conditions for reference in black
ax[0].plot(x, f1, 'k-')
ax[1].plot(x, f2, 'k-')

# Set up plot objects to be updated
p1, = ax[0].plot(x, f1, 'ro')
p2, = ax[1].plot(x, f2, 'ro')


# Fixed boundary conditions are satisfied by keeping first and last cells of f fixed (work with fj = f[1:Nx-1])
# fj+1 = f[2:] (Shift f by one cell to the right)
# fj-1 = f[:Nx-2] (Shift f by one cell to the left)
# Update f in each timestep
 
# Check for stability for Lax-Friedrichs : dt <= dx/u <= 1/0.1 <=10, here, dt = 1 <=10. Hence, stable for Lax-Friedrichs 

for tstep in range(Nt):
    
    # FTCS - update f1
    f1[1:Nx-1] = f1[1:Nx-1] - alpha * (f1[2:] - f1[:Nx-2])
    
    # Lax-Friedrichs - update f2
    f2[1:Nx-1] = 0.5 * (f2[2:] + f2[:Nx-2]) - alpha * (f2[2:] - f2[:Nx-2])
    
    # Update the plot

    p1.set_ydata(f1)
    p2.set_ydata(f2)
    
    fig.canvas.draw()
    plt.pause(0.0001)
        
    

