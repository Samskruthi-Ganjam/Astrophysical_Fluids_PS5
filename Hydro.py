"""
Astrophysical Fluids - PHYS 643: Problem Set 5

Section 5: 1D Hydro Solver

This code uses the donor cell convection scheme to follow the motion of sound waves in a uniform density gas, 
starting from a small Gaussian perturbation in density by solving hydro equations

@author: Samskruthi Ganjam

November 12th, 2020
"""


import numpy as np
import matplotlib.pyplot as plt

# Hydro solver function

def solve_hydro(rho, u, cs, x, Nx, Nt, dx, dt):
    
    # Define the functions and their initial conditions
    f1 = np.zeros(Nx)
    f2 = np.zeros(Nx)
    f1 = rho
    f2[:-1] = u*rho[:-1] #to the length of u
    
    
    # Set up the plot
    plt.ion()
    fig, ax = plt.subplots(1,1, figsize=(8,8))
    ax.set_title('Hydro Solver - Density Plot')
    ax.set_xlim(0, Nx)
    #ax.set_ylim(0.5, 0.7)
    ax.set_xlabel('x')
    ax.set_ylabel('rho(x)')
    
    # Plot initial perturbation for reference in black
    ax.plot(x, f1, 'k-')
    
    # Set up plot object to update 
    p1, = ax.plot(x, f1, 'r')
    
    # Define the flux arrays
    J1 = np.zeros(Nx-1)
    J2 = np.zeros(Nx-1)
    

    # Loop over timesteps
    for i in range(Nt):
        
        # First step is to find the velocity at the cell interface 

        # Since we cannot probe a fractional address j+1/2 in an array, we take u[j+1/2]=u[j]
        # Therefore, the equation for u gives us,         

        u = 0.5 * ((f2[0:Nx-1]/f1[0:Nx-1]) + (f2[1:Nx]/f1[1:Nx]))
                
        
        # Next step is to compute the flux terms

        J1[u>0] = u[u>0] * f1[0:Nx-1][u>0] # For indices where u>0
        J1[u<0] = u[u<0] * f1[1:Nx][u<0] # For indices where u<0
        
        J2[u>0] = u[u>0] * f2[0:Nx-1][u>0] # For indices where u>0
        J2[u<0] = u[u<0] * f2[1:Nx][u<0] # For indices where u<0
        
        
        # From Continuity Equation, we can now update f1
        f1[1:Nx-1] = f1[1:Nx-1] - (dt/dx)*(J1[1:Nx-1] - J1[:Nx-2])

        # From Euler Equation with zero pressure gradient, update f2

        f2[1:Nx-1] = f2[1:Nx-1] - (dt/dx)*(J2[1:Nx-1] - J2[:Nx-2])

        # Adding the source term for all cells except boundary

        #f1 = f1
        f2[1:Nx-1] = f2[1:Nx-1] - (dt * cs**2 / dx)*(f1[2:] - f1[:Nx-2]) # Update f2 from Euler equation considering the pressure gradient

        # Boundary conditions

        f1[0] = f1[0] - (dt/dx) * J1[0]
        f1[-1] = f1[-1] + (dt/dx) * J1[-1]

        f2[0] = f2[0] - (dt/dx) * J2[0]
        f2[-1] = f2[-1] + (dt/dx) * J2[-1]

        # Update the plot
        p1.set_ydata(f1)   
        fig.canvas.draw()
        plt.pause(0.0001)


#----------------------------------------------------------------------------------------------------------------------------------------------------

#Set up the grid

Nx = 100 #Space
Nt = 20000 #Time

dx = 0.05 #Space step
dt = 0.00003 #Time step

# Define space grid
x = np.arange(Nx)

# Define the initial gaussian perturbation

x0 = (x[0]+x[Nx-1])/2.0 # gaussian centered at
sig = 4.0 # sigma
amp = 0.05 # some small amplitude

# Increase amplitude from ~0.5 or ~0.6 onwards to see that the wave appears to shape into a shock

rho = 0.6 + amp * np.exp(-(x - x0)**2 / (2 * sig**2))

cs = 343
u = np.zeros(Nx-1)

# Call the hydro solver function
solve_hydro(rho, u, cs, x, Nx, Nt, dx, dt)


