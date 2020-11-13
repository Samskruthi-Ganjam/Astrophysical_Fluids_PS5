# Astrophysical_Fluids_PS5

 Numerical Assignment for Astrophysical Fluids (PHYS 643)

Name: Samskruthi Ganjam

Version: Python 3

List of python code files:

1. Advection.py : This code corresponds to section 3 of Problem Set 5. This is the solution to the Advection equation using FTCS and Lax-Friedrich methods.
2. AdvDiff.py : This code corresponds to section 4 of Problem Set 5. This is the solution to the Advection-Diffusion equation. 
3. Hydro.py : This code corresponds to section 5 of Problem Set 5. This is the code for the 1D Hydro Solver.

## 1. Advection Equation (refer to code Advection.py)

We see that the FCTS method fails to be stable. However, for appropriate time step and grid spacing (see comments in the code for stability criterion),
we see that the Lax-Friedrichs method gives a numerically stable solution.

## 2. Advection-Diffusion Equation (refer to code AdvDiff.py)

As we increase the diffusion coefficient (but within the stability range) the diffusion effect is more pronounced. 
The solution appears to be less steeper compared as one would expect from a pure advection scenario.


## 3. 1D Hydro solver (refer to code Hydro.py)

For a Gaussian perturbation in density with a small amplitude we see that the initial peak splits into two peaks of reduced amplitude, which travel outwards, 
reflect off the walls of the grid and move back towards the centre. As we increase the amplitude of the Gaussian perturbation, we see that 
the two travelling peaks steepen at one end which resembles a shock like feature. The width of this shock seems to be set by the strength of the perturbation 
i.e, the amplitude of the initial perturbation. 


