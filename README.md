# Riemann_Solver

A basic implementation of 1 dimensional riemann solver using SPH technique. 

Read about Riemann problems [here](https://en.wikipedia.org/wiki/Riemann_problem). 

As of now, it can simulate evolution of a field `F` governed by the following eqns
- dF/dx + dF/dt = 0     Given initial F 
- dF/dt + F\*dF/dx = 0    Given initial F

To add:

- XSPH Velocity correction for particles
- Variable particle velocities support : EGN as an evolving vector field 
- More equations for simulating cases such as [sod shocktube](https://en.wikipedia.org/wiki/Sod_shock_tube)
- Generalise periodicity from [-1,1] to any arbitary domain

