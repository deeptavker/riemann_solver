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

### Example

```Python

new_sim = RiemannSolver1D(particles = particles, phi_initial = phi, kernel = vec_kernel, EGN = 1)

new_sim.configure_solver(0.1, 40)

for i in range(int(new_sim.tf/new_sim.dt)):
    new_sim.update_rho()
    new_sim.update_field_euler()
    new_sim.update_position_euler(is_periodic=True, period=[-1,1])
    
plt.plot(new_sim.particles.x[3:new_sim.nopart-3], new_sim.particles.phi[3:new_sim.nopart-3])
```
