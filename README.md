# Riemann_Solver

A basic implementation of 1 dimensional riemann solver using SPH technique. 

Read about Riemann problems [here](https://en.wikipedia.org/wiki/Riemann_problem). 

As of now, it can simulate evolution of a field `F` governed by the following eqns
- dF/dx + dF/dt = 0     Given initial F , in this case EGN = 1
- dF/dt + F\*dF/dx = 0    Given initial F , in this case EGN = phi (field) #hardcoded 


### Example

Solving the following problem :
- `du/dt + du/dx = 0`
- `u(x, 0) = - sin(pi * x)`
- Periodic in `[-1, 1]`
- Find `u(x, 40)`

```Python

new_sim = RiemannSolver1D(particles = particles, phi_initial = phi, kernel = vec_kernel, EGN = 1)

new_sim.configure_solver(0.1, 40)

for i in range(int(new_sim.tf/new_sim.dt)):
    new_sim.update_rho()
    new_sim.update_field_euler()
    new_sim.update_position_euler(is_periodic=True, period=[-1,1])
    
plt.plot(new_sim.particles.x[3:new_sim.nopart-3], new_sim.particles.phi[3:new_sim.nopart-3])
```

### Sine wave captured while travelling. 
![sine](https://github.com/deeptavker/Riemann_Solver/blob/master/sine_9.png)


- `particles` are created using `ParticleArray` from `pysph.base.utils` , extra padding of
extra particles on either side of boundary is to be done while create the instance of solver
- `phi` here means the field, `F` in the equation. 
- `configure_solver()` is required for any simulation as it sets `dt` (time step) and `tf` (t_final)
- `update_rho()` uses the kernel approximation to update the density for each particle
- `update_field_euler()` uses the sph approximation for updating the field array

### Example

Solving the following problem :
- `du/dt + u * du/dx = 0`
- `u(x, 0) = 1` for `|x| < 1/3` and `0` otherwise
- Periodic in `[-1, 1]`
- Find `u(x, 0.6)`

```Python
new_sim = RiemannSolver1D(particles = particles, phi_initial = phi, kernel = vec_kernel, EGN = 1)

new_sim.configure_solver(0.01, 0.6)
for i in range(int(new_sim.tf/new_sim.dt)):
    new_sim.update_rho()
    new_sim.update_field_euler(field=True)
    new_sim.update_position_euler(is_periodic=True, period=[-1,1], field=True)
plt.scatter(new_sim.particles.x[3:new_sim.nopart-3], new_sim.particles.phi[3:new_sim.nopart-3])
```
