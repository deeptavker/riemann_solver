# Goal is to write a riemann solver

from pysph.base.utils import ParticleArray
import numpy as np

@np.vectorize
def vec_kernel(x, h, order):
    RPI = np.pi**0.5
    q = x/float(h)
    sigma = 1.0/(RPI*h)
    
    if order == 0:
        
        q = abs(q)
        
        if q <= 3.0:
            
            return sigma*np.exp(-q*q)
        
        else:
            
            return 0.0
        
    elif order == 1:
        
        if abs(q) <= 3:
            
            return sigma*(-2*q)*np.exp(-q*q)/h
        
        else:
            
            return 0.0
        
    else:
            pass


class RiemannSolver1D(object):
    
    def __init__(self, particles, phi_initial, kernel, EGN):
        '''
        Particles are assumed to have the following properties
        
        rho : Density
        phi : field vector
        x : Position vector
        tag : Ghost or real
        h : smoothing length
        name : a name
        
        EGN can be a scalar or a vector field. 
        '''
        self.particles = particles
        self.phi_initial = phi_initial
        self.kernel = kernel
        self.EGN = EGN
        self.nopart = len(self.particles.x)
        
    def configure_solver(self, dt, tf):
        self.dt = dt
        self.tf = tf
        
    def update_rho(self):
        '''
        This method updates the density
        The is for one timestep
        '''
        pos_arr = self.particles.x
        h = self.particles.h[0]
        for i in range(self.nopart):
            
            self.particles.rho[i] = np.sum(self.kernel(pos_arr[i] - pos_arr, h, 0))
            
    def update_field_euler(self, **kwargs):
        '''
        This method updates phi
        This is for one timestep
        '''        
        pos_arr = self.particles.x
        h = self.particles.h[0]
        phi = self.particles.phi
        tags = self.particles.tag
        rho = self.particles.rho
        phi_new = np.arange(self.nopart, dtype='float32')
        dphi = []
        
        if 'field' not in kwargs.keys():
        
            for i in range(self.nopart):
                
                temp = np.sum(self.kernel(pos_arr[i] - pos_arr, h, 1) * (phi - phi[i])) / rho[i]
                
                dphi.append( -(self.dt ) * temp)
                
            for i in range(self.nopart):

                phi_new[i] = phi[i] + dphi[i]

        else:

            self.EGN = self.particles.phi

            for i in range(self.nopart):
                
                temp = np.sum(self.kernel(pos_arr[i] - pos_arr, h, 1) * phi / rho[i])
                
                dphi.append( -(self.dt ) * temp * self.EGN[i])
                
            for i in range(self.nopart):

                phi_new[i] = phi[i] + dphi[i]
                
        phi = phi_new
        
    def compute_phi(self):
        sph_phi = []
        pos_arr = self.particles.x
        h = self.particles.h[0]
        phi = self.particles.phi
        tags = self.particles.tag
        rho = self.particles.rho
        
        for i in range(self.nopart):
            
            sph_phi.append(np.sum(self.kernel(pos_arr[i] - pos_arr, h, 0) * phi / rho[i]))
            
        self.particles.phi = sph_phi
            
    def update_position_euler(self, **kwargs):
        '''
        This method updates position
        This is for one timestep. 
        
        Kwargs : is_periodic, period, field
        
        Remember, with position update, there also needs to be a reordering of particles. 
        Use argsort for the same. 
        '''
        
        pos_arr = self.particles.x
        tags = self.particles.tag
        
        if 'is_periodic' in kwargs.keys():

            period = kwargs['period']

            dom_a = period[0]
            dom_b = period[1]
            diff = dom_b - dom_a
            
            if 'field' in kwargs.keys():
                
                #hardcoding velocity field for special case 
                self.EGN = self.particles.phi
                
                for i in range(self.nopart):
                    if tags[i] == 1:
                        new_pos = pos_arr[i] + self.EGN[i] * self.dt
                        if new_pos >= dom_b:
                            pos_arr[i] = new_pos - diff
                        elif new_pos <= dom_a:
                            pos_arr[i] = new_pos + diff
                        else:
                            pos_arr[i] = new_pos
                
            else:
                
                for i in range(self.nopart):
                    if tags[i] == 1:
                        new_pos = pos_arr[i] + self.EGN * self.dt
                        if new_pos >= dom_b:
                            pos_arr[i] = new_pos - diff
                        elif new_pos <= dom_a:
                            pos_arr[i] = new_pos + diff
                        else:
                            pos_arr[i] = new_pos

        else:


            if 'field' in kwargs.keys():
                
                #hardcoding velocity field for special case 
                self.EGN = self.particles.phi
                
                for i in range(self.nopart):
                    if tags[i] == 1:
                        new_pos = pos_arr[i] + self.EGN[i] * self.dt
                        pos_arr[i] = new_pos
                
            else:
                
                for i in range(self.nopart):
                    if tags[i] == 1:
                        new_pos = pos_arr[i] + self.EGN * self.dt
                        pos_arr[i] = new_pos



        self.sort_particles()
    
    def sort_particles(self):
        
        sortmask = np.argsort(self.particles.x)
        self.particles.x = self.particles.x[sortmask]
        self.particles.rho = self.particles.rho[sortmask]
        self.particles.phi = self.particles.phi[sortmask]
        self.particles.tag = self.particles.tag[sortmask]
        
        
        
            
    def xsph_correction(self):
        pass