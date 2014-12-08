# file: becke_grid.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it 

from __future__ import division
import numpy as np
from copy import deepcopy

from utilities import numerical_integration as nint

class BeckeGrid(object):
    """This class creates the Becke grid.
    see: Becke, A. D. A multicenter numerical integration scheme for polyatomic 
    molecules. J. Chem. Phys. 88, 2547 (1988).
    """
    def __init__(self, n_radial=15, n_angular=110, position=[0.,0.,0.]):
        super(BeckeGrid, self).__init__()

        self.data={}
        self.data['n_radial'] = n_radial
        self.data['n_angular'] = n_angular
        self.data['size'] = n_radial * n_angular

        self.data['x'] = np.zeros(self.data['size'], dtype=np.float64)
        self.data['y'] = np.zeros(self.data['size'], dtype=np.float64)
        self.data['z'] = np.zeros(self.data['size'], dtype=np.float64)
        self.data['w'] = np.zeros(self.data['size'], dtype=np.float64)

        # Radial distribution
        q_r, w_r = nint.chebgauss_t_radial_q_w(n_radial)

        #angular distribution
        x, y, z, w_a = nint.lebedev_q_w(n_angular)
    
        count = 0
        for i in xrange(n_radial):
            solid_angle = 8.0 * np.arccos(0.0) * q_r[i] * q_r[i]
            for j in xrange(n_angular):
                self.data['x'][count] = x[j] * q_r[i] 
                self.data['y'][count] = y[j] * q_r[i] 
                self.data['z'][count] = z[j] * q_r[i] 
                self.data['w'][count] = w_a[j] * w_r[i]  * solid_angle
                count += 1

        for i in xrange(self.get('size')):
            self.data['x'][i]+= position[0]
            self.data['y'][i]+= position[1]
            self.data['z'][i]+= position[2]


    def weights(self, particle_stack, particleID):
        """Computes the Becke weights as described in:
        Becke, A. D. A multicenter numerical integration scheme for polyatomic 
        molecules. J. Chem. Phys. 88, 2547 (1988).
        particle_stack: stack of particles. i.e. stack of e-
        """
        def cutoff_profile(mu):
            """Iterated cutoff profile. eq. 21
            """
            return 0.5 * mu * ( 3.0 - (mu * mu))

        distance_to_particle = np.zeros(particle_stack.size(), dtype=np.float64)
         
        for point in xrange(self.get('size')):
            cell_function = np.ones(particle_stack.size(), dtype=np.float64)
            #Distance from grid point i to particle (distance_to_particle) (eq. 11)
            for i in xrange(particle_stack.size()):
                distance_to_particle[i] = np.sqrt( 
                (self.get('x')[point] - particle_stack.get(i).get('position')[0])**2 +
                (self.get('y')[point] - particle_stack.get(i).get('position')[1])**2 +
                (self.get('z')[point] - particle_stack.get(i).get('position')[2])**2, dtype=np.float64 )

            for i in xrange(1, particle_stack.size()):
                for j in xrange(i):         
                    #Internuclear distance (R_ij eq. 11)
                    R_ij = np.sqrt(
                    (particle_stack.get(j).get('position')[0] - particle_stack.get(i).get('position')[0])**2 +
                    (particle_stack.get(j).get('position')[1] - particle_stack.get(i).get('position')[1])**2 +
                    (particle_stack.get(j).get('position')[2] - particle_stack.get(i).get('position')[2])**2 )
             
                    # \mu_ij eq. 11
                    mu_ij = (distance_to_particle[i] - distance_to_particle[j]) / R_ij                         
                    
                    #missing atomic size adjustment.
             
                    #f_k(\mu_ij) K = 3 eq. 20
                    f_ij = cutoff_profile(cutoff_profile(cutoff_profile(mu_ij)))
             
                    #Cutoff profile s_k(\mu_ij) eq. 21
                    s_ij = 0.5 * (1.0 - f_ij)

                    #Cell function P_i(r) eq. 13
                    cell_function[i] *=  s_ij            
       
            #cell_function sum denominator eq. 22
            cell_function_sum = 0.0
            for i in xrange(particle_stack.size()):
                cell_function_sum += cell_function[i]
       
            #Weight calculation for each point eq. 22
            self.get('w')[point] *= cell_function[particleID]/cell_function_sum

    def get(self, key):
        """Returns the value stored in key
        """
        assert type(key) == type('str')

        try:
            return self.data[key]
        except KeyError:
            raise

