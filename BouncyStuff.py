# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 14:28:47 2020

@author: Swoyam Panda
"""

import numpy as np
import pylab as pl
import math as mt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class Ball:
    """
    Contains all the attributes of a ball
    """
    def __init__(self, mass=0, radius=0, position=[0, 0], velocity=[0,0]):
        #everything is defined as float for consistency
        self.mass = float(mass)
        self.radius = float(radius)
        self._position = np.array(position, dtype=np.float64)
        self._velocity = np.array(velocity, dtype=np.float64)
        self.self_col_count = 0 
        
    
    def pos(self):
        return self._position
    
    def vel(self):
        return self._velocity
    
    def time_to_collision(self, other):
        #Enter the equation for delta t
        r = self.pos() - other.pos()
        v = self.vel() - other.vel()
        rv = np.dot(r, v)
        #_sqr variables contain vector magnitude squared
        v_sqr = np.dot(v, v)
        r_sqr = np.dot(r, r)
        """
        The radius of the container will be entered 
        as a negative number. This accounts for the 
        R1-R2 case. 
        """
        R = self.radius + other.radius
        
        #Check for valid values of the determinant
        if v_sqr == 0:
            determinant = float('nan')
        else:    
            determinant = (4*(rv**2)-4*v_sqr*(r_sqr-R**2))
            if determinant<0:
                determinant = float('nan')
            else:
                pos_dt = (-2*rv+(determinant)**0.5)/(2*v_sqr)
                neg_dt = (-2*rv-(determinant)**0.5)/(2*v_sqr)

        
        """
        After some diagramming, I realized that it is better to base the conditions 
        on the value of r dot v. The code ends up being more succinct. 
        """
        if mt.isnan(determinant):
            dt = float('inf')
        elif self.radius < 0 or other.radius<0: #accounts for ball-container collision
            dt = pos_dt  #dt is just going to be pos_dt
        elif rv<0:
            dt = min(pos_dt, neg_dt)
        else:
            dt = float('inf')
        
        return dt
        
    
    def move(self, dt):
        """
        Moves a ball the distance it would move in time dt 
        """
        self._position += self.vel()*dt
        
        

        
    def collide(self, other):
        """
        Collides two objects and changes their velocity accordingly. 
        """
        r = self._position - other.pos()
        mag_r = np.linalg.norm(r)
        normal_r = r/mag_r
        
        v1_par = np.dot(self.vel(), normal_r)*normal_r
        v2_par = np.dot(other.vel(), normal_r)*normal_r
        v1_perp = self.vel()- v1_par
        v2_perp = other.vel() - v2_par 
        
        m1 = self.mass
        m2 = other.mass
        
        #Post collision velocity
        vel1f = (m1-m2)*v1_par/(m2+m1) + 2*m2*v2_par/(m2+m1) # -v1_par + 2*v2_par # 
        vel2f =  2*m1*v1_par/(m2+m1) + (m2-m1)*v2_par/(m2+m1) #v2_par # 
        
        #Updated velocity
        v1_final = v1_perp + vel1f
        v2_final = v2_perp + vel2f
        
        #Change the value
        self._velocity = v1_final
        other._velocity = v2_final
        
    def patch_func(self, filling=False, fillc='b', edgec='b'):
        """
        Defines a patch. The patch represents a Ball object graphically. 
        """        
        #fill accepts a boolean value. 
        patch = pl.Circle(self.pos(), abs(self.radius), ec=edgec, fc=fillc, fill=filling, ls='solid')
        #ec is edge colour, fc is fill colour
        return patch

        
class Simulation(Ball):
    
    def __init__(self, mass, radius, position, velocity, ball_num):
        Ball.__init__(self, mass, radius, position, velocity)
        self.container = Ball(float(400000), -60, [0,0], [0,0]) #infinite mass container
        self.tot_time = 0 #total time elapsed
        self.del_p = 0 #just a scalar 
        self.collision_count = 0 #increments by 1 with each collision
        
        """
        These test balls were defined for debugging collisions
        """
        # self.ball_test1 = Ball(1, 1, [-3,0], [0,1])
        # self.ball_test2 = Ball(1, 1, [1,0], [0,1])
        
        # self._object_lst = [self.container, self.ball_test1, self.ball_test2] - testing with 2 balls
        # self._object_lst = [self.container, self.ball_test1] - use for testing with 1 ball
        """
        Comment them out after debugging
        """
        
        
        self._object_lst = [self.container]
        
        if ball_num == 1:
            self._object_lst.append(Ball(mass, radius, position, velocity))
        elif ball_num == 16:
            l = np.array([[2*i for i in range(-5, -1)], [2*i for i in range(-5, -1)]])
            l = l*4 #4 is an arbitrary constant to affect the spacing between the particles. 
            
            vel_mag = np.linalg.norm(velocity) #desired velocity magnitude
            for i in range(4):
                for j in range(4):
                    grid_pos = np.array([l[0][i], l[1][j]])
                    initial_vel = np.random.normal(scale = vel_mag, size=2)
                    self._object_lst.append(Ball(mass, radius, grid_pos, initial_vel))
        elif ball_num == 100:
            #Arrange 100 balls 
            l = np.array([[2*i for i in range(-5, 5)], [2*i for i in range(-5, 5)]])
            l = l*4 #4 is an arbitrary constant to affect the spacing between the particles. 
            
            vel_mag = np.linalg.norm(velocity) #desired velocity magnitude
            for i in range(10):
                for j in range(10):
                    grid_pos = np.array([l[0][i], l[1][j]])
                    initial_vel = np.random.normal(scale = vel_mag, size=2)
                    self._object_lst.append(Ball(mass, radius, grid_pos, initial_vel))
        else:
            raise Exception("This object can only be defined with 1, 16, or 100 particles")
            
    def col_count(self):
        """
        Returns the number of collision that the system has undergone. 
        """
        return self.collision_count
    
    def self_count(self):
        """
        Returns an array of the number of collisions each ball underwent.
        If the collisions work properly, the sum of all the count_arr elements should be double 
        the value returned by the col_count() method. 
        
        Another useful aspect of this function is to check how many of the collisions happen
        with the container vs. how many ball-ball collisions happen out of the total. 
        """
        count_arr = np.array([i.self_col_count for i in self._object_lst])
        return count_arr
    
    def particle_conserve(self):
        """
        Returns the number of particles still inside the container
        """
        count = 0
        for i in self._object_lst[1:]:
            #calculate the distance of a point on the ball furthest from the center of container
            #center of container is always at [0, 0]
            point_dist = np.sqrt(np.dot(i.pos(), i.pos())) + i.radius
            
            if point_dist <= abs(self.container.radius):
                count += 1
        
        if count == len(self._object_lst[1:]):
            print("Particle number is conserved")
        else:
            print("Particle number is not conserved")
        
        return count
    
    def short_time(self):
        """
        This method calculates the minimum time that 
        exists in a system of balls for a collision to 
        take place.
        """
        length = len(self._object_lst)
        # print(length)
        min_dt = float('inf')
        for i in range(length):
            if i<(length-1): #to not test for collision with itself
                for j in range(i,length-1):
                    temp1 = self._object_lst[i]
                    temp2 = self._object_lst[j+1]
                    dt = temp1.time_to_collision(temp2)
                    
                    if dt<min_dt:
                        min_dt = dt
                        index1 = i
                        index2 = j+1
                        col_data = [index1, index2, min_dt] #indices tell me which objects to collide
        return col_data
    
    def next_collision(self):         
        """
        Performs a collision. 
        Also calculates total time
        """
        collision_data = self.short_time()
        del_t = collision_data[2]
        index1 = collision_data[0]
        index2 = collision_data[1]
        
        v_1i = self._object_lst[index1].vel()
        v_2i = self._object_lst[index2].vel()
                
        if not mt.isnan(del_t): #this if statement is clearly no longer needed
            self.tot_time += del_t
            for i in self._object_lst: #a loop to move all the balls
                i.move(del_t)
                
        self._object_lst[index1].collide(self._object_lst[index2])
        self.collision_count += 1 #collision counter increments after collision
        
        v_1f = self._object_lst[index1].vel()
        v_2f = self._object_lst[index2].vel()
        
        vel_diff1 = v_1f - v_1i
        vel_diff2 = v_2f - v_2i
        
        #increment self_col_count
        
        self._object_lst[index1].self_col_count += 1
        self._object_lst[index2].self_col_count += 1
        
        
        #calculate delta momentum with container collisions
        #necessary for measuring calculating pressure
        if index2 == 0:
            vel_mag = np.sqrt(np.dot(vel_diff1, vel_diff1)) #magnitude of velocity change
            self.del_p += vel_mag*self._object_lst[index1].mass
        elif index1 == 0:
            vel_mag = np.sqrt(np.dot(vel_diff2, vel_diff2))
            self.del_p += vel_mag*self._object_lst[index2].mass
            
            
    def momentum(self):
        """
        Returns an array of the momentum of all the balls and the container
        """
        momentum = []
        xmom = []
        ymom = []
        for i in self._object_lst: #[1:] to avoid calculating container momentum
        
            xmom_temp = i.vel()[0]*i.mass 
            ymom_temp = i.vel()[1]*i.mass
            
            net_mom = (xmom_temp**2 + ymom_temp**2)**0.5
            momentum.append(net_mom)
            
            xmom.append(xmom_temp)
            ymom.append(ymom_temp)
            
            # vel_mag = np.sqrt(np.dot(i.vel(), i.vel()))
            # momentum_temporary = vel_mag*i.mass
            # momentum.append(momentum_temporary)
            
        return [np.array(xmom), np.array(ymom), np.array(momentum)]
    
    def kineticE(self):
        """
        Returns the kinetic energy of all the balls in an array
        Kinetic energy calculated using mean square velocity
        """
        ke = []
        for i in self._object_lst[1:]: #ignores container
            v = i.vel()
            v_sqr = np.dot(v,v)
            ke.append(0.5*i.mass*v_sqr)
        ke = np.array(ke)
        return ke
    
    def temperature(self):
        """
        Returns the temperature of the gas
        """
        avg_ke = np.mean(self.kineticE())
        kb = 1.38*10**(-23) #boltzmann constant
        T = avg_ke/kb #kinetic energy = (2/2)*kb*T, because 2 degrees of freedom in 2D case
        return T 
    
    def pressure(self):
        """
        Calculates pressure. 
        This is a 2D interpretation of pressure, with dimensions of force/length.
        Container circumference is the length used. 
        """
        t = self.tot_time
        d_p = self.del_p
        R = -self.container.radius
        if t == 0:
            return 0 
        else:
            return d_p/(t*2*np.pi*R)
    
    def volume(self):
        """
        Returns the volume of the container.
        Since this is 2D, the "volume" is the area of the container circle.
        """
        vol = np.pi*(-self.container.radius)**2
        return vol
        
        
    def run(self, num_frames, animate=False):
        """
        This code displays the container-gas system if animate is True, and it runs 
        collisions a num_frames number of times. 
        """
        if animate:
            pl.figure()
            cont_rad = -self._object_lst[0].radius #The container radius is negative
            ax = pl.axes(xlim=(-cont_rad, cont_rad), ylim=(-cont_rad, cont_rad))
            ax.set_xlabel("x direction (m)")
            ax.set_ylabel("y direction (m)")
            ax.set_aspect('equal', adjustable='box')
            ax.add_artist(self.container.patch_func())
            for i in range(1, len(self._object_lst)):
                ax.add_patch(self._object_lst[i].patch_func())
                ax.add_patch(self._object_lst[i].patch_func())
        
        for frame in range(num_frames):
            self.next_collision()
            if animate:
             
                # ax.add_patch(self.ball.patch_func())
                for i in range(1, len(self._object_lst)):
                    ax.add_patch(self._object_lst[i].patch_func())
                    ax.add_patch(self._object_lst[i].patch_func())
               
                pl.pause(0.02)
        if animate:
            pl.show()
        
        
        
    def histogram(self):
        """
        Displays histograms for ball distance from container center, and inter-ball distance.
        """
        #now, histogram of distance between each pair of balls
        #This method is used after the run() method
        final_pos = [i.pos() for i in self._object_lst]
        dist_mag = []
        length = len(final_pos)
        for i in range(length):
            if i<(length-1): #to not test for collision with itself
                for j in range(i,length-1):
                    temp1 = final_pos[i]
                    temp2 = final_pos[j+1]
                    dist = temp1 - temp2
                    dist_mag.append(np.sqrt(np.dot(dist, dist)))
                    
        plt.grid()
        plt.hist(dist_mag, 10)
        plt.xlabel('Inter Ball Distances (m)')
        plt.ylabel('Frequency')
        #plt.savefig("done 2-T-9.png", dpi=300)
        plt.show()
        
        
        #distance from center 
        center_disp = [np.sqrt(np.dot(i.pos(), i.pos())) for i in self._object_lst]
        
        plt.grid()
        plt.hist(center_disp, 8)
        plt.xlabel('Displacement from center ($m$)')
        plt.ylabel('Frequency')
        #plt.savefig("done 1-T-9.png", dpi=300)
        plt.show()
    
    def vel_hist(self):
        """
        Displays histogram of velocity magnitudes along with the theoretical Maxwell-Boltzmann fit.
        """
        
        def maxwell(v, data_max): #the maxwell distribution fitted on the data
            #np.array([np.dot(i.vel(), i.vel()) for i in self._object_lst])
            ke = 0.5*self._object_lst[1].mass*v**2 #mass of [1] is the mass of all balls
            kb = 1.38*10**(-23)
            temp = self.temperature()
            p_of_v = v*np.exp(-ke/(kb*temp))
            
            normalise = data_max/max(p_of_v)
            p_of_v = normalise*p_of_v
            
            return p_of_v
            
        graph_velocity = [np.sqrt(np.dot(i.vel(), i.vel())) for i in self._object_lst]
    
        plt.grid()
        data = plt.hist(graph_velocity, 10, label = "Data Histogram")
        max_freq = max(data[0])
        v = np.linspace(0, data[1][-1])
        distrib = maxwell(v, max_freq)
        plt.plot(v, distrib, '-', label = "Theoretical Maxwell-Boltzmann")
        plt.xlabel('Velocity ($ms^{-1}$)')
        plt.ylabel('Frequency')
        plt.legend()
        # plt.savefig('7-T-13.png', dpi = 300)
        plt.show()
               
            
    def temp_var(self):
        """
        A method to show the disparity between the temperatures generated using the mean velocity
        and the mean square velocity of the gas. 
        
        Function returns the temperature of the gascalculated using mean square velocity 
        (kinetic energy), and the temperature calculated using mean velocity
        """
        
        #consts
        kb = 1.38*10**(-23) #boltzmann constant
        ball_m = self._object_lst[1].mass #ball mass
        
        #temperature calculated by kinetic energy of gas, which is dependent on mean square velocity
        avg_sqr_vel_T = self.temperature()
        
        #temperature calculated using mean velocity
        graph_velocity = [np.sqrt(np.dot(i.vel(), i.vel())) for i in self._object_lst]
        mean_vel = np.mean(graph_velocity)
        avg_vel_T = (mean_vel**2)*(3*ball_m*np.pi)/(16*kb)
        
        return([avg_sqr_vel_T, avg_vel_T])
        
        
    def PVNKT(self):
        """
        This method is a convenient way to see the effects of changing the volume of the container
        and seeing the effects on the pressure.
        The container volume is expanded when run. This method is the equivalent of expanding a 
        container full of gas (abruptly, rather than continuously) while keeping temperature constant. 
        
        MUST REDEFINE THE OBJECT TO GET BACK THE ORIGINAL CONTAINER VOLUME. 
        
        Miscellaneous function. 
        """
        res = 100
        radii = np.linspace(60, 200, res)
        N = len(self._object_lst[1:])
        kb = 1.38*10**(-23)
        T = self.temperature()
        
        V_arr = []
        P_arr = []
        for i in range(res):
            self.container.radius = -radii[i]
            V = self.volume()
            
            self.run(200) #run a few frames to get a value for pressure
            P = self.pressure()
            
            V_arr.append(V)
            P_arr.append(P)
            
        V_arr = np.array(V_arr)
        P_arr = np.array(P_arr)
        
        P_theory = N*kb*T/V_arr
        
        plt.grid()
        plt.plot(V_arr, P_arr, 'o', label = "Simulation Data")
        plt.plot(V_arr, P_theory, label = "Ideal Gas")
        plt.xlabel("Volume ($m^2$)")
        plt.ylabel("Pressure (Pa)")
        plt.legend()
        # plt.savefig("P vs. V - .png", dpi = 300)
        plt.show()
        
        
        
        
            
            
            
            
        
     
        