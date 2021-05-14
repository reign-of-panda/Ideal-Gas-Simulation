# -*- coding: utf-8 -*-
#Author: Swoyam Panda
"""
TASK 1-3

These tasks mostly involved coding in the BouncyStuff.py file.
"""


#%%
"""
TASK 4
"""

import BouncyStuff as bs
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

#defining a Simulation class object for use later on.
s_b = bs.Simulation(1, 1, [-5, 0], [1,0], 1)

"""
The following two functions were used to run collisions frame by frame for two ball systems. They are 
now outdated as the classes in BouncyStuff have been edited extensively since the beginning of the 
project.  
"""

def col_output(obj): #runs a collision and displays the system of particles and container
    obj.next_collision()
    pl.figure()
    ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
    ax.set_xlabel("x direction (m)")
    ax.set_ylabel("y direction (m)")
    ax.set_aspect('equal', adjustable='box')
    ax.add_artist(obj.container.patch_func())
    ax.add_patch(obj._object_lst[1].patch_func())
    ax.add_patch(obj._object_lst[2].patch_func())    
    pl.show()
    
def only_disp(obj): #displays the system of particles and container
    pl.figure()
    ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
    ax.set_xlabel("x direction (m)")
    ax.set_ylabel("y direction (m)")
    ax.set_aspect('equal', adjustable='box')
    ax.add_artist(obj.container.patch_func())
    ax.add_patch(obj._object_lst[1].patch_func())
    ax.add_patch(obj._object_lst[2].patch_func())
    pl.show()
    

#%%
"""
This piece of code is to test the logic of my time_to_collision method.
Changing the values of r and v here to see how it affects my dt output. 
"""

r = np.array([-0.5,  0. ]) - np.array([1.5, 0. ])
v = np.array([-1, 0]) - np.array([1, 0])
rv = np.dot(r, v)
#_sqr variables contain vector magnitude squared
v_sqr = np.dot(v, v)
r_sqr = np.dot(r, r)
"""
The radius of the container will be entered 
as a negative number. This accounts for the 
R1-R2 case. 
"""
R = 1+1
        
#Use the positive value now
dt = (-2*rv+(4*(rv**2)-4*v_sqr*(r_sqr-R**2))**0.5)/(2*v_sqr)
print("dt", dt)

"""
WHOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
Problem discovered 
My delta r calculation is wonky for some reason. 

Update: Reason discovered and fixed. Error in quadractic equation entry. 
"""

#%%
"""
REPEAT
"""

r = np.array([-8,  0. ]) - np.array([9, 0. ])
v = np.array([1, 0]) - np.array([-1, 0])
rv = np.dot(r, v)
#_sqr variables contain vector magnitude squared
v_sqr = np.dot(v, v)
r_sqr = np.dot(r, r)
"""
The radius of the container will be entered 
as a negative number. This accounts for the 
R1-R2 case. 
"""
R = 1+1
        
#Use the positive value now
#ALTER THE SIGN
dt = (-2*rv-(4*(rv**2)-4*v_sqr*(r_sqr-R**2))**0.5)/(2*v_sqr)
print("dt", dt)


"""
CONCLUSION: It makes more sense to consider the polarity of r dot v when writing 
conditions for the time_to_collision method

Test cases examined:

1. Balls moving towards each other head on
2. Balls moving in the same direction with the same velocity
3. Balls moving in the same direction, with one ball moving faster than the other
3.1 Change which ball is moving faster
4. The two balls start out touching each other 
4.1 Initial velocity moving towards each other
4.2 Initial velocity moving away from each other

It's all working properly in class Ball now. 
"""

#%%
"""
TASK 5 - Testing to see if the animation works for one 

To view the animation, enter

%matplotlib auto

into the console.

"""
num_frames = 100
s_b = bs.Simulation(1, 1, [-5, 0], [1,0], 1)
s_b.run(num_frames, True)

#The animation works!

#%%
"""
TASK 6
"""
#This code is to calculate the kinetic energy

print("initial kinetic energy:", s_b.kineticE()) 
#NOTE: kinetic energy is returned as a list of inividual balls' kinetic energies. 

#%%
"""
TASK 7 

SIMULATION class can define objects with 1, 100, and, 16 balls. 
"""
#Edit the Simulation class to have collisions of more than one particle. 
#Currently set to have 16 particles
num_frames = 50
s_b2 = bs.Simulation(1, 1, [-5, 0], [1,0], 16)
s_b2.run(num_frames, True)


#%%
"""
A case where the ball doesn't collide and should circle the inside of the container. 

The ball doesn't move if the ball is touching the container and the velocity direction is tangent to 
the container's side. 

The particles are never touching the container initially, so this is not a problem one will run into.
"""

num_frames = 20
test = bs.Simulation(1, 1, [0, 59], [1,0], 1)
test.run(num_frames, True)

#%%
"""
Check if the correct number of collisions happen when the simulation is run. 
"""
intended_num = 30
counter = bs.Simulation(1, 1, [0, 0], [10, 0], 100)
counter.run(intended_num)
actual_num = counter.col_count() 
#col_count() is a method that keeps track of the number of colllisions

print("Number of collisions that happen:", actual_num)
print(actual_num == intended_num)
#The collisions are happening the correct number of times if the above statement prints true

#%%
"""
Define some constants
"""
mass_He = 6.64*10**(-27)
helium_r = 1.4*10**(-10)
kb = 1.38*10**(-23)

"""
CONSTANTS DEFINED

I have set the container radius to be 60m for most of the upcoming tasks. 
"""

#%%

"""
From this point onwards, the simulations are run with the visualiser switched off. 
This saves time when exploring the thermodynamic properties of this system. 
"""

#%%
"""
Try 100 balls and code the PV = NkbT relation.
"""
import BouncyStuff as bs
num_frames3 = 200
N = 100 #number of particles
kb = 1.38*10**(-23) #boltzmann constant
radius = 0.01 #vary to have a look at the effect on PV = NkbT
hundred = bs.Simulation(1, radius, [0, 0], [1,0], N)
hundred.run(num_frames3, False)

T_hundred = hundred.temperature()
col_pressure = hundred.pressure()
container_vol = hundred.volume()

PV = col_pressure*container_vol
RHS = N*kb*T_hundred

print("PV", PV)
print("NkbT", RHS)


#The values printed are fairly close
#THIS WORKS AND IS CONSISTENT WITH PHYSICS WHOOOOOOOOOO
#(for small radis)


#%%
"""
TASK 9 
Histograms

PLOTS 1. T-1 and 2. T-2

The histogram function plots the interball distances and the displacement from center

Use %matplotlib inline for these plots
"""
hist_exp = bs.Simulation(1, 1, [0, 0], [1000, 0], 100)
hist_exp.run(300, False)

hist_exp.histogram()


#%%
"""
Continuing with task 11 plots

KINETIC ENERGY CONSERVATION GRAPH

PLOT 3. T-11
"""
#system kinetic energy with time

#same thing with all the helium molecule variables
import matplotlib.pyplot as plt
import numpy as np

KE_experiment = bs.Simulation(mass_He, helium_r, [0, 0], [1000, 0], 100)

KE_measure = [sum(KE_experiment.kineticE())]
momentum = [0]
KE_time = [0]
#200 collisions
for i in range(45): #45 frames
    KE_experiment.run(1)
    #kineticE() returns a kinetic energy array
    KE_temporary = sum(KE_experiment.kineticE())
    KE_measure.append(KE_temporary)
    KE_time.append(KE_experiment.tot_time)

plt.grid()
plt.plot(KE_time, KE_measure, '-o', label = "Kinetic Energy in Simulation (45 frames)")
plt.xlabel("Time (s)")
plt.ylabel("Kinetic Energy (J)")
plt.legend()
# plt.savefig("3-T-11.png", dpi=300)
plt.show()

#%%
"""
Continuing with task 11 plots
MOMENTUM VS. TIME

PLOT 4. T-11
"""
import matplotlib.pyplot as plt
import BouncyStuff as bs
import numpy as np

momentum_exp = bs.Simulation(mass_He, helium_r, [0, 0], [1000, 0], 100)
momentum_measure = [sum(momentum_exp.momentum())]
momentum_time = [momentum_exp.tot_time]
momx = [sum(momentum_exp.momentum()[0])]
momy = [sum(momentum_exp.momentum()[1])]

for i in range(45): #45 frames
    momentum_exp.run(1)
    momentum_x = sum(momentum_exp.momentum()[0])
    momentum_y = sum(momentum_exp.momentum()[1])
        
    momx.append(momentum_x)
    momy.append(momentum_y)
    # momentum_measure.append(momentum_temporary)
    momentum_time.append(momentum_exp.tot_time)

the_sum = np.array(momx)+np.array(momy)
plt.grid()
plt.plot(momentum_time, momx, '-o', label = "x axis momentum")
plt.plot(momentum_time, momy, '-o', label = "y axis momentum")
# plt.plot(momentum_time, net_mom, '-o', label = "net momentum") This value is clearly not conserved
plt.plot
# plt.plot(momentum_time, momentum_measure, '-o', label= "Momentum in Simulation (45 frames)")
plt.xlabel("Time (s)")
plt.ylabel("Momentum ($kgms^{-1}$)")
plt.legend()
# plt.savefig("4-T-11.png", dpi=300)
plt.show()

#%%

"""
Onto TASK 11

PRESSURE VS. TEMPERATURE GRAPH

Graph generated: 5. T-11
"""
import matplotlib.pyplot as plt
import BouncyStuff as bs
#theoretical P-T graph (keep N and V constant)
V = np.pi*60**2 #This is an area as the simulation is 2D
N = 100 #number of particles
kb = 1.38*10**(-23) #boltzmann constant


#Now, measure pressure using the simulation - vary temperature by varying kinetic energy. 

#Generates 100 velocity distributions, and 100 temperatures
num_frames = 40
measurements = 20
mean_velocity = np.linspace(0.1, 1400, measurements) 
#Use the helium data to initialise the PT_experiment object. 
mass_He = 6.64*10**(-27)
helium_r = 1.4*10**(-10)

P_measure = []
T_measure = []
for i in range(measurements):
    PT_experiment = bs.Simulation(mass_He, helium_r, [0, 0], [mean_velocity[i],0], 100)
    PT_experiment.run(num_frames, False)
    P_temporary = PT_experiment.pressure()
    T_temporary = PT_experiment.temperature()
    P_measure.append(P_temporary)
    T_measure.append(T_temporary)
    
#convert into array because it's a convenient data type to work with
T_measure = np.array(T_measure)
P_measure = np.array(P_measure)

#some fitting:

#T and P below are for the theoretical P-T curve
T = np.linspace(min(T_measure), max(T_measure), 10000) 
P = N*kb*T/V


fitted, cov1 = np.polyfit(T_measure,P_measure,1, cov=True)
fit_func = np.poly1d(fitted)

plt.grid()
plt.plot(T_measure, P_measure, 'o', label = "Simulation Calculation")
plt.plot(T, fit_func(T), label = 'poly1d fit')
plt.plot(T, P, label = "Theoretical P-T graph")
#labels
plt.ylabel("Pressure (Pa)")
plt.xlabel("Temperature (K)")
plt.legend()
# plt.savefig("T-11 - radius=1.png", dpi=300)
plt.show()

#calculate uncertainty of fit

grad_uncertain = np.sqrt(cov1[0][0])
fit_gradient = fitted[0]

theoretical_m = N*kb/V

#print out the gradient and its uncertainty
print(fit_gradient, 'pm', grad_uncertain)

print(theoretical_m)


#%%
"""
Task 12 plot - EQUATION OF STATE - HOW IT VARIES WITH BALL RADIUS

PLOT 6. T-12

WARNING: Takes a reaaaaaaly long time to run. 
"""

import BouncyStuff as bs
num_frames4 = 600
N = 100 #number of particles
kb = 1.38*10**(-23) #boltzmann constant
radius = np.linspace(0, 1, 50)#vary to have a look at the effect on PV = NkbT

container_vol = hundred.volume() #this is just constant, no need to include in loop

del_value = [] #This contains an array of PV-NkbT values
for i in radius:
    hundred = bs.Simulation(1, i, [0, 0], [1,0], N)
    hundred.run(num_frames4, False)

    T_hundred = hundred.temperature()
    col_pressure = hundred.pressure()

    PV = col_pressure*container_vol
    RHS = N*kb*T_hundred
    del_value.append(PV-RHS)

#just placing the plotting part in a different cell

fit_r, cov_r = np.polyfit(radius, del_value, 1, cov=True)
plot_func = np.poly1d(fit_r)

plt.grid()
plt.plot(radius, del_value, 'o', label = '$PV-Nk_bT$ vs. Radius)')
plt.plot(radius, plot_func(radius), '-', label = 'Best Fit Line')
plt.xlabel('Radius (m)')
plt.ylabel('$PV-Nk_bT$')
plt.legend()
# plt.savefig("6-T-12 (attempt 2).png", dpi=300)
plt.show()

#600 collisions were run before taking a pressure and temperature "measurement" from the simulation.

#%%
"""
Task 13
Plot 7. T-13

Plotting the velocity histogram and comparing it to an ideal Maxwell-Boltzmann distribution
Not bad
"""
import BouncyStuff as bs

v_hist = bs.Simulation(mass_He, helium_r, [0, 0], [1000, 0], 100)
data = v_hist.vel_hist()

#%%
"""
A few tests to check for collision number conservation
"""
col_num = 45
import BouncyStuff as bs
hundred = bs.Simulation(mass_He, 1, [0, 0], [1000, 0], 100)
hundred.run(col_num) #Doing a few collisions here before running the upcoming lines

x1 = hundred.col_count()
x2 = hundred.self_count()

print("Total collisions in the system:", x1)
print("Sum of collision number of each object:", sum(x2)) #Should be twice the value of x1

#I was curious. I just wanted to know how many of the total collisions happen with the container. 
print("Number of collisions with the container:", x2[0])
#%%
"""
Check if particle number is conserved (ie. no particles escape the container)
"""
N_consv = hundred.particle_conserve() #this method returns the number of particles in the container
print(N_consv)


#%%
"""
Tests done.
Task 14 - Another PT graph

PLOT 8. T-14
"""
import matplotlib.pyplot as plt
import BouncyStuff as bs
#theoretical P-T graph (keep N and V constant)
V = np.pi*60**2 #This is an area as the simulation is 2D
N = 100 #number of particles
kb = 1.38*10**(-23) #boltzmann constant


#Now, measure pressure using the simulation - vary temperature by varying kinetic energy. 

#Generates 100 velocity distributions, and 100 temperatures
num_frames = 200
measurements = 20
mean_velocity = np.linspace(0.1, 1400, measurements) 
#Use the helium data to initialise the PT_experiment object. 
mass_He = 6.64*10**(-27)
helium_r = 1.4*10**(-10)

P_measure = []
T_measure = []
for i in range(measurements):
    PT_experiment = bs.Simulation(mass_He, 1, [0, 0], [mean_velocity[i],0], N)
    PT_experiment.run(num_frames, False)
    P_temporary = PT_experiment.pressure()
    T_temporary = PT_experiment.temperature()
    P_measure.append(P_temporary)
    T_measure.append(T_temporary)
    
#Fit a line to the data

#T and P below are for the theoretical P-T curve
T = np.linspace(min(T_measure), max(T_measure), 10000) 
P = N*kb*T/V

fitted, cov1 = np.polyfit(T_measure,P_measure,1, cov=True)
fit_func = np.poly1d(fitted)

plt.grid()
plt.plot(T_measure, P_measure, 'o', label = "Simulation Data")

plt.plot(T, fit_func(T), label = 'Best Fit (ball radius = 1m)')
plt.plot(T, P, label = "Ideal Gas P-T Curve")
#graph labels
plt.ylabel("Pressure (Pa)")
plt.xlabel("Temperature (K)")
plt.legend()
# plt.savefig("8-T-14 - radius=1.png", dpi=300)
plt.show()


#Calculate the value of b and the corresponding error. 
    
grad = fitted[0]
b = V/N - kb/grad #using the expression for the gradient of function P(T) for Van der Waal's law
grad_err = np.sqrt(cov1[0][0])

def b_err(grad = grad, grad_err=grad_err, N=N, kb=kb): #function to calculate error for b
    fraction_err = grad_err/grad
    b_abs_err = (kb/grad)*fraction_err
    return b_abs_err

b_error = b_err() #parameters already passed into function definition



print("Value of b:", b, "plus-minus", b_error)
print("Units: meters squared")

#for a helium radius measurement, b = 0.735617069133383 plus-minus 0.35638374264263173(maybe wrong?)

#Another measurement: Value of b: 0.647767450722796 plus-minus [0.23361882]
#also at helium radius


#%%
"""
EXTENSION - My code for task 14 was giving me negative values for b, so I wanted to see if it was 
because the number of collisions before each measurement was too low. 

Check the stability of pressure values.
Use 3 ball radii to see the impact of that.
"""
cols = 200 #change this to 1000 to get the same x axis as the diagram in my report
num = 100

pres_test = bs.Simulation(mass_He, 1, [0, 0], [1000,0], num)
pres_lst = [pres_test.pressure()]

pres_test2 = bs.Simulation(mass_He, 0.01, [0, 0], [1000,0], num)
pres_lst2 = [pres_test2.pressure()]

pres_test3 = bs.Simulation(mass_He, helium_r, [0, 0], [1000,0], num)
pres_lst3 = [pres_test3.pressure()]

frames = np.arange(0, cols+1, 1)
for i in range(cols):
    pres_test.run(1)
    pres_lst.append(pres_test.pressure())
    
    pres_test2.run(1)
    pres_lst2.append(pres_test2.pressure())
    
    pres_test3.run(1)
    pres_lst3.append(pres_test3.pressure())

pres_lst = np.array(pres_lst)
pres_lst2 = np.array(pres_lst2)
pres_lst3 = np.array(pres_lst3)

plt.grid()
plt.plot(frames, pres_lst, label = "Ball radius = 1m")
plt.plot(frames, pres_lst2, label = "Ball radius = 0.01m")
plt.plot(frames, pres_lst3, label = "Ball radius = $r_{He}$")


plt.xlabel("Number of Frames Run")
plt.ylabel("Pressure (Pa)")
plt.legend()
# plt.savefig("Pressure stability", dpi = 300)
plt.show()

#%%
"""
EXTENSION-2
Test the simulation to check if the relation between average velocity and temperature holds. 

The temp_var() method was defined in Class simulation specifically for this investigation.
"""
import BouncyStuff as bs
import matplotlib.pyplot as plt
import numpy as np

iterations = 30 #change this to 5000 to get similar graph as is on the report (takes a long time though)
mass_He = 6.64*10**(-27)
thermal_speed_He = 850   
helium_r = 1.4*10**(-10) #Van der Waal's radius


helium = bs.Simulation(mass_He, 1, [0, 9], [thermal_speed_He,0], 100)


x1 = [helium.temp_var()[0]] #gotta have the initial values
x2 = [helium.temp_var()[1]]
cols = np.array([n for n in range(iterations+1)]) #number of collision array
for i in range(iterations):
    helium.run(1, False)
    mean_sqr_vel = helium.temp_var()[0]
    mean_vel = helium.temp_var()[1]
    x1.append(mean_sqr_vel)
    x2.append(mean_vel)
x1 = np.array(x1)
x2 = np.array(x2)

plt.grid()
plt.ylim(min(x2)-10, max(x1)+15)
plt.plot(cols, x1, label = "T calculated using $\\langle v^2 \\rangle$")
plt.plot(cols, x2, '.', label = "T calculated using $\\langle v \\rangle$")
plt.xlabel("Number of Collisions")
plt.ylabel("Temperature (K)")
plt.legend()
#plt.savefig("Extension 1.2.png", dpi = 300)
plt.show()


