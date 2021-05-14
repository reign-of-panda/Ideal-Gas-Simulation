# Ideal-Gas-Simulation

This project aims to simulate an ideal gas and its thermodynamic properties using object oriented programming. 

The file BouncyStuff.py is a library containing the relevant classes for creating a container and gas particle objects. 
y2_thermo_sp1219.py imports the BouncyStuff library, and contains code to conduct thermodynamic experiments, and 
tests the library's methods. Open these two files as a single project and run code blocks in the y2_thermo_sp1219.py file 
to view the outputs that the simulation can generate. 

A particle-container system can be initialised using:

> pc_sys = bs.Simulation(mass, radius, position, velocity, particle_num)

The system can be set up with 1, 16, or 100 particles. For multiple particles, the position is pre-determined, hence an
arbitrary value can be entered for the parameter. 
To run the simulation, the following line is used:

> pc_sys.run(num_frames, animate)

where the number of frames determines how many particle collisions will occur.
When the animate parameter is set to True, a particle collision animation is displayed

In order to view the animation, the following line must be executed in the console (applicable for Spyder IDE):

> %matplotlib auto

Class Simulation (within BouncyStuff.py) contains various methods that are used to inspect thermodynamic properties. 

