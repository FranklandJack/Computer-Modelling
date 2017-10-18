# Computer-Modelling
Java programs to model the solar system
This was written by Jack Frankland, Micheal Spence and Martin Tasker as their second year Computer Modelling project at the University of Edinburgh.

Authors:
Jack Frankland s1404032	
Michael Spence s1419697
Martin Tasker s1414001

Instructions:

from the directory

run; make

run; java NBody input/ParticleInput.txt input/ParameterInput.txt output/Output.txt 



Where the first command line argument is ParticleInput.txt, the second ParameterInput.txt, and the third the output file to name in the command line. 

The NBody Program will automatically output a file called "Energy.txt" that outputs the energy of the system at each timestep to a file.
NBody will automatically output a file called "Orbits.txt" which prints out the number of orbits around the sun for the every particle in the simulation. 
NBody will also automatically output a file called "PerihelionAphelion" which gives the perhelion and aphelion of each planet. 

The ParticleInput.txt is formatted as:

number of particles
x1 y1 z1 
V1_x V1_y V1_z
mass1
label1
x2 y2 z2 
V2_x V2_y V2_z
mass2
label2
…
(repeat for number of particles)
…

position is given in Au
velocity is given in Au/s
mass is given in earth masses

ParameterInput.txt is formatted as:

timestep
number step 
gravitational constant (0.00000000000000000011905)

gravational constant is given in (Au)^3 (masses of the earth)^-1 (s)^-2.

Output.xyz is formatted as:

number of particles
point number1
x1 y1 z1 
x2 y2 z2 
x3 y3 z3
…
number of particles
point number2
x1 y1 z1 
x2 y2 z2 
x3 y3 z3
…
(repeated for number of steps)/10
…

position is given in Au.
