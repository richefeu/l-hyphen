Quick-start guid
================

What is L-hyphen ?
------------------

``l-hyphen`` is a DEM code written in C++11. blabla


Compilation
-----------

The compilation is managed by a makefile that should be similar to that one:

.. code-block:: makefile
   :caption: Makefile
   
   # The compiler to be used
   CXX = g++

   # The list of flags passed to the compiler
   CXXFLAGS = -Wall -Wextra -ansi -O3 -Wwrite-strings -Wstrict-prototypes -Wuninitialized \
   -Wunreachable-code -Wno-missing-braces  -Wno-missing-field-initializers \
   -std=c++0x -I ~/toofus

   # Notice that you need to install toofus before compiling

   APP = l-hyphen

   .PHONY: all clean

   all: $(APP)

   clean:
   	rm -f *.o $(APP)

   $(APP): $(APP).cpp $(APP).hpp
   	$(CXX) $(CXXFLAGS) $(APP).cpp -o $(APP)


To compile the program you just have to : 

.. code-block:: sh

   make

Sometimes, it necessar to delete the previous object and recompile so you can :

.. code-block:: sh

   make clean



Running a simulation
--------------------


To run a simulation, a configuration file has to be written. The format of such a file is described in the section Syntax for conf-files. We show here a simple example simulating a sphere bouncing on a plan.

.. code-block:: text
   :caption: input.txt
   
   # SIMULATION INPUTS ===================================================

   # time flow
   # t 0
   dt 1e-5 
   nstep  1000000

   # dissipation
   numericalDissipation 1e-4
   globalViscosity      0.0

   # volume forces
   gravity 0 0

   # proximity 
   distVerlet 0.01
   nstepPeriodVerlet 100

   # outputs
   isvg              0
   nstepPeriodSVG    5000
   iconf             0
   nstepPeriodConf   10000

   # contact adhésif frottant entre les cellules
   kn 1000000
   kt 1000000
   mu 0.3

   # SAMPLE ===================================================

   #             filename       barWidth  kn       kr    mz_max
   readNodeFile  nodeFile.txt  0.0045     1000000  100  100


   # IMPOSED CONTROLS
   setNodeControlInBox 0 1 -0.5 0.05 0 0.0 0 -1e-2
   setNodeControlInBox 0 1 0.95 1.5 0 0.0 0 0

   # PREPRO ===================================================
   # activation de la colle si la distance est inférieure à 0.02
   glue 0.005

   setGlueSameProperties 1e5 1e5 2000 2000 5

   findDisplayArea 1.0
   
To run the simulation with this input file, you have to write this command : 

.. code-block:: sh

   path/to/l-hyphen input.txt

