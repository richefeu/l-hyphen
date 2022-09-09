Quick-start guid
================

What is l-hyphen?
-----------------

``l-hyphen`` is a DEM code written in C++. blabla


Compilation
-----------

The compilation is managed by a ``makefile`` that should be similar to that one:

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


The applications run and see can be compiled with the following command:

.. code-block:: sh

   make

It is sometimes necessary to remove all object files (.o) together with the compiled applications. this can be made with:

.. code-block:: sh

   make clean


Running a simulation (brief overview)
-------------------------------------


To run a simulation, a configuration file has to be written. The format of such a file is described in the section "Format of input file" for conf-files. We show here a simple example simulating cell being crushed.

.. code-block:: text
   :caption: input.txt
   
   example truc (TODO)
   
   
   




