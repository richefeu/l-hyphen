
Format of input files 
=====================


The conf-files are the files that hold the whole configuration at a given time. They are used:

 1. for defining the initial configuration and parameters of a simulation, 
 2. for running some preprocessing commands,
 3. for saving periodically the history of a simulation. The keywords are defined in the following.

Header
------

A conf-file always starts with the header: ``Rockable dd-mm-yyyy`` (*e.g.*, ``Rockable 29-11-2018``). 
Each time a noticeable change is made in the format, the date of this change is also changed in the header of the file. 
It is used as the version of the format. In the source of the code it is defined in the preprocessor define of ``CONF_VERSION_DATE``

Timing
------

- ``t`` (*double*) **value**  
  Current time

- ``tmax`` (*double*) **value**  
  Maximum time (the simulation will end when time is ``tmax``)

- ``dt`` (*double*) **value**
  Time step increments



