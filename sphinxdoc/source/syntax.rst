
Format of input files 
=====================

The input files are the files that hold the whole configuration. 
They are used for defining the initial configuration and parameters of a simulation. 

Timing
------

- ``t`` (*double*) **value**  
  Current time. For a new simulation, this time is normally *0.0*

- ``dt`` (*double*) **value**
  Time step increments. It needs to be set thanks to similar condition to DEM (remember the masses are lumped to the cell-nodes)

- ``nstep`` (*double*) **value**
  Number of step increments. Final time is ``t + nstep * dt``

Dissipation
-----------

- ``numericalDissipation`` (*double*) **value**
  Value of energy dissipation (for example *1e-4*) 

- ``globalViscosity`` (*double*) **value**
  Add a viscosity on the sample to dissipate the energy

Volume Forces
-------------

- ``gravity`` (*double*) **value** (*double*) **value** 
  Apply a force (expressed by a 2D-vector) on nodes 


Build and rebuild of neighborhood
---------------------------------
 
- ``distVerlet`` (*double*) **value**
  Distance above which two elements are considered close enough to be part of the neighbor list.

- ``nstepPeriodVerlet`` (*double*) **value**
  Number of time-step between updates of the neighbor list.


Outputs
-------

- ``isvg`` (*integer*) **value**
  XXXXXX

- ``nstepPeriodSVG`` (*integer*) **value**
  XXXXXX
  
- ``findDisplayArea`` (*double*) **value**
  Find the limits of the display area in svg files. **value** is a scale factor.

- ``iconf`` (*integer*) **value**
  Number of configuration file (conf*x*.txt) that will be saved

- ``nstepPeriodConf`` (*integer*) **value**
  

Parameters for cell interactions
--------------------------------

- ``kn`` (*double*) **value**
  Normal contact-stiffness

- ``kt`` (*double*) **value**
  Tangential contact-stiffness 

- ``mu`` (*double*) **value**
  Friction coefficient (between cells)

- ``fadh`` (*double*) **value**
  Normal contact-adhesion force

Pre-processing (commands not saved in further conf-files)
---------------------------------------------------------

- ``readNodeFile`` (*const char*) **file.txt** (*double*) **barWidth** (*double*) **Kn** (*double*) **Kr** (*double*) **M_Y**
  Cette fonction permet de lire un fichier contenant une liste de positions x,y avec numéro de cellule. Peu importe les numéros tant qu'ils sont différents pour chaque cellule. **barWidth** is the thickness given to all bars. **Kn** and **Kr** sont les coefficients respectivement associés à la raideur axiale des barres et la raideur angulaire entre les barres adjacentes. **M_Y** est le moment seuil plastique.

Imposed Controls
----------------

- ``setNodeControl`` (*size_t*) **c** (*size_t*) **n** (*integer*) **xmode** (*double*) **xvalue** (*integer*) **ymode** (*double*) **yvalue**
  Ajoute un control au noeud **n** de la cellule **c**. **xmode** vaut 1 ou 0 respectivement pour le contrôle de la velocité ou le contrôle de la force. 

- ``setCellControl`` (*size_t*) **c** (*integer*) **xmode** (*double*) **xvalue** (*integer*) **ymode** (*double*) **yvalue**
  Ajoute un control à tous les noeuds de la cellule **c**. **xmode** vaut 1 ou 0 respectivement pour le contrôle de la velocité ou le contrôle de la force. 


- ``setNodeControlInBox`` (*double*) **xmin** (*double*) **xmax** (*double*) **ymin** (*double*) **ymax** (*integer*) **xmode** (*double*) **xvalue** (*integer*) **ymode** (*double*) **yvalue** 
  Définit un même control à tous les noeuds qui se trouvent dans une zone rectangulaire avec les marqueurs **xmin**, **xmax**, **ymin** et **ymax**. 
  
- ``glue`` (*double*) **value**
  Ajoute un lien de colle entre les barres de cellules différentes dont la distance est inférieure à **value**.
  
- ``setGlueSameProperties`` (*double*) **kn_coh** (*double*) **kt_coh** (*double*) **fn_coh_max** (*double*) **ft_coh_max** (*double*) **yieldPower**
  Doit être placée après la fonction ``glue``. **kn_coh** est la raideur normale de cohésion. **kt_coh** est la raideur tengentielle de cohésion. **fn_coh_max** est le seuil de la force normale de cohésion solide. **ft_coh_max** est le seuil de la force tangentielle de cohésion solide
  


