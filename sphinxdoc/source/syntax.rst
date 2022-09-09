
Format of input files 
=====================

The input files are the files that hold the whole configuration. They are used for defining the initial configuration and parameters of a simulation : 

Timing
------

- ``t`` (*double*) **value**  
  Current time

- ``tmax`` (*double*) **value**  
  Maximum time (the simulation will end when time is ``tmax``)

- ``dt`` (*double*) **value**
  Time step increments

- ``nstep`` (*double*) **value**
  Number of step increments

Dissipation
-----------

- ``numericalDissipation`` (*double*) **value**
  Value of energy dissipation (for example 1e-4) 

- ``globalViscosity`` (*double*) **value**
  Add a viscosity on the sample to dissipate the energy

Volume Forces
-------------

- ``gravity`` (*double*) **value** (*double*) **value** 
  Apply a force (expressed by a 2D-vector) on nodes 


Proximity
---------
 
- ``distVerlet`` (*double*) **value**
  

- ``nstepPeriodVerlet`` (*double*) **value**
  


Outputs
-------

- ``isvg`` (*integer*) **value**
  

- ``nstepPeriodSVG`` (*integer*) **value**
  

- ``iconf`` (*integer*) **value**
  

- ``nstepPeriodConf`` (*integer*) **value**
  

Adhesive Contact Rubbing Between the Cells
------------------------------------------

- ``kn`` (*double*) **value**
  raideur normale de contact  

- ``kt`` (*double*) **value**
  raideur tangentielle de contact

- ``mu`` (*double*) **value**
  coefficient de frottement (entre les cellules)

- ``fadh`` (*double*) **value**
  force normale d'adhésion au contact

Sample
------

- ``readNodeFile`` (*const char*) **file.txt** (*double*) **barWidth** (*double*) **Kn** (*double*) **Kr** (*double*) **Mz_max**
  Cette fonction permet de lire un fichier contenant une liste de positions x,y avec numéro de cellule. Peu importe les numéros tant qu'ils sont différents pour chaque cellule. **barWidth** est l'épaisseur de toutes les barres. **Kn** et **Kr** sont les coefficients respectivement associés à la raideur axiale des barres et la raideur angulaire entre les barres adjacentes. **Mz_max** est le moment seuil plastique.

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
  



Display
-------

- ``findDisplayArea`` (*double*) **value**
  trouver les limites de la zone dessinée dans les fichiers svg. **value** est un facteur multiplicateur
