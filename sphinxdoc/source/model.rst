Model
=====


Generalities about DEM methods
------------------------------

The model we have built is substantially similar to a form of DEM. The basic concept is to model a material as a collection of rigid elements that interact. The materials we are particularly interested in are cells of the plant tissue. So we agreed to create polygons to best represent our cells (Fig. :numref:`Cells`). The nodes/summits of the polygons would be our discrete elements and we connect them by bars. When the bars close the polygon, we call it a cell (Fig. :numref:`Cells`). 


.. _Cells:
.. figure:: images/cellModel.png
   
   Representation of a cell. Each node is numbered from 0 to 5 in this case.
   
   
Let us recall the DEM vu-quoc_3-d_2000. Let :math:`m` be the mass of the nodes, :math:`I` the inertia along the normal axis, :math:`\vec{\overset{\cdot \cdot}{x}}_i` the acceleration of the node, :math:`\vec{\overset{\cdot \cdot}{\theta}}_i` the angular acceleration, :math:`\vec{f^{j \rightarrow i}}` the force due to the action of the particle :math:`j` on the particle :math:`i`, :math:`\vec{b}` the lever-branch vector from particle center to the contact point and :math:`\vec{g}` the volume acceleration (like gravity in most cases). We then have :

.. math::
   \begin{split}
   m \vec{\overset{\cdot \cdot}{x}}_i &= \underset{j \neq i}{\sum} \vec{f^{j \rightarrow i}} + m \vec{g}\\
   I \vec{\overset{\cdot \cdot}{\theta}}_i &= \underset{j \neq i}{\sum} ||\vec{b} \times \vec{f^{j \rightarrow i}}|| +
   \underset{j \neq i}{\sum} C^{j \rightarrow i}
   \end{split}




