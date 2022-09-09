Model
=====


Generalities
------------

The model we have built is substantially similar to a form of DEM. The basic concept is to model a material as a collection of rigid elements that interact. The materials we are particularly interested in are cells of the plant tissue. So we agreed to create polygons to best represent our cells (Fig. :numref:`Cells`). The nodes of the polygons would be our discrete elements and we connect them by bars. When the bars close the polygon, we call it a cell (Fig. :numref:`Cells`). 


.. _Cells:
.. figure:: _static/img/cellModel.png
   
   Representation of a cell. Each node is numbered from 0 to 5 in this case.
   
   
Let us recall the DEM. Let :math:`m` be the mass of the nodes, :math:`I` the inertia along the normal axis, :math:`\vec{\overset{\cdot \cdot}{x}}_i` the acceleration of the node, :math:`\vec{\overset{\cdot \cdot}{\theta}}_i` the angular acceleration, :math:`\vec{f^{j \rightarrow i}}` the force due to the action of the particle :math:`j` on the particle :math:`i`, :math:`\vec{b}` the lever-branch vector from particle center to the contact point and :math:`\vec{g}` the volume acceleration (like gravity in most cases). We then have:

.. math::
   \begin{split}
   m \vec{\overset{\cdot \cdot}{x}}_i &= \underset{j \neq i}{\sum} \vec{f^{j \rightarrow i}} + m \vec{g}\\
   I \vec{\overset{\cdot \cdot}{\theta}}_i &= \underset{j \neq i}{\sum} ||\vec{b} \times \vec{f^{j \rightarrow i}}|| +
   \underset{j \neq i}{\sum} C^{j \rightarrow i}
   \end{split}


The DEM method essentially consists in integrating its equations to find a weak solution. The quadrature we use is the velocity-Verlet scheme: 

.. math::
   \begin{split}
   \vec{x}^{(t+\delta t)}_i &= \vec{x}^{(t)}_i + \delta t \vec{\overset{\cdot}{x}}^{(t)}_i + \frac{\delta t}{2}      
   \vec{\overset{\cdot \cdot}{x}}^{(t)}_i \\
   \vec{\overset{\cdot}{x}}^{(t+\delta t)}_i &= \vec{\overset{\cdot}{x}}^{(t)}_i + \delta t \frac{\vec{\overset{\cdot      \cdot}{x}}^{(t)}_i + \vec{\overset{\cdot \cdot}{x}}^{(t+ \delta t)}_i}{2} 
   \end{split}


We wish to study the rupture of a very heterogeneous structure such as that of a plant (Fig. :numref:`PlantStruc`). Naturally, the geometry of the object we wish to fragment allows us not to take into account its depth and thus to propose a 2D model. Indeed, a stem is essentially a tube whose section has a complex structure but a tube all the same. 


.. _PlantStruc:
.. figure:: _static/img/PlantStruc.png

   Micrograph of a plant structure. Some plant tissues like this one can be very regular with cells close to what a   
   Voronoi-Laguerre diagram would give us. We can observe (despite a great regularity), that some cells are broken or 
   rather crushed, giving them a more complex geometry}


Plant tissue can be seen as an aggregation of cells. A cell consists of an outer membrane that is more difficult to break. This observation leads us to construct cells from bars and nodes that are connected to create cells (Fig. :numref:`WhatACell`). The space between the cells is designated as the fragile zone. The bars that are parts of the membranes will be designated as the fragile areas of the model. We apply a glue bond. The rupture will therefore be created when the force that pulls two cells apart reaches a certain threshold.


.. _WhatACell:
.. figure:: _static/img/Cells.png

   Example of aggregation of identical, convex and regular cells.
   The cells must be joined by a glue relationship.
   
   
This model gives us a clear advantage over a classic DEM because it allows us to work on materials not only rigid.


External forces
---------------


Our model is close to a material point model for the nodes. The description of a node is reduced to the position of its center of gravity, its mass and its link (like the angle with the next bar) (Fig \ref{fig:NodesBarsOr}) with the second elements of our model: the bars (Fig . The bars work like springs which have a plasticity, a stiffness and a length. 

\begin{figure}[htbp]
  \centering
  \includegraphics[scale = 0.9]{Figures_histoDEM_orientation.pdf}
  \caption{Set of characteristics of the bar class taken into account by node}
  \label{fig:NodesBarsOr}
\end{figure}

In an implementation way, we do a test to know if there is contact or not (Fig \ref{fig:NodesBarsCon}). Our bars have a kind of contactbox that we use to determine if we should distribute a possible interaction to the nodes. The forces (interactions) linked to the bars are transferred to the concerned nodes 

\begin{figure}[htbp]
  \centering
  \includegraphics[scale = 1]{Figures_histoDEM_contact.pdf}
  \caption{Nodes and Bars interacting as contact. In this case, it is a node that "sees" a bar. $\vec{f}$ is the contact force and we decompose it into $\vec{n}$ and $\vec{t}$ the normal and the tangential component.}
  \label{fig:NodesBarsCon}
\end{figure}

In order to determine the forces and moments required for computing the accelerations, we need local computations able to calculate the forces or moments as a function of local parameters. 


These forces are those that are exerted between several physical objects defined in our model.

The first to describe is : The elastic repulsion (normal force) between object $i$ and $j$, $f^{e, i \rightarrow j}_n$ is calculated when the contactbox/control is activated. $f^{e, i \rightarrow j}_n$ is defined by : 

\begin{equation}
\begin{split}
    f^{e, i \rightarrow j} &= \mathcal{E}(v_i,v_j)
\end{split}
\end{equation}

where $\mathcal{E}$ is a function depending of the velocity and position of the objects $i$ and $j$. $v_i$ and $v_j$ are respectively the velocity of object $i$ and $j$.

In this model we also have a glue force between the cells (like in Fig \ref{fig:WhatACell}). Keeping the same notations, we will describe the glue force as $f^{g, i \rightarrow j}$ :

\begin{equation}
\begin{split}
    f^{g, i \rightarrow j} &= \mathcal{G}(v_i,v_j,k_n,k_t)
\end{split}
\end{equation}

where $\mathcal{G}$ is a function depending of the velocity and position of the objects $i$ and $j$ and two coefficient of cohesion $k_n$ and $k_t$. $v_i$ and $v_j$ are respectively the velocity of object $i$ and $j$.



Internal forces and moments
---------------------------


These are the forces that are exerted on the nodes directly without the need for redistribution to the nodes to obtain the accelerations. 

We use a global viscosity to dissipate energy physically. Let $N$ be a node, the viscosity on the node N is described by $f^{v}_N$ :
\begin{equation}
\begin{split}
    f^{v}_N &= -\nu \times Vel(N)     
\end{split}
\end{equation}
where Vel is the application that gives us the speed of a node. $\nu$ is a dissipation parameter that we fix at the beginning of the simulation.

The last step is to simulate the plasticity of the bars we transmit. The moments are directly managed in the Node class as shown in the Fig \ref{fig:NodesBarsOr}. Let $\Theta$ be a function depending of the velocity of the node N :
\begin{equation}
\begin{split}
    Mom_{t+1}(N) &= Mom_t(N) -k_r \times \Theta(Vel(N))
\end{split}
\end{equation}

where $Mom_t(N)$ is the function that determine the moment at the node $N$ by the incremental way. And for the plasticity, we apply a threshold like on Figure \ref{fig:Plasticity} : 

\begin{figure}[htbp]
  \centering
  \includegraphics[scale = 0.05]{Plasticity.png}
  \caption{Fonction de seuil pour la plasticité}
  \label{fig:Plasticity}
\end{figure}


\begin{remark}
It is interesting to think about the possibility of having cells that are not empty. In our case, they will be considered empty nevertheless, we could affect an internal pressure to our cells. 
\end{remark}

It is possible to impose a speed or a force on each node independently or by cell directly. 

\begin{figure}[htbp]
  \centering
  \includegraphics[height = 450pt, width = 170pt]{UML_L-hyphen.pdf}
  \caption{UML Diagram of the algorithm. The physical objects (nodes, bars, cells) are data boxes that evolve over time. The heart of L-hyphen lies in the sample and especially in the methods contactNodeForces and contactInteractionForces}
  \label{fig:UML}
\end{figure}








