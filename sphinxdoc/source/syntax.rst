
Format of input files
=====================

The input files hold the complete simulation configuration and parameters.
They are text files that define the initial geometry, material properties, boundary conditions,
and output settings for a simulation.

**File Format:**

- Lines starting with ``#``, ``/``, or ``!`` are comments
- Numeric values support arithmetic expressions: ``1.5 * 2 + 0.5`` → ``3.5``
- Constants can be defined with ``define`` and reused throughout the file
- Order of commands generally doesn't matter (except for geometric setup)

Timing and Integration
======================

These commands control the time stepping and integration scheme.

``t`` — Initial time
  | **Type:** double
  | **Default:** 0.0
  | **Description:** Starting time for the simulation. For new simulations, set to 0.0.
  | **Example:** ``t 0.0``

``dt`` — Time step
  | **Type:** double
  | **Units:** seconds
  | **Critical:** Must satisfy stability condition: ``dt < dt_crit`` (check ``diagnostic.txt``)
  | **Typical range:** 1e-6 to 1e-4 (depends on system stiffness and masses)
  | **Description:** Increment of time for each integration step. Smaller dt = more accurate but slower.
  | **Warning:** Too large dt causes instability and divergence!
  | **Example:** ``dt 1e-5``

``nstep`` — Number of time steps
  | **Type:** integer
  | **Description:** Total number of integration steps to perform.
  | **Total simulation time:** ``final_time = t + nstep * dt``
  | **Example:** ``nstep 100000`` (with dt=1e-5 → total time = 1 second)

``cyclicVelPeriod`` — Loading cycle period (optional)
  | **Type:** double
  | **Units:** seconds
  | **Description:** For cyclic loading: reverses loading direction every N seconds.
  | **Default:** Not used if omitted
  | **Example:** ``cyclicVelPeriod 2.0`` (loading reverses every 2 seconds)

Dissipation and Damping
=======================

Energy dissipation is crucial for stable simulations. Multiple damping mechanisms are available:

``numericalDissipation`` — Artificial velocity damping
  | **Type:** double (0.0 to 1.0)
  | **Default:** 0.0 (no numerical damping)
  | **Description:** Multiplies velocities by (1 - value) at each step. Pure numerical artifact for stability.
  | **Typical range:** 1e-5 to 1e-3 (small values to avoid unphysical behavior)
  | **Warning:** High values lead to unrealistic energy loss and stiffness
  | **Example:** ``numericalDissipation 1e-4``

``globalViscosity`` — Viscous damping
  | **Type:** double
  | **Default:** 0.0
  | **Description:** Adds viscous damping proportional to velocity: ``f_visc = -globalViscosity * v``
  | **Physical interpretation:** Models resistance like moving through a fluid
  | **Typical range:** 0.1 to 10.0 (depends on system characteristics)
  | **Example:** ``globalViscosity 2.0``

``setCellWallDampingRates`` — Proportional damping
  | **Type:** double double
  | **Parameters:**
  |   - alpha_s: damping rate for axial (stretching) modes
  |   - alpha_b: damping rate for bending modes
  | **Formula:** ``c = 2 * alpha * sqrt(m * k)`` (critical damping scaling)
  | **Typical range:** 0.01 to 0.1
  | **Description:** Sets damping based on mass and stiffness (more physical than global viscosity)
  | **Example:** ``setCellWallDampingRates 0.05 0.05``

``setCellWallDampings`` — Direct damping coefficients
  | **Type:** double double
  | **Parameters:**
  |   - nu_s: direct damping coefficient for stretching
  |   - nu_b: direct damping coefficient for bending
  | **Description:** Set damping coefficients directly (no scaling with mass/stiffness)
  | **Example:** ``setCellWallDampings 100 50``

**Comparison:** Use ``setCellWallDampingRates`` for a physically meaningful damping that scales with system properties.

Volume Forces
=============

``gravity`` — Gravitational acceleration
  | **Type:** double double (x-component, y-component)
  | **Units:** m/s²
  | **Default:** 0.0 0.0 (no gravity)
  | **Description:** Applies volumetric force to all nodes as ``F = m * g``
  | **Example:** ``gravity 0.0 -9.81`` (Earth gravity, downward) 


Neighbor Detection (Verlet Algorithm)
=====================================

The Verlet algorithm optimizes contact detection by maintaining an expanded neighbor list.

``distVerlet`` — Verlet skin distance
  | **Type:** double
  | **Units:** length units (same as geometry)
  | **Default:** 0.01 to 0.1 (10-100% of typical cell size)
  | **Description:** Thickness of the detection "skin" around cells. Two cells are neighbors if:
  |   ``distance < radius_i + radius_j + distVerlet``
  | **Physical interpretation:** Allows the list to remain valid between updates
  | **Trade-off:** Larger → fewer list updates but more contact checks; smaller → more updates
  | **Typical range:** 0.01 to 0.5
  | **Example:** ``distVerlet 0.05``

``nstepPeriodVerlet`` — Neighbor list update frequency
  | **Type:** integer
  | **Default:** 1 (update every step)
  | **Description:** Number of time-steps between neighbor list rebuilds.
  | **Constraint:** Must be chosen so max displacement < distVerlet between updates
  |   ``max_displacement = v_max * (nstepPeriodVerlet * dt) < distVerlet``
  | **Optimization:** Larger values = faster (fewer updates) but must stay within Verlet distance
  | **Typical range:** 1 to 1000 (depends on velocities)
  | **Warning:** Too large → missed contacts; too small → wasted computation
  | **Example:** ``nstepPeriodVerlet 100`` (update every 100 steps)


Outputs
-------

- ``isvg`` (*integer*) **value**
  Starting index for SVG file numbering (e.g., 0 generates sample0000.svg, sample0001.svg, ...)

- ``nstepPeriodSVG`` (*integer*) **value**
  Number of time-steps between SVG output saves. Set to 0 to disable SVG generation.

- ``findDisplayArea`` (*double*) **value**
  Find the limits of the display area in svg files. **value** is a scale factor (≥ 1.0).
  The displayed area is scaled by this factor from the system bounding box.

- ``iconf`` (*integer*) **value**
  Starting index for configuration file numbering (e.g., 0 generates conf0, conf1, ...)

- ``nstepPeriodConf`` (*integer*) **value**
  Number of time-steps between configuration file saves. Set to 0 to disable configuration saves.

- ``nstepPeriodRecord`` (*integer*) **value**
  Number of time-steps between scalar data recordings (if enabled).

- ``limits`` (*double*) **xmin** (*double*) **xmax** (*double*) **ymin** (*double*) **ymax**
  Manually set the display area limits (alternative to ``findDisplayArea``).


Contact Interactions
====================

These parameters control the forces between cells when they touch.

``kn`` — Normal contact stiffness
  | **Type:** double
  | **Units:** Force/length (N/m in SI)
  | **Stability impact:** CRITICAL - larger values require smaller dt
  | **Description:** Penalty stiffness for normal overlap: ``f_n = -kn * overlap``
  | **Effect:** Controls how "hard" cells are; higher kn → stiffer, less penetration
  | **Typical range:** 1e3 to 1e8 (depends on cell size and material)
  | **Rule of thumb:** Choose dt such that ``dt_crit / dt > 20`` (check diagnostic)
  | **Example:** ``kn 1000000``

``kt`` — Tangential contact stiffness
  | **Type:** double
  | **Units:** Force/length (N/m in SI)
  | **Description:** Penalty stiffness for tangential slip: ``f_t = -kt * slip_velocity * dt``
  | **Typical range:** 0.1 to 1.0 × kn (usually kt < kn)
  | **Example:** ``kt 800000``

``mu`` — Coulomb friction coefficient
  | **Type:** double (≥ 0)
  | **Default:** 0.0 (frictionless)
  | **Description:** Friction limiting tangential force: ``f_t_max = mu * |f_n|``
  | **Physical meaning:** Higher μ → more friction (stickier contacts)
  | **Typical range:** 0.0 to 1.0 (0.3 for many materials)
  | **Example:** ``mu 0.3``

``fadh`` — Normal adhesion force
  | **Type:** double (≥ 0)
  | **Default:** 0.0 (no adhesion)
  | **Units:** Force (N in SI)
  | **Description:** Attractive force when cells just separate: ``f_adhesion = -fadh``
  | **Physical interpretation:** Models weak surface attraction (van der Waals, etc.)
  | **Warning:** Only active when cells are NOT glued (see ``glue`` command)
  | **Typical range:** 0.0 to 1000
  | **Example:** ``fadh 10.0``

``adaptativeStiffness`` — Contact stiffness scaling
  | **Type:** integer (0 or 1)
  | **Default:** 0 (disabled)
  | **Description:** If enabled (1), scales contact stiffness by overlap:
  |   ``k_effective = kn * D / (D + overlap)``
  | **Effect:** Prevents cell-wall penetration, keeps forces reasonable
  | **Recommendation:** Enable (set to 1) for more realistic contact
  | **Example:** ``adaptativeStiffness 1``

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
  Defines the same control for all nodes within a rectangular region. Useful for boundary conditions.

Cohesion and Glue Bonds
======================

Glue bonds create permanent adhesive links between cell walls. Two rupture models are available:

**Force-based model:** Rupture when force threshold is exceeded (faster, simpler)
**Energy-based model (Gc):** Rupture when energy release rate exceeds critical value (more physical)

``glue`` — Create force-based glue bonds
  | **Type:** double
  | **Units:** length (same as geometry)
  | **Description:** Creates adhesive bonds between all cell pairs within the given distance.
  |   Two cells are glued if their distance < specified value.
  | **Must be followed by:** ``setGlueSameProperties``
  | **Effect:** Creates permanent bonds until rupture criterion is met
  | **Typical distance:** 0.005 to 0.05 (small fraction of cell size)
  | **Example:** ``glue 0.01``

``setGlueSameProperties`` — Configure force-based rupture criteria
  | **Type:** five doubles
  | **Parameters:**
  |   - **kn_coh** : normal stiffness of glue (N/m) — typically 1e4 to 1e6
  |   - **kt_coh** : tangential stiffness of glue (N/m) — typically 0.5 to 1.0 × kn_coh
  |   - **fn_coh_max** : maximum normal force before rupture (N) — typically 100 to 10000
  |   - **ft_coh_max** : maximum tangential force before rupture (N) — same order as fn_coh_max
  |   - **yieldPower** : exponent in rupture criterion (typically 1-5)
  | **Rupture criterion:** Bond breaks if:
  |   ``-fn_coh/fn_coh_max + (|ft_coh|/ft_coh_max)^yieldPower > 1``
  | **Interpretation:** Combines normal and tangential loading; higher yieldPower = more sensitive to shear
  | **Usage:**
  |   .. code-block:: text
  |
  |      glue 0.01
  |      setGlueSameProperties 100000 100000 2000 2000 2.0

``GcGlue`` — Create energy-based glue bonds
  | **Type:** double
  | **Units:** length (same as geometry)
  | **Description:** Creates adhesive bonds using fracture mechanics. Bonds break when strain energy release rate exceeds Gc.
  | **Advantage:** More physically founded; accounts for progressive damage
  | **Must be followed by:** ``setGcGlueSameProperties``
  | **Typical distance:** 0.005 to 0.05 (similar to force-based)
  | **Example:** ``GcGlue 0.01``

``setGcGlueSameProperties`` — Configure energy-based rupture criteria
  | **Type:** three doubles
  | **Parameters:**
  |   - **kn_coh** : normal stiffness of glue (N/m)
  |   - **kt_coh** : tangential stiffness of glue (N/m)
  |   - **Gc** : critical fracture energy per unit length (J/m)
  | **Physical basis:** Bond ruptures when elastic energy density > Gc
  | **Typical Gc values:** 0.1 to 100 (J/m) depending on material and interface
  | **Usage:**
  |   .. code-block:: text
  |
  |      GcGlue 0.01
  |      setGcGlueSameProperties 100000 100000 0.5

**When to use which model:**
  | **Force-based:** When you know the force threshold from experiments
  | **Energy-based:** When you have fracture energy data or want more physical behavior

Cell and Structural Properties
=============================

Mass Distribution
-----------------

``setCellMasses`` — Distribute mass per cell
  | **Type:** double
  | **Units:** kg (or system mass units)
  | **Description:** Gives each cell the specified total mass, distributed equally among its nodes.
  |   Each node gets: ``m_node = cell_mass / number_of_nodes``
  | **Use case:** When cell mass is known
  | **Example:** ``setCellMasses 0.1`` (each cell = 100g)

``setNodeMasses`` — Uniform node mass
  | **Type:** double
  | **Units:** kg
  | **Description:** All nodes get identical mass (simplest approach)
  | **Use case:** Initial tests, uniform materials
  | **Example:** ``setNodeMasses 0.01`` (each node = 10g)

``setCellWallDensities`` — Mass from geometry
  | **Type:** double double
  | **Parameters:** rho (kg/m³), thickness (m)
  | **Description:** Calculates mass from cell-wall geometry only.
  |   ``m_node = rho * bar_length * 2*radius * thickness / 2`` (half to node)
  | **Physical interpretation:** Models actual structural material mass
  | **Example:** ``setCellWallDensities 1400 0.001`` (tissue: 1.4 g/cm³, 1mm thick)

``setCellDensities`` — Total mass (structure + content)
  | **Type:** double double
  | **Parameters:** rho (kg/m³), thickness (m)
  | **Description:** Includes both cell-wall and interior mass distribution.
  | **Use case:** Realistic cells with internal fluid/content
  | **Example:** ``setCellDensities 1400 0.002``

Internal Pressure Models
------------------------

``cellContent`` — Pressure model type
  | **Type:** integer (0, 1, or 2)
  | **Options:**
  |   - ``0`` : CELL_EMPTY (no pressure, just structural)
  |   - ``1`` : CELL_ELASTIC_PV (elastic compression: p = K × Ω / Ω₀)
  |   - ``2`` : CELL_RIGID (constant volume, incompressible)
  | **Default:** 0 (empty)
  | **Example:** ``cellContent 1`` (elastic pressure model)

``compressFactor`` — Elastic pressure stiffness
  | **Type:** double
  | **Units:** Pressure × Volume / Pressure = Volume⁻¹ × Force
  | **Description:** Stiffness constant K for pressure: ``p = K * (Omega / Omega0)``
  | **Only used when:** cellContent = 1
  | **Typical range:** 1e4 to 1e8 (depends on desired pressure response)
  | **Example:** ``compressFactor 1000000`` (stiff elastic response)

``setCellInternalPressure`` — Initial pressure
  | **Type:** integer double
  | **Parameters:** cell_id (cell number), pressure (Pa or pressure units)
  | **Description:** Sets initial overpressure in a specific cell.
  | **Use case:** Pre-stressed cells (like turgid plant cells)
  | **Example:** ``setCellInternalPressure 0 100`` (100 Pa overpressure in cell 0)

Neighbor Detection
------------------

- ``linkCells`` (*double*) **lx** (*double*) **ly**
  Enable link-cells spatial hashing for faster neighbor detection (O(N) instead of O(N²)).
  **lx**, **ly** are the dimensions of the spatial cells. Recommended for large systems (N > 500).

- ``updateNeighbors``
  Manually trigger an update of the neighbor list (normally called automatically).

General Options
---------------

- ``nbThreads`` (*integer*) **n**
  Number of OpenMP threads for parallel execution (if compiled with OpenMP support).

- ``reorder`` (*integer*) **flag**
  Enable/disable (1/0) node reordering when reading configuration files.

- ``define`` (*string*) **NAME** (*double*) **value**
  Define a named constant for reuse in subsequent numeric expressions.

