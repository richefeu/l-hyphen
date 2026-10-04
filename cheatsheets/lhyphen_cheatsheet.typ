#import "@preview/boxed-sheet:0.1.0": *

#show: cheatsheet.with(
  title: [Configuration file for `lhyphen` (conf-files)],
  homepage: "https://github.com/richefeu/lhyphen",
  authors: "",
  write-title: false,
  title-align: left,
  title-number: true,
  title-delta: 2pt,
  scaling-size: true,
  font-size: 6.5pt,
  line-skip: 5.5pt,
  x-margin: 10pt,
  y-margin: 30pt,
  num-columns: 4,
  column-gutter: 9pt,
  numbered-units: false
)

#let command(body, fill: luma(90%)) = {
  set text(black, font:"Courier New", weight:"semibold")
  box(
    fill: fill,
    outset: 2pt,
    radius: 2pt,
    [#body]
  )
}

= Timing and simulation flow

#concept-block(body: [
  - #command("define <NAME> <value>") ~Define a named constant, usable in expressions written between two `$` (e.g. `dt $ T / 100 $`).
  - #command("t <value>") ~Current time.
  - #command("dt <value>") ~Time-step increment.
  - #command("nstep <value>") ~Total number of time-steps.
  - #command("cyclicVelPeriod <value>") ~Cyclic loading: imposed velocities change sign during the second half of each period (0 = off).
  - #command("nstepPeriodSVG <value>") ~Number of time-steps between SVG dumps.
  - #command("nstepPeriodCapture <value>") ~Number of time-steps between writes of the captureNodes files (0 = with each SVG dump).
  - #command("nstepPeriodRecord <value>") ~Number of time-steps between records in some files.
  - #command("nstepPeriodConf <value>") ~Number of time-steps between conf-file dumps.
  - #command("isvg <value>") ~Current SVG file id-number.
  - #command("iconf <value>") ~Current conf-file id-number.
])

= Global parameters

#concept-block(body: [
  - #command("nbThreads <value>") ~Number of OpenMP threads to be used.
  - #command("gravity <gx> <gy>") ~Gravity vector components.
  - #command("numericalDissipation <value>") ~A purely numerical dissipation. It consists to multiply the velocities by (`1 - value`) at each step.
  - #command("globalViscosity <value>") ~A dissipation that acts like a viscous fluid on the nodes (but this is not very physically sound).
  - #command("limits <xmin> <xmax> <ymin> <ymax>") ~The limits of the system (for display purpose).
  - #command("findDisplayArea <value>") ~Compute the limits with multiplying size-factor.
])

= Interaction parameters

#concept-block(body: [
  - #command("kn <value>") ~Normal contact stiffness (same value for all contacts).
  - #command("kt <value>") ~Tangential contact stiffness (same value for all contacts).
  - #command("adaptativeStiffness <0|1>") ~Make kn depend on overlap to avoid cell-wall penetration. Stiffness is multiplied by D / (D + d_n).
  - #command("mu <value>") ~Coulomb friction coefficient (same value for all contacts).
  - #command("viscnrate <value>") ~Normal contact viscosity, as a fraction of the critical damping: $c_n = "viscnrate" dot 2 sqrt(m_"eff" k_n)$.
  - #command("fadh <value>") ~Adhesion force (same value for all contacts). This adhesion force can act only for non-glued interactions.
])

= Glue parameters

#concept-block(body: [
  To set a glue parameter, first glue the adjacent cell-walls using this command:

  - #command("glue <distance_max>") ~For force-based rupture model.
  - #command("GcGlue <distance_max>") ~For energy-based rupture model (`distGcGlue` is a synonym).

  Whatever the rupture model, glue parameters are local to each interaction and not regularly refreshed. Once broken, they cannot be restored.

  - #command("setGlueSameProperties <kn_coh> <kt_coh> <fn_coh_max> <ft_coh_max> <power>") ~Sets the normal and tangential stiffnesses, thresholds, and the power used in the yield function.
  - #command("setGcGlueSameProperties <kn_coh> <kt_coh> <Gc>") ~Sets the normal and tangential stiffnesses, and the surface energy.
])

= Internal pressure

#concept-block(body: [
  - #command("cellContent <0|1|2>") ~Select a model for core-pressure (see _Cell content models_).
  - #command("compressFactor <value>") ~The elastic stiffness $K$ that links volume change to internal pressure: $p = -K (Omega - Omega_0) / Omega_0$.
  - #command("setCellInternalPressure <cellId> <pressure>") ~Set the internal pressure of a cell.
  - #command("setCellAsOpen <cellId>") ~Mark a cell as open (no closing bar). #command("setClose <cellId>") marks it as closed.
])

= cell and cell-wall parameters

#concept-block(body: [
  - #command("setCellMasses <cellMass>") ~Distribute mass equally among all nodes of each cell. Each node gets mass/n_nodes.
  - #command("setNodeMasses <nodeMass>") ~Set the same mass to all nodes in the system.
  - #command("setCellWallDensities <rho> <thickness>") ~Set cell-wall masses based on density and wall thickness. Mass is distributed at the nodes.
  - #command("setCellDensities <rho> <thickness>") ~Set masses for both cell-wall and interior. Combines wall density with interior volume distribution.
  - #command("setCellWallDampingRates <alpha_s> <alpha_b>") ~Set damping rates for stretching and bending.
  - #command("setCellWallDampings <nu_s> <nu_b>") ~Set the damping coefficients for stretching and bending directly.
])

= Cells

#concept-block(body: [
  - #command("cells <number>") ~Indicates the cell section of the given number of cells. Then for each cell:
    - #command("<radius> <nbNodes> <nbBars> <pressure> <surface> <surface0> <close>"), and for each node of the cell:
     - #command("<mass> <xpos> <ypos> <xvel> <yvel> <xforce> <yforce> <ictrl> <prevNode> <next> <kr> <mz> <mz_max>")

  `ictrl` is the control id-number, but for free node it is `x`. `prevNode` and `nextNode` are id-numbers of the previous and next node, respectively, in the cell. When a cell is not closed (it does not form a loop), a cell-wall extremity is indicated with `x` for `prevNode` or `nextNode`.
])

= Neighbors

#concept-block(body: [
  - #command("neighbors <number>") ~Indicates the neighbor section of the given number of neighbors. Then, for each neighbor:
   - #command("<icell> <jcell> <inode> <jnode> <nx> <ny> <contactState> <fn> <ft> <glueState>")
   Then, if `<glueState>` is `1`:
   - #command("<fn_coh> <ft_coh> <kn_coh> <kt_coh> <fn_coh_max> <ft_coh_max> <yieldPower>")
   and, if `<glueState>` is `2`:
   - #command("<fn_coh> <ft_coh> <kn_coh> <kt_coh> <Gc>")
])

= Neighbor-List of each cell

#concept-block(body: [
  - #command("linkCells") ~Use link-cells algorithm for neighbor search, O(N); the cell size is automatic (values left on the line by old files are ignored). Without it, brute-force O(N²) search is used.
  - #command("checkNeighbors") ~At start, check the link-cells neighbor list against brute-force search.
  - #command("distVerlet <value>") ~The Verlet skin distance. Two cells are neighbors if distance < contact radius + distVerlet. Increases list validity range.
  - #command("nstepPeriodVerlet <value>") ~Number of time-steps between neighbor list updates. Larger values = fewer updates but must stay within Verlet distance.
])

= Pre-processing

#concept-block(body: [
  - #command("addMultiLine <xo> <yo> <xe> <ye> <barWidth> <nbSegs> <Kn> <Kr> <Mz_max>") ~Create an open line with multiple segments from (xo,yo) to (xe,ye).
  - #command("addRegularPolygonalCell <nbFaces> <x> <y> <rot> <Rext> <barWidth> <Kn> <Kr> <Mz_max>") ~Add a regular polygonal cell (triangle, square, hexagon...).
  - #command("addSquareBrickWallCells <nx> <ny> <horizDist> <xleft> <ybottom> <barWidth> <Kn> <Kr> <Mz_max>") ~Create a brick-wall of square cells.
])

= Cells from a node-file

#concept-block(body: [
  - #command("readNodeFile <fileName> <barWidth> <Kn> <Kr> <Mz_max> <p_int>") ~Read cell nodes from a file. If barWidth < 0, it's auto-computed as half the min distance between different cells.

The node-file is a list of `<x> <y> <id>` values, with consecutive `id`-values belong to the same closed cell. The `id`-values are *not* the cell-id numbers.

  - #command("reorder <0|1>") ~Reorder the nodes when reading the node-file (default 1).
  - #command("cleanShortBars <ratio>") ~Merge bars shorter than ratio × mean bar length (e.g. 0.3), keeping neighbouring cells consistent. Use right after `readNodeFile`, before masses, dampings, controls and glue.
  - #command("momentForceMax <Fmax>") ~Only when a force transmitting a nodal moment (m / l) would exceed Fmax (short lever arm), the transmitted moment is reduced to ±Fmax · min(l_prev, l_next). Reversible (the elastic state mz is untouched), no spurious torque. Combine with alpha_b > 0. 0 = no limit.
])

= Node controls

#concept-block(body: [
  Control modes: 0 = VELOCITY_CONTROL, 1 = FORCE_CONTROL

  - #command("setNodeControl <cellId> <nodeId> <xmode> <xvalue> <ymode> <yvalue>") ~Apply force or velocity control to a single node.
  - #command("setCellControl <cellId> <xmode> <xvalue> <ymode> <yvalue>") ~Apply force or velocity control to all nodes of a cell.
  - #command("setNodeControlInBox <xmin> <xmax> <ymin> <ymax> <xmode> <xvalue> <ymode> <yvalue>") ~Apply control to all nodes within a rectangular region. The regions are saved in conf-files as a #command("controlBoxAreas <number>") section.
])

= Diagnostics and output

#concept-block(body: [
  - #command("captureNodes <file> <xmin> <xmax> <ymin> <ymax>") ~Record, in `file`, data of the nodes lying in the region when this line is read.
  - #command("followCell <cellId>") ~Follow a given cell (tracking data).

  A diagnostic report (geometry, time-step stability, Verlet parameters, search algorithm) is printed at each run and saved to `diagnostic.txt`.
])

= Events

#concept-block(body: [
  Checked at the start of each time step; each event fires its action once. Pending events are written in conf-files.

  - #command("event saveConfAtTime <t>") ~Save a conf-file when time reaches `t`.
  - #command("event saveConfAtBrokenLength <L>") ~Save a conf-file when the cumulated broken length exceeds `L` (`0` = first breakage).
  - #command("event stopAtBrokenLength <L>") ~Save a conf-file and stop the simulation when the cumulated broken length exceeds `L`.
  - #command("event stopAfterStressDrop <ictrl> <x|y> <drop%> <delay> <tStart> <tau>") ~Save a conf-file and stop the simulation `delay` after the reaction force on the nodes driven by control `ictrl` (0-based, in definition order) has dropped by `drop%` from its peak. The force is smoothed (exponential moving average, time constant `tau`, 0 = none) and the peak is tracked from `tStart`.
])

= Cell content models

#concept-block(body: [
  - `0` (`CELL_EMPTY`) ~Empty cell (no pressure model).
  - `1` (`CELL_CONSTANT_PV`) ~Gas: constant $p Omega$.
  - `2` (`CELL_ELASTIC_PV`) ~Liquid: $p = -K (Omega - Omega_0) / Omega_0$, with $K$ = `compressFactor`.
])
