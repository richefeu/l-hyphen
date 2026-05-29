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
  - #command("t <value>") ~Current time.
  - #command("dt <value>") ~Time-step increment.
  - #command("nstep <value>") ~Total number of time-steps.
  - #command("nstepPeriodSVG <value>") ~Number of time-steps between SVG dumps.
  - #command("nstepPeriodRecord <value>") ~Number of time-steps between records in some files.
  - #command("nstepPeriodConf <value>") ~Number of time-steps between conf-file dumps.
  - #command("isvg <value>") ~Current SVG file id-number.
  - #command("iconf <value>") ~Current conf-file id-number.
])

= Global parameters

#concept-block(body: [
  - #command("nbThreads <value>") ~Number of OpenMP threads to be used.
  - #command("gravity <gx> <gy>") ~Gravity vector components.
  - #command("numericalDissipation <value>") ~A purely numerical dissipation. It consists to multiply the velocities by (`1 - value`).
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
  - #command("fadh <value>") ~Adhesion force (same value for all contacts). This adhesion force can act only for non-glued interactions.
])

= Glue parameters

#concept-block(body: [
  To set a glue parameter, first glue the adjacent cell-walls using this command:

  - #command("glue <distance_max>") ~For force-based rupture model.
  - #command("GcGlue <distance_max>") ~For energy-based rupture model.

  Whatever the rupture model, glue parameters are local to each interaction and not regularly refreshed. Once broken, they cannot be restored.

  - #command("setGlueSameProperties <kn_coh> <kt_coh> <fn_coh_max> <ft_coh_max> <power>") ~Sets the normal and tangential stiffnesses, thresholds, and the power used in the yield function.
  - #command("setGcGlueSameProperties <kn_coh> <kt_coh> <Gc>") ~Sets the normal and tangential stiffnesses, and the surface energy.
])

= Internal pressure

#concept-block(body: [
  - #command("cellContent <value>") ~Select a model for core-pressure.
  - #command("compressFactor <value>") ~The elastic stiffness that links volume change to internal pressure (p = K × Ω / Ω₀).
])

= cell and cell-wall parameters

#concept-block(body: [
  - #command("setCellMasses <cellMass>") ~Distribute mass equally among all nodes of each cell. Each node gets mass/n_nodes.
  - #command("setNodeMasses <nodeMass>") ~Set the same mass to all nodes in the system.
  - #command("setCellWallDensities <rho> <thickness>") ~Set cell-wall masses based on density and wall thickness. Mass is distributed at the nodes.
  - #command("setCellDensities <rho> <thickness>") ~Set masses for both cell-wall and interior. Combines wall density with interior volume distribution.
  - #command("setCellWallDampingRates <alpha_s> <alpha_b>") ~Set damping coefficients for stretching and bending.
])

= Cells

#concept-block(body: [
  - #command("cells <number>") ~Indicates the cell section of the given number of cells. Then for each cell:
    - #command("<radius> <nbNodes> <nbBars> <pressure> <surface> <surface0> <close>"), and for each node of the cell:
     - #command("<mass> <xpos> <ypos> <xvel> <yvel> <xforce> <yforce> <ictrl> <prevNode> <next> <kr> <mz> <mz_max>")

  `ictrl` is the control id-number, but for free node it is `x`. `prevNode` and `nextNode` are id-numbers of the previous and next node, respectively, in the cell. When a cell is not closed (it does not form a loop), a cell-wall extremity is indicated with `x` for `prevNode` or `nextNode`.
])

= Neigbors

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
  - #command("linkCells <lx> <ly>") ~Use link-cells algorithm for neighbor search (lx, ly are cell sizes). If zero, brute-force O(N²) search is used.
  - #command("distVerlet <value>") ~The Verlet skin distance. Two cells are neighbors if distance < contact radius + distVerlet. Increases list validity range.
  - #command("nstepPeriodVerlet <value>") ~Number of time-steps between neighbor list updates. Larger values = fewer updates but must stay within Verlet distance.
  - #command("updateNeighbors") ~Manually update the neighbor list (automatically called at each Verlet period).
])

= Pre-processing

#concept-block(body: [
  - #command("addMultiLine <xo> <yo> <xe> <ye> <nbSegs> <barWidth> <Kn> <Kr> <Mz_max> <p_int>") ~Create a line with multiple segments from (xo,yo) to (xe,ye).
  - #command("addRegularPolygonalCell <nbFaces> <x> <y> <Rext> <barWidth> <rot> <Kn> <Kr> <Mz_max> <p_int>") ~Add a regular polygonal cell (triangle, square, hexagon...).
  - #command("addSquareBrickWallCells <nx> <ny> <horizDist> <xleft> <ybottom> <barWidth> <Kn> <Kr> <Mz_max> <p_int>") ~Create a brick-wall of square cells on triangular grid.
  - #command("addHoneycombCells <nx> <ny> <cellExtWidth> <xleft> <ybottom> <barWidth> <Kn> <Kr> <Mz_max> <p_int>") ~Create a honeycomb structure.
])

= Cells from a node-file

#concept-block(body: [
  - #command("readNodeFile <fileName> <barWidth> <Kn> <Kr> <Mz_max> <p_int>") ~Read cell nodes from a file. If barWidth < 0, it's auto-computed as half the min distance between different cells.

The node-file is a list of `<x> <y> <id>` values, with consecutive `id`-values belong to the same closed cell. The `id`-values are *not* the cell-id numbers.
])

= Node controls

#concept-block(body: [
  Control modes: 0 = FORCE_CONTROL, 1 = VELOCITY_CONTROL

  - #command("setNodeControl <cellId> <nodeId> <xmode> <xvalue> <ymode> <yvalue>") ~Apply force or velocity control to a single node.
  - #command("setCellControl <cellId> <xmode> <xvalue> <ymode> <yvalue>") ~Apply force or velocity control to all nodes of a cell.
  - #command("setNodeControlInBox <xmin> <xmax> <ymin> <ymax> <xmode> <xvalue> <ymode> <yvalue>") ~Apply control to all nodes within a rectangular region.
  - #command("setCellInternalPressure <cellId> <pressure>") ~Set the internal pressure of a cell.
])

= Diagnostics and output

#concept-block(body: [
  - #command("diagnostics") ~Print a diagnostic report of simulation parameters: geometry, time-step stability, Verlet parameters, search algorithm. Report is saved to `diagnostic.txt`.
  - #command("head") ~Display the L-HYPHEN ASCII art header.
  - #command("addRegularPolygonalCellsOnTriangularGrid <nx> <ny> <hDist> <vDist>") ~Add polygonal cells on a triangular grid.
])

= Cell content models

#concept-block(body: [
  - `CELL_EMPTY` ~Empty cell (no pressure model).
  - `CELL_ELASTIC_PV` ~Elastic volume with pressure-volume relation: $p = K_Omega Omega / Omega_0$.
  - `CELL_RIGID` ~Rigid cell (constant volume, not deformable).
  - `CELL_FLUID` ~Fluid pressure model.
])
