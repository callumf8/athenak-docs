Non-relativistic radiation transport
====================================

Introduction
------------

The ``<nr_radiation>`` module solves the **time-independent transfer equation on discrete
ordinates** by the method of **short characteristics**, in the **gray**, **LTE**,
**scattering-free** limit, and couples the result back to the gas as an operator-split
heating rate. The method follows `Davis, Stone & Jiang 2012
<https://ui.adsabs.harvard.edu/abs/2012ApJS..199....9D/abstract>`_.

The radiation field is quasi-static: the formal solution is recomputed once per cycle from
the current fluid state, with no time derivative of the intensity and hence no light-crossing
CFL condition.

.. note::

   This is a different module from :doc:`radiation`, and the two are not interchangeable.
   That one solves the general-relativistic, time-dependent, flux-conservative transport
   equation and requires ``<coord>/general_rel=true``. This one is valid only in flat
   Cartesian coordinates and aborts under GR, since the straight-ray short-characteristics
   solution does not hold on a curved background. They use separate input blocks, arrays and
   output variables; only one should appear in a given input file.

Governing equations
-------------------

Along each discrete ray direction :math:`\hat{n}`,

.. math::

   \hat{n} \cdot \nabla I(\hat{n}) = \chi \left[ S - I(\hat{n}) \right],

where :math:`\chi` is the extinction coefficient and :math:`S` the source function. With no
scattering :math:`\chi = \sigma_a`, and in LTE the source function is the Planck function of
the local gas temperature,

.. math::

   S = B = T^4, \qquad T = (\gamma - 1) e / \rho,

in code units. Those units are built on a reference density, temperature and length
:math:`\rho_0, T_0, L_0`, with :math:`\chi` measured per unit code length and the intensity
in units of the blackbody intensity :math:`B(T_0)`. It is that choice of intensity unit which
leaves the equation above free of physical constants and reduces :math:`B` to :math:`T^4`;
the constants reappear only in the heating rate (see **Coupling to the gas**). :math:`S`
therefore depends only on the gas state and is known before the solve begins.

Formal solution
~~~~~~~~~~~~~~~

A short characteristic runs from an upwind end point, through the cell centre, to a downwind
end point; superscripts :math:`-`, :math:`0` and :math:`+` denote values there, and
:math:`\Delta\tau^-`, :math:`\Delta\tau^+` the optical depths of the two halves. The
intensity at the centre is (Eq. 20 of Davis, Stone & Jiang 2012)

.. math::

   I^0 = e^{-\Delta\tau^-} I^- + \Psi^- S^- + \Psi^0 S^0 + \Psi^+ S^+ .

Representing :math:`S` as a parabola through the three points in optical depth and
integrating it against the attenuation gives the weights in closed form (`Kunasz & Auer 1988
<https://ui.adsabs.harvard.edu/abs/1988JQSRT..39...67K/abstract>`_): they depend on
:math:`\Delta\tau^-` and :math:`\Delta\tau^+` alone, not on :math:`S`, and sum to
:math:`1 - e^{-\Delta\tau^-}`. Both this parabola and the optical-depth integrals are
limited by the Bézier construction of `Auer 2003
<https://ui.adsabs.harvard.edu/abs/2003ASPC..288....3A/abstract>`_ and `Hayek et al. 2010
<https://ui.adsabs.harvard.edu/abs/2010A%26A...517A..49H/abstract>`_, which suppresses
overshoot at steep gradients and keeps :math:`I \ge 0` for :math:`S \ge 0`.

In multiple dimensions the end points do not fall on cell centres; :math:`\chi`, :math:`S`
and :math:`I` there are interpolated linearly (bilinear in 3D), making the scheme second
order along a ray and first order across it — the leading source of spatial error.

Angular quadrature and moments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Angles are the type-A set of `Bruls et al. 1999
<https://ui.adsabs.harvard.edu/abs/1999A%26A...348..233B/abstract>`_ (a Carlson symmetric
:math:`S_N` grid), built once at start-up and controlled by ``nmu`` alone. In 1D they are a
Gauss-Legendre set of :math:`2 n_\mu` directions; in 2D and 3D there are
:math:`n_\mu (n_\mu + 1)/2` rays per octant, so in 3D

.. math::

   n_{\mathrm{rays}} = 8, \; 24, \; 48, \; 80, \; 120, \; 168
   \quad \text{for } n_\mu = 1 \ldots 6,

and half as many in 2D. The moments are plain weighted sums over all rays,

.. math::
   :nowrap:

   \begin{align}
   J      &= \sum_a \omega_a \, I_a \\
   H_i    &= \sum_a \omega_a \, n_i \, I_a \\
   K_{ij} &= \sum_a \omega_a \, n_i \, n_j \, I_a
   \end{align}

with weights normalised so that :math:`\sum_a \omega_a = 1`, so :math:`J` is the mean
intensity and isotropic :math:`I` gives :math:`J = I`. The intensity, and hence each moment,
is measured in units of :math:`B(T_0) = a_R c T_0^4 / 4\pi`, with :math:`a_R = 4\sigma_B/c`
the radiation constant. From :math:`E_r = 4\pi J/c`,
:math:`F = 4\pi H` and :math:`P = 4\pi K/c` it follows that :math:`J` and :math:`K` are the
radiation energy density and pressure in units of :math:`a_R T_0^4`, and :math:`H` the flux
in units of :math:`c \, a_R T_0^4`. ``prat`` and ``crat`` convert these to the gas scales
(see **Inputs**).

Iteration
~~~~~~~~~

Since the medium is gray, in LTE and free of scattering, :math:`S = B` is known in advance
and set once per solve. There is consequently **no source-function iteration** — no
:math:`\Lambda`-iteration, no accelerated lambda iteration, no non-LTE update.

The loop that does exist carries intensity **across meshblock boundaries**. Within a
meshblock the cells are swept in upwind hyperplane order, so one sweep solves that block
exactly for the given ghost intensities; each sweep is followed by a boundary exchange, and
so advances the solution one meshblock further downstream along every ray. The sweep count
needed is therefore set by how many meshblocks a ray must cross, not by any convergence rate:
a single-meshblock problem converges in one sweep. The previous solve's :math:`J` is retained
as the starting reference, so a quasi-static field also exits after one sweep (two on the
first solve, where :math:`J` begins at zero).

Sweeping stops when :math:`J` stops changing,

.. math::

   \max_{\mathrm{cells}}
   \frac{\left| J^{(n)} - J^{(n-1)} \right|}{\max\!\left( J^{(n)}, J^{(n-1)} \right)}
   \;\le\; \epsilon_{\mathrm{tol}},

reduced over all MPI ranks. This tests the **change per sweep**, not the error: a faint
feature on a bright background can satisfy it before the information has crossed the mesh, so
``iter_tol`` should be tightened in proportion to the amplitude of interest. Exhausting
``iter_max`` warns but does not abort.

Inputs
------

All parameters below belong in the ``<nr_radiation>`` block.

**Angular grid**

- ``nmu`` (**mandatory**): number of polar levels of the quadrature, from ``1`` to ``6``.
  Larger values are rejected: the type-A construction determines the weights uniquely only up
  to six levels (:math:`S_{12}`, 168 rays; Sec. 3.2 of Davis, Stone & Jiang 2012).

**Opacity and emission**

- ``kappa_a`` (default ``0.0``): constant absorption opacity **per unit mass**, so
  :math:`\sigma_a = \kappa_a \rho`. Required when a fluid is present unless the problem
  generator enrolls its own opacity function; ignored when there is no fluid.

A problem generator may instead supply opacity and emission functions, called at the start of
every solve, which write :math:`\sigma_a` and :math:`B` cell by cell. That is how the bundled
tests impose a prescribed medium.

**Coupling to the gas**

The reference state fixes both of the constants below; neither is a free parameter. They
enter only through :math:`Q_{\mathrm{rad}}`, never the transfer solve, and there only as the
product :math:`c_{\mathrm{rat}} p_{\mathrm{rat}}`, which carries the radiation field from
its own unit :math:`B(T_0)` into the units of the gas energy equation.

- ``crat`` (**mandatory if** the radiation affects the fluid): :math:`c / a_0`, the speed of
  light over the isothermal sound speed :math:`a_0 = \sqrt{P_0 / \rho_0}` of the reference
  state. Its largeness is also what justifies dropping :math:`c^{-1} \partial_t I`.
- ``prat`` (**mandatory if** the radiation affects the fluid): :math:`a_R T_0^4 / P_0`, the
  radiation energy density at :math:`T_0` in units of the reference gas pressure
  :math:`P_0 = \rho_0 k_B T_0 / \mu m_H`.
- ``affect_fluid`` (default ``true`` when a fluid is present): whether the heating rate is
  added to the gas energy. Forced to ``false``, and not enableable, without a fluid.
- ``qrad_form`` (default ``hybrid``): how the heating rate is evaluated — ``integral``,
  ``divh`` or ``hybrid``; see **Coupling to the gas**.
- ``cfl_qrad`` (default ``0.25``, or ``0.0`` when ``qrad_form = integral``): fraction of a
  cell's internal energy the radiation may remove in one step. ``0`` disables the limit.

**Iteration**

These control how far intensity is carried across meshblock boundaries, not the convergence
of a source-function iteration, of which there is none in the present LTE solver.

- ``iter_max`` (default ``100``): maximum sweeps per solve. Must be large enough for a ray to
  cross the block decomposition.
- ``iter_min`` (default ``1``): minimum sweeps per solve.
- ``iter_tol`` (default ``1.0e-6``): threshold on the per-sweep relative change of :math:`J`.

**Sweep kernel**

These choose how the formal solution is swept through each meshblock in parallel. A ray
marches along the axis it crosses first (set by the cell shape, which may be anisotropic)
and reads the intensity only one step upwind along it, so a whole plane of cells is solved
at once. ``plane`` launches one plane at a time across the meshblock; ``tiled`` sweeps the
planes inside each tile and takes the tiles in diagonal wavefronts, as the footpoint
interpolation reaches into the transversely upwind tile. Both give bit-identical results:

.. image:: /_static/nr_radiation_sweep.png
   :alt: Cell planes swept across a whole meshblock, beside tiles taken in diagonal wavefronts with the same plane sweep inside each tile.
   :align: center
   :width: 85%

- ``sweep_kernel`` (default ``plane``): ``plane`` or ``tiled``; ``tiled`` is usually faster
  on GPUs when there are too few meshblocks to fill the device from one launch.
- ``tile_size`` (default ``0``): with ``tiled``, cells along a (cubic) tile edge; must divide
  the meshblock interior in every active dimension. ``0`` makes the meshblock one tile.
- ``team_size`` (default ``0``): with ``tiled``, threads per team, a multiple of 32; ``0``
  lets Kokkos choose, which is usually best.

**Boundary conditions**

The intensity accepts ``periodic``, ``outflow``, ``inflow``, ``reflect``, ``diode``,
``vacuum``, ``user`` and ``block``. Any other flag is rejected at start-up rather than
silently treated as a vacuum edge. The default ``inflow`` table is zero, i.e. a vacuum
boundary; a problem generator may overwrite it to inject radiation.

Outputs
-------

Request one of the following in an ``<output>/variable`` line (see :doc:`outputs` for output
file types):

- ``sc_rad``: all ten moments — ``rad_J``, ``rad_H1``, ``rad_H2``, ``rad_H3``, ``rad_K11``,
  ``rad_K22``, ``rad_K33``, ``rad_K12``, ``rad_K13``, ``rad_K23``. This does **not** include
  the source function or opacity, which must be requested by name.
- ``sc_rad_J``, ``sc_rad_H``, ``sc_rad_K``: the moments individually.
- ``sc_rad_src``: the source function, written as ``rad_S``.
- ``sc_rad_sigma``: the absorption coefficient, written as ``rad_sigma_a``.

Requesting any of these without an ``<nr_radiation>`` block is an error. Three points matter
when interpreting them:

- The moments are in units of :math:`B(T_0)` (see **Angular quadrature and moments**), so
  the dumped :math:`J`, :math:`H`, :math:`K` are :math:`E_r`, :math:`F`, :math:`P` in units
  of :math:`a_R T_0^4`, :math:`c \, a_R T_0^4` and :math:`a_R T_0^4`. In the gas units
  :math:`P_0` and :math:`P_0 a_0`: :math:`E_r = p_{\mathrm{rat}} J`,
  :math:`F = c_{\mathrm{rat}} p_{\mathrm{rat}} H`, :math:`P = p_{\mathrm{rat}} K`.
- In 2D the third direction cosine is unsigned, so ``rad_H3``, ``rad_K13`` and ``rad_K23``
  are not meaningful.
- **The** :math:`t = 0` **dump precedes the first solve**: ``rad_sigma_a`` is valid, but the
  moments and ``rad_S`` are zero. They are refreshed once per cycle, for the fluid state at
  the start of that cycle.

Neither the heating rate nor the specific intensity is available as an output variable. Two diagnostics help monitor a run. Every
``<time>/ndiag`` cycles the solver prints ``Iteration stops at niter: N relative error: R``,
where ``N`` is the number of sweeps that cycle's solve needed to carry the intensity across
the meshblocks (see **Iteration**) and ``R`` is the relative change of :math:`J` on the last
of them, which should be below ``iter_tol``; a warning is printed instead if ``iter_max`` was
exhausted. The ``log`` output adds a ``qrad_neg`` column: the number of cells per step in
which adding the heating rate would have left the internal energy non-positive. With the
depletion limit active (``cfl_qrad`` between 0 and 1) it stays at zero by construction, so
a non-zero count means that limit is off or too loose.

Coupling to the gas
-------------------

The coupling is operator split (Eq. 26 of Davis, Stone & Jiang 2012). One heating rate
:math:`Q_{\mathrm{rad}}` is built per cycle from the start-of-cycle fluid state, then held
fixed and added in each Runge-Kutta stage with that stage's weight,

.. math::

   u_{\mathrm{IEN}} \;\longrightarrow\; u_{\mathrm{IEN}}
       + \beta_{\mathrm{stage}} \, \Delta t \, Q_{\mathrm{rad}},

so the stages sum to :math:`\Delta t \, Q_{\mathrm{rad}}`. **Only the gas energy is
modified**; there is no radiation force on the momentum.

``qrad_form`` selects how the rate is evaluated, following Eqs. 27 and 28 of the same
paper. The integral form uses the moments directly and is accurate where a cell is
optically thin, losing precision where it is thick and near equilibrium, as it then
differences two nearly equal numbers against a large :math:`\sigma_a`. The differential
form takes the divergence of the face-centred radiation flux and is accurate when thick,
degrading as :math:`\sigma_a \to 0`:

.. math::

   Q_{\mathrm{int}} = c_{\mathrm{rat}} \, p_{\mathrm{rat}} \, \sigma_a (J - B),
   \qquad
   Q_{\mathrm{dif}} = -\, c_{\mathrm{rat}} \, p_{\mathrm{rat}} \, \nabla \cdot \mathbf{H}.

``integral`` and ``divh`` select these. The default, ``hybrid``, blends them on the optical
depth per cell :math:`\tau = \chi \min(\Delta x)` — integral below :math:`\tau = 0.3`,
differential above :math:`\tau = 3`, smoothly interpolated between:

.. math::

   Q_{\mathrm{rad}} = (1 - w) \, Q_{\mathrm{int}} + w \, Q_{\mathrm{dif}},
   \qquad w = t^2 (3 - 2t),
   \qquad t = \frac{\log_{10}(\tau / 0.3)}{\log_{10}(3 / 0.3)} .

On a refined mesh the face fluxes used by the differential forms are corrected across
coarse-fine boundaries, so both sides of a refinement jump difference the same flux.

Timestep
~~~~~~~~

The module contributes two limits and the smaller is taken. The first is the radiation
relaxation rate of Davis, Stone & Jiang 2012 (Eq. 40), with the optical-depth factor of their
Eq. 42 replaced by a rational form that matches it in the thin and thick limits and is
slightly more restrictive between; the second bounds the internal energy the radiation may
remove in one step, and is disabled by ``cfl_qrad = 0``:

.. math::

   \Delta t \le \frac{C}{\nu_{\mathrm{rad}}}, \qquad
   \nu_{\mathrm{rad}} = \frac{4 (\gamma - 1) \, T^3 \, (\sigma_a / \rho) \,
                              c_{\mathrm{rat}} \, p_{\mathrm{rat}}}
                             {1 + \frac{3}{\pi^2}
                              \left( \Delta x_{\min} \, \chi \right)^2},

.. math::

   \Delta t \le C \; \mathtt{cfl\_qrad} \;
   \min_{\mathrm{cells}} \frac{e_{\mathrm{int}}}{\left| Q_{\mathrm{rad}} \right|},

with :math:`C` the CFL number. The second is evaluated at the end of the transfer solve, so
it bounds the very step its heating rate will be applied to.

Verification: an emitting sphere
--------------------------------

A homogeneous sphere of radius :math:`R` with constant :math:`\sigma_a` and source
:math:`S = b` in vacuum has an exact interior solution. Along a ray at
:math:`\mu_r = \hat{n} \cdot \hat{r}` and radius :math:`r < R`,

.. math::

   I(r, \mu_r) = b \left[ 1 - e^{-\sigma_a s} \right], \qquad
   s = r \mu_r + \sqrt{R^2 - r^2 \left( 1 - \mu_r^2 \right)},

with :math:`s` the distance back along the ray to the surface. Taking the moments of this **on the same angular
quadrature the solver uses** cancels the angle integration exactly, so the comparison
isolates transport and interpolation error.

The test is ``sc_sphere`` (``tst/inputs/sc_sphere.athinput``), which uses no fluid and runs a
single cycle. Running it at three resolutions,

.. code-block:: bash

   for N in 32 64 128; do
     ./athena -i tst/inputs/sc_sphere.athinput \
              mesh/nx1=$N mesh/nx2=$N mesh/nx3=$N \
              meshblock/nx1=$N meshblock/nx2=$N meshblock/nx3=$N
   done

writes the relative errors of :math:`J`, :math:`H` and :math:`K` for each run to
``sc_sphere-errs.dat``:

.. image:: /_static/nr_radiation_sphere.png
   :alt: Emitting sphere test, mean intensity against the analytic profile and convergence
   :align: center
   :width: 90%

The mean intensity tracks the analytic profile to better than 0.2 per cent through the
interior, and the relative :math:`L_2` error of :math:`J` falls from
:math:`1.4 \times 10^{-2}` at :math:`N = 32` to :math:`2.9 \times 10^{-3}` at
:math:`N = 128`, an observed order of about 1.15. First order is expected: the residual panel
shows the error concentrated in the last few cells before the surface, where the sphere is a
staircase on a Cartesian grid, which is why the test compares only cells inside
:math:`0.6 R`. Repeating the finest case on eight meshblocks reproduces the same error row
exactly.

Two further tests are bundled: ``sc_beam``, crossing pencil beams through a statically
refined region, and ``sc_linwave``, a radiatively damped acoustic wave run under adaptive
refinement in 1D, 2D and 3D, whose fitted frequency and damping rate are compared with the
analytic dispersion relation (Eq. 38 of Davis, Stone & Jiang 2012).

Current limitations
-------------------

- **No variable Eddington tensor.** :math:`K_{ij}` is computed, but
  :math:`f_{ij} = K_{ij}/J` is never formed and no moment equations are integrated. Closing
  the moment equations with it, after `Jiang, Stone & Davis 2012
  <https://ui.adsabs.harvard.edu/abs/2012ApJS..199...14J/abstract>`_, is the intended
  successor to the present formal solver.
- **Gray only.** The intensity carries no frequency dimension. A multiband extension, with
  each group transported independently and summed in the coupling term, would not change the
  structure of the sweep.
- **LTE only.** With :math:`S = B` and no scattering there is no photon destruction
  probability and no coupling between rays through the source function, so non-LTE problems —
  including the uniform-temperature scattering atmosphere — need a source-function iteration
  that the gray LTE limit makes unnecessary.
- **First order in time**, from the operator splitting.
- ``nmu`` is limited to 6, i.e. 168 rays in 3D.
