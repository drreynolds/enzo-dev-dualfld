.. _radiative_transfer:

Radiative Transfer
==================
.. sectionauthor:: John Wise <jwise@physics.gatech.edu>
.. versionadded:: 2.0

For relevant parameters, please also see :ref:`radiative_transfer_ray_tracing` and :ref:`radiative_transfer_fld`.


Adaptive Ray Tracing
--------------------

Solving the radiative transfer equation can be computed with adaptive
ray tracing that is fully coupled with the hydrodynamics and energy /
rate solvers.  The adaptive ray tracing uses the algorithm of Abel &
Wandelt (2002) that is based on the HEALPix framework.

For the time being, a detailed description and test suite can be found
in the paper Wise & Abel (2011, MNRAS 414, 3458).


Flux Limited Diffusion
----------------------

More details can be found in the paper Reynolds et al. (2009, Journal
of Computational Physics 228, 6833).



.. _radiative_transfer_dualfld:

Dual Flux Limited Diffusion (UV + X-ray)
--------------------------------------------

The *DualFLD* Enzo solver module is based off of the previous "gray"
flux limited diffusion solver module *gFLDSplit*, that uses
operator-splitting for coupling the radiation with Enzo's chemistry
and gas energy fields.  The DualFLD solver is designed to be highly
scalable, using a field-based approach for propagating either or both
of X-ray and UV radiation within an Enzo simulation.  For the
remainder of this section, we will refer to these fields as
:math:`E_{Xr}` and :math:`E_{UV}`.  Both of these are assumed to be
scalar-valued fields defined throughout the computational domain, and
are evolved using a flux-limited diffusion (FLD) approximation for the
propagation, absorption and dilution of radiation in a
cosmologically-expanding volume.  Each radiation field is assumed to
have either a monochromatic SED at a specified input frequency, or is
treated as an integrated radiation energy density with an assumed
radiation spectrum.  Separate SED choices may be made for both
:math:`E_{Xr}` and :math:`E_{UV}`.  Coupling between these radiation
fields is assumed to occur only through interaction with the baryonic
chemistry and total energy fields, i.e. there is no direct X-ray to UV
radiation coupling.  The interactions with Enzo's chemistry and energy
fields occurs through four mechanisms:

(a) Absorption of radiation due to the monochromatic or
    frequency-integrated opacities valid over the frequency spectrum
    where each field is valid. 

(b) Photo-ionization of chemistry fields due to the radiation present
    within each spatial cell. 

(c) Photo-heating of the total energy field due to radiation present
    within each spatial cell. 

(d) Secondary ionizations of chemistry fields due to X-ray radiation
    at frequencies exceeding 100 eV. 

While the DualFLD solver handles the propagation of both radiation
fields, it does not handle steps (b)-(d) directly, and instead
computes photo-ionization and photo-heating rates that are 
placed inside auxiliary ``BaryonField`` arrays and passed to Enzo's
built-in chemical ionization and gas heating solvers.

The defining characteristics of this solver in comparison with the
gFLDSplit solver are three-fold.  First, similarly to gFLDSplit, the
DualFLD solver propagates radiation in a separate operator-split step
from the chemistry and heating, though in this solver the radiation
solve itself is split into separate UV and X-ray solves (since they
are assumed to be non-interacting), and the photo-heating and
photo-ionization rates are post-processed from the combined radiation
contributions.  Second, the DualFLD solver itself *cannot* perform
ionization or heating directly, instead relying on Enzo's existing
solvers for those physics.  Third, the DualFLD module *only* allows a
standard chemistry-dependent model (i.e. it does not allow the *local
thermodynamic equilibrium* model that is supported by gFLDSplit).

We describe the mathematical models and solution methods here.  A
discussion of all relevant solver parameters is provided in the
associated section on :ref:`DualFLD parameters <radiative_transfer_dualfld>`.



.. _rad_model:

Flux-limited diffusion radiation model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We consider the equation for flux-limited diffusion radiative
transfer in a cosmological medium (Reynolds et al., J. Comput. Phys. 2009),

.. math::
   :label: radiation_PDE

   \partial_{t} E_i + \frac1a \nabla\cdot\left(E_i{\bf v}_{b}\right) =
    \nabla\cdot\left(D_i\,\nabla E_i\right) - \frac{\dot{a}}{a} E_i - c\kappa_i E_i + \eta_i,

where here the comoving radiation energy density field :math:`E_i`,
emissivity :math:`\eta_i` and opacity :math:`\kappa_i` are functions
of space and time, and where :math:`i\in\{\text{Xr},\text{UV}\}`.  In
this equation, the frequency-dependence of the respective radiation
energy density field has been integrated away, under the premise of an
assumed radiation energy spectrum, 

.. math::
   :label: spectrumXr

   & E_{\nu}(\nu,{\bf x},t) = \tilde{E}_{Xr}({\bf x},t) \chi_{Xr}(\nu) +
      \tilde{E}_{UV}({\bf x},t) \chi_{UV}(\nu), \\ 
   \Rightarrow & \\
   & E_{Xr}({\bf x},t) = \int_{0}^{\infty} E_{\nu}(\nu,{\bf x},t)\,\mathrm{d}\nu 
     = \tilde{E}_{Xr}({\bf x},t) \int_{0}^{\infty}
     \chi_{Xr}(\nu)\,\mathrm{d}\nu, \\

and

.. math::
   :label: spectrumUV

   E_{UV}({\bf x},t) = \int_{0}^{\infty} E_{\nu}(\nu,{\bf x},t)\,\mathrm{d}\nu 
    = \tilde{E}({\bf x},t) \int_{0}^{\infty} \chi_{UV}(\nu)\,\mathrm{d}\nu,

where :math:`\tilde{E}_{Xr}` and :math:`\tilde{E}_{UV}` are
intermediate quantities (for analysis) that are never computed, and
where we have assumed that the two spectra, :math:`\chi_{Xr}(\nu)` and
:math:`\chi_{UV}(\nu)` do not overlap (i.e. :math:`\chi_{Xr}(\nu)`
disappears in the interval :math:`[0,\nu_{1})` and
:math:`\chi_{UV}(\nu)` disappears in the interval
:math:`[\nu_{1},\infty)`).  We note that if either assumed spectrum
is the Dirac delta function, :math:`\chi_i(\nu) =
\delta_{\nu_i}(\nu)`, then :math:`E_i` is a monochromatic radiation
energy density at the ionization threshold :math:`h\nu_i`, and the
:math:`-\frac{\dot{a}}{a}E` term in equation :eq:`radiation_PDE`,
obtained through integration by parts of the redshift term 
:math:`\frac{\dot{a}}{a}\partial_{\nu}E_{\nu}`, is omitted from
:eq:`radiation_PDE`.  Similarly, the emissivity functions
:math:`\eta_i({\bf x},t)` relate to the true emissivity
:math:`\eta_{\nu}(\nu,{\bf x},t)` through the formulas 

.. math::
   :label: emissivity

   \eta_{Xr}({\bf x},t) = \int_{0}^{\infty}\eta_{\nu}(\nu,{\bf x},t)\,\mathrm{d}\nu,
   \quad\text{and}\quad
   \eta_{UV}({\bf x},t) = \int_{0}^{\infty}\eta_{\nu}(\nu,{\bf x},t)\,\mathrm{d}\nu.

Within equation :eq:`radiation_PDE`, the function :math:`D_i` is the
*flux limiter* that depends on face-centered values of :math:`E_i`,
:math:`\nabla E_i` and the opacity :math:`\kappa_i`
(Morel, J. Quant. Spectrosc. Radiat. Transfer, 2000), 

.. math::

   D_i = \min\left\{c \left(9\kappa_{i,f}^2 + R^2\right)^{-1/2}, D_{max}\right\},\quad\mbox{and}\quad
   R = \max\left\{\frac{|\partial_{x} E_i|}{E_{i,f}},R_{min}\right\}.

Here the spatial derivative within :math:`R` is computed using
non-dimensional units at the computational face adjoining two
neighboring finite-volume cells, :math:`D_{max}=0.006\,c\,L_{unit}`
and :math:`R_{min}=10^{-20}/L_{unit}` with :math:`L_{unit}` the length
non-dimensionalization factor for the simulation, and the
face-centered radiation energy density and opacity are computed using
the arithmetic and harmonic means, respectively, 

.. math::

   E_{i,f} = \frac{E_{i,1} + E_{i,2}}{2}, \quad\text{and}\quad
   \kappa_{i,f} = \frac{2\kappa_{i,1} \kappa_{i,2}}{\kappa_{i,1} + \kappa_{i,2}},

where here :math:`E_{i,1}` and :math:`E_{i,2}` are the two values of
:math:`E_i` in the cells adjacent to the face.  Among the many
available limiter formulations we have tested (Hayes & Norman,
Ap. J. Supp., 2003; Morel, J. Quant. Spectrosc. Radiat. Transfer,
2000; Reynolds et al., J. Comput. Phys. 2009), this version performs
best at producing causal radiation propagation speeds in the
low-opacity limit typical of reionization simulations.



.. _chem_model:

Model couplings
^^^^^^^^^^^^^^^^^^^^^^^^^^


In general, radiation calculations in Enzo are used in simulations
where chemical ionization states are important.  For these situations, 
we couple the radiation equation :eq:`radiation_PDE` with equations
for both the conservation of gas energy and primordial chemistry
ionization/recombination,  

.. math::
   :label: cons_energy

   \partial_t e + \frac1a{\bf v}_{b}\cdot\nabla e &=
    - \frac{2\dot{a}}{a}e
    - \frac{1}{a\rho_b}\nabla\cdot\left(p{\bf v}_{b}\right) 
    - \frac1a{\bf v}_{b}\cdot\nabla\phi + G - \Lambda  + \dot{e}_{SF}, \\

and

.. math::
   :label: chemical_ionization

   \partial_t {\tt n}_j + \frac{1}{a}\nabla\cdot\left({\tt n}_j{\bf v}_{b}\right) &=
    \alpha_{j,k} {\tt n}_e {\tt n}_k - {\tt n}_j \Gamma_{j}^{ph}, \qquad
    j\in\{\text{HI, HII, HeI, HeII, HeIII}\}. 

Here, :math:`{\tt n}_{j}` is the comoving number density for each
chemical species, :math:`{\tt n}_k` corresponds to chemical species
that interact with species :math:`{\tt n}_j`, and :math:`{\tt n}_e` is
the electron number density.  In these equations, all terms are
evolved by Enzo's built-in chemistry and gas energy solvers, though
some of the relevant rates result from radiation-dependent
couplings. Specifically, the gas can be photo-heated by the radiation
through the term

.. math::
   :label: G

   G &= \frac{c\,E_{UV}\,\sum_{j} {\tt n}_j
     \int_{\nu_{j}}^{\infty} \sigma_{j}\, \chi_{UV}
     \left(1-\frac{\nu_{j}}{\nu}\right)\,
     d\nu}{\rho_b\,\int_{0}^{\infty} \chi_{UV} d\nu}
   + \frac{Y_{\Gamma}\,c\,E_{Xr}\,\sum_{j} {\tt n}_j
     \int_{\nu_{j}}^{\infty} \sigma_{j}\, \chi_{Xr}
     \left(1-\frac{\nu_{j}}{\nu}\right)\,
     d\nu}{\rho_b\,\int_{0}^{\infty} \chi_{Xr} d\nu},

for :math:`j\in\{\text{HI,HeI,HeII}\}`, where the X-ray secondary
photo-heating coefficient :math:`Y_{\Gamma}` depends on the electron
fraction :math:`\xi` in a cell via the formula  

.. math::
   :label: Ygamma

   Y_{\Gamma} = 0.9971 \left[1 - \left(1-\xi^{0.2663}\right)^{1.3163}\right].

Within the Enzo code base, :math:`G` is stored in the baryon field
``PhotoGamma``, for communication between DualFLD and Enzo's
heating/cooling solvers.

Additionally, the photo-ionization rates :math:`\Gamma_{j}^{ph}`
within equation :eq:`chemical_ionization` depend on the X-ray and UV
radiation fields via the formulas 

.. math::
   :label: Gamma

   \Gamma_j^{ph} &= \frac{Y_{j}\, c\, E_{Xr}\, \int_{\nu_j}^{\infty}
     \frac{\sigma_j(\nu) \chi_{Xr}(\nu)}{\nu}\,\mathrm d\nu}{h\,
     \int_{0}^{\infty} \chi_{Xr}(\nu)\,\mathrm d\nu} 
   + \frac{c\, E_{UV}\, \int_{\nu_j}^{\infty}
     \frac{\sigma_j(\nu) \chi_{UV}(\nu)}{\nu}\,\mathrm d\nu}{h\,
     \int_{0}^{\infty} \chi_{UV}(\nu)\,\mathrm d\nu}.

In this formula, we employ the X-ray photo-ionization coefficients

.. math::
   :label: YH

   Y_{HI} &= 0.3908 \left(1 - \xi^{0.4092}\right)^{1.7592}, \\
   Y_{HeI} &= 0.0554 \left(1 - \xi^{0.4614}\right)^{1.666}, \\
   Y_{HeII} &= 0.

Within the Enzo code base, the rates :math:`\Gamma_{HI}^{ph}`,
:math:`\Gamma_{HeI}^{ph}` and :math:`\Gamma_{HeII}^{ph}` are held in
the baryon fields ``kphHI``, ``kphHeI`` and ``kphHeII`` for
communication between DualFLD and Enzo's chemistry solvers.

Lastly, the frequency-integrated opacities depend on the chemical
state at each spatial location, i.e.

.. math::
   :label: opacityXr

   \kappa_{Xr} \ = \ \frac{\sum_{j} 
   {\tt n}_{j} \int_{\nu_{j}}^{\infty} \chi_{Xr}\,\sigma_{j}\,d\nu}{
   \int_{0}^{\infty} \chi_{Xr}\,d\nu}, \quad j\in\{\text{HI,HeI,HeII}\}

and

.. math::
   :label: opacityUV

   \kappa_{UV} \ = \ \frac{\sum_{j} 
   {\tt n}_{j} \int_{\nu_{j}}^{\infty} \chi_{UV}\,\sigma_{j}\,d\nu}
   {\int_{0}^{\infty} \chi_{UV}\,d\nu}, \quad j\in\{\text{HI,HeI,HeII}\},

where these integrals with the assumed radiation spectra
:math:`\chi_{Xr}(\nu)` and :math:`\chi_{UV}(\nu)` handle the change
from the original frequency-dependent radiation equation to the
integrated gray radiation equations.


Within the DualFLD module, the baryon field ``kdissH2I`` is always set
to 0.




.. _solution_approach:

Numerical solution approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We evolve these models in an operator-split fashion, wherein we solve
the radiation equations :eq:`radiation_PDE` separately from the gas
energy correction and chemistry equations :eq:`cons_energy` and
:eq:`chemical_ionization`, which are evolved together.  These solves
are coupled to Enzo's existing operator-split solver framework in the
following manner: 

a. Evolve the radiation fields implicitly in time using an
   up-to-second-order accurate method. 

b. Project the dark matter particles onto the finite-volume mesh to
   generate a dark-matter density field; 

c. Solve for the gravitational potential and compute the
   gravitational acceleration field; 

d. Advect the dark matter particles with the Particle-Mesh method;

e. Evolve the hydrodynamics equations using an up-to-second-order
   accurate explicit method, and have the velocity :math:`{\bf v}_{b}` 
   advect both the X-ray and UV radiation fields, :math:`E_{Xr}`
   and :math:`E_{UV}`; 

f. Evolve the coupled gas energy correction and chemistry evolution
   equations;  

The implicit solution approach for the radiation step is similar to
the one from (Norman et al., Ap. J. Supp., 2015); in this
documentation we describe only enough to highlight the available user
parameters, and more fully describe some additional options available
in the solver.

In solving the radiation equations we use a method of lines approach
for the space-time discretization of :eq:`radiation_PDE`, in that we
first discretize the equations in space using a second-order-accurate,
uniform-grid, finite volume discretization, and then evolve the
resulting system of ODEs in time.




.. _rad_solve:

Radiation subsystem
"""""""""""""""""""""""

Assuming that all spatial derivatives are treated using standard
second-order centered difference approximations on our finite-volume
grid, we need only discuss the time-discretization of our radiation
equation :eq:`radiation_PDE`.  Our approach follows a standard
two-level :math:`\theta`-method, 

.. math::
   :label: radiation_PDE_theta

   E_i^n - E_i^{n-1} &- \theta\Delta t\left(\nabla\cdot\left(D_i^{n-1}\,\nabla E_i^n\right)
     - \frac{\dot{a}}{a} E_i^n - c\kappa_i^n E_i^n + \eta_i^n\right) \\ 

   & - (1-\theta)\Delta t\left(\nabla\cdot\left(D_i^{n-1}\,\nabla E_i^{n-1}\right) -
     \frac{\dot{a}}{a} E_i^{n-1} - c\kappa_i^{n-1} E_i^{n-1} +
     \eta_i^{n-1}\right) = 0, 

where the parameter :math:`0\le\theta\le 1` defines the
time-discretization, and where we have assumed that the advective
portions of :eq:`radiation_PDE` have already been taken care of
through Enzo's hydrodynamics solvers. Recommended values of
:math:`\theta` are 1 (backwards Euler) and 0.5 (trapezoidal,
a.k.a.~Crank-Nicolson).

Whichever :math:`\theta` value we use (as long as it is nonzero), the
equation :eq:`radiation_PDE_theta` is linearly-implicit in the
time-evolved radiation energy density :math:`E_i^n`.  We solve this in
predictor-corrector form (for ease of boundary condition
implementation).  Writing the linear system and corrector update as

.. math::
   :label: linear_system

   J s = b, \qquad E_i^n = E_i^{n-1} + s,

we approximately solve this linear equation for :math:`s` ,to a tolerance :math:`\delta`, 

.. math::
   :label: linear_system_approx

   \| J s - b \|_2 \le \delta,

using using a multigrid-preconditioned Krylov linear solver.
Different Krylov solvers (e.g. CG, BiCGStab, GMRES) may be used for
each field based on input parameters.

Both the X-ray and UV fields are solved using the same
time-discretization parameter :math:`\theta`, but may use different
sets of internal multigrid solver parameters.




.. _dt_selection:

Time-step selection
""""""""""""""""""""""""""""

Time steps are chosen adaptively in an attempt to control error in the
calculated solution.  To this end, we first define an heuristic
measure of the time accuracy error in a radiation field :math:`E_i` as 

.. math::
   :label: time_error
 
   err = \left(\frac1N \sum_{j=1}^N
    \left(\frac{E_{i,j}^{n}-E_{i,j}^{n-1}}{\omega_j}\right)^p\right)^{1/p}, 

where the weighting vector :math:`\omega` is given by 

.. math:: 
   :label: time_weighting

   \omega_j &= \sqrt{E_{i,j}^n E_{i,j}^{n-1}} + 10^{-3}, \quad j=1,\ldots,N.

i.e.~we scale the radiation change by the geometric mean of the old
and new states, adding on a floor value of 1e-3 in case any of the
states are too close to zero.  This approach works well when the
internal solution variables are unit-normalized, or at least close to
unit-normalized, since the difference between the old and new
solutions, divided by this weighting factor :math:`\omega`, should
give a reasonable estimate of the number of significant digits that
are correct in the solution. 

With these error estimates :eq:`time_error` for both :math:`E_{Xr}`
and :math:`E_{UV}`, we set the new time step size for each subsystem
based on the previous time step size and a user-input tolerance
:math:`\tau_{\text{tol}}` as 

.. math::
   :label: time_estimate

   \Delta t^{n} = \frac{\tau_{\text{tol}} \Delta t^{n-1}}{err}.

Since :math:`E_{Xr}` and :math:`E_{UV}` are evolved separately, we
allow either solver to subcycle at a faster rate if necessary to allow
convergence of the underlying linear solver.  However, in general we
enforce that both fields utilize the same step size,

.. math::
   :label: FLD_time_estimate

   \Delta t^{n} &= \min\{\Delta t_{Xr}^{n},\Delta t_{UV}^{n},\Delta t_{Enzo}^{n}\},

where :math:`\Delta t_{\text{Enzo}}` is the time step size that Enzo's
other routines (e.g.~hydrodynamics) would normally take.  We further
note that the DualFLD solver module will force Enzo to similarly take
this more conservative time step size, due to the tight physical
coupling between radiation transport and chemical ionization.

Additionally, a user may override these adaptive time step controls
with the input parameters :math:`\Delta t_{\text{max}}` and
:math:`\Delta t_{\text{min}}`.  However, even with such controls in
place the overall time step will still be selected to adhere to the
bound required by Enzo's other physical modules, i.e.

.. math::

   \Delta t^{n} &= \min\{\Delta t_{\text{min}}^{n},\Delta t_{Enzo}^{n}\}.





.. _variable_rescaling:

Variable rescaling
""""""""""""""""""""""""""""

In case Enzo's standard unit non-dimensionalization using 
``DensityUnits``, ``LengthUnits`` and ``TimeUnits`` is insufficient to
render the resulting solver values :math:`E_{Xr}` and :math:`E_{UV}`
to have nearly unit magnitude, the user may input additional variable
scaling factors to be used inside the DualFLD module.  The basic
variable non-dimensionalization of these fields is to create a
non-dimensional radiation field value by dividing the physical value
(in ergs/cm\ :sup:`3`) by the factor ``DensityUnits`` *
``LengthUnits``\ :sup:`2` * ``TimeUnits``\ :sup:`-2`.  As would be
expected, the values of :math:`E_{Xr}` and :math:`E_{UV}` may differ
by orders of magnitude, so it is natural that they should be
non-dimensionalized differently.

To this end, if we denote these user-input values as :math:`s_{Xr}`,
and :math:`s_{UV}`, then the DualFLD module defines the rescaled variables 

.. math::
   :label: variable_rescaling

   \tilde{E}_{Xr} = E_{Xr} / s_{Xr}, \qquad \tilde{E}_{UV} = E_{UV} / s_{UV},

and then uses the rescaled variables :math:`\tilde{E}_{Xr}` and
:math:`\tilde{E}_{UV}` in its internal routines instead of Enzo's
"non-dimensionalized" internal variables :math:`E_{Xr}` and
:math:`E_{UV}`.  If the user does not know appropriate values for
these scaling factors *a-priori*, a generally-applicable rule of thumb
is to first run their simulation for a small number of time steps and
investigate Enzo's HDF5 output files to see the magnitude of the
values stored internally by Enzo; if these are far from
unit-magnitude, appropriate scaling factors :math:`s_{Xr}` and
:math:`s_{UV}` should be supplied in the DualFLD parameter input file.




.. _boundary_conditions:

Boundary conditions
""""""""""""""""""""""""""""

As the radiation equation :eq:`radiation_PDE` is parabolic, boundary
conditions must be supplied on the radiation field :math:`E_i`.  The
DualFLD module allows three types of boundary conditions to be placed
on the radiation field: 

0. Periodic,
1. Dirichlet, i.e. :math:`E_i(x,t) = g(x), \; x\in\partial\Omega`, and
2. Neumann, i.e. :math:`\nabla E_i(x,t)\cdot n = g(x), \;
   x\in\partial\Omega`. 

In most cases, the boundary condition types (and values of :math:`g`) are
problem-dependent.  When adding new problem types, these conditions
should be set near the bottom of the file ``DualFLD\_Initialize.C``, 
otherwise these will default to either (a) periodic, or (b) will use
:math:`g=0`, depending on the user input boundary condition type.

