.. _Physics:

=================
Physics in ASCOT5
=================

Coordinate system
=================

.. default-role:: math

The coordinate system ASCOT5 uses is `COCOS3 <https://www.sciencedirect.com/science/article/abs/pii/S0010465512002962>`_ which is elaborated here.

The cylindrical coordinates that are used are right-handed meaning:

  - The vertical z-axis runs through the geometrical center of the torus and the conventional azimuthal angle, here called the toroidal angle `\phi`, increases counter-clockwise. The radial coordinate corresponds to the major radius of the device, R 
  - Toroidal field and plasma current are positive when they point in the same direction as `\hat{\phi}`.
  - In the azimuthal plane, the Cartesian coordinate axes `(x,y)` are aligned so that the toroidal angle is measured from the `x`-axis, i.e., the x-axis is chosen so that `\phi=0` corresponds to `y=0`.

If we slice the torus vertically, two cross sections become visible in the `(R,z)` plane. These are called poloidal cross sections and, according to our coordinate system, looking at the cross section on the right:

  - we are looking at the direction of (positive) `\hat{\phi}`.
  - for a positive plasma current, the poloidal field points clockwise.

.. list-table:: Sign of toroidal field, plasma current, and poloidal field direction in various machines in ASCOT5 convention
   :widths: 10 5 5 5
   :header-rows: 1

   * - Machine
     - `B_\mathrm{phi}`
     - `I_p`
     - `B_\mathrm{pol}`
   * - ITER
     - `-`
     - `+`
     - CW
   * - ASDEX Upgrade (AUG)
     - `-`
     - `+`
     - CW
   * - JET (RIP)
     - `-`
     - `-`
     - CCW
   * - JT60-SA
     - `-`
     - `-`
     - CCW

An observant reader notices that all tokamaks have negative toroidal field. The underlying reason is that they all have their X-point below the device and, empirically, the L-H transition is easier to achieve with the gradient drift pointing towards the X-point.

Another fact related to the directions is that the default direction of the beam injection is the direction of the plasma current. It is only with this configuration that the beam ions are generated on well-confined orbits even at the edge.

To describe the features inside the plasma, instead of the cylindrical coordinate system, the flux coordinate system is adopted:

  - The radial coordinate is rho, `\rho = \sqrt{(\Psi-\Psi_0) / (\Psi_1 - \Psi0)}`, where `\Psi` is the poloidal magnetic flux and `\Psi_0` and `\Psi_1` are the values on the axis and at the separatrix, respectively.
  - Therefore, \rho=0 corresponds to the magnetic axis and \rho=1 to the separatrix.
  - The poloidal angle `\theta` is the geometrical one, with positive `\hat{\theta}` pointing counter-clockwise when looking at the direction of positive `\hat{\phi}`.
  - Therefore, for positive `I_p`, `\hat{\theta}` and `\hat{B}_\mathrm{pol}` point in opposite directions.
  - `\theta = 0` at the outer mid-plane.
  - Note that in the code *theta* refers to the Boozer poloidal angle and *pol* refers to the geometrical angle.
  - For stellarators the poloidal angle we use is ill-defined since the magnetic axis is not a straight loop.

In the velocity space, the dynamic state of the particle is defined by its energy, `E`, and its pitch, `\xi`. For the particle pitch there are different conventions, but we use `\xi=v_\parallel/v`.
This means that `\xi=(-)1` corresponds to strongly (counter-)passing particle and `\xi=0` is a deeply trapped particle.

A note on the nomenclature: “marker” refers to an object that is simulated, “particle” to a physical particle (a marker with a  weight factor assigned to it), and “test particle” to an approximation where the particle under consideration does not affect the surrounding system.

The Big Picture
===============

ASCOT5 is a tool that solves *the Fokker-Planck equation* for a test particle species `s`:

.. math::
   :name: fokker-planck

   \frac{\partial }{\partial t}f(\mathbf{z},t)_s =
   -\frac{\partial }{\partial \mathbf{z}}\cdot \left[\mathbf{K}(\mathbf{z},t)f(\mathbf{z},t)_s\right]
   +\frac{\partial^2}{\partial\mathbf{z}\partial\mathbf{z}}:\left[\mathbf{D}(\mathbf{z},t)f(\mathbf{z},t)_s\right].

Here `f(\mathbf{z}, t)_s` is the distribution function to be solved, `\mathbf{z}` are the phase-space coordinates and `\mathbf{K}` and `\mathbf{D}` are the Fokker-Planck advection and diffusion coefficients.

The Fokker-Planck equation is not solved explicitly because we are apes and apes use brute force: the corresponding *Langevin equation*

.. math::
   :name: langevin

   d\mathbf{z} = \boldsymbol{\mu}(\mathbf{z},t)dt+ \boldsymbol{\sigma}(\mathbf{z},t)\cdot d\mathbf{W}

is solved for a large number of markers instead.
Here `\boldsymbol{\mu}=\mathbf{K}`, `(1/2)\boldsymbol{\sigma}\boldsymbol{\sigma}^\intercal=\mathbf{D}`, and `\mathbf{W}` is the so-called *Wiener process*.
(Those familiar with stochastic-differential equations may note that Eq. :math:numref:`langevin` is in Itô form which comes handy when implementing the collision operator).

Eqs. :math:numref:`fokker-planck` and :math:numref:`langevin` are connected so that when a marker motion is governed by Eq. :math:numref:`langevin` then the distribution of those markers

.. math::
   :name: marker-distribution

   f(\mathbf{z},t)_s\approx \sum_i^{N_\mathrm{markers}}w_i\delta(\mathbf{z}-\mathbf{z}_i(t)),

where `w_i` are individual marker weights, is a solution to Eq. :math:numref:`fokker-planck`.

This is in essence what ASCOT5 does.
The code can output marker *ini- and endstates*, which are just `f(\mathbf{z},t=t_0)` and `f(\mathbf{z},t=t_f)`, *orbit trajectories*, which are solutions to Eq. :math:numref:`langevin`, and *distributions* that approximate the distribution function with discrete histograms as

.. math::
   :name: marker-histogram

   f(z,t)\approx \sum_i\sum_{\alpha,\beta}\frac{f_{i,\alpha\beta}}{z_{\alpha+1}-z_\alpha}
   \boldsymbol{1}_{[z_\alpha\leq z < z_{\alpha+1}]}\boldsymbol{1}_{[t_\beta\leq t < t_{\beta+1}]},


where `\boldsymbol{1}` is the indicator function, and `\alpha` and `\beta` are the histogram bin labels.
Individual marker contribution to a histogram bin is calculated at each time-step as

.. math::
   :name: marker-contribution

   f_{i,\alpha\beta,k+1} = f_{i,\alpha\beta,k} + w_i(t_{k+1}-t_k) \boldsymbol{1}_{[z_\alpha\leq z_i(t_{k+1}) < z_{\alpha+1}]}\boldsymbol{1}_{[t_\beta\leq t_{k+1}< t_{\beta+1}]}.

ASCOT5 has two main limitations:

  1. It uses test-particle approximation.
  2. It assumes that the test particle population is produced by a constant source.

The test-particle approximation means that there is no feedback from test particle population to the background plasma nor there are interactions between the test particles in simulation time.
One can relax this approximation by running short simulations repeatedly and adjusting the background quantities between the simulations based on how the test particle population evolved.
Note that this doesn't mean that ASCOT5 is cabable of simulating minority species only: bulk plasma species can be simulated as long as one keeps these limitations in mind (e.g. estimating transport coefficients in steady-state plasma is fine).

The other approximation affects how one should interpret the quantities where marker weights are involved, i.e. wall loads and distributions.
The weight is not actually "how many physical particles this marker represents" but it is a particle flux and has units "particles/s".
This means that the wall loads are not in units of Joule but in units of Watts.
Again for stead-state plasmas this works perfectly fine but one must be careful when studying transient phenomena.
The distributions are steady-state distributions since every time a distribution is updated in a simulation, we place "weight * dt" in a bin corresponding to marker's current position.
This means that the resulting histogram has units of "particles".
Therefore one must be careful when interpreting distributions (or wall loads) in a simulation with an existing particle population that is not given by a constant source, e.g. runaway electrons in a disruption.

Orbit-following
===============

Markers can be traced using one of the following three schemes:

Field-line-tracing
******************

Marker is assumed to have no mass and travelling at the speed of light along the magnetic field lines.
The equation of motion is

.. math::
   :name: fieldline-equationsofmotion

   \dot{\mathbf{x}} = c\hat{\mathbf{b}},

which is solved with `the Cash-Karp method <https://doi.org/10.1145/79505.79507>`_ that uses an adaptive time-step.

Gyro-orbit a.k.a. particle
**************************

Marker is a physical particle and its whole gyro-motion is solved.
The Hamiltonian of a charged particle in an electromagnetic field is

.. math::
   :name: gyro-hamiltonian

   \mathcal{H}_\mathrm{prt} \equiv \gamma mc^2 +q \Phi,

where `\Phi` is electric potential and `\gamma = 1/\sqrt{1-v^2/c^2}` or, equivalently `\gamma= \sqrt{1+(p/mc)^2}`, is *the Lorentz factor*, which relates particle kinetic energy to its rest mass as `\gamma=1+E_\mathrm{kin}/mc^2`.
Hamiltonian dynamics yield the particle equations of motion:

.. math::
   :name: gyro-equationsofmotion

    \dot{\mathbf{x}} &= \frac{1}{\gamma m} \mathbf{p}\\
    \dot{\mathbf{p}} &= q\left(\mathbf{E}+\dot{\mathbf{x}}\times\mathbf{B}\right).

The numerical scheme used to solve these equations is `the Volume-Preserving Algorithm (VPA) <https://doi.org/10.1063/1.4916570>`_ which can be though as a relativistic variant of *the Boris scheme* since it preserves marker energy.
Usually the time-step is small (a fraction of gyro time) when using this scheme, so take care not to set it too small as then the limited machine precision starts to accumulate error.
This probably happens somewhere below `1\times10^{-12}` s.

Note that this scheme is valid also when the marker charge is zero, and therefore it is used in the code when tracing neutrals.

Guiding-center
**************

The gyro-orbit effects can be ignored to obtain faster simulations if:

  - Collecting guiding center distribution is sufficient.
  - Wall loads doesn't have to be exact.
  - Magnetic field doesn't vary *much* in time and space during a single gyro-orbit and therefore the magnetic moment is an adiabatic invariant.

In this case one can use the guiding-center approximation.
It is advised to approach a new study by first running both gyro-orbit and guiding-center simulations with limited number of markers to see if the guiding-center approximation is valid.
Usually it is unless the machine is small or a spherical tokamak with strong magnetic field gradients.

For the guiding center dynamics we employ non-canonical coordinates: guiding center position, `\mathbf{X}`, momentum component parallel to the magnetic field, `p_\parallel`, magnetic moment, `\mu`, and gyroangle, `\zeta`.
The so-called `guiding center transformation <https://doi.org/10.1017/S0022377815000744>`_, which is a coordinate transformation from particle phase space, `\mathbf{z}=(\mathbf{x},\mathbf{p})`, to guiding center phase space, `\mathbf{Z}=(\mathbf{X},p_\parallel,\mu,\zeta)`, is a near-identity transformation,

.. math::
   :name: gc-transformation

   \mathbf{X}  &= \mathbf{X}_0 + \epsilon\mathbf{X}_1 + \epsilon^2\mathbf{X}_2 + \ldots, \\
   p_\parallel &= p_{\parallel,0} + \epsilon p_{\parallel,1} + \epsilon^2p_{\parallel,2} + \ldots,\\
   \mu         &= \mu_0 + \epsilon\mu_1 + \epsilon^2\mu_2 + \ldots, \\
   \zeta       &= \zeta_0 + \epsilon\zeta_1 + \epsilon^2\zeta_2 + \ldots,

where `\epsilon` is a dimensionless ordering parameter which is used to group terms of similar size.
The transformation from particle to guiding center coordinates is performed to the first order in ASCOT5.
This can be adjusted from options, where the first order terms can be dropped, but this serves mainly one's curiosity and not practical applications.

The zeroth order terms in the transformation are

.. math::
   :name: gc-transformation0th

   \mathbf{X}_0    &= \mathbf{x},\\
   vp{\parallel,0} &= \mathbf{p}\cdot\hat{\mathbf{b}},\\
   \mu_0           &= \frac{p_\perp^2}{2mB},\\
   \zeta_0         &= \arctan2(-\hat{\boldsymbol{\rho}}\cdot\hat{\mathbf{e}}_2, \hat{\boldsymbol{\rho}}\cdot\hat{\mathbf{e}}_1).

The zeroth order term of the gyroangle is somewhat arbitrary as it is defined by basis vectors `\hat{\mathbf{e}}_1` and `\hat{\mathbf{e}}_2`:

.. math::
   :name: gc-basisvectors

   \hat{\boldsymbol{\rho}}  &=  \cos\zeta_0 \hat{\mathbf{e}}_1 - \sin\zeta_0 \hat{\mathbf{e}}_2\\
   \hat{\boldsymbol{\perp}} &= -\sin\zeta_0 \hat{\mathbf{e}}_1 - \cos\zeta_0 \hat{\mathbf{e}}_2.

These vectors can be chosen arbitrarily as long as `(\hat{\mathbf{e}}_1,\;\hat{\mathbf{e}}_2,\;\hat{\mathbf{b}})` form an orthogonal right-handed system.
Since `\hat{\mathbf{b}}` is fixed, we are free to choose `\hat{\mathbf{e}}_1`.
For cylindrical coordinates in tokamaks, a suitable choice is `\hat{\mathbf{e}}_1 = \hat{\mathbf{b}}\times\hat{\mathbf{z}}` because there is always a toroidal field present.

As for the first order terms, the first-order position-term is the gyro-vector

.. math::
   :name: gc-transformation1stpos

   \mathbf{X}_1=\boldsymbol{\rho}_g
   \equiv \frac{1}{q}\sqrt{\frac{2m\mu_0}{B}}\hat{\mathbf{b}}\times\hat{\mathbf{v}},

which is quite intuitive.
The first-order momentum space terms are less so:

.. math::
   :name: gc-transformation1stmom

   p_{\parallel,1} &= -p_{\parallel,0}\boldsymbol{\rho}_g\cdot\boldsymbol{\kappa}+\frac{m\mu_0}{q}\left( \tau_B+ \mathbf{a}_1:\nabla\hat{\mathbf{b}}\right),\\
   \mu_1           &= \boldsymbol{\rho}_g\cdot \left( \mu_0\nabla\ln B + \frac{p_{\parallel,0}}{mB}\boldsymbol{\kappa} \right)
   -\frac{\mu_0p_{\parallel,0}}{qB}\left( \tau_B + \mathbf{a}_1:\nabla\hat{\mathbf{b}} \right).

Here the dyadic is, `\mathbf{a}_1\equiv -\frac{1}{2}\left(\hat{\boldsymbol{\rho}}\hat{\boldsymbol{\perp}}+\hat{\boldsymbol{\perp}}\hat{\boldsymbol{\rho}}\right)`,
where the vectors `\hat{\boldsymbol{\rho}}` and `\hat{\boldsymbol{\perp}}` form an orthogonal right-handed basis
`(\hat{\boldsymbol{\rho}},\hat{\boldsymbol{\perp}},\hat{\mathbf{b}})` and
`\hat{\boldsymbol{\rho}}=\hat{\mathbf{b}}\times\hat{\mathbf{v}}`.
The magnetic field torsion, `\tau_B= \hat{\mathbf{b}}\cdot \nabla\times\hat{\mathbf{b}}`,
and the magnetic field twist, `\boldsymbol{\kappa} = \hat{\mathbf{b}}\cdot\nabla\hat{\mathbf{b}}`,
are related by the relation, `\nabla\times\hat{\mathbf{b}} = \tau_B\hat{\mathbf{b}} + \hat{\mathbf{b}}\times\boldsymbol{\kappa}`.
Finally, the first order gyroangle term is

.. math::
   :name: gc-transformation1stang

   \zeta_1 = -\boldsymbol{\rho}_g\cdot\mathbf{R} + \frac{p_{\parallel,0}}{qB} \left(\mathbf{a}_2:\nabla\hat{\mathbf{b}}\right) 
   + \frac{\rho_g}{B}\hat{\boldsymbol{\perp}}\cdot\left(\nabla B + \frac{p_{\parallel,0}^2}{2m\mu_0}\boldsymbol{\kappa}\right),

where `\mathbf{R}=\nabla\hat{\mathbf{e}}_1\cdot\hat{\mathbf{e}}_2` is the *Littlejohn's gyrogauge vector* and 

.. math::
   :name: gc-a2

   \mathbf{a}_2\equiv \frac{1}{4}\left(\hat{\boldsymbol{\perp}}\hat{\boldsymbol{\perp}}-\hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}\right).

Once the particle Hamiltonian has undergone the guiding-center transformation, it becomes the guiding center Hamiltonian

.. math::
   :name: gc-hamiltonian

   \mathcal{H}_\mathrm{gc} \equiv \gamma mc^2 +q \Phi(\mathbf{X},t),

where the Lorentz factor in the new coordinates is

.. math::
   :name: gc-gamma

   \gamma = \sqrt{1 + (2/mc^2)\mu B(\mathbf{X},t) + (p_\parallel/mc)^2}.

Note that the Hamiltonian does not depend on `\zeta`, which is as expected since the basis of the guiding center formalism is the decoupling of the gyro-motion, meaning guiding center dynamics must be independent of `\zeta`.
However, the gyroangle can be included as one of the phase space coordinates, for which Hamiltonian dynamics give the following (first-order) `equations of motion <http://dx.doi.org/10.1063/1.2773702>`_

.. math::
   :name: gc-equationsofmotion

   \dot{\mathbf{X}}         &= \frac{p_\parallel}{\gamma m} \frac{\mathbf{B}^*}{B_\parallel^*} + \mathbf{E}^*\times\frac{\hat{\mathbf{b}}}{B_\parallel^*},\\
   \dot{p}_\parallel        &= q\mathbf{E}^*\cdot\frac{\mathbf{B}^*}{B_\parallel^*},\\
   \dot{\boldsymbol{\mu}}   &= 0,\\
   \dot{\boldsymbol{\zeta}} &= \frac{qB}{\gamma m} + \dot{\mathbf{X}}\cdot\left(\mathbf{R} +\frac{\tau_B}{2}\hat{\mathbf{b}}\right),

with the effective fields being defined as

.. math::
   :name: gc-effbande

   \mathbf{B}^*&= \mathbf{B} + \frac{p_\parallel}{q}\nabla\times\hat{\mathbf{b}},\\
   \mathbf{E}^*&= \mathbf{E} -\frac{1}{q}\left( \frac{mc^2\mu}{\gamma}\nabla B -p_\parallel\frac{\partial \hat{\mathbf{b}}}{\partial t} \right),

and `B^*_\parallel=\hat{\mathbf{b}}\cdot\mathbf{B}^*`.

In the code, the guiding-center equations of motion can be solved with either RK4 (fixed time-step) or Cash-Karp (adaptive time-step).
These methods don't preserve the marker energy, but one can choose the time-step to bee small enough so that the resulting error is insignificant.

Note that when tracing guiding centers, the gyro angle is not solved so this information is lost and the transformation back to particle coordinates effectively uses a random gyro angle.

Collisions
==========

The test particle collision operator is based on *the Landau collision operator*, which can be expressed in the form of a Fokker-Planck equation, and the corresponding Langevin equation is

.. math::
   :name: collision-particle

   d\mathbf{p} = -K\mathbf{p}dt + \left(\sqrt{2D_\parallel}\hat{\mathbf{p}}\hat{\mathbf{p}} 
   + \sqrt{2D_\perp}\left(\mathbf{I}-\hat{\mathbf{p}}\hat{\mathbf{p}}\right)\right)\cdot d\mathbf{W},

where we have assumed that the background plasma is isotropic.
By further assuming that the plasma is Maxwellian, the advection coefficient, parallel diffusion coefficient and perpendicular diffusion coefficient have the explicit forms

.. math::
   :name: collision-coefficients

   K(v)           &= \sum_b\left(1+\frac{m}{m_b}\right)\frac{2\Gamma_{b}}{m^2 v_b^2}\frac{G(v/v_b)}{v},\\
   D_\parallel(v) &= \sum_b\frac{\Gamma_{b}}{v}G(v/v_b),\\
   D_\perp(v)     &= \sum_b\frac{1}{2}\frac{\Gamma_{b}}{v}\left(\mathrm{erf}(v/v_b)-G(v/v_b)\right),

where `\Gamma_{b}=n_bq^2q_b^2\ln\Lambda/4\pi\epsilon_0^2` (where `q` is charge, `\ln\Lambda` is the Coulomb logarithm and `\epsilon_0` is the vacuum permittivity) and the special functions `\mathrm{erf}(x)` and `G(x)` are the *error function*,

.. math::
   :name: errorfun

   \mathrm{erf}(x) \equiv \frac{2}{\sqrt{\pi}}\int_0^x e^{-s^2} ds,

and the *Chandrasekhar function*,

.. math::
   :name: chandrasekhar

   G(x) = \frac{\mathrm{erf}(x)-\mathrm{erf}'(x)}{2x^2} = \frac{\mathrm{erf}(x)-\frac{2x}{\sqrt{\pi}} e^{-x^2}}{2x^2}.

Note that the collision operator here is non-relativistic.
A relativistic variant exists but hasn't been implemented yet due to lack of interest.

The collision operator above is used in the gyro-orbit simulations.
The guiding-center test particle collision operator is given by three equations

.. math::
   :name: collision-gc

   dp          &= Q dt + \sqrt{2 D_\parallel} dW_p,\\
   d\xi        &= -\xi \nu dt + \sqrt{(1-\xi^2)\nu}dW_\xi, \\
   d\mathbf{X} &= \sqrt{2 D_\mathbf{X}}(\mathbf{I}-\hat{\mathbf{b}}\hat{\mathbf{b}})\cdot d\mathbf{W}_\mathbf{X},

where `W_p`, `W_\xi`, and `\mathbf{W}_\mathbf{X}` are independent Wiener processes.
The pitch collision frequency `\nu=\frac{2D_\perp}{p^2}` must be much smaller than the gyro-motion or otherwise the guiding-center approximation is broken by the collisions.
The drag term is given by

.. math::
   :name: collision-q

   Q=-Fp+\frac{\partial D_\parallel}{\partial p} +2\frac{D_\parallel}{p},

where

.. math::
   :name: collision-f

   F=\sum_b\frac{2\Gamma_{b}}{m^2 v_b^2}\frac{G(v/v_b)}{v}.

The spatial diffusion coefficient corresponding to the classical diffusion is

.. math::
   :name: collision-dx

   D_\mathbf{X} = \left[(D_\parallel - D_\perp)\frac{1-\xi^2}{2} + D_\perp\right]\frac{c^2}{\omega_g}.

The guiding center collision operator is obtained by transforming the particle Fokker-Planck equation to the guiding-center phase-space and performing gyro-averaging for the result.
Therefore collisions make the gyro angle intractable in guiding center simulations.
Furthermore, the guiding center collision operator uses a different set of momentum coordinates, `(p,\xi)`, than what is used in the orbit-following, `(p_\parallel,\mu)`.
This is because the collision operator cannot be diagonalized (which is desired for the numerical implementation) in the latter set of coordinates.
Note that in the code it is possible to toggle individual components in the guiding-center collision operator, which is useful e.g. if one wishes to isolate the effects of the pitch angle scattering on particle transport.

The collisions are applied separately in the code right after the orbit-step has been taken.
Physics-wise these should be evaluated simultaneously with the deterministic motion due to the Lorentz force residing inside the advection coefficient `\mathrm{K}`, but for practical applications this is not feasible.
This is because the orbit-integration requires a high-order numerical scheme whereas those schemes does not exist or they are overly complicated in case of stochastic differential equations.

The collisions in the particle picture are solved with the Euler-Maruyama method (a SDE version of the Euler method),

.. math::
   :name: eulermaruyama

   \mathbf{p}^{k+1} = \mathbf{p}^k - K\mathbf{p}^k\Delta t
   + \sqrt{2D_\parallel\Delta t}(\hat{\mathbf{p}}\cdot\boldsymbol{\beta})\hat{\mathbf{p}},
   + \sqrt{2D_\perp\Delta t}(\boldsymbol{\beta}-(\hat{\mathbf{p}}\cdot\boldsymbol{\beta})\hat{\mathbf{p}}).

Here `\boldsymbol{\beta}` is a random vector whose values are sampled from the Normal distribution.
In the guiding center picture, the collisions are resolved with the Euler-Maruyama methdon when a fixed time-step is used.
For the adaptive step it is necessary to use a higher-order method, for which we have chosen the Milstein method,

.. math::
   :name: milstein

   p^{k+1} &= p^k + Q \Delta t + \sqrt{2 D_\parallel } \Delta W_p
   + \frac{1}{2} \frac{\partial D_\parallel}{\partial p} \left((\Delta W_p)^2 - \Delta t\right), \\
   \xi^{k+1} &= \xi^k - \xi\nu\Delta t + \sqrt{(1-\xi^2)\nu}W_\xi
   - \frac{1}{2} \xi \nu \left((\Delta W_\xi)^2 - \Delta t\right).

The spatial component is still resolved with the Euler-Maruyama method.

If a time-step is rejected, it is not enough to simply repeat the time-step with a smaller `\Delta t`.
This is because we have already realized a value for the Wiener process `W(t)`.
Those values must be stored because if we have an interval `[t0, t1]` where Wiener processes has been realized on both ends, the values on the interval are not normally distributed with zero mean and variance $\Delta t$.
Instead the mean and the variance are given by the so-called Brownian bridge:

.. math::
   :name: brownianbridge

   E[W(t)]   &= W(t_0)+(W(t_1)-W(t_0))\frac{t-t_0}{t_1-t_0}, \\
   Var[W(t)] &= \frac{(t-t_0)(t_1-t)}{t_1-t_0}.

In other words, whenever a time-step is rejected, the generated Wiener process is stored and used to generate new values until the simulation passes that moment.

Finally, the collision operator have to deal with pitch being limited to range `[-1,1]` and the guiding center collision coefficients diverging at `p=0`.
These are dealt with by using reflecting boundary conditions for pitch and for momentum at some small value of `p` (a fraction of thermal momentum).


Wall model
==========

Wall model is either 2D contour or 3D mesh consisting of triangles.
If a straight line from marker initial position (at the beginning of the time-step) to its final position intersects the contour or one of the wall elements, a collision with the wall is recorded and simulation for that marker is terminated.

The collision algorithm in 2D is straight-forward: the wall is assumed to form a closed loop (and this is enforced by the code) and on each time-step `a winding number <https://www.engr.colostate.edu/~dga/documents/papers/point_in_polygon.pdf>`_ is calculated to determine if the marker is inside the wall polygon or not.
If the marker is outside, an algorithm is used to find which wall element the marker intersected.

In 3D, the bounding box of the wall model is divided along the axes into eight identical boxes which are then successively divided into smaller and smaller boxes for a fixed number of times.
This so-called `octree <https://en.wikipedia.org/wiki/Octree>`_ structure is used in the simulation to perform collision checks only with the elements that are in the same box as the marker.
The intersection between the line segment and the wall triangle is found with `the Möller–Trumbore algorithm <https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm>`_.

Atomic reactions
================

Currently atomic reactions are only available when using the gyro-orbit simulation mode.

TBD

Neutral beam injection
======================

Neutral markers are generated from the injector geometry using the beamlet-based model.
Ballistic trajectories of the neutral markers are then traced until i) the marker is ionized ii) the marker has intersected the wall.

WIP

Fusion source
=============

TBD

Magnetic field interpolation
============================

ASCOT5 uses modular inputs meaning that there are no specific way that inputs are interpolated during the simulation and new schemes can be included with tolerable effort.
However, the magnetic field data has huge impact on how accurate are the results are and how fast are the the simulations.
Therefore we review the magnetic field interpolation schemes here.
See here for details on other inputs.

One of the input types is the analytical representation of a tokamak field, which is fast and super-accurate, but rarely useful.
More commonly used input is the axisymmetric tokamak field, where the field is interpolated in two parts.
First the equilibrium component is evaluated from the poloidal flux `psi`,

.. math::
   :name: b2ds

   B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z},\\
   B_z &=  \frac{1}{R}\frac{\partial\psi}{\partial R},

and then we include `B_\mathrm{phi}` by interpolating values tabulated in `(R,z)` grid with cubic splines.
It is also possible to include tabulated values of `B_R` and `B_z` and sum those with Eq. :math:numref:`b2ds`, but this is rarely used as usually the poloidal field is completely defined by `psi`.
One possible use case is when `psi` is of poor quality and it is scaled so that it doesn't contribute to `\mathbf{B}_\mathrm{pol}` (but it can still be used to evaluate `\rho`), and the field is completely interpolated from the tabulated values of `\mathbf{B}`.

In 3D, the magnetic field evaluation works in a similar fashion except now `\mathbf{B}` is tabulated in `(R,\phi,z)` grid and `B_R` and `B_z` are usually non-zero as they contain the perturbation components.


Interaction with MHD modes
==========================

In the simulation it is possible to introduce EM-perturbations `\tilde{\mathbf{A}} = \alpha\mathbf{B}` and `\tilde{\Phi}` of the form

.. math::
   :name: mhd-alphaphidefinition

   \alpha       &= \sum_{nm} \lambda_{nm} \alpha_{nm}(\rho, t) \cos\left(n\zeta-m\theta-\omega_{nm}t\right),\\
   \tilde{\Phi} &= \sum_{nm} \lambda_{nm} \Phi_{nm}(\rho, t)   \cos\left(n\zeta-m\theta-\omega_{nm}t\right),

where `n` is toroidal mode number, `m` is poloidal mode number, `\omega_{nm}` is mode frequency, and `\alpha_{nm}` and `\Phi_{nm}` are mode eigenfunctions that may or may not depend on time.
The mode amplitude `\lambda_{nm}` is a scaling factor used to adjust `\tilde{B}/B` to a desired value.
The perturbations are evaluated in straight-field line coordinates `(\psi(\rho),\theta,\zeta)`, which are discussed separately below.

Mapping to straight-field-line coordinates
******************************************

During the simulation, the marker cylindrical coordinates are mapped to straight-field line coordinates if the MHD perturbations are enabled.
This mapping is implemented only for stationary tokamak fields and we further assume that the field is axisymmetric (but these can be used in non-axisymmetric fields as well).
Our choice of the coordinate system are the Boozer coordinates `(\psi,\theta,\zeta)`, where `\psi` is the poloidal flux, `\theta` is the Boozer poloidal angle which points in same direction as the geometrical poloidal angle `\theta_\mathrm{geo}` (counter-clockwise when looking at the same direction as positive `\hat{\phi}`, and `\zeta = \phi - \nu` is the Boozer toroidal angle with the same positive direction as the cylindrical toroidal angle.
Both Boozer angular coordinates have the periodicity of `2\pi`.
To faciliate the mapping in run-time, we precalculate `\theta(\rho,\theta_\mathrm{geo})` and `\nu(\rho,\theta)` in an uniform grid and use the tabulated values together with the cubic-spline interpolation to perform the mapping.

The computation of `\theta` and `\nu` is based on `these notes <https://youjunhu.github.io/research_notes/tokamak_equilibrium.pdf>`_ and is performed as follows.
First we find an equicontour of `\psi` on the `(R,z)` plane.
This process might fail near the axis or very close to the separatrix, which is why it is possible to set limits `[\rho_\mathrm{min},\rho_\mathrm{max}]` where the Boozer coordinates are defined.

The coordinate transformation requires the calculation of the Boozer Jacobian,

.. math::
   :name: boozer-jacobian

   J = \frac{I+qg}{B^2},

where `q(\psi)` is the safety factor, `g=RB_\mathrm{phi}`, and `I(\psi)` is toroidal current function that is related to enclosed plasma current as `I_p(\psi) = (2\pi/\mu_0)I(\psi)`, where `\mu_0` is the magnetic constant.
The safetu factor and the toroidal current function are evaluated using line integrals, where the integration starts from the outer mid-plane and proceeds in the same direction as `\hat{\theta}_\mathrm{geo}`:

.. math::
   :name: boozer-Iqg

   q &= \frac{1}{2\pi}\oint  \mathbf{B}_\mathrm{pol}\cdot d\mathbf{l},\\
   I &= \frac{1}{2\pi}\oint \frac{g}{R^2B^2_\mathrm{pol}} \mathbf{B}_\mathrm{pol}\cdot d\mathbf{l}.

Now the Boozer poloidal angle can be evaluated as

.. math::
   :name: boozer-theta

   \theta(\psi,\theta_\mathrm{geo}) = \frac{1}{2\pi}\int_0^{\theta_\mathrm{geo}} \frac{1}{JB^2_\mathrm{pol}}  \mathbf{B}_\mathrm{pol}\cdot d\mathbf{l},

and the Boozer toroidal angle with

.. math::
   :name: boozer-nu

   \nu(\psi,\theta) = -\frac{1}{2\pi}\int_0^\theta \frac{g}{R^2B^2_\mathrm{pol}} \mathbf{B}_\mathrm{pol} + q\theta

where q is the local safety factor.

The accuracy of the transformation can be assessed by verifying that `JB^2` is a flux surface function and `\mathbf{B}` is correct when these are evaluated using the Boozer coordinates:

.. math::
   :name: boozer-bvec

   \mathbf{B} &= q\nabla\theta\times\nabla\psi + \nabla\psi\times\nabla\zeta, \\
   J^{-1}     &= \nabla\theta\times\nabla\zeta\cdot\nabla\psi.

Including MHD in orbit-following
********************************

The perturbation is included in simulations when calculating the orbit-following part.
In gyro-orbit and field lines simulations, the pertubation components `\tilde{\mathbf{B}}` and `\tilde{\mathbf{E}}` are computed and directly added to the background field when solving the equations of motion.
Note that in the field-line simulations, the time is frozen so that the modes are not rotating and constructing field-line Poincaré plots show a snapshot of the field structure.

For the guiding center simulations, the perturbation is included by modifying the effective potentials, Eq. :math:numref:`gc-effbande`:

.. math::
   :name: mhd-effbande

   \mathbf{B}^{**}&= \mathbf{B}^{ *} + \nabla\times(\alpha\mathbf{B}), \\
   \mathbf{E}^{**}&= \mathbf{E}^{ *} - \frac{\partial \alpha\mathbf{B}}{\partial t} - \nabla\tilde{\Phi}

If the modes are rapidly rotating so that the electrons are able to balance any electric field parallel to the field lines, we have a condition `E_\parallel=0` which makes the magnetic and electric perturbations co-dependent:

.. math::
   :name: mhd-alphafromphi

   \omega_{nm}\alpha_{nm} = \frac{nq - m}{I+gq}\Phi_{nm}.

This is not enforced in the code so ensuring it is user's responsibility.

Another useful property is the conservation of

.. math::
   :name: mhd-h

   K = H - \omega_n P / n,

where `H` is the Hamiltonian and `P` canonical angular toroidal momentum.

Backward Monte-Carlo
====================

This section is not done yet, but you can find the reference `here <https://iopscience.iop.org/article/10.1088/1741-4326/ac3a1b>`_.
