=============
Citing ASCOT5
=============

It is sufficient to acknowledge ASCOT5 usage with something like:

  *Simulations were carried out using ASCOT5 orbit-following code [REF].*

Where the refence is either the official (unpublished) ASCOT5 reference paper :cite:labelpar:`varje2019high`, Jari's thesis :cite:labelpar:`Varjephd`, Konsta's thesis :cite:labelpar:`Sarkimakiphd`, or ASCOT4 reference paper :cite:labelpar:`hirvijoki2014ascot`.

If the work contains considerable ASCOT5 work, the referees might insist on a more detailed description.
At minimum, describe i) how markers were generated and how many were used, ii) what physics were included or excluded, iii) when markers where terminated, and iv) whether the simulations were done in guiding-center or gyro-orbit mode.

  *Markers representing fusion-born alpha particles were sampled from a fusion source distribution.*
  *Total of 10^6 markers were traced until they made contact with the wall or cooled below 2 x Te.*
  *Markers were traced using the guiding-center approximation and the simulations included Coulomb collisions.*

It is also a good practice to refer as *markers* the points whose trajectories ASCOT5 simulates, and *particles* when discussing physical particles (markers multiplied by their weight).
Furthermore, *full orbit* can refer to both *poloidal orbit* and *gyro-orbit* meaning it is better to write "gyro-orbit" explicitly (e.g. when comparing the results to guiding-center simulations).
Also, pitch can be defined either as the angle between velocity and magnetic field vectors but it can also be the cosine of that angle (the definition what ASCOT5 uses) so please be explicit when writing that "pitch is zero".

=======
License
=======

This software is distributed under the `LGPLv3 <https://www.gnu.org/licenses/lgpl-3.0.html>`_ license.

=======
Gallery
=======

.. figure:: ../figures/sparcre.png
   :class: with-border
   :align: center
   :width: 300px

   Runaway electron transport coefficients calculated with SPARC.
   Ref. :cite:t:`tinguely2021modeling` studied an external coil for runaway electron mitigation.
   Markers were traced with ASCOT5 using magnetic fields computed with NIMROD, and the resulting coefficients allowed DREAM to take losses into account when calculating runaway electron beam evolution.

========================
Publications using ASCOT
========================

These include ones that have used ASCOT5 or its predecessors ASCOT4 and ASCOT3.

.. bibliography::
   :all:
   :style: plain
