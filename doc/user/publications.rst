Publications
============

Citing ASCOT5
*************

It is sufficient to acknowledge ASCOT usage with something like:

  *Simulations were carried out using the ASCOT5 orbit-following code [REF].*

where the refence is the "official" ASCOT5 reference paper :cite:labelpar:`Varje2019ASCOT5`.

If the work contains considerable ASCOT work, the referees might insist on a more detailed description.
At minimum, describe **i)** how markers were generated and how many were used, **ii)** what physics were included or excluded, **iii)** when markers where terminated, and **iv)** whether the simulations were done in guiding-center or gyro-orbit mode.

  *Markers representing fusion-born alpha particles were sampled from a fusion source distribution.*
  *Total of* :math:`10^6` *markers were traced until they made contact with the wall or cooled below* :math:`2\times T_i`.
  *Markers were traced using the guiding-center approximation and the simulations included Coulomb collisions.*

.. note::

    It is also a good practice to refer as *markers* the point-like objects whose trajectories ASCOT simulates, and *particles* when discussing physical particles (markers multiplied by their weight).

    Furthermore, *full orbit* can refer to both *poloidal orbit* and *gyro orbit*, meaning it is better to write "gyro-orbit" explicitly (e.g. when comparing the results to guiding-center simulations).

    Also, pitch can be defined either as the angle between velocity and magnetic field vectors but it can also be the cosine of that angle (the definition what ASCOT uses).
    Please be explicit in what definition is used.

Gallery of selected works
*************************

.. card-carousel:: 1

    .. card:: Power loads enhanced by CX in W7-X
        :width: 100%

        .. figure:: ../_static/gallery_Kiviniemi2025W7X_SOL.png

            Power loads in the lower target region as predicted by ASCOT.
            The result shows increased losses due to fast ion interaction with the background neutrals.
            Image reproduced and enhanced from Ref. :cite:labelpar:`Kiviniemi2025W7X_SOL`.

    .. card:: Alpha particle confinement in Infinity Two stellarator
        :width: 100%

        .. figure:: ../_static/gallery_Carbajal2025InfinityTwo.png

            The birth locations of alpha particles that were lost in the ASCOT simulation :cite:labelpar:`Carbajal2025InfinityTwo`.

    .. card:: Runaway electron transport in SPARC
        :width: 100%

        .. figure:: ../_static/gallery_Tinguely2021Modeling.png

            Runaway electron transport coefficients calculated by ASCOT in post-disruption field with passive runaway electron mitigation coil present :cite:labelpar:`Tinguely2021Modeling`.

    .. card:: Losses in phase space in MAST-U
        :width: 100%

        .. figure:: ../_static/gallery_Ollus2024MASTU_CX.png

            Simulated (ASCOT+FIDASIM) and measured FIDA spectra :cite:labelpar:`Ollus2024MASTU_CX`.

    .. card:: Thea Energy: Alpha particle losses in quasi-axisymmetrix stellarator
        :width: 100%

        .. figure:: ../_static/gallery_Bonofiglo2025QAStellarator.png

            Investigation of fast-ion loss channels with ASCOT.
            Figure shows the origin of the losses :cite:labelpar:`Bonofiglo2025QAStellarator`.

    .. card:: Post-ELM fast ion losses in TCV
        :width: 100%

        .. figure:: ../_static/gallery_Poley2025First.png

            Synthetic FILD signal calculated with ASCOT and FILDSIM compared to the experimental result :cite:labelpar:`Poley2025First`.

    .. card:: CX induced beam ion losses in EAST
        :width: 100%

        .. figure:: ../_static/gallery_Zhang2026Investigation.png

            A snapshot of surviving beam-ion population after 0.5 ms calculated with ASCOT :cite:labelpar:`Zhang2026Investigation`.
            An illustration of a trapped particle orbit is shown in red and green.

    .. card:: Fast-ion redistribution in DIII-D
        :width: 100%

        .. figure:: ../_static/gallery_Gonzalez2022Flow.png

            Flow of fast ions simulated with NOVA-K, ASCOT and INPASIM at different time instances :cite:labelpar:`Gonzalez2022Flow`.

    .. card:: Synthetic ITER FILD signal
        :width: 100%

        .. figure:: ../_static/gallery_Reyner2026Optimized.png

            Velocity-space distributions of the lost alpha particles reaching the FILD pinhole for different radial positions, computed with ASCOT :cite:labelpar:`Reyner2026Optimized`.

    .. card:: Beam-driven chirping modes in ST40
        :width: 100%

        .. figure:: ../_static/gallery_Bland2023Interplay.png

            The resonance condition for beam ions calculated with ASCOT :cite:labelpar:`Bland2023Interplay`.
            The mode frequencies are labeled with red arrows.

.. rubric:: Other works

ASCOT5 was featured as one of the contest applications in the `ISC Student Cluster Competition 2026 <https://www.hpcadvisorycouncil.com/events/2026/student-cluster-competition/benchmarking.php>`_.

Studies using ASCOT5
********************

.. bibliography::
    :all:
