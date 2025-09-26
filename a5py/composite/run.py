from typing import Optional

import unyt
import numpy as np

from a5py.data.access import OutputVariant
from a5py.data.marker.state import MarkerState

class Run(OutputVariant):

    def _setup(self, mrk):
        self._diagnostics["endstate"] = MarkerState.from_params(mrk)

    def _load(self, file):
        super()._load(file)
        for name in self._file.children:
            if name == "endstate":
                self._diagnostics[name] = MarkerState()

    def save(self):
        super().save()
        for name, diag in self._diagnostics.items():
            if name == "endstate":
                file = self._file.access_data("endstate")
                diag.save(file)

    def getstate(
            self, *qnt: str,
            mode: str="gc",
            state: str="ini",
            ids: Optional[list]=None,
            endcond: Optional[str]=None,
            ) -> np.ndarray | unyt.unyt_array:
        """Evaluate quantities from markers' ini- and endstate.

        Inistate is marker's phase-space position right at the start of the
        simulation and endstate is the position at the end of the simulation.

        This function not only returns the marker phase space coordinates but
        also other quantities that can be inferred from it. For a complete list
        of available quantities, see. The returned arrays are sorted by marker
        ID.

        All simulations store both particle and guiding center phase-space
        coordinates. *By default, the returned quantity is in guiding-center
        picture*.

        Parameters
        ----------
        *qnt : str
            Names of the quantities.
        mode : {"gc", "prt"}, *optional*
            Choose whether to return guiding-center or particle quantity.
        state : {"ini", "end"}, *optional*
            Choose whether to return inistate or endstate quantity.
        ids : array_like, *optional*
            Filter markers by their IDs.

            If given, the values in the returned arrays in the same order as
            the given IDs. Cannot be used with ``endcond``.
        endcond : str or list [str], optional
            Filter markers by their end conditions.

            See for a list of all possible end conditions or to list end
            conditions that are currently present in the data.

            Markers may have multiple end conditions active simultaneously. If
            just the name of the end condition e.g. "POLMAX" is passed, then all
            markers that have (at least) the ``POLMAX`` end condition are
            returned.

            If the end condition is preceded by "NOT", e.g. "NOT POLMAX", then
            markers that don't have that end condition are returned.

            Passing multiple end conditions in a single string returns markers
            that have all listed end conditions active, e.g. "MAXPOL MAXTOR"
            returns markers that have both ``POLMAX`` and ``TORMAX`` active
            simultaneously.

            Passing end condition strings as separate list items acts as
            a logical OR, e.g. ["POLMAX", "TORMAX"] returns markers that have
            either ``POLMAX`` or ``TORMAX`` active.

        Returns
        -------
        *val : array_like or tuple(array_like)
            The evaluated quantity sorted by marker ID.

            If multiple quantities are queried, they are returned as a tuple in
            the order they were listed in ``*qnt``.

        Raises
        ------
        ValueError
            Raised when the queried quantity could not be interpreted.
        AscotDataException
            Raised when data required for the operation is not present.
        """
        if state not in ["ini", "end"]:
            raise ValueError("Unrecognized state: " + state)
        self._require("marker")
        if endcond is not None or state == "end":
            self._require("endstate")
        if state == "end":
            self._require("_endstate")

        # Get or evaluate the quantity
        data = getattr(self, "_" + state + "state").get(*qnt, mode=mode)

        # Parse endcond
        idx = np.ones(data[0].shape, dtype=bool)
        if endcond is not None:
            if not isinstance(endcond, list):
                endcond = [endcond]

            # Go through each unique end cond and mark that end cond valid or
            # not. This can then be used to make udix as boolean mask array.
            uecs, uidx = np.unique(self._endstate.get("endcond")[0],
                                   return_inverse=True)
            mask = np.zeros(uecs.shape, dtype=bool)
            for i, uec in enumerate(uecs):
                for ec in endcond:
                    ec = ec.replace(" and ", " ")
                    accept  = State.endcond_check(uec, ec)
                    mask[i] = mask[i] or accept

            idx = mask[uidx]

