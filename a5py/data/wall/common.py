"""Common tools shared by wall instances."""
import copy

import numpy as np

from a5py.data.access import InputVariant

class WallVariant(InputVariant):
    """Common methods shared by wall instances."""


    _labels: dict[str, int] | None
    """Human readable labels for the flag values.

    These are not needed and, hence, not stored in libascot.so. Therefore we
    store them here when this instance is staged.
    """

    @property
    def n(self) -> int:
        """Number of wall vertices."""
        raise NotImplementedError

    @property
    def flag(self) -> np.ndarray:
        r"""Integer label for grouping wall elements."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("flag", shape=(self.n,))
        assert self._file is not None
        return self._file.read("flag")

    @property
    def labels(self) -> dict[str,int]:
        r"""Human readable labels for the flag values."""
        if self._cdata is not None:
            assert self._labels is not None
            return copy.deepcopy(self._labels)
        assert self._file is not None
        return {
            label.decode("utf-8"): int(flag) for label, flag in
            zip( self._file.read("labelkeys"), self._file.read("labelvalues") )
            }

    def _save_flags_and_labels(self) -> None:
        labels = np.char.encode(list( self.labels.keys() ), "utf-8")
        flags = np.array( list(self.labels.values()) )
        assert self._file is not None
        self._file.write("flag", self.flag)
        self._file.write("labelkeys", labels)
        self._file.write("labelvalues", flags)

    @staticmethod
    def validate_labels(
        labels: dict[str, int] | None,
        flag: np.ndarray,
        ) -> dict[str, int]:
        """Validate labels.

        Helper method to be used when the labels are provided by user.
        """
        if labels is None:
            labels = {"wall": 0}

        for label, label_flag in labels.items():
            if not isinstance(label, str):
                raise ValueError("Labels must be strings.")
            if not isinstance(label_flag, int):
                raise ValueError("Flags in 'labels' must be integers.")

        assert isinstance(flag, np.ndarray)
        unlabeled = [f for f in flag if f not in labels.values()]
        if unlabeled:
            raise ValueError(
                f"Flags `{unlabeled}` do not have a corresponding label."
                )

        return labels
