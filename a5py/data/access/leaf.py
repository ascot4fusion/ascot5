"""Constants, methods, and classes which define and manipulate the metadata."""
from __future__ import annotations

from enum import Enum
import random
import datetime
from typing import Any, NamedTuple, Union, Optional, TYPE_CHECKING

import numpy as np

from a5py import utils
from a5py.exceptions import AscotIOException
from .dataholder import Status, DataStruct

if TYPE_CHECKING:
    from .tree import TreeManager
    from hdf5 import HDF5MiniManager


DEFAULT_TAG: str = "TAG"
"""Default tag and note."""

QIDLEN: int = 10
"""Number of characters (numbers) in quasi-unique identifier."""


class MetaData(NamedTuple):
    """Tuple for compact transfer of the meta data.

    Attributes
    ----------
    qid : str
        Unique identifier.
    date : str
        Date when the data was created.
    note : str
        Short note for the user to document the data.
    variant : str
        What data variant the data represents.
    """
    qid: str
    date: str
    note: str
    variant: str


class Leaf():
    """Leaf of a tree which has no further children managed by the tree.

    Holds metadata and acts as a leaf for the tree data structure. The metadata
    is immutable with the exception of `note`.

    Classes that have access to the actual data should inherit from this class.

    Attributes
    ----------
    qid : str
        Unique identifier for this data.
    date : str
        Date when this data was created.
    note : str
        Short note for the user to document this data.
    variant : str
        What data variant this instance represents.
    _treemanager: TreeManager
        Manager of the tree this leaf belongs to.
    """

    _registry = {}
    """Mapping between variant name and it's class.

    This registry is updated by decorating a class that inherits from Leaf with
    Leaf.register. This registry is used to read instances from HDF5 file in
    a way that the correct class is created.
    """

    def __init__(
            self,
            qid: str,
            date: str,
            note: str,
            variant: str,
            **kwargs: Any,
            ) -> None:
        """Initialize object which does not belong to a tree initially."""
        super().__init__(**kwargs)
        self._qid: str = qid
        self._date: str = date
        self._note: str = note
        self._variant: str = variant
        self._file: HDF5MiniManager | None = None
        self._treemanager: Optional[TreeManager] = None

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(qid={self._qid}, date={self._date}, "
            f"variant={self._variant})>"
            )

    def __setattr__(self, name: str, value: Any) -> None:
        """Change the value of the attribute unless it is read-only.

        Parameters
        ----------
        name : str
            Name of the attribute.
        value : any
            Value to set.

        Raises
        ------
        AscotIOException
            If the attribute is read-only.
        """
        if name in ["qid", "qqid", "date", "variant", "name"]:
            raise AscotIOException(
                f"Please don't do that ಠ_ಠ (Attribute {name} is read only)."
                )
        super().__setattr__(name, value)

    def _extract_tag(self) -> str:
        """Extracts a tag from note.

        Note is converted to a tag like this:

        1. The first word in the note is chosen, i.e., everything before
            the first whitespace.
        2. All special characters are removed from the first word.
        3. The word is converted to uppercase which becomes the tag.
        4. If the tag is invalid (empty string) or it starts with a number,
            the default tag is returned instead.

        Note that the returned value is not the actual tag used in the tree if
        there are other leafs with identical tags.

        Returns
        -------
        tag : str
            The tag.
        """
        tag_candidate = self.note.split(" ")[0]
        tag_candidate = "".join(ch for ch in tag_candidate if ch.isalnum())
        tag_candidate = tag_candidate.upper()
        if not tag_candidate or tag_candidate[0] in "1234567890":
            return DEFAULT_TAG
        return tag_candidate

    @property
    def _metadata(self) -> MetaData:
        """Metadata collected in a single object."""
        return MetaData(
            qid=self._qid, date=self._date, variant=self._variant,
            note=self._note,
        )

    @property
    def qid(self) -> str:
        """Unique identifier for this data."""
        return self._qid

    @property
    def qqid(self) -> str:
        """Unique identifier for this data with preceding 'q'."""
        return "q" + self._qid

    @property
    def date(self) -> str:
        """Date when this data was created."""
        return self._date

    @property
    def note(self) -> str:
        """Short note for the user to document this data."""
        return self._note

    @note.setter
    def note(self, note: str) -> None:
        """Set the note.

        Changing note may change the tag as well.
        """
        self._note = note
        if self._treemanager:
            self._treemanager.note_changed(self)

    @property
    def variant(self) -> str:
        """What variant of data this object represents."""
        return self._variant

    @property
    def name(self) -> str:
        """Full name of the data."""
        return f"{self.variant}_{self.qid}"

    def activate(self) -> None:
        """Set this data as active.

        Active inputs are used when the simulation is run. Active data variants
        are also used during post-processing by default unless otherwise
        specified.
        """
        if self._treemanager:
            self._treemanager.activate_leaf(self)

    def destroy(self, *, repack: bool=True) -> None:
        """Remove this data permanently.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file to reduce it's disk size.

            Deleting data in an HDF5 file with h5py only removes its references,
            not the actual data. Repacking copies the data to a new file and
            replaces the original, freeing up space. May take a moment for large
            files.
        """
        if self._treemanager:
            self._treemanager.destroy_leaf(self, repack=repack)


    def save(self):
        """Store data to disk."""
        if self.status is Status.SAVED:
            raise AscotIOException("Data is already saved to disk.")
        self._treemanager.save(self)
        for node, group in self._treemanager.nodes.items():
            if self in group:
                break
        self._file = self._treemanager.hdf5manager.get_minimanager(
            node, self.variant, self.qid
            )

    def status(self):
        """
        """
        if not self._file is None:
            return Status.SAVED
        return Status.STAGED


    @classmethod
    def register(cls, variant: str) -> function:
        """Add variant class to registry.

        Use this as a decorator for subclasses.
        """
        def decorator(subclass):
            cls._registry[variant] = subclass
            print(variant, subclass)
            return subclass
        return decorator

    @classmethod
    def create_leaf(cls, meta: MetaData, **kwargs) -> Leaf:
        """Initialize a Leaf (subclass) based on the name of the variant."""
        leaf_cls = cls._registry.get(meta.variant, cls)
        return leaf_cls(qid=meta.qid, date=meta.date, note=meta.note,
                        variant=meta.variant, **kwargs)


@Leaf.register("input")
class InputLeaf(Leaf):

    def __init__(self, qid, date, note, variant):
        super().__init__(qid=qid, date=date, note=note, variant=variant)
        self._cdata: DataStruct | None = None

    def export(self) -> Dict[str, np.ndarray]:
        """Return a dictionary with sufficient data to duplicate this instance.

        Returns
        -------
        data : dict[str, np.ndarray or unyt.unyt_array]
            Data that can be passed to the create method to duplicate this
            instance.
        """
        return {}

    def stage(self):
        """Make this data ready for simulation or evaluation.

        The memory consumption may increase significantly, so remember to
        unstage afterwards with :meth:`unstage`.
        """
        if self.status & Status.STAGED:
            raise AscotIOException("This instance is already staged.")


    def unstage(self):
        """Undo the effect of :meth:`stage` and free consumed memory."""
        if self.status is Status.STAGED:
            raise AscotIOException(
                "Cannot unstage instance which has not been saved beforehand."
                )
        if self.status is Status.SAVED:
            raise AscotIOException("This instance is not staged.")


def get_qid(
        dataset: Union[Leaf, str], with_prefix: bool=False,
        ) -> str:
    """Returns the QID from a given object or string.

    Parameters
    ----------
    dataset : `Leaf` or str
        Object that has QID, name of the object, or it's QID with or without
        the preceding 'q'.
    with_prefix : bool, optional
        Return the QID with the preceding 'q'.

    Returns
    -------
    qid : str
        A string with 10 digits.
    """
    if isinstance(dataset, str):
        if len(dataset) >= QIDLEN and dataset[-QIDLEN:].isdigit():
            qid = dataset[-QIDLEN:]
        else:
            raise ValueError(
                f"This doesn't appear to be a valid QID: {dataset}"
                )
    else:
        qid = dataset.qid

    if with_prefix:
        qid = f"q{qid}"
    return qid


def generate_metadata() -> tuple[str, str, str]:
    """Generate QID, date and default note/tag.

    Calls random number generator to create 32 bit string which is then
    converted as a QID string using left-padding with zeroes if necessary.

    Returns
    -------
    qid : str
        (Quasi-)unique identifier for the data.
    data : str
        Date right now.
    note : str
        Default note describing the data which the user can overwrite.
    """
    qid = str( np.uint32( random.getrandbits(32) ) ).rjust(QIDLEN, "0")
    date = utils.format2universaldate(datetime.datetime.now())
    note = DEFAULT_TAG
    return qid, date, note
