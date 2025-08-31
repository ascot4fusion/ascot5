"""Constants, methods, and classes which define and manipulate the metadata."""
from __future__ import annotations

import random
import datetime
import importlib
from typing import Tuple, Any, NamedTuple, Dict, Union
from importlib.metadata import version as importlib_version

import numpy as np

from ... import utils
from ...exceptions import AscotIOException

VERSION: str = importlib_version("a5py")
"""Current version of the code."""

DEFAULT_TAG: str = "TAG"
"""Default tag and note."""

QIDLEN: int = 10
"""Number of characters (numbers) in quasi-unique identifier."""

input_categories: Tuple[str, ...] = (
    "options", "bfield", "efield", "marker", "plasma", "neutral", "wall",
    "boozer", "mhd", "asigma", "nbi",
    )
"""All input categories."""

run_variants: Tuple[str, ...] = ("run", "afsi", "bbnbi")

simulation_diagnostics: Tuple[str, ...] = (
    "inistate", "endstate", "state", "orbits",
)
"""All simulation diagnostics."""


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


class MetaDataHolder():
    """Holds metadata and acts as a leaf for the tree data structure.

    Classes that have access to the actual data should inherit from this class.

    The metadata is immutable with the exception of `note`.

    Attributes
    ----------
    qid : str
        Unique identifier for this data.
    date : str
        Date when this data was created.
    note : str
        Short note for the user to document this data.
    variant : str
        What is the variant of the data this instance represents.
    """

    def __init__(
            self,
            qid: str,
            date: str,
            note: str,
            variant: str,
            **kwargs: Any,
            ) -> None:
        """Initialize object which does not belong to a tree initially.
        """
        super().__init__(**kwargs)

        self._qid: str = qid
        self._date: str = date
        self._note: str = note
        self._variant: str = variant

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
        """Set a short note to document this data."""
        self._note = note

    @property
    def variant(self) -> str:
        """What variant of data this object represents."""
        return self._variant

    @property
    def name(self) -> str:
        """Full name of the data."""
        return f"{self.variant}_{self.qid}"


def get_input_category(variant: str) -> str:
    """Return the input category corresponding to the given variant.

    Parameters
    ----------
    variant : str
        Variant of the data.

    Returns
    -------
    category : str
        Input category corresponding to the given variant.

    Raises
    ------
    ValueError
        If the variant is unknown.
    """
    for category in input_categories:
        if hasattr(importlib.import_module(f"a5py.data.{category}"), variant):
            return category
    raise ValueError(f"Unknown variant: {variant}")


def get_qid(
        dataset: Union[MetaDataHolder, str], with_prefix: bool = False,
        ) -> str:
    """Returns the QID from a given object or string.

    Parameters
    ----------
    dataset : `MetaDataHolder` or str
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


def generate_metadata() -> Tuple[str, str, str]:
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
