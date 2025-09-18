from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Generator

class ImmutableStorage():
    """Object which supports dictionary-like assignment and which can be made
    immutable.

    Attributes
    ----------
    frozen : bool
        Indicates whether the node is frozen, preventing attribute modification.
    """

    def __init__(self, **kwargs: Any) -> None:
        """Initialize an empty node which is unfrozen.

        Parameters
        ----------
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        self._frozen: bool = False
        super().__init__(**kwargs)

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return f"<{self.__class__.__name__}(frozen={self._frozen})>"

    def __setitem__(self, key: str, value: Any) -> None:
        """Add a new attribute this node in dictionary style.

        Parameters
        ----------
        key : str
            Name of the attribute.
        value
            Value of the attribute.

        Raises
        ------
        AttributeError
            Raised if the node is frozen.
        """
        if self._frozen:
            raise AttributeError(
                "The attributes of this class are immutable."
                )
        setattr(self, key, value)

    def __setattr__(self, key: str, value: Any) -> None:
        """Add a new attribute this node.

        Parameters
        ----------
        key: str
            Name of the attribute.
        value:
            Value of the attribute.

        Raises
        ------
        AttributeError
            Raised if the node is frozen.
        """
        if key != "_frozen" and self._frozen:
            raise AttributeError(
                "The attributes of this class are immutable."
                )
        super().__setattr__(key, value)

    def __getitem__(self, key: str) -> Any:
        """Retrieve attribute in dictionary-like manner.

        Parameters
        ----------
        key : str
            Name of the attribute.

        Returns
        -------
        value
            Value of the attribute.
        """
        return getattr(self, key)

    def _freeze(self) -> None:
        """Make this node immutable."""
        self._frozen = True

    def _unfreeze(self) -> None:
        """Make this node mutable."""
        self._frozen = False

    @contextmanager
    def _modify_attributes(self) -> Generator[ImmutableStorage, None, None]:
        """Open a context where attributes can be modified."""
        self._unfreeze()
        try:
            yield self
        finally:
            self._freeze()
