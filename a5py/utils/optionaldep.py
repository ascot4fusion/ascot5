"""Contains tools to import optional modules without causing import errors if
the module is missing.
"""
import importlib
from types import ModuleType


#pylint: disable=too-few-public-methods
class OptionalDependency:
    """A lazy, optional import wrapper for third-party packages.

    This class allows you to declare an optional dependency without causing
    an ImportError at module import time. Instead, the dependency is loaded
    only when first accessed. If the dependency is not available in the
    environment, an informative ImportError is raised *at the point of use*.

    Example
    -------
    Instead of::

        import optionalpackage as optional

    use::

        optional = OptionalDependency("optionalpackage")

    Then, when calling::

        optional.do_stuff()

    an ImportError will be raised if ``optionalpackage`` is missing.
    If the package is present, it will be imported on demand and cached for
    subsequent use.

    To satisfy type checkers and linters, one can::

        from typing import TYPE_CHECKING

        if TYPE_CHECKING:
            import optionalpackage

    Notes
    -----
    - **Graceful degradation**: code can be imported and executed even if
      optional dependencies are not installed, as long as they are not used.
    - **Lazy loading**: dependencies are imported only when needed, which
      can improve startup time.
    - **Caching**: once imported, the module object is cached, so repeated
      accesses incur no extra overhead.
    """

    _name: str
    """Name of the module (i.e. the optional package)."""

    _loaded: bool
    """Whether the module has been loaded or not."""

    _module: ModuleType
    """The module object."""

    def __init__(self, name: str) -> None:
        self._name = name
        self._loaded = False

    def _load(self) -> ModuleType:
        if not self._loaded:
            try:
                self._module = importlib.import_module(self._name)
            except ImportError:
                raise ImportError(
                    f"Optional dependency '{self._name}' is not installed. "
                    f"Please install it to use this feature."
                    ) from None
            self._loaded = True
        return self._module

    def __getattr__(self, item: str) -> ModuleType:
        return getattr(self._load(), item)
