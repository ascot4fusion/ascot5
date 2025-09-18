import importlib

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

    Benefits
    --------
    - **Graceful degradation**: code can be imported and executed even if
      optional dependencies are not installed, as long as they are not used.
    - **Lazy loading**: dependencies are imported only when needed, which
      can improve startup time.
    - **Caching**: once imported, the module object is cached, so repeated
      accesses incur no extra overhead.
    """

    def __init__(self, name: str):
        self._name = name
        self._module = None
        self._loaded = False

    def _load(self):
        if not self._loaded:
            try:
                self._module = importlib.import_module(self._name)
            except ImportError:
                raise ImportError(
                    f"Optional dependency '{self._name}' is not installed. "
                    f"Please install it to use this feature."
                )
            self._loaded = True
        return self._module

    def __getattr__(self, item):
        return getattr(self._load(), item)
