Testing
*******

All tests are located in `tests` folder.
Tests are implemented with pytest and they can be executed with:

.. code-block:: bash

    make tests

Three types of test cases are omitted by default and they need to be run separately.

#.  MPI routines are tested with (see ``tests/test_engine/test_solve.py`` for examples)

    .. code-block:: bash

        mpirun -n 2 python -m pytest --with-mpi -k mpi

#.  Tests for the plotting routines (see ``tests/test_plotting/test.py`` for examples)

    .. code-block:: bash

        pytest --mpl --mpl-baseline-path=baseline tests/test_plotting/test.py

#.  Physics tests are executed with

    .. code-block:: bash

        pytest -s -v tests/test_physics

    These tests take considerably longer to execute and, hence, they are run separately.



Individual test files can be run explicitly

or even individual tests (note the brackets that specify the )

.. note::

    When testing C-functions directly use ``ctypes`` to make wrappers, see for an example.

Linters
*******

The C code is linted and type-checked by the compiler, but for Python we use ``pylint`` for linting and ``mypy`` for type checking.

These can be run with

.. code-block:: bash

    make lint

or explicitly for individual files:

.. code-block:: bash

    pylint --disable=similarities --generated-members=unyt.*
	mypy --follow-imports=silent --disable-error-code=import-untyped --disallow-untyped-defs --disallow-incomplete-defs

Remember that ``pylint`` is not the objective truth and sometimes exceptions are justified.
Still, we aim for 10.0/10.0 code quality just because that makes it easier to spot any new complains that ``pylint`` might have.
To suppress individual warnings, use

.. code-block:: bash



However, do not do this before thinking hard if there is a smarter way to write the code.

Coding style
************

-   Python: `black <https://pypi.org/project/black/>`_.
-   C: ClangFormat with following custom settings (can also be found in ``.clangformat``):

    .. code-block::

        AlignAfterOpenBracket: AlwaysBreak
        IndentWidth: 4
        UseTab: Never
        BreakBeforeBraces: Allman

-   Documentation: Write every sentence to a new line to make it easier to see changes in ``git diff``.

Debug mode
**********

The code can be run in debug mode to spot memory leaks in the engine.

.. code-block:: bash

    make clean
    make ascot DEBUG=1
    export LD_PRELOAD=$(gcc -print-file-name=libasan.so)
    python runascot.py
    export LD_PRELOAD= # Once you are done debugging

Regression tests
****************

In addition to unit tests (test functionality) and physics test (code validation) we also have regression tests for code verification.
These are tests that rerun archived simulations where ASCOT has been benchmarked against other codes or verified against experiments.
Since these are essentially "real" simulations they require proper computational resources and even a single test can take several hours to complete.

Making a release
****************

Tutorials
*********



Building the documentation
**************************

The documentation is built with Doxygen ➔ Breathe ➔ Sphinx toolchain.
You can build just the user/developer documentation with ``make doc-user`` / ``make doc-dev``.
Note that doxygen generated files are located in ``doc/_static/doxygen`` and these only removed with ``make cleandoc``.
