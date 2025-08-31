.. _Codingstyle:

============
Contributing
============

How to contribute
*****************

.. tab-set::

   .. tab-item:: Anyone

      .. card::

         Fork the repository and submit your code changes for review via a pull request.

	      Forking can be conveniently done through `GitHub <https://github.com/ascot4fusion/ascot5/fork>`_.
	      A forked repository is an independent copy of the original, allowing you to freely push your changes without affecting the original code.

	      Pull requests (PRs) are special requests that include the modifications you wish to incorporate into the main repository.
	      They are reviewed by the maintainers, who may request modifications before accepting them.
	      You can even submit a PR when your work is not yet complete to ask for help.

	      To submit a PR for a forked repository, simply click *Contribute* on the main page of your fork.

   .. tab-item:: Developers

      .. card::

         Create your own feature branch and merge it via a `pull request <https://github.com/ascot4fusion/ascot5/pulls>`_ after approval.

         If you're solely improving the documentation, you can directly commit and push your changes to the ``docs`` branch.

         To become a developer, reach out to the maintainers for further information.

   .. tab-item:: Users

      .. card::

         Provide the maintainers with the HDF5 file whenever you conduct a study where ASCOT5 is benchmarked against other similar codes or when numerical results are compared to experiments.

         Our goal is to compile a record of simulations that *validate* the code.
         The data will be hosted on `Zenodo <https://zenodo.org/communities/ascot4fusion/records>`_ (with closed access if necessary) and will be used for regression testing during code development.

ASCOT5 is being actively developed, so we welcome new contributions.
If you have a post-processing script or plot that you believe could benefit others, please submit it via a pull request, even as a standalone script.
We can then collaborate to integrate it into the code.
In particular, we're missing many tools that serve as interfaces between ASCOT5 and other codes and data formats.

Git usage
*********

The code is versioned using Git so knowledge of basic commands is necessary.

.. dropdown:: Retrieving repository and switching to an existing branch:

   .. code-block:: bash

      git clone git@github.com:ascot4fusion/ascot5.git
      cd ascot5
      git checkout <branch name>

   Please note that cloning via SSH requires generating and setting SSH keys in your GitHub profile.
   These keys are specific to your machine, so you'll need to repeat this process when transitioning to another platform.

   The repository has three special branches:

   - ``main``: Contains stable code intended for production use.
   - ``develop``: Collects upcoming features.
   - ``docs``: For updating documentation without creating a new release.

.. dropdown:: Creating a new branch:

   .. code-block:: bash

      git checkout -b <branch name>
      git push origin

   Alternatively, you can create branches on GitHub and link them to issues (this option is available in the sidebar when viewing an issue).

   All new branch names must start with either ``hotfix/`` or ``feature/``, followed by the issue number and a descriptive name.
   For instance, ``hotfix/95-smoke-from-cpu`` and ``feature/22-adas-interface`` are valid branch names.

   Hotfix branches should branch from ``main`` since, by definition, they fix a bug in ``main``.
   Otherwise branch from ``develop``.
   If you accidentally branch from ``main`` you can always rebase your branch onto ``develop``.
   Keep in mind that Git branches are lightweight, so don't hesitate to further branch your feature branch if needed.

.. dropdown:: Keeping your branch synchronized with main repository.

   .. code-block:: bash

      git fetch
      git stash
      git pull --rebase origin/develop
      git stash apply
      git push origin --force

   Note the use of the rebase option, which is necessary because we maintain a linear history on the main branch.
   Consider a scenario where you've created your feature branch, made several commits, and someone pushes a hotfix to the ``main`` branch (or a new feature to ``develop``).
   To maintain a linear history, this hotfix must appear in the commit log before your commits.
   Rebase moves the base of your branch to the tip of the ``develop`` branch and then reapplies your commits.
   Forgetting to use ``--rebase`` can lead to complications, and following Git's "helpful tips" in such cases can exacerbate the situation.
   In such scenarios, it's best to use ``git merge --abort`` and start afresh.

   Stashing stores any local changes you have made but not yet committed.

   .. warning::

      Force-pushing overwrites the branch in the main repository with your local version.
      After rebasing, you've essentially rewritten the history of your branch so that the commits on the ``main`` or ``develop`` branch appear before your own.
      Consequently, Git will complain that your branch isn't up to date with the remote, necessitating the use of ``--force`` when pushing.

   There is no need to manually rebase ``develop`` or ``docs`` to ``main`` as this is handled by a bot.

.. dropdown:: Making a record of your changes

   .. code-block:: bash

      git add <new files or files you have modified>
      git commit -m 'Your descriptive commit message'
      git push origin

   This makes a commit in the version control system and publishes it in the main repository.
   To see overview of what modifications you have made, use ``git status`` and ``git diff <file>`` for details.

   Don't push binaries, pictures, or anything else but text files.
   If you accidentally pushed some large file, notify the maintainers who will remove it from the Git history.

   If you committed something accidentally, you can revert the previous commit with ``git reset HEAD~``.
   If you also pushed it in your feature branch, you can force push your local branch after committing.
   If you pushed it into ``develop``, that commit is going to stay there forever so you need to push another that undo the changes.

Testing and actions
*******************

Some tests are always run every time the repository is updated.
You can view ongoing and finished tests on the *Actions* tab in GitHub.

Testing consists of several layers:

.. card:: Build

   Tests that the code compiles and the Python package can be installed.
   The compilation is done within the Conda environment in ``environment-dev.yaml``.
   (For the MPI build we additionally use ``conda install openmpi``.)
   This test is always run and other tests are not even attempted if this one fails.

.. card:: Unit tests

   Tests that verify that individual parts of the code or some specific algorithms are working properly, e.g. when the user does X check that the code does Y.
   These tests don't verify the physics.
   Unit tests are always run and you can run them locally via (these tests create a temporary ``unittest.h5`` file).
   To run a specific test case:

   .. code-block:: bash

      cd a5py/testascot
      python unittest.py

   .. note::

      There are also unit tests for the C kernel but those are not up to date and thus are not currently included in testing.

.. card:: Physics tests

   Tests that verify that ASCOT5 models correctly the physics that it is expected to model.
   For example that the neoclassical transport in tokamaks is modelled properly.

   These test take about an hour to complete, so they are run only when develop is updated or pull request is made to main.
   If any of the tests fail, the simulation file ``testascot.h5`` is uploaded as an artifact and it can be downloaded from the workflow run in GitHub for the next 24 hours.
   This file can be used locally to plot the results:

   .. code-block:: bash

      cd a5py/testascot
      python physicstest.py

   Executing ``python`` when ``testascot.h5`` is not present runs the tests locally.

   .. note::

      Some physics tests fail occasionally when run in GitHub.
      Particularly the neoclassical transport test is prone to do that, and the printed values for the transport coefficients make no sense.
      Even downloading the file and plotting the results locally makes the test pass.
      The cause of this behaviour is unknown.

.. card:: Regression tests

   Tests that run complete simulations and verify that the results have not changed between the current and the previous versions.

   These tests require HPC resources so they are run on demand by the maintainers.

In addition to building and testing, there is a workflow that builds the documentation, and also publishes it if the commit was made to the ``docs`` branch.
Building the documentation involves running the tutorials, which may also act as tests, so the documentation is built also when ``develop`` is updated or pull request is made to ``main``.

Finally there is a workflow that rebases ``develop`` and ``docs`` to ``main`` whenever ``main`` is updated.
If there is a merge conflict, the maintainers have to do the rebasing manually.

Coding style
************

- Indentation is 4 spaces.
  - OpenMP pragmas indented same as normal code
  - Compiler pragmas are not indented.
- Use spaces; no tabs (except in Makefile).
- Maximum of 80 characters per line.
  Only exception are the ``*.rst`` files (where you should have a line per sentence).

.. dropdown:: Emacs settings to enforce the rules

   Copy-paste these to your ``~/.emacs`` file:

   .. code-block::

      ;; Standard style for C
      (setq c-default-style "linux" c-basic-offset 4)
      ;; Indent compiler pragmas
      (c-set-offset (quote cpp-macro) 0 nil)
      ;; Add column marker at line 80
      (require 'whitespace)
      (setq whitespace-style '(face empty tabs lines-tail trailing))
      (global-whitespace-mode t)
      ;; Indent switch-case statements
      (c-set-offset 'case-label '+)
      ;; Kill all tabs
      (setq-default indent-tabs-mode nil)
      ;; Bonus: makes emacs save all backup files (*~) in this directory
      (setq backup-directory-alist `(("." . "~/.emacssaves")))'

.. warning::

   If the aforementioned rules are violated, a sinkhole opens beneath you and you drop in a rat-filled pit with no escape.

In addition to the mandatory rules, we aim to uphold `PEP8 <https://peps.python.org/pep-0008/>`_ when writing Python.
For both Python and C, we have our own standard practices that are listed here.

1. Align variables and comments:

   .. code-block:: python

      real a =  0; # Setting a to zero
      int b  = -1; # Setting b to minus one
                   # because it is not zero

2. Leave empty space between brackets and math operations if it increases readibility:

   .. code-block:: python

      a = ( x + 3*y ) / np.pi # Correct
      a=(x+3*y)/3.141         # Wrong

3. No space between keyword and bracket; opening curl in the same line:

   .. code-block:: C

      /* Correct */
      if(a == 0) {
          // Do X
      }
      else {
          // Do Y
      }

      /* Wrong */
      if (a == 0)
      {
          // Do X
      } else
      {
          // Do Y
      }

When in doubt, look the surrounding code and mimic that.
But then again, these are minor issues that can be fixed when the pull request is submitted.

.. note::

   These custom Python conventions might be applied in the future, but they should then be applied on a file level at the minimum:

   - Type hints in function arguments and returned values.
   - Using double quotes ``"`` for strings (or anything that could contain natural language) and single quotes ``'`` for string literals e.g. ``color='red'``.

Documentation
*************

To build the documentation, make sure you are in a Conda environment specified by ``environment-dev.yaml`` and that you have installed ``a5py`` with optional ``doc`` packages as per the developer installation instructions.
Then the documentation (the one that you are now reading) is built with

.. code-block:: bash

   make doc

The output files are in ``build/doc`` and the main page is ``build/doc/index.html``.
Don't commit any output files since the documentation on the web is automatically generated with a GitHub action

.. rubric:: Comment blocks in the source files

The Python docstrings must conform to the `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`_ guidelines.
`This <https://developer.lsst.io/python/numpydoc.html>`_ is a very nice guide with examples of how the docstrings should look.

.. dropdown:: Example of Python docstring

   .. code-block:: python

      def sum(values):
          """Sum numbers in an array.

          Essentially a duplicate of :obj:`~np.sum`.

          Parameters
          ----------
          values : [float]
              Array whose values are summed.

          Returns
          -------
          sum : float
              Sum of ``values``.
          """
          return np.sum(values)

The documentation for the C code is generated using Doxygen, which `specifies how the comments should look <https://www.doxygen.nl/manual/docblocks.html>`_ (we are using the Javadoc style).

.. dropdown:: Example of C documentation

   .. code-block:: C

      /**
      * @file thisfile.c
      * @brief This file contains an example documentation.
      *
      * Longer description.
      *
      * Doxygen has a Latex support if needed:
      * \f$ \gamma = \sqrt{\frac{1}{1-v^2/c^2}}\f$
      *
      */

      int variable /**< Variable documentation */

      /**
      * @brief Short description of the function.
      *
      * @param arg Description of the parameter and also units if applicable.
      *
      * @return Description of the return value.
      */
      int myfun(arg) {
          ...
      }

.. rubric:: The documentation source

The documentation is built from the code source and the documentation source located in the ``doc`` folder using Sphinx.
Sphinx enables one to embed the comment blocks from the code to the actual documentation written in RST.
See the documentation Python docstrings are referenced via `autodoc <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_:

.. dropdown:: Example of embedding Python docstring

   Source:

   .. code-block:: rst

      .. automethod:: a5py.Ascot.input_eval

   Output:

   .. automethod:: a5py.Ascot.input_eval
      :noindex:

The Doxygen output, consisting of XML files, is linked to Sphinx via `Breathe <https://breathe.readthedocs.io/en/latest/>`_.
This allows one to reference C documentation as (see the Breathe manual for all directives and their options):

.. dropdown:: Example of embedding C comment block

   Source:

   .. code-block:: rst

      .. doxygenfunction:: simulate_gc_fixed

   Output:

   .. doxygenfunction:: simulate_gc_fixed

Code structure and development
******************************

For IDE, we recommend VScode which is free and with extensions work well with a multilingual code like ASCOT5.
Furthermore it can be used to run the Jupyter notebook tutorials more conveniently than managing the notebook server yourself.
For Windows machines you can use the VScode from Windows while ASCOT5 is being run in WSL (Windows Subsystem for Linux).

.. rubric:: Compiling in debug mode

The debug mode makes it easier to catch segmentation faults, memory leaks and similar that are usually caught with valgrind.
To run the code in debug mode:

.. code-block:: bash

   conda activate ascot-dev
   make clean
   make libascot DEBUG=1
   export LD_PRELOAD=$(gcc -print-file-name=libasan.so)
   python runascot.py
   export LD_PRELOAD= # Once you are done debugging

The debug mode uses AddressSanitizer which is more convenient than valgrind in hybrid Python-C codes such as ASCOT5.
Change ``gcc`` in this example to a different compiler depending on what you are using.

.. rubric:: How the code is structured

.. card:: General layout

   All pre- and post-processing is done by a single :class:`~.Ascot`.
   This class (as well as the executables) are for the most part the only thing that an user interacts with.
   Even GUI is built so that it only provides a canvas and plotting is done by the :class:`~.Ascot` object.

   Data is passed between ``a5py`` and the C-kernel either passively via the HDF5 file or actively via ``ascotpy`` + ``libascot.so``.

   .. mermaid::

      block-beta
          columns 3
	  User:1 space:2
	  space:3
	  Ascot["Ascot"]:1 space:1 GUI
	  space:3
	  block:width:3
              HDF5:1 space:1 libascot:1
	  end
	  space:3
	  kernel["The C-kernel"]:3

	  User-->Ascot
      GUI-->Ascot
      Ascot-->HDF5
      Ascot-->libascot
      HDF5-->kernel
      libascot-->kernel

.. card:: Input generation and evaluation

   This is how input is written in python

   .. mermaid::

      flowchart TB

      classDef routine font-family:Monospace
      classDef empty width:0px,height:0px;

      %% Define cells
      ID1[/"Input data"/]
      ID2[/"Input data (generated)"/]
      IP1[/"Template parameters"/]
      CI1[" "]:::empty
      CI2[" "]:::empty
      IG1("Input DataGroup\n write_hdf5\n ascot5io/&lt;input&gt;.py"):::routine
      Template("templates"):::routine
      HDF5[("HDF5")]
      HI("hdf5_interface.c\nhdf5io/*")
      SIM1(["Simulation"])
      INIH5("hdf5_&lt;input&gt;_init")
      INIOFF("&lt;input&gt;_init_offload")
      INI("&lt;input&gt;_init")
      A5PYINI1[" "]:::empty
      A5PYINI2[" "]:::empty
      coreio("AscotIO\n coreio/treeview.py\n coreio/fileapi.py"):::routine

      %% Subgraphs
      subgraph common1[" "]
      CI1
      common1name["AscotIO\ncreate_input"]:::routine
      CI2
      end

      subgraph common2[" "]
      A5PYINI1
      common2name("ascotpy\ninit_input"):::routine
      A5PYINI2
      end

      subgraph RunGroup
      OutputGroup
      end

      subgraph ascot5io[" "]
      coreio
      InputGroup
      RunGroup
      end

      %% Link
      HDF5-->coreio-.-InputGroup & OutputGroup
      ID1---A5PYINI2-->INIH5
      HDF5-..-A5PYINI1
      A5PYINI1-->INIOFF

      ID1---CI1---->IG1-->HDF5
      IP1---CI2-->Template-->ID2-->IG1
      SIM1-->HI-->HDF5
      HDF5-->INIH5-->INIOFF-->INI

      %% Tooltips
      click IP1 callback "Template-specific parameters which
      can be something as simple as filename."
      click HDF5 callback "ascot.h5
      ├─ bfield/efield/plasma/neutral/ #60;input parent#62;
      │    wall/boozer/mhd/asigma/nbi
      │    ├─ active (active group qid) #60;attribute#62;
      │    └─ #60;input type#62;_#60;qid#62;
      │         ├─ date #60;attribute#62;
      │         ├─ desc #60;attribute#62;
      │         └─ ...  #60;attribute#62;
      └─ results
           ├─ active (active group qid)
           └─ run/bbnbi/afsi
                ├─ date
                ├─ desc
                └─ inistate/endstate/orbit/ #60;diagnostics#62;
                     dist#60;dist#62;/transcoef
                      └─ ... #60;output data#62;"

.. rubric:: Roadmap

Here are some long-term goals on how ASCOT5 should look like in the future:

- No executables.
  Only ``libascot.so`` remains and all simulations are run via Python.

- HDF5 interface relocated from C to Python.

- Possibility to import inputs and run simulations via GUI.

- Code is packaged via Conda with a minimal loss in performance.
