.. _Codingstyle:

==========================
Coding style and practices
==========================

- Indentation is 4 spaces.
  - OpenMP pragmas indented same as normal code
  - Compiler pragmas are not indented.
- Don't use tabs.
- Stick to the 80-column rule.

Copy paste these settings to your ``~/.emacs`` file to enforce the rules:

.. code-block::

   ;; Standard style for C
   (setq c-default-style "linux" c-basic-offset 4)
   ;; Makes emacs save all backup files in this directory
   (setq backup-directory-alist `(("." . "~/.emacssaves")))
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

Git usage
=========

`Git usage <https://about.gitlab.com/images/press/git-cheat-sheet.pdf>`_

- For all features and bugs, set up an issue on the issue tracker.
- Create a new branch for the issue

  - For a feature, name the branch as "feature/<issue number>_<issue name>" and branch from ``develop``.
  - For a bug, name the branch as "hotfix/<issue number>_<issue name>" and branch from ``master`` (unless of course the bug is present only in ``develop``).

- Automatic tests are run after each commit.

  - See Testing
  - You can review the process in GitLab at ``CI/CD->Pipelines``.
  - There is a "Download artifacts" button where you can access the output produced by the tests.

- Only permitted developers can merge to ``master`` as each merge corresponds to a new version.

Documentation
=============

The C code is documented with the help of Doxygen.
Here's a sample of how the code blocks should look:

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

The Python docstrings must conform to the `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`_ guidelines.
`This <https://developer.lsst.io/python/numpydoc.html>`_ is a very nice guide with examples of how the docstrings should look, here's an example:

.. code-block:: python

   def sum(values):
       """Sum numbers in an array.

       Parameters
       ----------
       values : iterable
           Python iterable whose values are summed.

       Returns
       -------
       sum : `float`
           Sum of ``values``.
       """
       pass

Generating documentation is a two step process.

1. In ascot5 folder, the Doxygen documentation from C source files is compiled as

   .. code-block:: bash

      doxygen Doxyfile

   The output consists of xml files stored in ``docs/_static/xml`` and equivalent html files in ``docs/build/capi``.

2. The actual documentation (the one you are reading) is compiled with Sphinx.

   .. code-block:: bash

      cd docs
      make clean html
      firefox build/index.html

Note however, that the documentation is generated automatically when the code is pushed to ``master`` and uploaded to GitLab pages.
Therefore you only have to compile the documentation to check that your modifications look as you expect them to look.
Only commit source files to the repository.

The Doxygen output is linked to Sphinx via `Breathe <https://breathe.readthedocs.io/en/latest/>`_.
The xml files allow one to reference C documentation as (see the Breathe manual for all directives and their options):

.. code-block:: rst

   .. doxygenfunction:: B_field_eval_B

which produces this output:

.. doxygenfunction:: B_field_eval_B

However, Breathe cannot automatically construct the whole C API (as Doxygen does), especially since we have documentation in both ``.c`` and ``.h`` files, which is why we must link to the Doxygen html files directly in C API.

Python docstrings are referenced via `autodoc <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_:

.. code-block:: rst

   .. automethod:: a5py.Ascot.input_eval
      :noindex:

which produces this output:

.. automethod:: a5py.Ascot.input_eval
   :noindex:

===============
Parallelization
===============

TODO

=======
Testing
=======

TODO

=====================================
`C API <_static/doxygen/index.html>`_
=====================================
