ASCOT5 Python library
=====================

`a5py` is a python library for pre- and postprocessing `ASCOT5` relevant
data. The library contains tools to access and modify `HDF5` files used for
`ASCOT5` input and output. Some higher level scripts for running `ASCOT5` are
also provided.

This package is compatible with Python 3.X. Required libraries are `numpy`,
`scipy`, `h5py`. Additionally, `matplotlib` is required for plotting routines
and `doxypypy` is needed for generating Doxygen documentation.

Websites
--------

* Source code: https://version.aalto.fi/gitlab/ascot/python/tree/master/a5py
* ASCOT5: https://version.aalto.fi/gitlab/ascot/ascot5

Installation
------------

Make sure you have `pip` and `virtualenv` installed first (and `python3`).
Begin by creating your own virtual environment:::

    virtualenv -p python3 --system-site-packages ascot
    source ascot/bin/activate

If you can't use/install virtualenv, you can use `venv` instead:::

    python3 -m venv --system-site-packages ascot
    source ascot/bin/activate

Flag --system-site-packages makes your virtual environment inherit packages
already present in your system. This is necessary if you wish to do any plotting
as installing `matplotlib` is not possible with just `pip`.

Now you have your very own virtual environment. See what packages are present
with:::

    pip list

Then install `a5py` and the requirements will be picked up from
`requirements.txt`:::

    pip install -e /path/to/a5py

If not, install the required packages:::

    pip install -r requirements.txt

The `-e` in the command above makes the package "editable", meaning links to
the source files are made in the site-packages folder rather than copies. This
means you can change and update the source files, and this will be immediately
reflected in the `a5py` Python package.

You can now import `a5py` as any library. Note that `a5py` contains scripts that
can be used directly from command line after the installation e.g.:::

    a5gui ascot.h5
    (Open GUI on ascot.h5)

To exit your virtual environment, type::

    deactivate

Updating a5py
--------------

If you did not use the `-e` flag to install `a5py` in your virtual environment,
then everytime the source code is modified, e.g. by pulling a new version from
git::

    git pull --rebase

or if you modify the source locally, `a5py` needs to be reinstalled. You can do
this with::

    pip install --ignore-installed path/to/a5py

or::

    pip install -I path/to/a5py

Note that if a file was removed from a5py, it still lingers in virtualenv. This
can be fixed by creating a new virtualenv.

Subpackages
-----------

- ascot5io        -- Read and write `ASCOT5` `HDF5` file. Modify existing
                     HDF5 files (copy input fields etc.) Use this module for
                     generating input and accessing output.

- marker          -- Module for evaluating marker specific quantities and
                     plotting markers. This module is mostly used by State
                     and Orbit objects from Ascot output.

- dist            -- Module for evaluating and plotting distributions. Used
                     by Ascot distribution output.

- gui             -- Graphical user interface for accessing Ascot5 HDF5 file.
                     In GUI, one can e.g. set descriptions and which inputs are
                     active. Also some tools for plotting input and output are
                     available.

- preprocessing   -- Generate input data from other formats such as `EQDSK`.

- postprocessing  -- Higher level data-analysis.

- ascot4interface -- Produce `ASCOT5` input from `ASCOT4` data.

Scripts
-------

Higher level scripts for running and preparing `ASCOT5` are found in a5py/bin.
To add a new script, place it in the folder and add the script name to setup.py
inside scripts=[] brackets.

By convention, each script should begin with an "a5".
