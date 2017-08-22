

ASCOT5 Python library
=====================

`a5py` is a python library for pre- and postprocessing `ASCOT5` relevant 
data. The library contains tools to access and modify `HDF5` files used for 
`ASCOT5` input and output. Some higher level scripts for running `ASCOT5` are 
also provided.

This package is compatible with Python 3.X. Required libraries are `numpy`, `scipy`, 
`h5py`.

Websites
--------

* Source code: https://version.aalto.fi/gitlab/ascot/python/a5py
* ASCOT5: https://version.aalto.fi/gitlab/ascot/ascot5

Installation
------------

Make sure you have `pip` and `virtualenv` installed first (and `python3`). 
Begin by creating your own virtual environment:

>> virtualenv -p python3 --system-site-packages ascot
>> source ascot/bin/activate

Flag --system-site-packages makes your virtual environment inherit packages already
present in your system. This is necessary if you wish to do any plotting as installing
`matplotlib` is painful with pip.

Now you have your very own virtual environment. See what packages are present with:

>> pip list

If necessary, install the required packages:

>> pip install numpy
>> pip install scipy
>> pip install h5py

and then you can install `a5py`:

>> pip install /path/to/a5py

You can now import `a5py` as any library. Note that `a5py` contains scripts that can
be used directly from command line after the installation e.g.:

>> a5removeruns ascot.h5
   (removes all results from ascot.h5)

To exit your virtualenvironment, type

>> deactivate

For developers
--------------

Everytime you modify the source code, you need to re-install `a5py`. You can 
do this with

>> pip install --ignore-installed path/to/a5py

or

>> pip install -I path/to/a5py

Subpackages
-----------

- preprocessing   -- Generate input data from other formats
                     such as `EQDSK`.

- ascot5io        -- Read and write `ASCOT5` `HDF5` file. Modify existing 
                     HDF5 files (copy input fields etc.)

- posprocessing   -- Generate plots and analyze data

- ascot4interface -- Produce `ASCOT5` input from `ASCOT4` data

Scripts
-------

Higher level scripts for running and preparing `ASCOT5` are found in a5py/bin.
To add a new script, place it in the folder and add the script name setup.py
inside scripts=[] brackets.

By convention, each script should begin with a "a5".
