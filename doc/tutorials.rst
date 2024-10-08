.. _Tutorial:

=========
Tutorials
=========

These tutorials are implemented as Jupyter notebooks meaning that, while they appear as ordinary HTML documents here, you can also run them interactively.
The notebooks are found in ``doc/tutorials`` folder in your cloned repository.
To run them interactively, you need to have Jupyter installed and then configure it to use the virtual environment that you've created when installing ASCOT5:

.. code-block::

   pip install jupyter
   ipython kernel install --user --name=ascotenv

Now you can run the examples interactively with e.g. VS Code, or even with your browser (remember to select the ``ascotenv`` kernel),

.. code-block::

   jupyter-notebook

or you can convert these examples to ordinary Python scripts and execute them (executing them in ipython terminal is advised as otherwise the plotted figures close immediately):

.. code-block::

   jupyter-nbconvert --to python <name of the notebook e.g. tutorial.ipynb>
   python tutorial.py

.. _Examples:

.. rubric:: Basics

.. nbgallery::
   :hidden:

   tutorials/introduction.ipynb
   tutorials/orbits.ipynb
   tutorials/distributions.ipynb

.. rubric:: Marker generation

.. nbgallery::
   :hidden:

   tutorials/afsi.ipynb
   tutorials/bbnbi.ipynb
   tutorials/markergen.ipynb

.. rubric:: Physics

.. nbgallery::
   :hidden:

   tutorials/slowingdown.ipynb
   tutorials/wallload.ipynb
   tutorials/mhd.ipynb
   tutorials/atomic.ipynb
   tutorials/reversetime.ipynb

.. rubric:: Tools

.. nbgallery::
   :hidden:

   tutorials/poincare.ipynb
   tutorials/orbitanalysis.ipynb
   tutorials/biosaw.ipynb
