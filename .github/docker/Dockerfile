# This container should include all the necessary libraries and tools
# to compile and use ascot5. It can be used e.g. as the image for CI testing.
FROM ubuntu:20.04

COPY ../../requirements.txt requirements.txt

RUN apt-get update && apt-get install -y \
    libhdf5-openmpi-dev libhdf5-dev h5utils emacs make gcc libclang-12-dev \
    git doxygen pandoc pip

RUN pip install -y ctypeslib2 clang==12.0.1
RUN pip install -y -r requirements.txt
RUN pip install -y ipykernel sphinx numpydoc nbsphinx breathe \
    sphinxcontrib-bibtex sphinx-rtd-theme sphinx-gallery

RUN pip uninstall vtk -y
RUN install vtk-osmesa --extra-index-url https://wheels.vtk.org