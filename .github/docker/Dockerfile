# This container should include all the necessary libraries and tools
# to compile and use ascot5. It can be used e.g. as the image for CI testing.

FROM ubuntu:20.04

ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN yes | apt-get update
RUN yes | apt-get install libhdf5-openmpi-dev
RUN yes | apt-get install libhdf5-dev
RUN yes | apt-get install h5utils
RUN yes | apt-get install emacs
RUN yes | apt-get install make
RUN yes | apt-get install gcc
RUN yes | apt-get install libclang-12-dev
RUN yes | apt-get install pip
RUN yes | pip install ctypeslib2 clang==12.0.1
RUN yes | apt-get install git
RUN yes | apt-get install -y doxygen
RUN yes | apt-get install pandoc
COPY requirements.txt requirements.txt
RUN yes | pip install -r requirements.txt
RUN yes | pip install ipykernel
RUN yes | pip install sphinx
RUN yes | pip install numpydoc
RUN yes | pip install nbsphinx
RUN yes | pip install breathe
RUN yes | pip install sphinxcontrib-bibtex
RUN yes | pip install sphinx-rtd-theme
RUN yes | pip install sphinx-gallery
RUN yes | pip uninstall vtk -y
RUN yes | pip install vtk-osmesa --extra-index-url https://wheels.vtk.org