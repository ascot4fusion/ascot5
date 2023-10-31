# This container should include all the necessary libraries and tools
# to compile ascot5. It can be used e.g. as the image for CI testing.

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

# Copies your code file from your action repository to the filesystem
# path `/` of the container
COPY entrypoint.sh /entrypoint.sh

# Code file to execute when the docker container starts up (`entrypoint.sh`)
ENTRYPOINT ["/entrypoint.sh"]