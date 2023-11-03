# This container should include all the necessary libraries and tools
# to compile and use ascot5. It can be used e.g. as the image for CI testing.

FROM ubuntu:20.04

ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    libhdf5-openmpi-dev libhdf5-dev h5utils emacs make gcc libclang-12-dev git \
    libgl1 python3-pip
RUN pip3 install ctypeslib2 clang==12.0.1 --no-cache-dir
COPY requirements.txt .
RUN pip3 install -r requirements.txt
