#!/bin/bash

if [[ "$mpi" == "nompi" ]]; then
  USEMPI=0
else
  USEMPI=1
fi

make libascot MPI=$(USEMPI) -j
cp build/libascot.so $PREFIX/lib/

make ascot5_main MPI=$(USEMPI) -j
make bbnbi5 MPI=$(USEMPI) -j
cp build/ascot5_main $PREFIX/bin/
cp build/bbnbi5 $PREFIX/bin/

$PYTHON -m pip install . --prefix=$PREFIX -vv