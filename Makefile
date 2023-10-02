.PHONY: doc

ascot5_main:
	$(MAKE) -C src ascot5_main
	mkdir -p build
	mv src/ascot5_main build/ascot5_main

libascot:
	$(MAKE) -C src libascot
	mkdir -p build
	mv src/libascot.so build/libascot.so

bbnbi5:
	$(MAKE) -C src bbnbi5
	mkdir -p build
	mv src/bbnbi5 build/bbnbi5

ascot2py.py:
	$(MAKE) -C src ascot2py.py
	python .setcdllascot2py.py
	mv src/ascot2py.py a5py/ascotpy/ascot2py.py

doc:
	$(MAKE) -C src doc
	$(MAKE) -C doc

clean:
	$(MAKE) -C src clean
	rm -rf build
