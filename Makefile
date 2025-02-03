.PHONY: doc, tests

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
	python3 .setcdllascot2py.py
	mv src/ascot2py.py a5py/ascotpy/ascot2py.py

ascot5: libascot
	pip install -e .

lint:
	pylint --generated-members=unyt.*
	mypy --follow-imports=silent --disable-error-code=import-untyped

tests:
	pytest -s -v

doc:
	$(MAKE) -C src doc
	$(MAKE) -C doc

clean:
	$(MAKE) -C src clean
	rm -rf build
