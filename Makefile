.PHONY: doc, tests

libascot:
	$(MAKE) -C src libascot
	mkdir -p build
	mv src/libascot.so build/libascot.so

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
	rm -rf build/libascot.so

cleanall: clean
	rm -rf build
