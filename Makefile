.PHONY: doc, tests

libascot:
	$(MAKE) -C src libascot

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

cleanall: clean
	rm -rf build


.PHONY: doc
