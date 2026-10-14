ascot:
	$(MAKE) -C src ascot

lint:
	pylint --disable=similarities --generated-members=unyt.*
	mypy --follow-imports=silent --disable-error-code=import-untyped \
	--disallow-untyped-defs --disallow-incomplete-defs

tests:
	pytest -s -v -W "error::UserWarning" --ignore="tests/test_physics"

tutorial:
	python doc/user/tutorials/execute.py

doc-user:
	$(MAKE) -C doc PROJECT=user doc

doc-dev:
	$(MAKE) -C doc PROJECT=dev doc

clean:
	$(MAKE) -C src clean

cleandoc:
	$(MAKE) -C doc cleanall

cleanall: clean
	$(MAKE) -C doc cleanall
	rm -rf build

.PHONY: ascot lint tests tutorial doc-user doc-dev doc clean cleandoc cleanall
