ascot:
	$(MAKE) -C src ascot

#ascot: libascot
#	pip install -e .

lint:
	pylint --disable=similarities --generated-members=unyt.*
	mypy --follow-imports=silent --disable-error-code=import-untyped \
	--disallow-untyped-defs --disallow-incomplete-defs

tests:
	pytest -s -v -W "error::UserWarning"

tutorial:
	@echo TODO

doc-user:
	$(MAKE) -C doc PROJECT=user doc

doc-dev:
	$(MAKE) -C doc PROJECT=dev doc

doc: doc-user doc-dev

clean:
	$(MAKE) -C src clean

clean-doc:
	$(MAKE) -C doc cleanall

cleanall: clean
	rm -rf build

.PHONY: libascot ascot5 lint tests tutorial doc-user doc-dev doc clean \
		clean-doc cleanall
