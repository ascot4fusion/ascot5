ascot:
	$(MAKE) -C src ascot

ascot-debug:
	$(MAKE) -C src ascot DEBUG=1

# https://stackoverflow.com/questions/75847374/using-memory-sanitizer-asan-on-c-c-library-loaded-to-python-with-ctypes
# Run export LD_PRELOAD=$(cc -print-file-name=libasan.so) beforehand
# when using AddressSanitizer as python does not compile with it
ascot-sanitizer:
	$(MAKE) -C src ascot SANITIZER=1

lint:
	pylint --disable=similarities --generated-members=unyt.*
	mypy --follow-imports=silent --disable-error-code=import-untyped \
	--disallow-untyped-defs --disallow-incomplete-defs

tests:
	python -m pytest -s -v -W "error::UserWarning" --ignore="tests/test_physics"

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
