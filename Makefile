.PHONY: doc

ascot5_main:
	$(MAKE) -C src ascot5_main
	mkdir -p build
	mv src/ascot5_main build/ascot5_main

libascot:
	$(MAKE) -C src libascot
	mkdir -p build
	mv src/libascot.so build/libascot.so

doc:
	$(MAKE) -C src doc
	$(MAKE) -C doc

clean:
	$(MAKE) -C src clean
	rm -rf build
