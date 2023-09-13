ascot5_main:
	$(MAKE) -C src ascot5_main

libascot:
	$(MAKE) -C src libascot

doc:
	$(MAKE) -C src doc
	$(MAKE) -C doc

clean:
	$(MAKE) -C src clean
