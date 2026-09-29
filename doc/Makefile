PROJECT ?= user

BASEBUILDDIR = ../build/doc/
SOURCEDIR = .
BUILDDIR  = $(BASEBUILDDIR)$(PROJECT)
DOXYGENOUTPUT = _static/doxygen/

.PHONY: doc clean

doc:
	#doxygen Doxyfile
	#sphinx-build --nitpicky -b html $(SOURCEDIR) $(BUILDDIR)
	sphinx-build -E -a --nitpicky -b html -v $(SOURCEDIR) $(BUILDDIR)

pdf:
	doxygen Doxyfile
	sphinx-build --nitpicky -b latex $(SOURCEDIR) $(BUILDDIR)/latex
	$(MAKE) -C $(BUILDDIR)/latex

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(DOXYGENOUTPUT)

cleanall:
	rm -rf $(BASEBUILDDIR)
	rm -rf $(DOXYGENOUTPUT)
