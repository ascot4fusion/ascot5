PROJECT ?= user

BASEBUILDDIR = ../build/doc/
SOURCEDIR = .
BUILDDIR  = $(BASEBUILDDIR)$(PROJECT)

.PHONY: doc clean

doc:
	doxygen Doxyfile
	sphinx-build --nitpicky -b html $(SOURCEDIR) $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR)

cleanall:
	rm -rf $(BASEBUILDDIR)
