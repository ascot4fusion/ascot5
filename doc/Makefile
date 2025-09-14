PROJECT ?= user

BASEBUILDDIR = ../build/doc/
SOURCEDIR = .
BUILDDIR  = $(BASEBUILDDIR)$(PROJECT)

.PHONY: doc clean

doc:
	sphinx-build -b html $(SOURCEDIR) $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR)

cleanall:
	rm -rf $(BASEBUILDDIR)
