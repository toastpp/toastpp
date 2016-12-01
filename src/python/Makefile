include $(TOASTDIR)/Makefile.incl
SHELL = /bin/csh
VPATH = $(TSRC)/python

all:: python doc

python::
	python setup.py build
	python setup.py install --prefix=$(TOASTVER)

doc::
	cd doc; \
	make html
