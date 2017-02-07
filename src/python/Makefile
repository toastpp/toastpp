include $(TOASTDIR)/Makefile.incl
SHELL = /bin/csh
VPATH = $(TSRC)/python

all:: python doc

python::
	python setup.py build --build-base=$(TOASTVER)/python
	python setup.py install --prefix=$(TOASTVER)
#	python setup.py install

doc::
	cd doc; \
	make html
