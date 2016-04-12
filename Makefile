# file: Makefile
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

TOPDIR=.
include $(TOPDIR)/config.make

SUBDIRS = napmo src

default: SERIAL

$(BUILD):
	cd src && $(MAKE) $@
	python setup.py install --record files.txt --user
	python2 setup.py install --record files.txt --user

clean:
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done
	cat files.txt | xargs rm -rf
	rm -rf __pycache__ build dist napmo.egg-info files.txt

.PHONY: clean default
