# file: Makefile
# nAPMO package 
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

TOPDIR=.
include $(TOPDIR)/config.make

SUBDIRS = napmo src

default: OMP

$(BUILD):
	/usr/bin/env python3 setup.py install --record files.txt --user

clean:
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done
	cat files.txt | xargs rm -rf
	rm -rf __pycache__ build dist napmo.egg-info files.txt

.PHONY: clean default
