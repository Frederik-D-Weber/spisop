#!/usr/bin/env make -f
# XMLTree Makefile called by {SPM}/src/Makefile
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 3261 2011-03-31 15:06:56Z roboos $

include ../../src/Makefile.var

SPMMEX = xml_findstr.$(SUF)

all: $(SPMMEX)

clean:

distclean: clean
	$(DEL) $(SPMMEX)

install:

tarball: all
	$(TAR) cf spm_mex.tar $(SPMMEX)

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)
