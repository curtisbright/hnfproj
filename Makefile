.PHONY: default all clean lib

include defs.mk

default: examples lib

all: lib shared examples tuning version

##################

lib:
	$(MAKE) $(MAKEOPTS) --directory=src lib

shared:
	$(MAKE) $(MAKEOPTS) --directory=src shared

examples: lib shared
	$(MAKE) $(MAKEOPTS) --directory=examples

tests: lib
	$(MAKE) $(MAKEOPTS) --directory=tests tests

version: lib
	$(MAKE) $(MAKEOPTS) --directory=tests version

tuning: lib
	$(MAKE) $(MAKEOPTS) --directory=tests tuning

clean:
	$(MAKE) $(MAKEOPTS) --directory=src clean
	rm -f objs/* lib/* bin/*; true;

