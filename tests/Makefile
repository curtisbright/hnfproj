
.PHONY: tests run tuning version

include ../defs.mk

TEST_SRCS := $(wildcard test_src/*.c)

CFLAGS += -I$(BASEDIR)include
CFLAGS += -I$(SRCDIR)core
CFLAGS += -I$(SRCDIR)lift
CFLAGS += -I$(SRCDIR)linsys
CFLAGS += -I$(SRCDIR)spinv
CFLAGS += -I$(SRCDIR)iherm

tests: run
run : $(TEST_SRCS) $(LIFTLIB)
	$(CC)  $? $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o run
	./run

tuning:
	$(MAKE) $(MAKEOPTS) --directory=tuning

VERSION_EXE := $(BASEDIR)bin/version
version: $(VERSION_EXE)
$(VERSION_EXE): version.c
	$(CC) $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm run

