
.PHONY: all 

UNICERT_EXE := $(BASEDIR)bin/unicert_main
IHERM_EXE := $(BASEDIR)bin/iherm_main
HOR_EXE := $(BASEDIR)bin/hor_main
SIMPLE_IHERM_EXE := $(BASEDIR)bin/simple_iherm

CFLAGS += -I$(BASEDIR)include

all: $(UNICERT_EXE) $(IHERM_EXE) $(SIMPLE_IHERM_EXE) $(HOR_EXE)

$(UNICERT_EXE): unicert_main.c $(LIFTLIB)
	$(CC) $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS) -I$(SRCDIR)core -I$(SRCDIR)spinv -o $@

$(HOR_EXE): hor_main.c $(LIFTLIB)
	$(CC) $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS) -I$(SRCDIR)core -I$(SRCDIR)spinv -o $@

$(SIMPLE_IHERM_EXE): iherm_simple.c $(SHAREDLIB)
	$(CC) $< $(CFLAGS) $(LDFLAGS) -L$(BASEDIR)lib -lhnfproj $(LDLIBS) -o $@

$(IHERM_EXE): iherm_main.c $(LIFTLIB)
	$(CC) $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS) -I$(SRCDIR)core -I$(SRCDIR)linsys -I$(SRCDIR)spinv -I$(SRCDIR)iherm -o $@

clean:
	rm $(UNICERT_EXE) $(IHERM_EXE) $(SIMPLE_IHERM_EXE)
