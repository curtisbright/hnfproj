.PHONY: core lift maple spinv iherm clean

lib: $(LIFTLIB)

shared: $(SHAREDLIB)

$(LIFTLIB): $(OBJDIR)*.o
	$(AR) -r -c -s $@ $^

$(SHAREDLIB): $(OBJDIR)*.o
	$(CC) $^ $(CFLAGS) $(MAPLEOPTS) -shared $(LDFLAGS) $(LDLIBS) -o $@

$(OBJDIR)*.o: iml core lift spinv iherm
ifdef MAPLEDIR
  $(OBJDIR)*.o: maple
endif

CFLAGS += -I$(SRCDIR)iml-1.0.3/include

iml:
	$(MAKE) $(MAKEOPTS) --directory=iml-1.0.3

core:
	$(MAKE) $(MAKEOPTS) --directory=core

lift:
	$(MAKE) $(MAKEOPTS) --directory=lift

spinv:
	$(MAKE) $(MAKEOPTS) --directory=spinv

iherm:
	$(MAKE) $(MAKEOPTS) --directory=iherm

maple:
	$(MAKE) $(MAKEOPTS) --directory=maple

clean:
	$(MAKE) $(MAKEOPTS) --directory=iml-1.0.3 clean
