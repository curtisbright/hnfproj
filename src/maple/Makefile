.PHONY: maple

OBJS    := #empty
OBJS    += $(OBJDIR)maple_call.o

CFLAGS += $(MAPLEOPTS)
CFLAGS += -I$(BASEDIR)include
CFLAGS += -I$(SRCDIR)core
CFLAGS += -I$(SRCDIR)linsys
CFLAGS += -I$(SRCDIR)spinv
CFLAGS += -I$(SRCDIR)iherm

maple: $(OBJS)

$(OBJDIR)%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

