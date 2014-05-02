
.PHONY: spinv

#OBJDIR  := objs/
OBJS    := #empty
OBJS    += $(OBJDIR)highorder.o
OBJS    += $(OBJDIR)spinv.o
OBJS    += $(OBJDIR)unicert.o
OBJS    += $(OBJDIR)unicert_utils.o

CFLAGS += -I$(BASEDIR)include
CFLAGS += -I$(SRCDIR)core
CFLAGS += -I$(SRCDIR)lift
CFLAGS += -I$(SRCDIR)linsys

spinv: $(OBJS)

$(OBJDIR)%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm $(OBJS)
