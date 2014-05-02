
.PHONY: lift

#OBJDIR  := objs/
OBJS    := #empty
OBJS    += $(OBJDIR)imlsolve.o
OBJS    += $(OBJDIR)lift.o
OBJS    += $(OBJDIR)rns_conversion.o

CFLAGS += -I$(SRCDIR)core

lift: $(OBJS)

$(OBJDIR)%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm $(OBJS)
