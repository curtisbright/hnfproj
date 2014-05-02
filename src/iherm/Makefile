
.PHONY: iherm

#OBJDIR  := objs/
OBJS    := #empty
OBJS    += $(OBJDIR)iherm.o
OBJS    += $(OBJDIR)iherm_rect.o
OBJS    += $(OBJDIR)iherm_utils.o
OBJS    += $(OBJDIR)hnfproj.o
OBJS    += $(OBJDIR)pk_matrix.o
OBJS    += $(OBJDIR)uni_inv.o

CFLAGS += -I$(SRCDIR)core
CFLAGS += -I$(SRCDIR)lift
CFLAGS += -I$(SRCDIR)spinv
CFLAGS += -I$(SRCDIR)linsys
CFLAGS += -I$(BASEDIR)include

iherm: $(OBJS)

$(OBJDIR)%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm $(OBJS)
