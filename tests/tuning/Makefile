
.PHONY: all clean

CFLAGS += -I$(BASEDIR)include
CFLAGS += -I$(SRCDIR)core
CFLAGS += -I$(SRCDIR)lift
CFLAGS += -I$(SRCDIR)linsys
CFLAGS += -I$(SRCDIR)spinv
CFLAGS += -I$(SRCDIR)iherm

TGTS := #empty
#TGTS += iherm
TGTS += gemm
TGTS += rns_conv
TGTS += pk_gemm

all: $(TGTS)

% : %.c
	$(CC) $? $(LIFTLIB) $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o $(BASEDIR)/bin/$@

clean:
	rm $(TGTS)

