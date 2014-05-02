export BASEDIR := $(dir $(CURDIR)/$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))

export CC      := gcc
#export LDLIBS  := -liml -lopenblas -lgmp -lm
export LDLIBS  := -lcblas -latlas -lgmp -lm

export LDFLAGS := #empty
export CFLAGS  := #empty

export OBJDIR := $(BASEDIR)objs/
export SRCDIR := $(BASEDIR)src/

export LIFTLIB := $(BASEDIR)lib/libhnfproj.a
export SHAREDLIB := $(BASEDIR)lib/libhnfproj.so

ifdef ATLAS_LIB_DIR
  LDFLAGS += -L$(ATLAS_LIB_DIR)
endif
ifdef ATLAS_INCLUDE_DIR
  CFLAGS += -I$(ATLAS_INCLUDE_DIR)
endif
ifdef GMP_LIB_DIR
  LDFLAGS += -L$(GMP_LIB_DIR)
endif
ifdef GMP_INCLUDE_DIR
  CFLAGS += -I$(GMP_INCLUDE_DIR)
endif


CFLAGS  += -fPIC
CFLAGS  += -pedantic
CFLAGS  += -Wall
CFLAGS  += -Wextra
CFLAGS  += -Wshadow
CFLAGS  += -Wpointer-arith
CFLAGS  += -Wcast-align
CFLAGS  += -Wstrict-prototypes
CFLAGS  += -Wmissing-prototypes
CFLAGS  += -Wno-long-long
CFLAGS  += -Wno-variadic-macros

ifdef NOTIMER
  CFLAGS  += -DNOTIMER
endif
ifdef NOPRINT
  CFLAGS  += -DNOPRINT
endif

ifdef DEBUG
  CFLAGS  += -O0
  CFLAGS  += -g
  CFLAGS  += -DDEBUG
else
  CFLAGS  += -O3
  CFLAGS  += -DNDEBUG
endif

ifdef THREAD
  CFLAGS += -DTHREAD
  LDLIBS += -lpthread
endif

ifeq ($(CC),icc)
  CFLAGS += -wd1782 #pragma once is okay
endif

export MAKEOPTS := #empty
MAKEOPTS += --no-print-directory
