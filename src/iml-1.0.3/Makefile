.PHONY: clean default

SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c, $(OBJDIR)%.o, $(SRCS))

default: $(OBJS)

$(OBJDIR)%.o: src/%.c
	$(CC) $< -c $(CFLAGS) -o $@

clean:
	rm -f lib/libiml.a objs/*.o; true

