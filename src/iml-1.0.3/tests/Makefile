.PHONY: all clean

TGTS := #empty
TGTS += test-largeentry
TGTS += test-smallentry


all: $(TGTS)

%: %.c
	$(CC) $^ -L../lib -liml -lcblas -latlas -lgmp -o $@

clean:
	rm $(TGTS)
