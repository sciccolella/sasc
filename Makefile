UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CFLAGS = -std=gnu99 -DNDEBUG -O3
endif
ifeq ($(UNAME), Darwin)
	CFLAGS = -std=c99 -DNDEBUG -O3
endif
CC = /usr/bin/gcc

.PHONY: all

all: sasc

sasc: sasc.o mt19937ar.o sastep.o tree.o utils.o vector.o
	@echo "* Linking SASC"
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	@echo '* Compiling $<'
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o sasc
