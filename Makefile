UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC = /usr/bin/gcc
	CFLAGS = -std=gnu99 -DNDEBUG -O3 -fopenmp
endif
ifeq ($(UNAME), Darwin)
	CC = /usr/local/bin/gcc-8
	CFLAGS = -std=c99 -DNDEBUG -O3 -fopenmp
endif

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
