
CC = gcc
CFLAGS = -O3 -msse4.2 -mpopcnt -march=corei7

all: read

read: read.o
	$(CC) $(LDFLAGS) -o $@ $^

read.o: read.c
	$(CC) $(CFLAGS) -c -o $@ $<