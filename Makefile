
CC = gcc
CFLAGS = -g -msse4.2 -mpopcnt -march=corei7 -Wall -pedantic -std=c99

all: read

read: read.o
	$(CC) $(LDFLAGS) -o $@ $^

read.o: read.c
	$(CC) $(CFLAGS) -c -o $@ $<