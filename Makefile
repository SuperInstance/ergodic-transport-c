CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2
LDFLAGS = -lm

SRC = src/markov.c src/ergodic.c src/budget.c src/wasserstein_control.c
OBJ = $(SRC:.c=.o)

.PHONY: all test clean

all: libergodictransport.a test_ergodic_transport

%.o: %.c include/ergodic_transport.h
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@

libergodictransport.a: $(OBJ)
	ar rcs $@ $^

test_ergodic_transport: tests/test_ergodic_transport.c libergodictransport.a
	$(CC) $(CFLAGS) -Iinclude $< -L. -lergodictransport $(LDFLAGS) -o $@

test: test_ergodic_transport
	./test_ergodic_transport

clean:
	rm -f $(OBJ) libergodictransport.a test_ergodic_transport
