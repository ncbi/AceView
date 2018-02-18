
CFLAGS=-g -pedantic -Wall

all:	test_sleep test_pi test_pi2

wego.o: wego.c wego.h

test_sleep:	test_sleep.o wego.o
	cc -g -o test_sleep test_sleep.o wego.o -lpthread

test_pi: test_pi.o wego.o
	cc -g -o test_pi test_pi.o wego.o -lpthread -lm

test_pi2: test_pi2.o wego.o
	cc -g -o test_pi2 test_pi2.o wego.o -lpthread -lm

test_sleep.o:	test_sleep.c wego.h

clean:
	rm -f *.o test_sleep test_pi test_pi2


