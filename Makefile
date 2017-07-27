SRC = .
FLINT_DIR = "/usr/local/include/flint"
CFLAGS = -W -Wall -I$(FLINT_DIR) -Ofast -march=native
LDFLAGS = -lflint -lmpfr -lgmp -lpthread
CC = gcc

OBJ1 = o/infmax.o o/aux.o
OBJ2 = o/vacc.o o/aux.o
OBJ3 = o/senti.o o/aux.o

all : infmax vacc senti

infmax: $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

vacc: $(OBJ2)
	$(CC) -o $@ $^ $(LDFLAGS)

senti: $(OBJ3)
	$(CC) -o $@ $^ $(LDFLAGS)

o/infmax.o : $(SRC)/infmax.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/infmax.c -o $@

o/vacc.o : $(SRC)/vacc.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/vacc.c -o $@

o/senti.o : $(SRC)/senti.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/senti.c -o $@

o/aux.o : $(SRC)/aux.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/aux.c -o $@
