SRC = .
FLINT_DIR = "/usr/local/include/flint"
IGRAPH_DIR = "/usr/local/include/igraph"
CFLAGS = -W -Wall -I$(FLINT_DIR) -I$(IGRAPH_DIR) -Ofast -march=native -fomit-frame-pointer
LDFLAGS = -lflint -lmpfr -lgmp -lpthread -ligraph
CC = gcc

OBJ1 = o/infmax.o o/aux.o o/stepdown.o o/get_inf_rec.o
OBJ2 = o/vacc.o o/aux.o o/stepdown.o o/get_inf_rec.o
OBJ3 = o/senti.o o/aux.o o/get_inf_rec.o

all : infmax vacc senti

infmax : $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

vacc : $(OBJ2)
	$(CC) -o $@ $^ $(LDFLAGS)

senti : $(OBJ3)
	$(CC) -o $@ $^ $(LDFLAGS)

o/infmax.o : $(SRC)/infmax.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/infmax.c -o $@

o/vacc.o : $(SRC)/vacc.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/vacc.c -o $@

o/senti.o : $(SRC)/senti.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/senti.c -o $@

o/aux.o : $(SRC)/aux.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/aux.c -o $@

o/stepdown.o : $(SRC)/stepdown.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/stepdown.c -o $@

o/get_inf_rec.o : $(SRC)/get_inf_rec.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/get_inf_rec.c -o $@
