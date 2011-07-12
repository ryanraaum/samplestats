CC=gcc
CFLAGS=-O2
LFLAGS=-lm
OBJECTS=sample_stats3.o tajd.o fs.o r2.o simple_getopt.o
EXECUTABLE=sample_stats3

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LFLAGS) -o $@ $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o

clobber: clean
	rm -f $(EXECUTABLE)
