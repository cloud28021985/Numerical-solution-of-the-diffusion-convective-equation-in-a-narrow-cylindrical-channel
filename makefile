# nonlinear theory of macroscopic flow induced in a drop of ferrofluid
# makefile; instruction for compiling, linking and running the program


СС = cc
CFLAGS = -c -Ofast -Wall
SRCPATH = src/
OBJPATH = obj/
DATAPATH = data/
FIGPATH = figs/
SOURCES = $(wildcard $(SRCPATH)*.c)
OBJECTS = $(patsubst $(SRCPATH)%.c, $(OBJPATH)%.o, $(SOURCES))
EXECUTABLE = $(OBJPATH)prog


all:
	$(EXECUTABLE)
	python3 plot.py


all: $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -lm -o $@


$(OBJPATH)%.o: $(SRCPATH)%.c
	$(CC) $(CFLAGS) $< -o $@


$(OBJECTS): $(SRCPATH)header.h makefile


# command $make clean
clean:
	rm -rf $(OBJPATH)
	rm -rf $(DATAPATH)
	rm -rf $(FIGPATH)
	mkdir $(OBJPATH)
	mkdir $(DATAPATH)
	mkdir $(DATAPATH)vector_field
	mkdir $(FIGPATH)
