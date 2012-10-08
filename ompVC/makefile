# Generic makefile for flowVC
#
#
ifeq ($(mode),debug)
   CFLAGS = -g -Wall -O0 "-fopenmp" 
else
   mode = release
   CFLAGS =  -O3 "-fopenmp" 
endif
LFLAGS =  -lm "-fopenmp" 
CC = gcc
SRC = exposuretime.c ftle.c globals.c integration.c io.c main.c memory.c mesh.c mymath.c parameters.c \
	residencetime.c strainrate.c tracers.c velocity.c velout.c
OBJ = $(SRC:.c=.o)
EXE = ompflowVC

.PHONY: all
all: info $(EXE)

info:
ifneq ($(mode),release)
ifneq ($(mode),debug)
    @echo "Invalid build mode." 
    @echo "Please use 'make mode=release' or 'make mode=debug'"
    @exit 1
endif
endif
	@echo "Building in "$(mode)" mode..."

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -vf *.o $(EXE)


