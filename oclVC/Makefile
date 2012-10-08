##############################################
# 	Make file for OpenCL code of flowVC		##
##############################################

ifeq ($(mode),debug)
   CFLAGS = -g -Wall -O0 -std=gnu99
else
   mode = release
   CFLAGS = -Wall -O3 -std=gnu99 

endif



LFLAGS = -lOpenCL -lm
CC = gcc

SRC =  globals.c openclfuncs.c main.c error.c parameters.c mesh.c allocate.c clean.c memcopy.c tracers.c fileoutput.c GlobalSearch.c \
		velocity.c index.c mymath.c integration3D.c
OBJ = $(SRC:.c=.o)
EXE = OpenVC

INC = -I ~/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc/ 
LIBS = -L ~/NVIDIA_GPU_Computing_SDK/OpenCL/common/lib/Linux64/

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
	$(CC) $(OBJ) $(LIBS) $(LFLAGS)  -o $@

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: clean
clean:
	rm -vf *.o $(EXE)



