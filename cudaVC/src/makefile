# Generic Make file for cudaVC

ifeq ($(mode),debug)
   CFLAGS = -g -G -Xptxas -v -arch sm_30
else
   mode = release
   CFLAGS = -arch sm_30

endif

# Build tools
NVCC = nvcc

# Source File
SRC = main.cu

# objects
GPUOBJS = main.o 

# Executable
EXE = cudaVC

#CFLAGS = -arch sm_30
LFLAGS = -lm

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
	
	
# make and compile
$(EXE): $(GPUOBJS)
	$(NVCC) -o $(EXE) $(GPUOBJS) $(LFLAGS)

main.o: $(SRC)
	$(NVCC) -c $(SRC) $(CFLAGS)

.PHONY: clean
clean:
	rm -vf *.o $(EXE)
