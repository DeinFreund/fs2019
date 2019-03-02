#CFLAGS = -O2
CFLAGS = -O3 -Wall -ftree-vectorize -fopt-info-all -march=native -msse2 -D NOALIAS
#-ffast-math
#-D NOFUNCCALL -D NOALIAS
#CFLAGS = -O0 -g
LFLAGS = 

BINARIES = heat2d

.SECONDARY: 
.PHONY: all 
all: $(BINARIES)

heat2d: auxiliar/auxiliar.o heat2d.o 
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $^

auxiliar/auxiliar.o:
	$(MAKE) -C auxiliar

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< 

.PHONY: clean
clean:
	$(MAKE) -C auxiliar clean
	$(RM) $(BINARIES) *.o *.ti *.optrpt 


