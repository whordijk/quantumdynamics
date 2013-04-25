FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O3 -ffast-math
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LDLFLAGS =
LIBS = -llapack
LIBS += -lfftw3
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS = 
OBJS += cranknicolson.o
OBJS += splitoperator.o
OBJS += wavepacket.o
OBJS += QuantumDynamics.o

all: QuantumDynamics

QuantumDynamics: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)
%.o: %.f90
	$(COMPILE) -o $@ -c $<
clean:
	$(RM) $(OBJS) *.mod
