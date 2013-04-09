FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O3
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LDLFLAGS =
LIBS = -llapack
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS = 
OBJS += BiCGSTAB.o
OBJS += wave1d.o
OBJS += QuantumDynamics.o

all: QuantumDynamics

QuantumDynamics: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)
%.o: %.f90
	$(COMPILE) -o $@ -c $<
clean:
	$(RM) $(OBJS) *.mod
