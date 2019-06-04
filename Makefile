METHOD ?= BAOAB
NN ?= 1
NRESPA ?= 10

FORT = /usr/bin/gfortran
FOPTS = -Ofast -march=native -ffast-math -cpp -DMETHOD=$(METHOD) -DNN=$(NN) -DNRESPA=$(NRESPA)

EXEC = harmonic-$(METHOD)-n$(NN)-respa$(NRESPA)

.PHONY: all clean

all: $(EXEC)

clean:
	rm -f $(EXEC)

clean-all:
	rm -f $(EXEC) *.o

$(EXEC): harmonic.o mRandom.o mThermostat.o
	$(FORT) $(FOPTS) -o $@ $^

harmonic.o: harmonic.f90 mRandom.o mThermostat.o
	$(FORT) $(FOPTS) -c -o $@ $<

mRandom.o: mRandom.f90
	$(FORT) $(FOPTS) -c -o $@ $<

mThermostat.o: mThermostat.f90
	$(FORT) $(FOPTS) -c -o $@ $<
