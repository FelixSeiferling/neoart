include Makefile-neoart.include

neoart: $(NEOART_OBJECTS) MAIN.o
	$(FC) $(FFLAGS) -o ../neoart.x $^

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f90 
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o ../neoart.x *.mod

