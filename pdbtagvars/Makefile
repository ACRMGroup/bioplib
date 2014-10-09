CFLAGS = -g -ansi -pedantic -Wall -I$(HOME)/include
CFLAGS := $(CFLAGS) $(shell xml2-config --cflags)
LFLAGS = -L$(HOME)/lib -lbiop -lgen -lxml2
OFILES = main.o pdbtagvars.o 

mytest : $(OFILES)
	cc $(CFLAGS) -o $@ $(OFILES) $(LFLAGS)

clean :
	\rm -f $(OFILES)
