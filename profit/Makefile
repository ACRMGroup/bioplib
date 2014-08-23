# You are advised to change this to provide warning messages for your
# complier - the sample line is for GCC
#CC    = cc -O3
#CC    = cc -O3 -ansi -Wall -pedantic 
CC = gcc


# To allow use of the GNU Readline library, uncomment the following two
# lines - You may need to install the GNU readline development libraries
# first!
#READLINE    = -DREADLINE_SUPPORT
#READLINELIB = -lreadline -lcurses

# Uncomment if you want to use the rotate and refit code for avoiding 
# local minimum problem. If compiled with gcc -O3 this problem seems to
# go away in any case and only affects fitting of identical structures
# A fix will soon be provided in fit.c instead
#ROTATEREFIT = -DROTATE_REFIT

# Compile options
COPT = 
COPT := $(COPT) -O3
#COPT := $(COPT)  -ansi
COPT := $(COPT) -fPIC -Wall -pedantic -g

# Decompression
COPT  := $(COPT)  -DGUNZIP_SUPPORT

# Link to BiopLib
COPT    := $(COPT) -I ${HOME}/include -L ${HOME}/lib
BIOPLIB  = -lbiop -lgen -lxml2



# Uncomment these and fix the paths if you have the XMAS library and
# want XMAS support
# XMAS  = -DUSE_XMAS  -I/home/amartin/inpharmatica/cvs/software/include
# XMASLIB = -L/home/amartin/inpharmatica/cvs/software/lib -lxmas -lacrm

LINK1 = 
LINK2 =
ANSI  = ansi -p 

OFILES = main.o todo.o fitting.o NWAlign.o
PFILES = main.p todo.p fitting.p NWAlign.p
IFILES = ProFit.h protos.h

profit : $(OFILES)
	$(CC) $(COPT) -o profit $(OFILES) $(XMAS) -lm $(LINK2) $(XMASLIB) $(READLINELIB) $(BIOPLIB)

.c.o : $(IFILES)
	$(CC) $(COPT) $(READLINE) $(ROTATEREFIT) -o $@ -c $<

clean :
	/bin/rm $(OFILES) profit

protos : $(PFILES)

.c.p   :
	$(ANSI) $< $@
