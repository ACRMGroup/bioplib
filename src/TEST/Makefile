# Simple makefile for building Bioplib tests

# Define C compiler
CC = gcc

# Options for the C compiler
COPT = -fPIC -ansi -Wall -pedantic -g

# Add /usr/local to search paths (not required for unix systems)
COPT := $(COPT) -I /usr/local/include -L /usr/local/lib

# Link to libxml2 library.
# 
# These option are required for xml support and must be included if
# BiopLib has been compiled with the '-D XML_SUPPORT' option.
# These options must be omitted if libxml2 is not installed.
#
# Note:
#  The program, xml2-config, is installed with libxml2 and will provide
#  -I -L and -l options for the compiler. Libxml2 does not require
#  additional libraries however, xml2-config will provide links for some
#  additional libraries if they are present.
#  Use -lxml2 to link solely to libxml2.

# Libxml2 include file path (-I option):
XML_OPT = $(shell xml2-config --cflags)

# Link to libxml2 only:
#XML_LIB = -lxml2

# Link to libxml2 and additional libraries (-L and -l options):
XML_LIB = $(shell xml2-config --libs)


# Test source code
TEST_SRC = src/*.c

# Bioplib object files
BIOP_OBJ = ../*.o


# Compile tests
tests : 
	$(CC) $(COPT) -o run_tests $(TEST_SRC) $(BIOP_OBJ) -lcheck $(XML_OPT) $(XML_LIB)
