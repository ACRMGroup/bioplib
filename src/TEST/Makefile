# Simple makefile for building Bioplib tests

# Define C compiler
CC = gcc

# Test code
TEST_SRC = src/*.c
#TEST_OBJ = src/*.o

# Bioplib object files
BIOP_OBJ = ../*.o

# Include libxml2
# Note: xml2-config is installed with libxml2
#       Link to libxml2 with -lxml2
XML_OPT := $(shell xml2-config --cflags)


# Suppress error messages for deprecated functions.
# Uncomment compile options below to either display the error message or 
# silence the error message:
#DEP_OPT = -D BIOPLIB_DEPRECATED_CHECK
DEP_OPT = -D BIOPLIB_DEPRECATED_QUIET


# Compile tests            
tests : 
	$(CC) -o run_tests $(TEST_SRC) $(BIOP_OBJ) $(XML_OPT) $(DEP_OPT) -lcheck -lxml2
#clean : 
#	rm -f run_tests $(TEST_OBJ)
