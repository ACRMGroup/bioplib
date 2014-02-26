# Simple makefile for building Bioplib tests

# Define C compiler
CC = gcc

# Test code
TEST_SRC = src/*.c
#TEST_OBJ = src/*.o

# Bioplib object files
BIOP_OBJ = ../*.o

# Compile tests            
tests : 
	$(CC) -o run_tests $(TEST_SRC) $(BIOP_OBJ) -lcheck
#clean : 
#	rm -f run_tests $(TEST_OBJ)
