# Makefile for qfs

# Compiler
CC=gfortran
# Object flags
OFLAGS=-o 
# Other flags
FLAGS=-g -Wall -fdefault-real-8
# FLAGS=-g -Wall
# Executable
EXECUTABLE=raveshift
# Source 
SOURCES=ravenshift4.f
# Objects
OBJECTS=ravenshift4.o

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(FLAGS) $(SOURCES) $(OFLAGS) $(EXECUTABLE)

clean: 
	rm -rf $(EXECUTABLE) fort.*  
