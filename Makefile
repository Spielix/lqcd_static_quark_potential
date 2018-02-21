CC = gcc
CFLAGS  = -g -O4 -Wall -pedantic 
LOADLIBES = -lm -lgsl -lgslcblas
VPATH = ./source

#INCLUDES = ranlxs.h

OBJS = static_quark_potential.o

static_quark_potential: ${OBJS}

${OBJS}: ${INCLUDES} Makefile
