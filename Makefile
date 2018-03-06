CC = gcc
CFLAGS  = -g -O3 -Wall -pedantic # -fsanitize=address
LOADLIBES = -lm -lgsl -lgslcblas # -lasan <-- has to be the first flag
VPATH = ./source

#INCLUDES = ranlxs.h

OBJS = static_quark_potential.o rand_su2.o

static_quark_potential: ${OBJS}

${OBJS}: ${INCLUDES} Makefile
