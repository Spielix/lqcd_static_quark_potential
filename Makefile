CC = gcc
CFLAGS  = -O3 -Wall -pedantic # -g -fsanitize=address
LOADLIBES = -lm -lgsl -lgslcblas # -lasan <-- has to be the first flag
VPATH = ./source

#INCLUDES = ranlxs.h

OBJS = static_quark_potential.o rand_su2.o data_connect.o

static_quark_potential: static_quark_potential.o rand_su2.o # ${OBJS}

data_connect: data_connect.o

${OBJS}: ${INCLUDES} Makefile
