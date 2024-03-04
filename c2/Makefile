CC = gcc
CFLAGS = -O2
LDFLAGS =
LIBS = -lm -lz

PROGS = dist_best classify_v1 classify_v2 classify_info classify_rseq dist_bipart dist_matrix
PROGS_SRCS = $(PROGS:=.c)

ROUTINES_SRCS = routines_taxonomy.c routines_sequence.c routines_model.c routines_cli.c
ROUTINES_OBJS = $(ROUTINES_SRCS:.c=.o)

DEFS = defs.h

all: $(PROGS)

objs: $(PROGS_SRCS:.c=.o) $(ROUTINES_OBJS)

%.o: %.c $(DEFS)
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS) : %: %.o $(ROUTINES_OBJS) $(DEFS)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ $(ROUTINES_OBJS) $(LIBS)

