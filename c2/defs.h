#ifndef _DEFS_
#define _DEFS_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define NUCLEOTIDES_IN_WORD 16
#define MAXLINE 10240
#define UNKNAME "unk"

typedef struct {
  double rth;
  int len, n_rseq, n_iseq, min_len;
} InputOptions;

typedef struct {
  int nid, pid, level;
  char *name, *mxcode;
  int num_cnodes;
  int *cnode_index;
  int num_rseqs;
  int *rseq_index;
  int isunk, no_rseqs;
  double prior, prob, sumcprob_no_rseqs;
  int ind1, ind2;
  double dist1, dist2;
} TaxonomyNode;

typedef struct {
  int num_seqs, alen;
  char **seq, **id;
  int *start, *end;
} SequenceSet;

typedef struct {
  int num_seqs, alen, ulen, mulen;
  char **id;
  long unsigned int **b, **m;
  int *start, *end;
} SequenceSetB;

typedef struct {
  int num_levels, dim;
  double **params;
} Model;

/* routines_taxonomy.c */

TaxonomyNode *read_taxonomy(char *filename, int *num_nodes);
int add_rseq2taxonomy(char *filename, TaxonomyNode *node);
int print_taxonomy(TaxonomyNode *node, int num_nodes);

/* routines_cli.c */

InputOptions get_input_options(const int argc, char * const * argv);
InputOptions get_input_options_custom(const int argc, char * const * argv, const char * options);

/* routines_sequence.c */

void scan_aligned_sequences(const char *filename, int *len, int *num_seqs);

SequenceSet *read_aligned_sequences(const char *filename, const int len, const int num_seqs);
void read_sequence_sets(InputOptions iopt, const char * rfile, const char * ifile,
                         SequenceSet **rseq, SequenceSet **iseq);
double pdist(const char *a, const char *b, const int start, const int end);
int compute_distances(const SequenceSet *a, const char *seq,
                      const int start, const int end,
                      const int min_len,
                      double *pdistances);

int nucleotide2binary(const char *s, const int n, long unsigned int *b, long unsigned int *m, int *start, int *end);
SequenceSetB *read_aligned_sequencesB(const char *filename, const int len, const int num_seqs);
void read_sequence_setsB(InputOptions iopt, const char * rfile, const char * ifile,
                         SequenceSetB **rseq, SequenceSetB **iseq);
double pdistB(const long unsigned int *a, const long unsigned int *ma,
              const long unsigned int *b, const long unsigned int *mb,
              const int start, const int end);
int compute_distancesB(const SequenceSetB *a,
                       const long unsigned int *b, const long unsigned int *m,
                       const int start, const int end,
                       const int min_len,
                       double *pdistances);

/* routines_model.c */

Model *read_model(char *filename);
double **read_level_scalings(char *filename, int *num_levels);
int compute_cnode_probs_best2(const char *qid, TaxonomyNode *node, int nid, double prevprob, const Model *m, const double **scs, double pth, double rth, const double *pdistances);
int print_model(Model *m);


#endif
