#ifndef _DEFS_
#define _DEFS_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define NUCLEOTIDES_IN_WORD 16
#define MAXLINE 1024
#define UNKNAME "unk"

typedef struct {
  double rth;
  int len, n_rseq, n_iseq;
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
} SequenceSet;

typedef struct {
  int num_seqs, alen, ulen, mulen;
  char **id;
  long unsigned int **b, **m;
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
double pdist(const char *a, const char *b, const int len);
int compute_distances(const SequenceSet *a, const char *seq, double *pdistances);

int nucleotide2binary(const char *s, const int n, long unsigned int *b, long unsigned int *m);
SequenceSetB *read_aligned_sequencesB(const char *filename, const int len, const int num_seqs);
void read_sequence_setsB(InputOptions iopt, const char * rfile, const char * ifile,
                         SequenceSetB **rseq, SequenceSetB **iseq);
double pdistB(const long unsigned int *a, const long unsigned int *ma,
              const long unsigned int *b, const long unsigned int *mb,
              const int n, const int n2);
int compute_distancesB(const SequenceSetB *a,
                       const long unsigned int *b, const long unsigned int *m,
                       double *pdistances);

/* routines_model.c */

Model *read_model(char *filename);
double **read_level_scalings(char *filename, int *num_levels);
int compute_cnode_probs_best2(const char *qid, TaxonomyNode *node, int nid, double prevprob, const Model *m, const double **scs, double pth, double rth, const double *pdistances);
int print_model(Model *m);


#endif
#ifndef _DEFS_
#define _DEFS_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NUCLEOTIDES_IN_WORD 16
#define MAXLINE 1024
#define UNKNAME "unk"

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
} SequenceSet;

typedef struct {
  int num_seqs, alen, ulen, mulen;
  char **id;
  long unsigned int **b, **m;
} SequenceSetB;

typedef struct {
  int num_levels, dim;
  double **params;
} Model;

/* routines_taxonomy.c */

TaxonomyNode *read_taxonomy(char *filename, int *num_nodes);
int add_rseq2taxonomy(char *filename, TaxonomyNode *node);
int print_taxonomy(TaxonomyNode *node, int num_nodes);

/* routines_sequence.c */

SequenceSet *read_aligned_sequences(char *filename);
double pdist(char *a, char *b, int len);
int compute_distances(SequenceSet *a, char *seq);
int nucleotide2binary(char *s, int n, long unsigned int *b, long unsigned int *m);
SequenceSetB *read_aligned_sequencesB(char *filename);
double pdistB(long unsigned int *a, long unsigned int *ma, long unsigned int *b, long unsigned int *mb, int n, int n2);
int compute_distancesB(SequenceSetB *a, long unsigned int *b, long unsigned int *m, double *pdistances);

/* routines_model.c */

Model *read_model(char *filename);
double **read_level_scalings(char *filename, int *num_levels);
int compute_cnode_probs_best2(TaxonomyNode *node, int nid, double prevprob, Model *m, double **scs, double pth, double *pdistances);

int print_model(Model *m);


#endif
