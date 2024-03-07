#include <time.h>
#include "defs.h"

int get_mini(double *pdistances, int n) {
  int i, mini;
  double best;

  mini = 0;
  best = pdistances[mini];

  for (i=0; i<n; i++) {
    if (pdistances[i] < best) {
      mini = i;
      best = pdistances[i];
    }
  }

  return (mini);
}

int main (int argc, char **argv) {
  InputOptions iopt;
  int i, mini;
  SequenceSetB *rseq,*iseq;
  char *rfile, *ifile;
  double *pdistances;
  clock_t start_time, now_time;

  setvbuf(stdout, NULL, _IOLBF, 0);

  iopt = get_input_options(argc, argv);

  if (argc - optind != 2) {
    fprintf(stderr,"dist_best: calculate all pdistances between input and reference sequences and print the best for each input sequence\n");
    fprintf(stderr,"usage: dist_best [-l len] [-r n_rseq] [-i n_iseq] rseqFASTA inputFASTA\n");
    exit(0);
  }

  rfile = argv[optind++];
  ifile = argv[optind++];

  read_sequence_setsB(iopt, rfile, ifile, &rseq, &iseq);

  if ((pdistances = (double *) malloc(rseq->num_seqs * sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR: cannot malloc %d doubles for pdistances.\n",rseq->num_seqs);
    perror(""); exit(-1);
  }

  start_time = clock();
  for (i=0; i<iseq->num_seqs; i++) {
    /* compute all distances between input sequence i and reference sequences */
    compute_distancesB(rseq, iseq->b[i], iseq->m[i], iseq->start[i], iseq->end[i], pdistances);

    /* simple application: get smallest distance (neglecting the possible ties..) */
    mini = get_mini(pdistances, rseq->num_seqs);
    printf("%s %s %f\n",iseq->id[i], rseq->id[mini], pdistances[mini]);
  }
  now_time = clock();
  fprintf(stderr,"timing: %f seconds\n",(double) (now_time - start_time) / CLOCKS_PER_SEC);

  return(0);
}
