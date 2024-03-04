#include "defs.h"
#include <time.h>

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
  int i, mini, rlen, ilen, n_rseq, n_iseq;
  SequenceSetB *rseq,*iseq;
  double *pdistances;
  clock_t start_time, now_time;
  
  if (argc < 3) {
    fprintf(stderr,"dist_test: calculate all pdistances between input and reference sequences and print the best for each input sequence\n");
    fprintf(stderr,"usage: dist_test rseqFASTA inputFASTA\n");
    exit(0);	    
  }
  scan_aligned_sequences(argv[1], &rlen, &n_rseq);
  rseq = read_aligned_sequencesB(argv[1], rlen, n_rseq);
  scan_aligned_sequences(argv[2], &ilen, &n_iseq);
  iseq = read_aligned_sequencesB(argv[2], ilen, n_iseq);
  
  if (rseq->alen != iseq->alen) {
    fprintf(stderr,"ERROR: sequence lengths different in two files (%d,%d), files '%s','%s'.\n",rseq->alen,iseq->alen,argv[1],argv[2]);
    exit(0);
  }

  if ((pdistances = (double *) malloc(rseq->num_seqs * sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR: cannot maloc %d doubles for pdistances.\n",rseq->num_seqs);
    perror(""); exit(-1);
  }

  start_time = clock();
  for (i=0; i<iseq->num_seqs; i++) {
    /* compute all distances between input sequence i and reference sequences */
    compute_distancesB(rseq, iseq->b[i], iseq->m[i], pdistances);

    /* simple application: get smallest distance (neglecting the possible ties..) */
    mini = get_mini(pdistances, rseq->num_seqs);
    printf("%s %s %f\n",iseq->id[i], rseq->id[mini], pdistances[mini]);

  }
  now_time = clock();
  fprintf(stderr,"timing: %f seconds\n",(double) (now_time - start_time) / CLOCKS_PER_SEC);

  
  return(0);
}
