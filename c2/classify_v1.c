#include <time.h>
#include "defs.h"

int main (int argc, char **argv) {
  int i,j, num_tnodes, num_sclevels;
  SequenceSet *rseq,*iseq;
  TaxonomyNode *taxonomy;
  Model *model;
  double pth, **scs;
  double *pdistances;
  clock_t start_time, now_time;
  
  if (argc < 7) {
    fprintf(stderr,"classify sequences and measure timing (sequences represented as character strings)\n");
    fprintf(stderr,"usage: classify taxonomy rseqFASTA taxonomy2rseq modelparameters scalingfile probability_threshold inputFASTA\n");
    exit(0);	    
  }

  taxonomy = read_taxonomy(argv[1], &num_tnodes);
  rseq = read_aligned_sequences(argv[2]);
  add_rseq2taxonomy(argv[3], taxonomy);
  model = read_model(argv[4]);
  scs=read_level_scalings(argv[5], &num_sclevels);

  if (model->num_levels != num_sclevels) {
    fprintf(stderr,"ERROR: %d model levels but %d scaling levels, files '%s' and '%s'.\n", model->num_levels, num_sclevels, argv[4], argv[5]);
    exit(0);
  }
  
  pth = atof(argv[6]);  
  iseq = read_aligned_sequences(argv[7]);
  
  if (rseq->alen != iseq->alen) {
    fprintf(stderr,"ERROR: sequence lengths different in two files (%d,%d), files '%s','%s'.\n",rseq->alen,iseq->alen,argv[2],argv[6]);
    exit(0);
  }

  /*
    print_taxonomy(taxonomy, num_tnodes);
    print_model(model);
  */

  if ((pdistances = (double *) malloc(rseq->num_seqs * sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR: cannot maloc %d doubles for pdistances.\n",rseq->num_seqs);
    perror(""); exit(-1);
  }
  
  for (i=0; i<iseq->num_seqs; i++) {
    start_time = clock();
    compute_distances(rseq, iseq->seq[i], pdistances);
    now_time = clock();
    printf("distances: %f seconds\n",(double) (now_time - start_time) / CLOCKS_PER_SEC);
    start_time = clock();
    printf("%s",iseq->id[i]);
    compute_cnode_probs_best2(taxonomy, 0, 1.0, model, scs, pth, pdistances);
    printf("\n");
    now_time = clock();
    printf("classification: %f seconds\n",(double) (now_time - start_time) / CLOCKS_PER_SEC);
  }
  
  return(0);
}
