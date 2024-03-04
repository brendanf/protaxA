#include <time.h>
#include "defs.h"

int main (int argc, char **argv) {
  InputOptions iopt;
  char *rfile, *ifile;
  int i,j, num_tnodes, num_sclevels;
  SequenceSet *rseq,*iseq;
  TaxonomyNode *taxonomy;
  Model *model;
  double pth, **scs;
  double *pdistances;
  clock_t start_time, now_time;

  iopt = get_input_options(argc, argv);

  if (argc - optind != 7) {
    fprintf(stderr,"classify sequences and measure timing (sequences represented as character strings)\n");
    fprintf(stderr,"usage: classify_v1 [-l len] [-r n_rseq] [-i n_iseq] taxonomy rseqFASTA taxonomy2rseq modelparameters scalingfile probability_threshold inputFASTA\n");
    exit(0);
  }

  taxonomy = read_taxonomy(argv[optind++], &num_tnodes);
  rfile = argv[optind++];
  add_rseq2taxonomy(argv[optind++], taxonomy);
  model = read_model(argv[optind++]);
  scs=read_level_scalings(argv[optind++], &num_sclevels);

  if (model->num_levels != num_sclevels) {
    fprintf(stderr,"ERROR: %d model levels but %d scaling levels, files '%s' and '%s'.\n", model->num_levels, num_sclevels, argv[4], argv[5]);
    exit(0);
  }

  pth = atof(argv[optind++]);
  ifile = argv[optind++];

  read_sequence_sets(iopt, rfile, ifile, &rseq, &iseq);

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
    fprintf(stderr, "distances: %f seconds\n",(double) (now_time - start_time) / CLOCKS_PER_SEC);
    start_time = clock();
    compute_cnode_probs_best2(iseq->id[i], taxonomy, 0, 1.0, model, (const double **)scs, pth, pdistances);
    now_time = clock();
    fprintf(stderr, "classification: %f seconds\n",(double) (now_time - start_time) / CLOCKS_PER_SEC);
  }

  return(0);
}
