#include "defs.h"

int main (int argc, char **argv) {
  InputOptions iopt;
  char *thresh_param, *thresh_result, *ifile;
  int i, j;
  SequenceSetB *iseq;
  double thresh, d;

  iopt = get_input_options_custom(argc, argv, ":l:i:");

  if (argc - optind != 2) {
    fprintf(stderr,"dist_matrix: calculate all pdistances between input sequences\n");
    fprintf(stderr,"usage: dist_matrix [-l len] [-i n_iseq] threshold inputFASTA\n");
    exit(-1);
  }

  thresh_param = argv[optind++];
  thresh = strtod(thresh_param, &thresh_result);
  if (thresh_param == thresh_result) {
    fprintf(stderr, "ERROR: could not interpret argument '%s' as a number.\n", thresh_param);
    exit(-1);
  }
  if (thresh < 0.0 || thresh > 1.0) {
    fprintf(stderr, "ERROR: threshold must be in the range [0,1]. value: %f\n", thresh);
    exit(-1);
  }

  ifile = argv[optind++];

  if (iopt.len * iopt.n_iseq == 0) {
    scan_aligned_sequences(ifile, &iopt.len, &iopt.n_iseq);
  }

  iseq = read_aligned_sequencesB(ifile, iopt.len, iopt.n_iseq);

  for (i=1; i<iseq->num_seqs; i++) {
    for (j=0; j<i; j++) {
      d = pdistB(iseq->b[i], iseq->m[i], iseq->b[j], iseq->m[j], iseq->ulen, iseq->mulen);
      if (d <= thresh) printf("%s %s %f\n",iseq->id[i], iseq->id[j], d);
    }
  }

  return(0);
}
