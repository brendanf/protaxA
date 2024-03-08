#include "defs.h"

int main (int argc, char **argv) {
  InputOptions iopt;
  char *thresh_param, *thresh_result, *ifile, *rfile;
  int i, j, start, end;
  SequenceSetB *rseq,*iseq;
  double thresh, d;

  setvbuf(stdout, NULL, _IOLBF, 0);

  iopt = get_input_options(argc, argv);

  if (argc - optind != 3) {
    fprintf(stderr,"dist_bipart: calculate all pdistances between input and reference sequences\n");
    fprintf(stderr,"usage: dist_matrix [-l len] [-r n_rseq] [-i n_iseq] [-m min_overlap] threshold refFASTA inputFASTA\n");
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

  rfile = argv[optind++];
  ifile = argv[optind++];

  read_sequence_setsB(iopt, rfile, ifile, &rseq, &iseq);

  for (i=0; i<rseq->num_seqs; i++) {
    if (rseq->end[i] - rseq->start[i] < iopt.min_len) continue;
    for (j=0; j<iseq->num_seqs; j++) {
      start = rseq->start[i] > iseq->start[j] ? rseq->start[i] : iseq->start[j];
      end = rseq->end[i] < iseq->end[j] ? rseq->end[i] : iseq->end[j];
      if (end - start < iopt.min_len) continue;
      d = pdistB(rseq->b[i], rseq->m[i], iseq->b[j], iseq->m[j], start, end);
      if (d <= thresh) printf("%s %s %f\n",rseq->id[i], iseq->id[j], d);
    }
  }

  return(0);
}
