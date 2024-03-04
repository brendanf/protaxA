#include "defs.h"

InputOptions get_input_options(const int argc, char * const * argv) {
  return get_input_options_custom(argc, argv, ":l:r:i:");
}

InputOptions get_input_options_custom(const int argc, char * const * argv,
                                      const char * options) {
  int opt;
  InputOptions iopt = {0, 0, 0};

  while ((opt = getopt(argc, argv, options)) > 0) {
    switch(opt) {
    case 'l':
      iopt.len = atoi(optarg);
      if (iopt.len <= 0) {
        fprintf(stderr,"ERROR: sequence length must be a positive integer: %s", optarg);
        exit(-1);
      }
      fprintf(stderr, "Using user supplied sequence length: %d\n", iopt.len);
      break;
    case 'r':
      iopt.n_rseq = atoi(optarg);
      if (iopt.n_rseq <= 0) {
        fprintf(stderr,"ERROR: count of reference sequences must be a positive integer: %s", optarg);
        exit(-1);
      }
      fprintf(stderr, "Using user supplied count of reference sequences: %d\n", iopt.n_rseq);
      break;
    case 'i':
      iopt.n_iseq = atoi(optarg);
      if (iopt.n_iseq <= 0) {
        fprintf(stderr,"ERROR: count of input sequences must be a positive integer: %s", optarg);
        exit(-1);
      }
      fprintf(stderr, "Using user supplied count of input sequences: %d\n", iopt.n_iseq);
      break;
    case ':':
      fprintf(stderr, "Option -%c requires an operand\n", optopt);
      exit(-1);
      break;
    default:
      fprintf(stderr,"ERROR: unknown option: -%c\n", optopt);
      exit(-1);
    }
  }

  return iopt;
}
